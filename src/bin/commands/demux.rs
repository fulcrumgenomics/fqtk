use crate::commands::command::Command;
use anyhow::{anyhow, Result};
use clap::Parser;
use fgoxide::io::Io;
use fqtk_lib::barcode_matching::BarcodeMatcher;
use fqtk_lib::samples::SampleGroup;
use log::info;
use proglog::ProgLogBuilder;
use read_structure::ReadStructure;
use read_structure::ReadStructureError;
use read_structure::SegmentType;
use seq_io::fastq::write_to;
use seq_io::fastq::OwnedRecord;
use seq_io::fastq::Reader as FastqReader;
use seq_io::fastq::Record;
use seq_io::fastq::RecordsIntoIter;
use std::collections::HashSet;
use std::io::{BufWriter, Read, Write};
use std::{
    fs,
    io::BufReader,
    path::{Path, PathBuf},
};

type VecOfReaders = Vec<BufReader<Box<dyn Read>>>;

/// The bases and qualities associated with a segment of a FASTQ record.
#[derive(Debug, Clone)]
struct FastqSegment {
    /// bases of the FASTQ subsection
    seq: Vec<u8>,
    /// qualities of the FASTQ subsection
    quals: Vec<u8>,
}

/// One unit of FASTQ records separated into their component read segments.
#[derive(Debug, Clone)]
struct ReadSet {
    /// Name of the FASTQ record
    name: Vec<u8>,
    /// The full sample barcode for this record (concatenated bases from individual sample barcode
    /// segments)
    sample_barcode_sequence: Vec<u8>,
    /// The template segments for this set of reads.
    template_segments: Option<Vec<FastqSegment>>,
    /// The sample barcode segments ofr this set of reads (if being output, otherwise None).
    sample_barcode_segments: Option<Vec<FastqSegment>>,
    /// The molecular barcode segments of this set of reads (if being ouput, otherwise None).
    molecular_barcode_segments: Option<Vec<FastqSegment>>,
}

/// A struct for iterating over the records in multiple FASTQ files simultaneously, destructuring
/// them according to the provided read structures, yielding ``ReadSet`` objects on each iteration.
struct ReadSetIterator {
    /// Read structure objects describing the structure of the reads to be demultiplexed, one per
    /// input file.
    read_structures: Vec<ReadStructure>,
    /// number of template segments per set of FASTQ records
    num_templates: usize,
    /// number of molecular barcode segments per set of FASTQ records
    num_molecular_barcodes: usize,
    /// number of sample barcode segments per set of FASTQ records
    num_sample_barcodes: usize,
    /// Iterators over the files containing FASTQ records, one per input file.
    sources: Vec<RecordsIntoIter<BufReader<Box<dyn Read>>>>,
    /// If true this iterator will populate the template segments field of the ReadSet object it
    /// returns
    generate_template_segments: bool,
    /// If true this iterator will populate the sample barcode segments field of the ReadSet
    /// objects it returns
    generate_sample_barcode_segments: bool,
    /// If true this iterator will populate the molecular barcode segments field of the ReadSet
    /// objects it returns
    generate_molecular_barcode_segments: bool,
}

impl Iterator for ReadSetIterator {
    type Item = ReadSet;

    fn next(&mut self) -> Option<Self::Item> {
        let mut rec_results = Vec::with_capacity(self.sources.len());
        for source in &mut self.sources {
            if let Some(rec) = source.next() {
                rec_results.push(rec);
            }
        }

        if rec_results.is_empty() {
            return None;
        }
        assert!(
            rec_results.len() == self.sources.len(),
            "FASTQ sources out of sync at records: {:?}",
            rec_results
        );

        let next_sources = rec_results
            .into_iter()
            .collect::<Result<Vec<OwnedRecord>, seq_io::fastq::Error>>()
            .expect("Unexpected error parsing FASTQs");

        let read_name = next_sources[0].head.clone();

        let mut sample_barcode = Vec::new();
        let mut template_segments = if self.generate_template_segments {
            Some(Vec::with_capacity(self.num_templates))
        } else {
            None
        };
        let mut sample_barcode_segments = if self.generate_sample_barcode_segments {
            Some(Vec::with_capacity(self.num_sample_barcodes))
        } else {
            None
        };
        let mut molecular_barcode_segments = if self.generate_molecular_barcode_segments {
            Some(Vec::with_capacity(self.num_molecular_barcodes))
        } else {
            None
        };

        for (read_structure, fastq_record) in self.read_structures.iter().zip(next_sources.iter()) {
            for read_segment in read_structure.iter() {
                let (seq, quals) = read_segment
                    .extract_bases_and_quals(fastq_record.seq(), fastq_record.qual())
                    .unwrap_or_else(|e| {
                        panic!(
                            "Error extracting read segment bases or quals from FASTQ record {}; {}",
                            String::from_utf8(fastq_record.head().to_vec()).unwrap(),
                            e
                        )
                    });
                match read_segment.kind {
                    SegmentType::Template => {
                        if let Some(ref mut templates) = template_segments {
                            templates
                                .push(FastqSegment { seq: seq.to_vec(), quals: quals.to_vec() });
                        }
                    }
                    SegmentType::SampleBarcode => {
                        if let Some(ref mut sample_segments) = sample_barcode_segments {
                            sample_segments
                                .push(FastqSegment { seq: seq.to_vec(), quals: quals.to_vec() });
                        }
                        sample_barcode.extend(seq.iter());
                    }
                    SegmentType::MolecularBarcode => {
                        if let Some(ref mut mol_segments) = molecular_barcode_segments {
                            mol_segments
                                .push(FastqSegment { seq: seq.to_vec(), quals: quals.to_vec() });
                        }
                    }
                    _ => continue,
                }
            }
        }

        Some(ReadSet {
            name: read_name,
            sample_barcode_sequence: sample_barcode,
            template_segments,
            sample_barcode_segments,
            molecular_barcode_segments,
        })
    }
}

struct SampleWriters {
    /// Vec of the writers for template read segments.
    template_writers: Option<Vec<BufWriter<Box<dyn Write>>>>,
    /// Vec of the writers for the sample barcode read segments.
    sample_barcode_writers: Option<Vec<BufWriter<Box<dyn Write>>>>,
    /// Vec of the writers for the molecular barcode read segments.
    molecular_barcode_writers: Option<Vec<BufWriter<Box<dyn Write>>>>,
}

impl SampleWriters {
    fn write(&mut self, read_set: &ReadSet) -> Result<()> {
        for (writers_opt, segments_opt) in [
            (&mut self.template_writers, &read_set.template_segments),
            (&mut self.sample_barcode_writers, &read_set.sample_barcode_segments),
            (&mut self.molecular_barcode_writers, &read_set.molecular_barcode_segments),
        ] {
            if let (Some(writers), Some(segments)) = (writers_opt, segments_opt) {
                for (writer, segment) in writers.iter_mut().zip(segments.iter()) {
                    write_to(writer, &read_set.name, &segment.seq, &segment.quals)?;
                }
            }
        }
        Ok(())
    }

    fn num_template_writers(&self) -> usize {
        if let Some(v) = &self.template_writers {
            v.len()
        } else {
            0
        }
    }
    fn num_sample_writers(&self) -> usize {
        if let Some(v) = &self.sample_barcode_writers {
            v.len()
        } else {
            0
        }
    }
    fn num_molecular_writers(&self) -> usize {
        if let Some(v) = &self.molecular_barcode_writers {
            v.len()
        } else {
            0
        }
    }
}

/// Performs sample demultiplexing on FASTQs.
///
/// The sample barcode for each sample in the metadata TSV will be compared against the sample
/// barcode bases extracted from the FASTQs, to assign each read to a sample.  Reads that do not
/// match any sample within the given error tolerance will be placed in the ``unmatched_prefix``
/// file.
///
/// FASTQs and associated read structures for each sub-read should be given:
///
/// - a single fragment read (with inline index) should have one FASTQ and one read structure
/// - paired end reads should have two FASTQs and two read structures
/// - a dual-index sample with paired end reads should have four FASTQs and four read structures
///   given: two for the two index reads, and two for the template reads.
///
/// If multiple FASTQs are present for each sub-read, then the FASTQs for each sub-read should be
/// concatenated together prior to running this tool
/// (e.g. `zcat s_R1_L001.fq.gz s_R1_L002.fq.gz | bgzip -c > s_R1.fq.gz`).
///
/// (Read structures)[<https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures>] are made up of
/// `<number><operator>` pairs much like the `CIGAR` string in BAM files.
/// Four kinds of operators are recognized:
///
/// 1. `T` identifies a template read
/// 2. `B` identifies a sample barcode read
/// 3. `M` identifies a unique molecular index read
/// 4. `S` identifies a set of bases that should be skipped or ignored
///
/// The last `<number><operator>` pair may be specified using a `+` sign instead of number to
/// denote "all remaining bases". This is useful if, e.g., fastqs have been trimmed and contain
/// reads of varying length. Both reads must have template bases.  Any molecular identifiers will
/// be concatenated using the `-` delimiter and placed in the given SAM record tag (`RX` by
/// default).  Similarly, the sample barcode bases from the given read will be placed in the `BC`
/// tag.
///
/// Metadata about the samples should be given as a headered metadata TSV file with two columns
/// 1. name - the name of the sample.
/// 2. barcode - the expected barcode sequence for that sample.
///
/// The read structures will be used to extract the observed sample barcode, template bases, and
/// molecular identifiers from each read.  The observed sample barcode will be matched to the
/// sample barcodes extracted from the bases in the sample metadata and associated read structures.
/// ## Example Command Line
///
/// As an example, if the sequencing run was 2x100bp (paired end) with two 8bp index reads both
/// reading a sample barcode, as well as an in-line 8bp sample barcode in read one, the command
/// line would be
///
/// ```
/// fqtk demux \
///     --inputs r1.fq i1.fq i2.fq r2.fq \
///     --read-structures 8B92T 8B 8B 100T \
///     --sample-metadata metadata.tsv \
///     --output output_folder
/// ```
///
#[derive(Parser, Debug)]
pub(crate) struct Demux {
    /// One or more input fastq files each corresponding to a sequencing (e.g. R1, I1).
    #[clap(long, short = 'i', required = true, num_args = 1..)]
    inputs: Vec<PathBuf>,

    /// The read structures, one per input FASTQ in the same order.
    #[clap(long, short = 'r' , required = true, num_args = 1..)]
    read_structures: Vec<ReadStructure>,

    /// The read structure types to write to their own files (Must be one of T, B, or M for
    /// template reads, sample barcode reads, and molecular barcode reads).
    #[clap(long, short='b', default_value="T", num_args = 1.. )]
    output_types: Vec<char>,

    /// A file containing the metadata about the samples.
    #[clap(long, short = 's', required = true)]
    sample_metadata: PathBuf,

    /// The output directory into which to write per-sample FASTQs.
    #[clap(long, short = 'o', required = true)]
    output: PathBuf,

    /// Output prefix for FASTQ file(s) for reads that cannot be matched to a sample.
    #[clap(long, short = 'u', default_value = "unmatched")]
    unmatched_prefix: String,

    /// Maximum mismatches for a barcode to be considered a match.
    #[clap(long, default_value = "1")]
    max_mismatches: usize,

    /// Minimum difference between number of mismatches in the best and second best barcodes for a
    /// barcode to be considered a match.
    #[clap(long, short = 'd', default_value = "2")]
    min_mismatch_delta: usize,

    /// The number of threads to use. Cannot be less than 3.
    #[clap(long, short = 't', default_value = "3")]
    threads: usize,
}

impl Demux {
    /// Creates one writer per read segment in the read structures on this object, restricted by
    /// requested output type.
    /// # Errors:
    ///     - Will error if opening the output fails for any reason.
    fn create_sample_writers(
        read_structures: &[ReadStructure],
        prefix: &str,
        output_types: &HashSet<SegmentType>,
        output_dir: &Path,
    ) -> Result<SampleWriters> {
        let fg_io = Io::default();
        let mut template_writers = None;
        let mut sample_barcode_writers = None;
        let mut molecular_barcode_writers = None;
        for output_type in output_types {
            let mut output_type_writers = Vec::new();

            let file_type_code = match output_type {
                SegmentType::Template => 'R',
                SegmentType::SampleBarcode => 'I',
                SegmentType::MolecularBarcode => 'U',
                _ => 'S',
            };

            let segment_count: usize =
                read_structures.iter().map(|s| s.segments_by_type(*output_type).count()).sum();

            for idx in 1..=segment_count {
                output_type_writers.push(fg_io.new_writer(
                    &output_dir.join(format!("{}.{}{}.fq.gz", prefix, file_type_code, idx)),
                )?);
            }
            match output_type {
                SegmentType::Template => template_writers = Some(output_type_writers),
                SegmentType::SampleBarcode => {
                    sample_barcode_writers = Some(output_type_writers);
                }
                SegmentType::MolecularBarcode => {
                    molecular_barcode_writers = Some(output_type_writers);
                }
                _ => {}
            }
        }

        Ok(SampleWriters { template_writers, sample_barcode_writers, molecular_barcode_writers })
    }

    /// Creates one writer per sample per read segment in the provided read structures for
    /// requested output type.
    /// # Errors:
    ///     - Will error if opening the output fails for any reason.
    fn create_writers(
        read_structures: &[ReadStructure],
        sample_group: &SampleGroup,
        output_types: &HashSet<SegmentType>,
        output_dir: &Path,
    ) -> Result<Vec<SampleWriters>> {
        let mut samples_writers = Vec::with_capacity(sample_group.samples.len());
        for sample in &sample_group.samples {
            samples_writers.push(Demux::create_sample_writers(
                read_structures,
                &sample.name,
                output_types,
                output_dir,
            )?);
        }
        Ok(samples_writers)
    }
    /// Checks that inputs to demux are valid and returns open file handles for the inputs.
    /// Checks:
    ///     - That the number of input files and number of read structs provided are the same
    ///     - That the output directory is not read-only
    ///     - That the input files exist
    ///     - That the input files have read permissions.
    fn validate_and_prepare_inputs(&self) -> Result<(VecOfReaders, HashSet<SegmentType>)> {
        let mut constraint_errors = vec![];

        if self.inputs.len() != self.read_structures.len() {
            let preamble = "The same number of read structures should be given as FASTQs";
            let specifics = format!(
                "{} read-structures provided for {} FASTQs",
                self.read_structures.len(),
                self.inputs.len()
            );
            constraint_errors.push(format!("{preamble} {specifics}"));
        }

        if !self.output.exists() {
            info!("Output directory {:#?} didn't exist, creating it.", self.output);
            fs::create_dir_all(&self.output)?;
        }

        if self.output.metadata()?.permissions().readonly() {
            constraint_errors
                .push(format!("Ouput directory {:#?} cannot be read-only", self.output));
        }

        let output_segment_types_result = self
            .output_types
            .iter()
            .map(|&c| SegmentType::try_from(c))
            .collect::<Result<HashSet<_>, ReadStructureError>>();
        if let Err(e) = &output_segment_types_result {
            constraint_errors.push(format!("Error parsing segment types to report: {}", e));
        }

        for input in &self.inputs {
            if !input.exists() {
                constraint_errors.push(format!("Provided input file {:#?} doesn't exist", input));
            }
        }
        // Attempt to open the files for reading.
        let io = Io::default();
        let fq_readers_result = self
            .inputs
            .iter()
            .map(|p| io.new_reader(p))
            .collect::<Result<Vec<_>, fgoxide::FgError>>();
        if let Err(e) = &fq_readers_result {
            constraint_errors.push(format!("Error opening input files for reading: {}", e));
        }

        if self.threads < 3 {
            constraint_errors
                .push(format!("Threads provided {} was too low! Must be 3 or more.", self.threads));
        }

        if constraint_errors.is_empty() {
            let output_segment_types = output_segment_types_result?;
            if output_segment_types.is_empty() {
                constraint_errors.push(
                    "No output types requested, must request at least one output segment type."
                        .to_owned(),
                );
            } else {
                return Ok((fq_readers_result?, output_segment_types));
            }
        }
        let mut details = "Inputs failed validation!\n".to_owned();
        for error_reason in constraint_errors {
            details.push_str(&format!("    - {}\n", error_reason));
        }
        Err(anyhow!("The following errors with the input(s) were detected:\n{}", details))
    }
}

impl Command for Demux {
    /// Executes the demux command
    fn execute(&self) -> Result<()> {
        let (fq_readers, output_segment_types) = self.validate_and_prepare_inputs()?;

        let sample_group = SampleGroup::from_file(&self.sample_metadata)?;
        info!(
            "{} samples loaded from file {:?}",
            sample_group.samples.len(),
            &self.sample_metadata
        );
        let fq_sources = fq_readers
            .into_iter()
            .map(|fq| FastqReader::new(fq).into_records())
            .collect::<Vec<_>>();

        let mut unassigned_writers = Demux::create_sample_writers(
            &self.read_structures,
            &self.unmatched_prefix,
            &output_segment_types,
            &self.output,
        )?;
        let mut sample_writers = Demux::create_writers(
            &self.read_structures,
            &sample_group,
            &output_segment_types,
            &self.output,
        )?;
        info!("Created sample (and unassigned) writers");

        let barcode_matcher = BarcodeMatcher::new(
            &sample_group.samples.iter().map(|s| s.barcode.as_str()).collect::<Vec<_>>(),
            u8::try_from(self.max_mismatches)?,
            u8::try_from(self.min_mismatch_delta)?,
        );

        let fq_iterator = ReadSetIterator {
            read_structures: self.read_structures.clone(),
            sources: fq_sources,
            num_templates: unassigned_writers.num_template_writers(),
            num_sample_barcodes: unassigned_writers.num_sample_writers(),
            num_molecular_barcodes: unassigned_writers.num_molecular_writers(),
            generate_template_segments: output_segment_types.contains(&SegmentType::Template),
            generate_sample_barcode_segments: output_segment_types
                .contains(&SegmentType::SampleBarcode),
            generate_molecular_barcode_segments: output_segment_types
                .contains(&SegmentType::MolecularBarcode),
        };

        let logger = ProgLogBuilder::new()
            .name("fqtk")
            .noun("records")
            .verb("demultiplexed")
            .unit(1_000_000)
            .level(log::Level::Info)
            .build();

        for read_set in fq_iterator {
            if let Some(barcode_match) = barcode_matcher.assign(&read_set.sample_barcode_sequence) {
                sample_writers[barcode_match.best_match].write(&read_set)?;
            } else {
                unassigned_writers.write(&read_set)?;
            }
            logger.record();
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fqtk_lib::samples::Sample;
    use rstest::rstest;
    use std::str;
    use std::str::FromStr;
    use tempfile::TempDir;

    const SAMPLE1_BARCODE: &str = "GATTGGG";

    /// Given a record name prefix and a slice of bases for a set of records, returns the contents
    /// of a FASTQ file as a vec of Strings, one string per line of the FASTQ.
    fn fq_lines_from_bases(prefix: &str, records_bases: &[&str]) -> Vec<String> {
        let mut result = Vec::with_capacity(records_bases.len() * 4);
        for (i, &bases) in records_bases.iter().enumerate() {
            result.push(format!("@{}_{}", prefix, i));
            result.push(bases.to_owned());
            result.push("+".to_owned());
            result.push(";".repeat(bases.len()));
        }
        result
    }

    /// Generates a FASTQ file in the tmpdir with filename "{prefix}.fastq" from the record bases
    /// specified and returns the path to the FASTQ file.
    fn fastq_file(
        tmpdir: &TempDir,
        filename_prefix: &str,
        read_prefix: &str,
        records_bases: &[&str],
    ) -> PathBuf {
        let io = Io::default();

        let path = tmpdir.path().join(format!("{filename_prefix}.fastq"));
        let fastq_lines = fq_lines_from_bases(read_prefix, records_bases);
        io.write_lines(&path, fastq_lines).unwrap();

        path
    }

    fn metadata_lines_from_barcodes(barcodes: &[&str]) -> Vec<String> {
        let mut result = Vec::with_capacity(barcodes.len() + 1);
        result.push(Sample::deserialize_header_line());
        for (i, &barcode) in barcodes.iter().enumerate() {
            result.push(format!("Sample{:04}\t{}", i, barcode));
        }
        result
    }

    fn metadata_file(tmpdir: &TempDir, barcodes: &[&str]) -> PathBuf {
        let io = Io::default();

        let path = tmpdir.path().join("metadata.tsv");
        let metadata_lines = metadata_lines_from_barcodes(barcodes);
        io.write_lines(&path, metadata_lines).unwrap();

        path
    }

    fn metadata(tmpdir: &TempDir) -> PathBuf {
        metadata_file(tmpdir, &[SAMPLE1_BARCODE])
    }

    fn read_fastq(file_path: &PathBuf) -> Vec<OwnedRecord> {
        let fg_io = Io::default();

        FastqReader::new(fg_io.new_reader(file_path).unwrap())
            .into_records()
            .collect::<Result<Vec<_>, seq_io::fastq::Error>>()
            .unwrap()
    }

    // ############################################################################################
    // Test that ``Demux:: execute`` can succeed.
    // ############################################################################################
    #[test]
    fn validate_inputs_can_succeed() {
        let tmp = TempDir::new().unwrap();
        let input_files = vec![
            fastq_file(&tmp, "read1", "ex", &["GATTACA"]),
            fastq_file(&tmp, "read2", "ex", &["TAGGATTA"]),
            fastq_file(&tmp, "index1", "ex", &[&SAMPLE1_BARCODE[0..3]]),
            fastq_file(&tmp, "index2", "ex", &[&SAMPLE1_BARCODE[3..]]),
        ];
        let sample_metadata = metadata(&tmp);

        let read_structures = vec![
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
        ];

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output_types: vec!['T'],
            output: tmp.path().to_path_buf(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 3,
        };
        demux_inputs.execute().unwrap();
    }

    #[rstest]
    #[should_panic(expected = "The same number of read structures should be given as FASTQs")]
    #[case(vec![
        ReadStructure::from_str("+T").unwrap(),
        ReadStructure::from_str("+T").unwrap(),
        ReadStructure::from_str("+B").unwrap(),
        ])]
    #[should_panic(expected = "The same number of read structures should be given as FASTQs")]
    #[case(vec![
        ReadStructure::from_str("+T").unwrap(),
        ReadStructure::from_str("+T").unwrap(),
        ReadStructure::from_str("+B").unwrap(),
        ReadStructure::from_str("+B").unwrap(),
        ReadStructure::from_str("+B").unwrap(),
    ])]
    fn test_different_number_of_read_structs_and_inputs_fails(
        #[case] read_structures: Vec<ReadStructure>,
    ) {
        let tmp = TempDir::new().unwrap();
        let input_files = vec![
            fastq_file(&tmp, "read1", "ex", &["GATTACA"]),
            fastq_file(&tmp, "read2", "ex", &["TAGGATTA"]),
            fastq_file(&tmp, "index1", "ex", &[&SAMPLE1_BARCODE[0..3]]),
            fastq_file(&tmp, "index2", "ex", &[&SAMPLE1_BARCODE[3..]]),
        ];
        let sample_metadata = metadata(&tmp);

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output_types: vec!['T'],
            output: tmp.path().to_path_buf(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 3,
        };
        demux_inputs.execute().unwrap();
    }

    #[test]
    #[should_panic(expected = "cannot be read-only")]
    fn test_read_only_output_dir_fails() {
        let tmp = TempDir::new().unwrap();
        let read_structures = vec![
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
        ];

        let input_files = vec![
            fastq_file(&tmp, "read1", "ex", &["GATTACA"]),
            fastq_file(&tmp, "read2", "ex", &["TAGGATTA"]),
            fastq_file(&tmp, "index1", "ex", &[&SAMPLE1_BARCODE[0..3]]),
            fastq_file(&tmp, "index2", "ex", &[&SAMPLE1_BARCODE[3..]]),
        ];
        let sample_metadata = metadata(&tmp);
        let mut permissions = tmp.path().metadata().unwrap().permissions();
        permissions.set_readonly(true);
        fs::set_permissions(tmp.path(), permissions.clone()).unwrap();

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output_types: vec!['T'],
            output: tmp.path().to_path_buf(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 3,
        };
        let demux_result = demux_inputs.execute();
        permissions.set_readonly(false);
        fs::set_permissions(tmp.path(), permissions).unwrap();
        demux_result.unwrap();
    }

    #[test]
    #[should_panic(expected = "doesn't exist")]
    fn test_inputs_doesnt_exist_fails() {
        let tmp = TempDir::new().unwrap();
        let read_structures = vec![
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
        ];

        let input_files = vec![
            tmp.path().join("this_file_does_not_exist.fq"),
            fastq_file(&tmp, "read2", "ex", &["TAGGATTA"]),
            fastq_file(&tmp, "index1", "ex", &[&SAMPLE1_BARCODE[0..3]]),
            fastq_file(&tmp, "index2", "ex", &[&SAMPLE1_BARCODE[3..]]),
        ];
        let sample_metadata = metadata(&tmp);

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output_types: vec!['T'],
            output: tmp.path().to_path_buf(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 2,
        };
        demux_inputs.execute().unwrap();
    }

    #[test]
    #[should_panic(expected = "Threads provided 2 was too low!")]
    fn test_too_few_threads_fails() {
        let tmp = TempDir::new().unwrap();
        let read_structures = vec![
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
        ];

        let input_files = vec![
            fastq_file(&tmp, "read1", "ex", &["GATTACA"]),
            fastq_file(&tmp, "read2", "ex", &["TAGGATTA"]),
            fastq_file(&tmp, "index1", "ex", &[&SAMPLE1_BARCODE[0..3]]),
            fastq_file(&tmp, "index2", "ex", &[&SAMPLE1_BARCODE[3..]]),
        ];
        let sample_metadata = metadata(&tmp);

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output_types: vec!['T'],
            output: tmp.path().to_path_buf(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 2,
        };
        demux_inputs.execute().unwrap();
    }

    #[test]
    fn test_demux_fragment_reads() {
        let tmp = TempDir::new().unwrap();
        let read_structures = vec![ReadStructure::from_str("17B100T").unwrap()];
        let s1_barcode = "AAAAAAAAGATTACAGA";
        let sample_metadata = metadata_file(
            &tmp,
            &[s1_barcode, "CCCCCCCCGATTACAGA", "GGGGGGGGGATTACAGA", "GGGGGGTTGATTACAGA"],
        );
        let input_files =
            vec![fastq_file(&tmp, "ex", "ex", &[&(s1_barcode.to_owned() + &"A".repeat(100))])];

        let output_dir = tmp.path().to_path_buf().join("output");

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output_types: vec!['T'],
            output: output_dir.clone(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 3,
        };
        demux_inputs.execute().unwrap();

        let expected_output = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads = read_fastq(&expected_output);

        assert_eq!(fq_reads.len(), 1);
        assert_eq!(
            fq_reads[0],
            OwnedRecord {
                head: b"ex_0".to_vec(),
                seq: "A".repeat(100).as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            }
        );
    }

    #[test]
    fn test_demux_paired_reads_with_in_line_sample_barcodes() {
        let tmp = TempDir::new().unwrap();
        let read_structures = vec![
            ReadStructure::from_str("8B100T").unwrap(),
            ReadStructure::from_str("9B100T").unwrap(),
        ];
        let s1_barcode = "AAAAAAAAGATTACAGA";
        let sample_metadata = metadata_file(
            &tmp,
            &[s1_barcode, "CCCCCCCCGATTACAGA", "GGGGGGGGGATTACAGA", "GGGGGGTTGATTACAGA"],
        );
        let input_files = vec![
            fastq_file(&tmp, "ex_R1", "ex", &[&(s1_barcode[..8].to_owned() + &"A".repeat(100))]),
            fastq_file(&tmp, "ex_R2", "ex", &[&(s1_barcode[8..].to_owned() + &"T".repeat(100))]),
        ];

        let output_dir = tmp.path().to_path_buf().join("output");

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output_types: vec!['T'],
            output: output_dir.clone(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 3,
        };
        demux_inputs.execute().unwrap();

        let expected_output1 = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads1 = read_fastq(&expected_output1);

        assert_eq!(fq_reads1.len(), 1);
        assert_eq!(
            fq_reads1[0],
            OwnedRecord {
                head: b"ex_0".to_vec(),
                seq: "A".repeat(100).as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            }
        );
        let expected_output2 = output_dir.join("Sample0000.R2.fq.gz");
        let fq_reads2 = read_fastq(&expected_output2);

        assert_eq!(fq_reads2.len(), 1);
        assert_eq!(
            fq_reads2[0],
            OwnedRecord {
                head: b"ex_0".to_vec(),
                seq: "T".repeat(100).as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            }
        );
    }

    #[test]
    fn test_demux_dual_indexed_paired_end_reads() {
        let tmp = TempDir::new().unwrap();
        let read_structures = vec![
            ReadStructure::from_str("8B").unwrap(),
            ReadStructure::from_str("100T").unwrap(),
            ReadStructure::from_str("100T").unwrap(),
            ReadStructure::from_str("9B").unwrap(),
        ];
        let s1_barcode = "AAAAAAAAGATTACAGA";
        let sample_metadata = metadata_file(
            &tmp,
            &[s1_barcode, "CCCCCCCCGATTACAGA", "GGGGGGGGGATTACAGA", "GGGGGGTTGATTACAGA"],
        );
        let input_files = vec![
            fastq_file(&tmp, "ex_I1", "ex", &[&s1_barcode[..8]]),
            fastq_file(&tmp, "ex_R1", "ex", &[&"A".repeat(100)]),
            fastq_file(&tmp, "ex_R2", "ex", &[&"T".repeat(100)]),
            fastq_file(&tmp, "ex_I2", "ex", &[&s1_barcode[8..]]),
        ];

        let output_dir = tmp.path().to_path_buf().join("output");

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output_types: vec!['T'],
            output: output_dir.clone(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 3,
        };
        demux_inputs.execute().unwrap();

        let expected_output1 = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads1 = read_fastq(&expected_output1);

        assert_eq!(fq_reads1.len(), 1);
        assert_eq!(
            fq_reads1[0],
            OwnedRecord {
                head: b"ex_0".to_vec(),
                seq: "A".repeat(100).as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            }
        );
        let expected_output2 = output_dir.join("Sample0000.R2.fq.gz");
        let fq_reads2 = read_fastq(&expected_output2);

        assert_eq!(fq_reads2.len(), 1);
        assert_eq!(
            fq_reads2[0],
            OwnedRecord {
                head: b"ex_0".to_vec(),
                seq: "T".repeat(100).as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            }
        );
    }

    // TODO - expand this test to test molecular barcode once we add that info a la fgbio version
    // of this test.
    #[test]
    fn test_demux_a_wierd_set_of_reads() {
        let tmp = TempDir::new().unwrap();
        let read_structures = vec![
            ReadStructure::from_str("4B4M8S").unwrap(),
            ReadStructure::from_str("4B100T").unwrap(),
            ReadStructure::from_str("100S3B").unwrap(),
            ReadStructure::from_str("6B1S1M1T").unwrap(),
        ];
        let sample1_barcode = "AAAAAAAAGATTACAGA";
        let sample_metadata = metadata_file(
            &tmp,
            &[sample1_barcode, "CCCCCCCCGATTACAGA", "GGGGGGGGGATTACAGA", "GGGGGGTTGATTACAGA"],
        );
        let input_files = vec![
            fastq_file(&tmp, "example_1", "ex", &["AAAACCCCGGGGTTTT"]),
            fastq_file(&tmp, "example_2", "ex", &[&"A".repeat(104)]),
            fastq_file(&tmp, "example_3", "ex", &[&("T".repeat(100) + "GAT")]),
            fastq_file(&tmp, "example_4", "ex", &["TACAGAAAT"]),
        ];

        let output_dir = tmp.path().to_path_buf().join("output");

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output_types: vec!['T'],
            output: output_dir.clone(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 3,
        };
        demux_inputs.execute().unwrap();

        let expected_output1 = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads1 = read_fastq(&expected_output1);

        assert_eq!(fq_reads1.len(), 1);
        assert_eq!(
            fq_reads1[0],
            OwnedRecord {
                head: b"ex_0".to_vec(),
                seq: "A".repeat(100).as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            }
        );
        let expected_output2 = output_dir.join("Sample0000.R2.fq.gz");
        let fq_reads2 = read_fastq(&expected_output2);

        assert_eq!(fq_reads2.len(), 1);
        assert_eq!(
            fq_reads2[0],
            OwnedRecord {
                head: b"ex_0".to_vec(),
                seq: "T".as_bytes().to_vec(),
                qual: ";".as_bytes().to_vec(),
            }
        );
    }

    #[test]
    fn test_demux_a_read_structure_with_multiple_templates_in_one_read() {
        let tmp = TempDir::new().unwrap();
        let read_structures = vec![ReadStructure::from_str("17B20T20S20T20S20T").unwrap()];
        let s1_barcode = "AAAAAAAAGATTACAGA";
        let sample_metadata = metadata_file(
            &tmp,
            &[s1_barcode, "CCCCCCCCGATTACAGA", "GGGGGGGGGATTACAGA", "GGGGGGTTGATTACAGA"],
        );
        let input_files = vec![fastq_file(
            &tmp,
            "ex",
            "ex",
            &[&(s1_barcode.to_owned()
                + &"A".repeat(20)
                + &"C".repeat(20)
                + &"T".repeat(20)
                + &"C".repeat(20)
                + &"G".repeat(20))],
        )];

        let output_dir = tmp.path().to_path_buf().join("output");

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output_types: vec!['T'],
            output: output_dir.clone(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 3,
        };
        demux_inputs.execute().unwrap();

        let expected_output1 = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads1 = read_fastq(&expected_output1);

        assert_eq!(fq_reads1.len(), 1);
        assert_eq!(
            fq_reads1[0],
            OwnedRecord {
                head: b"ex_0".to_vec(),
                seq: "A".repeat(20).as_bytes().to_vec(),
                qual: ";".repeat(20).as_bytes().to_vec(),
            }
        );
        let expected_output2 = output_dir.join("Sample0000.R2.fq.gz");
        let fq_reads2 = read_fastq(&expected_output2);

        assert_eq!(fq_reads2.len(), 1);
        assert_eq!(
            fq_reads2[0],
            OwnedRecord {
                head: b"ex_0".to_vec(),
                seq: "T".repeat(20).as_bytes().to_vec(),
                qual: ";".repeat(20).as_bytes().to_vec(),
            }
        );
        let expected_output2 = output_dir.join("Sample0000.R3.fq.gz");
        let fq_reads2 = read_fastq(&expected_output2);

        assert_eq!(fq_reads2.len(), 1);
        assert_eq!(
            fq_reads2[0],
            OwnedRecord {
                head: b"ex_0".to_vec(),
                seq: "G".repeat(20).as_bytes().to_vec(),
                qual: ";".repeat(20).as_bytes().to_vec(),
            }
        );
    }

    #[test]
    #[should_panic(
        expected = "No output types requested, must request at least one output segment type."
    )]
    fn test_fails_if_zero_read_structures_have_template_bases() {
        let tmp = TempDir::new().unwrap();
        let input_files = vec![
            fastq_file(&tmp, "read1", "ex", &["GATTACA"]),
            fastq_file(&tmp, "read2", "ex", &["TAGGATTA"]),
            fastq_file(&tmp, "index1", "ex", &[&SAMPLE1_BARCODE[0..3]]),
            fastq_file(&tmp, "index2", "ex", &[&SAMPLE1_BARCODE[3..]]),
        ];
        let sample_metadata = metadata(&tmp);

        let read_structures = vec![
            ReadStructure::from_str("+M").unwrap(),
            ReadStructure::from_str("+M").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
        ];

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output_types: vec![],
            output: tmp.path().to_path_buf(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 3,
        };
        demux_inputs.execute().unwrap();
    }

    #[test]
    #[should_panic(expected = "The same number of read structures should be given as FASTQs")]
    fn test_fails_if_not_enough_fastq_records_are_passed() {
        let tmp = TempDir::new().unwrap();
        let input_files = vec![
            fastq_file(&tmp, "read1", "ex", &["GATTACA"]),
            fastq_file(&tmp, "index1", "ex", &[&SAMPLE1_BARCODE[0..3]]),
            fastq_file(&tmp, "index2", "ex", &[&SAMPLE1_BARCODE[3..]]),
        ];
        let sample_metadata = metadata(&tmp);

        let read_structures = vec![
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
            ReadStructure::from_str("+T").unwrap(),
        ];

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output_types: vec!['T'],
            output: tmp.path().to_path_buf(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 3,
        };
        demux_inputs.execute().unwrap();
    }

    #[test]
    #[should_panic(expected = "The same number of read structures should be given as FASTQs")]
    fn test_fails_if_too_many_fastq_records_are_passed() {
        let tmp = TempDir::new().unwrap();
        let input_files = vec![
            fastq_file(&tmp, "read1", "ex", &["GATTACA"]),
            fastq_file(&tmp, "index1", "ex", &[&SAMPLE1_BARCODE[0..3]]),
            fastq_file(&tmp, "index2", "ex", &[&SAMPLE1_BARCODE[3..]]),
            fastq_file(&tmp, "read2", "ex", &["TAGGATTA"]),
        ];
        let sample_metadata = metadata(&tmp);

        let read_structures = vec![
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
        ];

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output_types: vec!['T'],
            output: tmp.path().to_path_buf(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 3,
        };
        demux_inputs.execute().unwrap();
    }
}
