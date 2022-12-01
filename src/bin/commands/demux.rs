use crate::commands::command::Command;
use anyhow::{anyhow, Result};
use clap::Parser;
use fgoxide::io::Io;
use fqtk_lib::samples::SampleGroup;
use log::info;
use read_structure::ReadStructure;
use std::io::Read;
use std::{fs, io::BufReader, path::PathBuf};

/// Performs sample demultiplexing on FASTQs.
///
/// The sample barcode for each sample in the metadata TSV will be compared against the sample
/// barcode bases extracted from the FASTQs, to assign each read to a sample.  Reads that do not
/// match any sample within the given error tolerance will be placed in the ``unmatched_prefix``
/// file.
///
/// FASTQs and associated read structures for each sub-read should be given:
///
/// - a single fragment read should have one FASTQ and one read structure
/// - paired end reads should have two FASTQs and two read structures
/// - a dual-index sample with paired end reads should have four FASTQs and four read structures
///   given: two for the two index reads, and two for the template reads.
///
/// If multiple FASTQs are present for each sub-read, then the FASTQs for each sub-read should be
/// concatenated together prior to running this tool
/// (ex. `cat s_R1_L001.fq.gz s_R1_L002.fq.gz > s_R1.fq.gz`).
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
/// Metadata about the samples should be given as a headered metadata CSV file with two columns
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
    #[clap(long, short = 's', default_value = "1")]
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
    /// Checks that inputs to demux are valid and returns open file handles for the inputs.
    /// Checks:
    ///     - That the number of input files and number of read structs provided are the same
    ///     - That the output directory is not read-only
    ///     - That the input files exist
    ///     - That the input files have read permissions.
    fn validate_inputs(&self) -> Result<Vec<BufReader<Box<dyn Read>>>> {
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
            Ok(fq_readers_result?)
        } else {
            let mut details = "Inputs failed validation!\n".to_owned();
            for error_reason in constraint_errors {
                details.push_str(&format!("    - {}\n", error_reason));
            }
            Err(anyhow!("The following errors with the input(s) were detected:\n{}", details))
        }
    }
}

impl Command for Demux {
    // (For now all this does is check the inputs to the demux command, but will eventually execute
    // TODO - remove this disclaimer once working core is implemented.)
    /// Executes the demux command
    fn execute(&self) -> Result<()> {
        self.validate_inputs()?;
        let _sample_group = SampleGroup::from_file(&self.sample_metadata)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fqtk_lib::samples::Sample;
    use rstest::rstest;
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
    fn fastq_file(tmpdir: &TempDir, prefix: &str, records_bases: &[&str]) -> PathBuf {
        let io = Io::default();

        let path = tmpdir.path().join(format!("{prefix}.fastq"));
        let fastq_lines = fq_lines_from_bases(prefix, records_bases);
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

    // ############################################################################################
    // Test that ``Demux:: execute`` can succeed.
    // ############################################################################################
    #[test]
    fn validate_inputs_can_succeed() {
        let tmpdir = TempDir::new().unwrap();
        let input_files = vec![
            fastq_file(&tmpdir, "read1", &["GATTACA"]),
            fastq_file(&tmpdir, "read2", &["TAGGATTA"]),
            fastq_file(&tmpdir, "index1", &[&SAMPLE1_BARCODE[0..3]]),
            fastq_file(&tmpdir, "index2", &[&SAMPLE1_BARCODE[3..]]),
        ];
        let sample_metadata = metadata(&tmpdir);

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
            output: tmpdir.path().to_path_buf(),
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
        let tmpdir = TempDir::new().unwrap();
        let input_files = vec![
            fastq_file(&tmpdir, "read1", &["GATTACA"]),
            fastq_file(&tmpdir, "read2", &["TAGGATTA"]),
            fastq_file(&tmpdir, "index1", &[&SAMPLE1_BARCODE[0..3]]),
            fastq_file(&tmpdir, "index2", &[&SAMPLE1_BARCODE[3..]]),
        ];
        let sample_metadata = metadata(&tmpdir);

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output: tmpdir.path().to_path_buf(),
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
        let tmpdir = TempDir::new().unwrap();
        let read_structures = vec![
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
        ];

        let input_files = vec![
            fastq_file(&tmpdir, "read1", &["GATTACA"]),
            fastq_file(&tmpdir, "read2", &["TAGGATTA"]),
            fastq_file(&tmpdir, "index1", &[&SAMPLE1_BARCODE[0..3]]),
            fastq_file(&tmpdir, "index2", &[&SAMPLE1_BARCODE[3..]]),
        ];
        let sample_metadata = metadata(&tmpdir);
        let mut permissions = tmpdir.path().metadata().unwrap().permissions();
        permissions.set_readonly(true);
        fs::set_permissions(tmpdir.path(), permissions.clone()).unwrap();

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output: tmpdir.path().to_path_buf(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 3,
        };
        let demux_result = demux_inputs.execute();
        permissions.set_readonly(false);
        fs::set_permissions(tmpdir.path(), permissions).unwrap();
        demux_result.unwrap();
    }

    #[test]
    #[should_panic(expected = "doesn't exist")]
    fn test_inputs_doesnt_exist_fails() {
        let tmpdir = TempDir::new().unwrap();
        let read_structures = vec![
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
        ];

        let input_files = vec![
            tmpdir.path().join("this_file_does_not_exist.fq"),
            fastq_file(&tmpdir, "read2", &["TAGGATTA"]),
            fastq_file(&tmpdir, "index1", &[&SAMPLE1_BARCODE[0..3]]),
            fastq_file(&tmpdir, "index2", &[&SAMPLE1_BARCODE[3..]]),
        ];
        let sample_metadata = metadata(&tmpdir);

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output: tmpdir.path().to_path_buf(),
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
        let tmpdir = TempDir::new().unwrap();
        let read_structures = vec![
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
        ];

        let input_files = vec![
            fastq_file(&tmpdir, "read1", &["GATTACA"]),
            fastq_file(&tmpdir, "read2", &["TAGGATTA"]),
            fastq_file(&tmpdir, "index1", &[&SAMPLE1_BARCODE[0..3]]),
            fastq_file(&tmpdir, "index2", &[&SAMPLE1_BARCODE[3..]]),
        ];
        let sample_metadata = metadata(&tmpdir);

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output: tmpdir.path().to_path_buf(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 2,
        };
        demux_inputs.execute().unwrap();
    }
}
