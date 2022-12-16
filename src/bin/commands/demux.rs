#![warn(missing_docs)]
#![warn(clippy::missing_docs_in_private_items)]

use crate::commands::command::Command;
use anyhow::{anyhow, Result};
use clap::Parser;
use fgoxide::io::{DelimFile, Io};
use fqtk_lib::barcode_matching::BarcodeMatcher;
use fqtk_lib::samples::Sample;
use fqtk_lib::samples::SampleGroup;
use gzp::BUFSIZE;
use log::info;
use pooled_writer::{bgzf::BgzfCompressor, Pool, PoolBuilder, PooledWriter};
use proglog::ProgLogBuilder;
use read_structure::ReadSegment;
use read_structure::ReadStructure;
use read_structure::ReadStructureError;
use read_structure::SegmentType;
use seq_io::fastq::Reader as FastqReader;
use seq_io::fastq::Record;
use serde::Serialize;
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufWriter, Read, Write};
use std::str;
use std::{
    fs,
    path::{Path, PathBuf},
};

/// Number of fields in the FASTQ header before the UMI / molecular barcode field
const NUMBER_BEFORE_MOL_BARCODE: usize = 7;
/// Number of colons in the FASTQ header before the colon immediately preceeding the UMI / molecular barcode field
const COLONS_BEFORE_MOL_BARCODE: usize = 6;

/// type alias for the purpose of getting clippy to not complain about complex types further down
/// in the code.
type VecOfReaders = Vec<Box<dyn Read>>;

#[allow(clippy::too_many_arguments)]
/// Writes components of a FASTQ record to a provided writer
fn write_fastq_record_to_writer<W: Write>(
    writer: &mut W,
    before_mol_barcode_header_chunk: &[u8],
    mol_barcode_header_chunk: &[u8],
    record_number_header_chunk: &[u8],
    after_record_number_before_sample_barcode_header_chunk: &[u8],
    sample_barcode_header_chunk: &[u8],
    sequence: &[u8],
    quals: &[u8],
) -> Result<()> {
    writer.write_all(b"@")?;
    writer.write_all(before_mol_barcode_header_chunk)?;
    writer.write_all(mol_barcode_header_chunk)?;
    writer.write_all(record_number_header_chunk)?;
    writer.write_all(after_record_number_before_sample_barcode_header_chunk)?;
    writer.write_all(sample_barcode_header_chunk)?;
    writer.write_all(b"\n")?;
    writer.write_all(sequence)?;
    writer.write_all(b"\n+\n")?;
    writer.write_all(quals)?;
    Ok(())
}

#[derive(Debug, Clone, Serialize)]
struct DemuxMetrics {
    sample_name: String,
    sample_barcode: String,
    number_reads: usize,
}

impl DemuxMetrics {
    pub fn get_metrics(
        number_reads: &[usize],
        samples: &[Sample],
        unmatched_name: String,
    ) -> Vec<Self> {
        let mut result = Vec::with_capacity(number_reads.len());
        assert_eq!(
            number_reads.len(),
            samples.len() + 1,
            "Number of samples with metrics being tracked differed from expected!"
        );
        for (sample, &number_sample_reads) in samples.iter().zip(number_reads.iter()) {
            result.push(Self {
                sample_name: sample.name.clone(),
                sample_barcode: sample.barcode.clone(),
                number_reads: number_sample_reads,
            });
        }
        result.push(Self {
            sample_name: unmatched_name,
            sample_barcode: "N/A".to_string(),
            number_reads: number_reads[number_reads.len() - 1],
        });
        result
    }
}

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
    /// Portion of the FASTQ record header before the molecular barcode.
    before_mol_barcode_header_chunk: Vec<u8>,
    /// Portion of the FASTQ record header that occurs between the record number and the
    /// sample barcode
    after_record_number_before_sample_barcode_header_chunk: Vec<u8>,
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

impl ReadSet {
    /// The full molecular barcode for this record (concatenated bases from individual molecular
    /// barcode segments, separated by '+' char)
    fn molecular_barcode_header_component(&self) -> Vec<u8> {
        if let Some(mol_barcodes) = &self.molecular_barcode_segments {
            if mol_barcodes.is_empty() {
                Vec::new()
            } else {
                let segments = mol_barcodes
                    .iter()
                    .map(|m| {
                        str::from_utf8(&m.seq).expect("unexpected error parsing molecular barcodes")
                    })
                    .collect::<Vec<_>>();
                let num_segments = segments.len();
                // If it were just the barcode it would be num_segments - 1 but we're accounting
                // for the : that needs to be appended if it's there.
                let size = segments.iter().map(|&s| s.len()).sum::<usize>() + num_segments;
                let mut result = Vec::with_capacity(size);
                result.push(b':');
                for (i, &segment) in segments.iter().enumerate() {
                    result.extend_from_slice(segment.as_bytes());
                    if i < num_segments - 1 {
                        result.push(b'+');
                    }
                }
                result
            }
        } else {
            Vec::new()
        }
    }

    /// Adds the number of ':' characters necessary to reach the desired level of padding for this
    /// header.
    /// If not padding for this header: allocates a new vector for the beginning portion of the
    /// existing FASTQ header and clones it into it
    fn create_padded_header_contents(
        existing_header: &[u8],
        before_space_ending_offset: usize,
        num_colons_before_space: usize,
        pad_for_mol_barcode: bool,
    ) -> Vec<u8> {
        // Determine the approx number of characters we'll need in the new header and use that to
        // instantiate a new vec with known capacity.
        let before_colon_pad =
            if !pad_for_mol_barcode || num_colons_before_space >= COLONS_BEFORE_MOL_BARCODE {
                0
            } else {
                COLONS_BEFORE_MOL_BARCODE - num_colons_before_space
            };

        let mut before_mol_barcode_chunk =
            Vec::with_capacity(before_space_ending_offset + before_colon_pad);
        before_mol_barcode_chunk.extend_from_slice(&existing_header[..before_space_ending_offset]);
        // Pad with colons to make it semi-valid Illumina format (if adding molecular barcode and
        // the colons didn't already exist)
        if before_colon_pad != 0 {
            before_mol_barcode_chunk.extend(vec![b':'; before_colon_pad]);
        };
        before_mol_barcode_chunk
    }

    /// Parses the existing FASTQ header and determines the relevant boundary points for
    /// deconstructing it and reconstructing it with new barcodes.
    /// If adding molecular barcode to header, will trim any existing readname starting at the
    /// point of the 7th ':' character before a space in the header.
    fn determine_boundaries_in_existing_header(
        existing_header: &[u8],
        add_molecular_barcode_to_header: bool,
    ) -> (usize, usize, Option<(usize, usize)>) {
        let mut space_offset = None;
        let mut contents_after_space = false;

        let mut num_colons_before_space = 0usize;
        let mut before_space_ending_offset = 0usize;

        let mut num_colons_after_space = 0usize;
        let mut after_space_starting_offset = None;
        let mut after_space_ending_offset = None;
        for (i, &c) in existing_header.iter().enumerate() {
            if c == b' ' {
                space_offset = Some(i);
            } else {
                if space_offset.is_some() {
                    contents_after_space = true;
                }
                if c == b':' {
                    if space_offset.is_some() {
                        num_colons_after_space += 1;
                        if num_colons_after_space == 1 {
                            after_space_starting_offset = Some(i);
                        } else if num_colons_after_space == 3 {
                            after_space_ending_offset = Some(i + 1); // slices are half open
                        }
                    } else {
                        num_colons_before_space += 1;

                        if num_colons_before_space == NUMBER_BEFORE_MOL_BARCODE {
                            // Don't include the 7th : as the mol barcode will include it if it has
                            // something to add.
                            before_space_ending_offset = i;
                        }
                    }
                } else if num_colons_before_space <= COLONS_BEFORE_MOL_BARCODE
                    && space_offset.is_none()
                {
                    before_space_ending_offset = i + 1;
                }
            }
        }
        if num_colons_before_space < NUMBER_BEFORE_MOL_BARCODE || !add_molecular_barcode_to_header {
            before_space_ending_offset = if let Some(space_offset_inner) = space_offset {
                space_offset_inner
            } else {
                existing_header.len()
            };
        }

        if contents_after_space {
            // I can see an argument for not doing this, but for now I'm enforcing it
            assert!(num_colons_after_space != 0, "Invalid FASTQ contents after space, must be either no contents after space or illumina format (https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm).");
        }

        assert!(num_colons_after_space == 0 || num_colons_after_space == 3, "Invalid FASTQ contents after space, must be either no contents after space or illumina format (https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm).");
        (
            before_space_ending_offset,
            num_colons_before_space,
            after_space_starting_offset.zip(after_space_ending_offset),
        )
    }

    /// Modifies the read header field (``name``) and ``read_number_position`` for outputting
    /// read headers:
    ///     1. removes the read number from the header
    ///     2. pads the read name with an appropriate number of ':' characters to make the
    ///        resulting FASTQ headers semi-valid FASTQ in the event they are not valid already.
    ///     3.
    fn process_existing_read_header(
        existing_header: &[u8],
        pad_for_mol_barcode: bool,
        add_molecular_barcode_to_header: bool,
        add_sample_barcode_to_header: bool,
    ) -> (Vec<u8>, Vec<u8>) {
        let (before_space_ending_offset, num_colons_before_space, after_space_offsets) =
            Self::determine_boundaries_in_existing_header(
                existing_header,
                add_molecular_barcode_to_header,
            );

        let before_mol_barcode_chunk = Self::create_padded_header_contents(
            existing_header,
            before_space_ending_offset,
            num_colons_before_space,
            pad_for_mol_barcode && add_molecular_barcode_to_header,
        );

        let after_record_chunk =
            if let Some((after_space_start, after_space_end)) = after_space_offsets {
                if add_sample_barcode_to_header {
                    existing_header[after_space_start..after_space_end].to_vec()
                } else {
                    existing_header[after_space_start..].to_vec()
                }
            } else if add_sample_barcode_to_header {
                b":N:0:".to_vec()
            } else {
                Vec::new()
            };
        (before_mol_barcode_chunk, after_record_chunk)
    }
}

/// A struct for iterating over the records in multiple FASTQ files simultaneously, destructuring
/// them according to the provided read structures, yielding ``ReadSet`` objects on each iteration.
#[allow(clippy::struct_excessive_bools)]
struct ReadSetIterator {
    /// Read structure objects describing the structure of the reads to be demultiplexed, one per
    /// input file.
    read_structures: Vec<ReadStructure>,
    /// expected length of the sample barcode given the read structures
    sample_barcode_length: usize,
    /// number of template segments per set of FASTQ records
    num_templates: usize,
    /// number of molecular barcode segments per set of FASTQ records
    num_molecular_barcodes: usize,
    /// number of sample barcode segments per set of FASTQ records
    num_sample_barcodes: usize,
    /// Iterators over the files containing FASTQ records, one per input file.
    sources: Vec<FastqReader<Box<dyn Read>>>,
    /// If true this iterator will populate the template segments field of the ReadSet object it
    /// returns
    generate_template_segments: bool,
    /// If true this iterator will populate the sample barcode segments field of the ReadSet
    /// objects it returns
    generate_sample_barcode_segments: bool,
    /// If true this iterator will populate the molecular barcode segments field of the ReadSet
    /// objects it returns
    generate_molecular_barcode_segments: bool,
    /// If true this iterator will add molecular barcodes to the header of the resulting reads.
    add_molecular_barcode_to_header: bool,
    /// If true this iterator will add molecular barcodes to the header of the resulting reads.
    add_sample_barcode_to_header: bool,
}

impl Iterator for ReadSetIterator {
    type Item = ReadSet;

    fn next(&mut self) -> Option<Self::Item> {
        let source_len = self.sources.len();
        let mut next_fq_recs = Vec::with_capacity(source_len);
        for source in &mut self.sources {
            if let Some(rec) = source.next() {
                next_fq_recs.push(rec.expect("Unexpected error parsing FASTQs"));
            }
        }

        if next_fq_recs.is_empty() {
            return None;
        }
        assert_eq!(
            next_fq_recs.len(),
            source_len,
            "FASTQ sources out of sync at records: {:?}",
            next_fq_recs
        );

        let mut sample_barcode = Vec::with_capacity(self.sample_barcode_length);
        let mut template_segments = if self.generate_template_segments {
            Some(Vec::with_capacity(self.num_templates))
        } else {
            None
        };
        let mut sample_barcode_segments =
            if self.generate_sample_barcode_segments || self.add_sample_barcode_to_header {
                Some(Vec::with_capacity(self.num_sample_barcodes))
            } else {
                None
            };
        let mut molecular_barcode_segments =
            if self.generate_molecular_barcode_segments || self.add_molecular_barcode_to_header {
                Some(Vec::with_capacity(self.num_molecular_barcodes))
            } else {
                None
            };

        for (read_structure, fastq_record) in self.read_structures.iter().zip(next_fq_recs.iter()) {
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

        let (before_chunk, after_chunk) = ReadSet::process_existing_read_header(
            next_fq_recs[0].head(),
            self.num_molecular_barcodes != 0
                && molecular_barcode_segments.is_some()
                && self.add_molecular_barcode_to_header,
            self.add_molecular_barcode_to_header,
            self.add_sample_barcode_to_header,
        );

        Some(ReadSet {
            before_mol_barcode_header_chunk: before_chunk,
            after_record_number_before_sample_barcode_header_chunk: after_chunk,
            sample_barcode_sequence: sample_barcode,
            template_segments,
            sample_barcode_segments,
            molecular_barcode_segments,
        })
    }
}

impl ReadSetIterator {
    /// Instantiates a new iterator over the read sets for a set of FASTQs with defined read
    /// structures
    pub fn new(
        read_structures: Vec<ReadStructure>,
        fq_sources: Vec<FastqReader<Box<dyn Read>>>,
        output_segment_types: &HashSet<SegmentType>,
        add_molecular_barcode_to_header: bool,
        add_sample_barcode_to_header: bool,
    ) -> Self {
        let sample_barcode_length: usize = read_structures
            .iter()
            .map(|rs| {
                rs.segments_by_type(SegmentType::SampleBarcode)
                    .filter_map(ReadSegment::length)
                    .sum::<usize>()
            })
            .sum();
        let num_templates: usize = read_structures
            .iter()
            .map(|rs| rs.segments_by_type(SegmentType::Template).count())
            .sum();
        let num_sample_barcodes: usize = read_structures
            .iter()
            .map(|rs| rs.segments_by_type(SegmentType::SampleBarcode).count())
            .sum();
        let num_molecular_barcodes: usize = read_structures
            .iter()
            .map(|rs| rs.segments_by_type(SegmentType::MolecularBarcode).count())
            .sum();

        Self {
            read_structures,
            sources: fq_sources,
            sample_barcode_length,
            num_templates,
            num_sample_barcodes,
            num_molecular_barcodes,
            generate_template_segments: output_segment_types.contains(&SegmentType::Template),
            generate_sample_barcode_segments: output_segment_types
                .contains(&SegmentType::SampleBarcode),
            generate_molecular_barcode_segments: output_segment_types
                .contains(&SegmentType::MolecularBarcode),
            add_molecular_barcode_to_header,
            add_sample_barcode_to_header,
        }
    }
}

/// Stores the writers for a single sample in demultiplexing. Fields can be None if that type of
/// ``ReadSegment`` is not being written.
struct SampleWriters<W: Write> {
    /// Name of the sample this set of writers is for
    name: String,
    /// Vec of the writers for template read segments.
    template_writers: Option<Vec<W>>,
    /// Vec of the writers for the sample barcode read segments.
    sample_barcode_writers: Option<Vec<W>>,
    /// Vec of the writers for the molecular barcode read segments.
    molecular_barcode_writers: Option<Vec<W>>,
}

impl<W: Write> SampleWriters<W> {
    /// Destroys this struct and decomposes it into its component types. Used when swapping
    /// writers for pooled writers.
    #[allow(clippy::type_complexity)]
    fn into_parts(self) -> (String, Option<Vec<W>>, Option<Vec<W>>, Option<Vec<W>>) {
        (
            self.name,
            self.template_writers,
            self.sample_barcode_writers,
            self.molecular_barcode_writers,
        )
    }

    /// Writes a set of reads (defined as a ``ReadSet``) to the appropriate writers on this
    /// ``Self`` struct.
    /// Reads in the read set should be 1:1 with writers in the writer set however this is not
    /// checked at runtime as doing so substantially slows demulitplexing.
    fn write(
        &mut self,
        read_set: &ReadSet,
        write_mol_barcode: bool,
        write_sample_barcode: bool,
    ) -> Result<()> {
        let empty_vec = Vec::new();

        let mol_barcode_chunk = if write_mol_barcode {
            read_set.molecular_barcode_header_component()
        } else {
            Vec::new()
        };

        let sample_barcode_chunk =
            if write_sample_barcode { &read_set.sample_barcode_sequence } else { &empty_vec };

        let write_read_no = read_set.after_record_number_before_sample_barcode_header_chunk.len()
            + sample_barcode_chunk.len()
            > 0;

        for (writers_opt, segments_opt) in [
            (&mut self.template_writers, &read_set.template_segments),
            (&mut self.sample_barcode_writers, &read_set.sample_barcode_segments),
            (&mut self.molecular_barcode_writers, &read_set.molecular_barcode_segments),
        ] {
            if let (Some(writers), Some(segments)) = (writers_opt, segments_opt) {
                for (i, (writer, segment)) in writers.iter_mut().zip(segments.iter()).enumerate() {
                    write_fastq_record_to_writer(
                        writer,
                        &read_set.before_mol_barcode_header_chunk,
                        &mol_barcode_chunk,
                        (if write_read_no {
                            // encode the space here because otherwise we'll include it before the
                            // newline even in the case where no after-space segment exists
                            format!(" {}", i + 1)
                        } else {
                            String::new()
                        })
                        .as_bytes(),
                        &read_set.after_record_number_before_sample_barcode_header_chunk,
                        sample_barcode_chunk,
                        &segment.seq,
                        &segment.quals,
                    )?;
                }
            }
        }
        Ok(())
    }
}

impl SampleWriters<PooledWriter> {
    /// Attempts to gracefully shutdown the writers in this struct, consuming the struct in the
    /// process
    /// # Errors
    ///     - Will error if closing of the ``PooledWriter``s fails for any reason
    fn close(self) -> Result<()> {
        for writers in
            [self.template_writers, self.sample_barcode_writers, self.molecular_barcode_writers]
                .into_iter()
                .flatten()
        {
            writers.into_iter().try_for_each(PooledWriter::close)?;
        }

        Ok(())
    }
}

/// Constructs a pooled writer from individual sample writers and a number of threads, and converts
/// the writers to pooled writer objects.
/// Returns the pooled writers organized into ``SampleWriters`` structs and the pool
/// struct.
/// # Errors:
///     - Should never error but will if there is an internal logic error that results in
///         zero template read writers being produced during conversion.
/// # Panics
///     - Should never panic but will if there is an internal logic issue that results in a None
///         type being unwrapped when convering writers.
fn build_pool(
    sample_writers: Vec<SampleWriters<BufWriter<File>>>,
    compression_level: usize,
    threads: usize,
) -> Result<(Vec<SampleWriters<PooledWriter>>, Pool)> {
    let mut new_sample_writers = Vec::with_capacity(sample_writers.len());
    let mut pool_builder = PoolBuilder::<_, BgzfCompressor>::new(threads * 2, threads)
        .compression_level(u8::try_from(compression_level)?)?;
    for sample in sample_writers {
        let (name, template_writers, barcode_writers, mol_writers) = sample.into_parts();
        let mut new_template_writers = None;
        let mut new_sample_barcode_writers = None;
        let mut new_molecular_barcode_writers = None;

        for (optional_ws, target) in vec![
            (template_writers, &mut new_template_writers),
            (barcode_writers, &mut new_sample_barcode_writers),
            (mol_writers, &mut new_molecular_barcode_writers),
        ] {
            if let Some(ws) = optional_ws {
                let mut new_writers = Vec::with_capacity(ws.len());
                for writer in ws {
                    new_writers.push(pool_builder.exchange(writer));
                }
                _ = target.insert(new_writers);
            }
        }
        new_sample_writers.push(SampleWriters {
            name,
            template_writers: new_template_writers,
            sample_barcode_writers: new_sample_barcode_writers,
            molecular_barcode_writers: new_molecular_barcode_writers,
        });
    }
    let pool = pool_builder.build()?;
    Ok((new_sample_writers, pool))
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

    /// The level of compression to use to compress outputs.
    #[clap(long, short = 'c', default_value = "5")]
    compression_level: usize,

    /// The filename and extension to use for the generated metrics file (will be placed in
    /// ``output`` directory)
    #[clap(long, default_value = "fqtk-demux-metrics.tsv")]
    metrics_filename: String,

    /// If true, will not output metrics file
    #[clap(long)]
    no_metrics: bool,

    /// If true, caching of barcodes will be disabled
    #[clap(long, hide = true)]
    no_cache: bool,

    /// If true, will not attempt to add molecular barcodes to the FASTQ headers
    #[clap(long)]
    no_molecular_barcode_in_header: bool,

    /// If true, will not attempt to add sample barcodes to the FASTQ headers
    #[clap(long)]
    no_sample_barcode_in_header: bool,
}

impl Demux {
    /// Creates one writer per read segment in the read structures on this object, restricted by
    /// requested output type.
    /// # Errors:
    ///     - Will error if opening the output fails for any reason.
    ///     - Will error if no template segments were found in the read structures on this object.
    fn create_sample_writers(
        read_structures: &[ReadStructure],
        prefix: &str,
        output_types: &HashSet<SegmentType>,
        output_dir: &Path,
    ) -> Result<SampleWriters<BufWriter<File>>> {
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
                output_type_writers.push(BufWriter::new(File::create(
                    output_dir.join(format!("{}.{}{}.fq.gz", prefix, file_type_code, idx)),
                )?));
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

        Ok(SampleWriters {
            name: prefix.to_string(),
            template_writers,
            sample_barcode_writers,
            molecular_barcode_writers,
        })
    }

    /// Creates one writer per sample per read segment in the provided read structures for
    /// requested output type.
    /// # Errors:
    ///     - Will error if opening the output fails for any reason.
    fn create_writers(
        read_structures: &[ReadStructure],
        sample_group: &SampleGroup,
        output_types: &HashSet<SegmentType>,
        unmatched_prefix: &str,
        output_dir: &Path,
    ) -> Result<Vec<SampleWriters<BufWriter<File>>>> {
        let mut samples_writers = Vec::with_capacity(sample_group.samples.len());
        for sample in &sample_group.samples {
            samples_writers.push(Self::create_sample_writers(
                read_structures,
                &sample.name,
                output_types,
                output_dir,
            )?);
        }
        // Add the unmatched 'sample'
        samples_writers.push(Self::create_sample_writers(
            read_structures,
            unmatched_prefix,
            output_types,
            output_dir,
        )?);
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
            .map(|p| {
                let inner: Box<dyn Read> = Box::new(io.new_reader(p)?);
                Ok(inner)
            })
            .collect::<Result<VecOfReaders, fgoxide::FgError>>();
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
            .map(|fq| FastqReader::with_capacity(fq, BUFSIZE))
            .collect::<Vec<_>>();

        let (mut sample_writers, mut pool) = build_pool(
            Self::create_writers(
                &self.read_structures,
                &sample_group,
                &output_segment_types,
                &self.unmatched_prefix,
                &self.output,
            )?,
            self.compression_level,
            self.threads - 1,
        )?;
        let unmatched_index = sample_writers.len() - 1;
        info!("Created sample (and unassigned) writers");

        let mut barcode_matcher = BarcodeMatcher::new(
            &sample_group.samples.iter().map(|s| s.barcode.as_str()).collect::<Vec<_>>(),
            u8::try_from(self.max_mismatches)?,
            u8::try_from(self.min_mismatch_delta)?,
            !self.no_cache,
        );

        let fq_iterator = ReadSetIterator::new(
            self.read_structures.clone(),
            fq_sources,
            &output_segment_types,
            !self.no_molecular_barcode_in_header,
            !self.no_sample_barcode_in_header,
        );

        let logger = ProgLogBuilder::new()
            .name("fqtk")
            .noun("records")
            .verb("demultiplexed")
            .unit(1_000_000)
            .level(log::Level::Info)
            .build();
        let mut number_reads: Vec<usize> = vec![0; sample_writers.len()];
        for read_set in fq_iterator {
            if let Some(barcode_match) = barcode_matcher.assign(&read_set.sample_barcode_sequence) {
                sample_writers[barcode_match.best_match].write(
                    &read_set,
                    !self.no_molecular_barcode_in_header,
                    !self.no_sample_barcode_in_header,
                )?;
                number_reads[barcode_match.best_match] += 1;
            } else {
                sample_writers[unmatched_index].write(
                    &read_set,
                    !self.no_molecular_barcode_in_header,
                    !self.no_sample_barcode_in_header,
                )?;
                number_reads[unmatched_index] += 1;
            }
            logger.record();
        }
        info!("Finished reading input FASTQs.");
        sample_writers.into_iter().map(SampleWriters::close).collect::<Result<Vec<_>>>()?;
        pool.stop_pool()?;
        info!("Output FASTQ writing complete.");
        if !self.no_metrics {
            let delim_file = DelimFile::default();
            let metrics_output = self.output.join(&self.metrics_filename);
            delim_file.write_tsv(
                &metrics_output,
                DemuxMetrics::get_metrics(
                    &number_reads,
                    &sample_group.samples,
                    self.unmatched_prefix.clone(),
                ),
            )?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bytecount;
    use fqtk_lib::samples::Sample;
    use rstest::rstest;
    use seq_io::fastq::OwnedRecord;
    use std::str;
    use std::str::FromStr;
    use tempfile::TempDir;

    const SAMPLE1_BARCODE: &str = "GATTGGG";

    fn print_bytes(bytes: &[u8]) {
        println!("{}", str::from_utf8(bytes).unwrap());
    }

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
            compression_level: 5,
            no_metrics: true,
            metrics_filename: "fqtk-demux-metrics.tsv".to_string(),
            no_cache: false,
            no_molecular_barcode_in_header: false,
            no_sample_barcode_in_header: false,
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
            compression_level: 5,
            no_metrics: true,
            metrics_filename: "fqtk-demux-metrics.tsv".to_string(),
            no_cache: false,
            no_molecular_barcode_in_header: false,
            no_sample_barcode_in_header: false,
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
            compression_level: 5,
            no_metrics: true,
            metrics_filename: "fqtk-demux-metrics.tsv".to_string(),
            no_cache: false,
            no_molecular_barcode_in_header: false,
            no_sample_barcode_in_header: false,
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
            compression_level: 5,
            no_metrics: true,
            metrics_filename: "fqtk-demux-metrics.tsv".to_string(),
            no_cache: false,
            no_molecular_barcode_in_header: false,
            no_sample_barcode_in_header: false,
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
            compression_level: 5,
            no_metrics: true,
            metrics_filename: "fqtk-demux-metrics.tsv".to_string(),
            no_cache: false,
            no_molecular_barcode_in_header: false,
            no_sample_barcode_in_header: false,
        };
        demux_inputs.execute().unwrap();
    }

    #[rstest]
    #[case(true)]
    #[case(false)]
    fn test_demux_fragment_reads(#[case] no_cache: bool) {
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
            compression_level: 5,
            no_metrics: true,
            metrics_filename: "fqtk-demux-metrics.tsv".to_string(),
            no_cache,
            no_molecular_barcode_in_header: false,
            no_sample_barcode_in_header: false,
        };
        demux_inputs.execute().unwrap();

        let expected_output = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads = read_fastq(&expected_output);

        assert_eq!(fq_reads.len(), 1);
        assert_eq!(
            fq_reads[0],
            OwnedRecord {
                head: b"ex_0 1:N:0:AAAAAAAAGATTACAGA".to_vec(),
                seq: "A".repeat(100).as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            }
        );
    }

    #[rstest]
    #[case(true)]
    #[case(false)]
    fn test_demux_paired_reads_with_in_line_sample_barcodes(#[case] no_cache: bool) {
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
            compression_level: 5,
            no_metrics: true,
            metrics_filename: "fqtk-demux-metrics.tsv".to_string(),
            no_cache,
            no_molecular_barcode_in_header: false,
            no_sample_barcode_in_header: false,
        };
        demux_inputs.execute().unwrap();

        let expected_output1 = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads1 = read_fastq(&expected_output1);

        assert_eq!(fq_reads1.len(), 1);
        assert_eq!(
            fq_reads1[0],
            OwnedRecord {
                head: b"ex_0 1:N:0:AAAAAAAAGATTACAGA".to_vec(),
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
                head: b"ex_0 2:N:0:AAAAAAAAGATTACAGA".to_vec(),
                seq: "T".repeat(100).as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            }
        );
    }

    #[rstest]
    #[case(true)]
    #[case(false)]
    fn test_demux_dual_indexed_paired_end_reads(#[case] no_cache: bool) {
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
            compression_level: 5,
            no_metrics: true,
            metrics_filename: "fqtk-demux-metrics.tsv".to_string(),
            no_cache,
            no_molecular_barcode_in_header: false,
            no_sample_barcode_in_header: false,
        };
        demux_inputs.execute().unwrap();

        let expected_output1 = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads1 = read_fastq(&expected_output1);

        print_bytes(fq_reads1[0].head());

        assert_eq!(fq_reads1.len(), 1);
        assert_eq!(
            fq_reads1[0],
            OwnedRecord {
                head: b"ex_0 1:N:0:AAAAAAAAGATTACAGA".to_vec(),
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
                head: b"ex_0 2:N:0:AAAAAAAAGATTACAGA".to_vec(),
                seq: "T".repeat(100).as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            }
        );
    }

    // TODO - expand this test to test molecular barcode once we add that info a la fgbio version
    // of this test.
    #[rstest]
    #[case(true)]
    #[case(false)]
    fn test_demux_a_wierd_set_of_reads(#[case] no_cache: bool) {
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
            compression_level: 5,
            no_metrics: true,
            metrics_filename: "fqtk-demux-metrics.tsv".to_string(),
            no_cache,
            no_molecular_barcode_in_header: false,
            no_sample_barcode_in_header: false,
        };
        demux_inputs.execute().unwrap();

        let expected_output1 = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads1 = read_fastq(&expected_output1);

        print_bytes(fq_reads1[0].head());

        assert_eq!(fq_reads1.len(), 1);
        assert_eq!(
            fq_reads1[0],
            OwnedRecord {
                head: b"ex_0:::::::CCCC+A 1:N:0:AAAAAAAAGATTACAGA".to_vec(),
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
                head: b"ex_0:::::::CCCC+A 2:N:0:AAAAAAAAGATTACAGA".to_vec(),
                seq: "T".as_bytes().to_vec(),
                qual: ";".as_bytes().to_vec(),
            }
        );
    }

    #[rstest]
    #[case(true)]
    #[case(false)]
    fn test_demux_a_read_structure_with_multiple_templates_in_one_read(#[case] no_cache: bool) {
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
            compression_level: 5,
            no_metrics: true,
            metrics_filename: "fqtk-demux-metrics.tsv".to_string(),
            no_cache,
            no_molecular_barcode_in_header: false,
            no_sample_barcode_in_header: false,
        };
        demux_inputs.execute().unwrap();

        let expected_output1 = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads1 = read_fastq(&expected_output1);

        print_bytes(fq_reads1[0].head());

        assert_eq!(fq_reads1.len(), 1);
        assert_eq!(
            fq_reads1[0],
            OwnedRecord {
                head: b"ex_0 1:N:0:AAAAAAAAGATTACAGA".to_vec(),
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
                head: b"ex_0 2:N:0:AAAAAAAAGATTACAGA".to_vec(),
                seq: "T".repeat(20).as_bytes().to_vec(),
                qual: ";".repeat(20).as_bytes().to_vec(),
            }
        );
        let expected_output3 = output_dir.join("Sample0000.R3.fq.gz");
        let fq_reads3 = read_fastq(&expected_output3);

        assert_eq!(fq_reads3.len(), 1);
        assert_eq!(
            fq_reads3[0],
            OwnedRecord {
                head: b"ex_0 3:N:0:AAAAAAAAGATTACAGA".to_vec(),
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
            compression_level: 5,
            no_metrics: true,
            metrics_filename: "fqtk-demux-metrics.tsv".to_string(),
            no_cache: false,
            no_molecular_barcode_in_header: false,
            no_sample_barcode_in_header: false,
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
            compression_level: 5,
            no_metrics: true,
            metrics_filename: "fqtk-demux-metrics.tsv".to_string(),
            no_cache: false,
            no_molecular_barcode_in_header: false,
            no_sample_barcode_in_header: false,
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
            compression_level: 5,
            no_metrics: true,
            metrics_filename: "fqtk-demux-metrics.tsv".to_string(),
            no_cache: false,
            no_molecular_barcode_in_header: false,
            no_sample_barcode_in_header: false,
        };
        demux_inputs.execute().unwrap();
    }

    #[test]
    fn test_molecular_barcode_header_component() {
        let mol_set = ReadSet {
            before_mol_barcode_header_chunk: b"doesnt matter".to_vec(),
            after_record_number_before_sample_barcode_header_chunk: b"doesnt matter".to_vec(),
            sample_barcode_sequence: b"doesnt matter".to_vec(),
            template_segments: None,
            sample_barcode_segments: None,
            molecular_barcode_segments: None,
        };

        assert_eq!(mol_set.molecular_barcode_header_component(), Vec::<u8>::new());

        let mol_set = ReadSet {
            before_mol_barcode_header_chunk: b"doesnt matter".to_vec(),
            after_record_number_before_sample_barcode_header_chunk: b"doesnt matter".to_vec(),
            sample_barcode_sequence: b"doesnt matter".to_vec(),
            template_segments: None,
            sample_barcode_segments: None,
            molecular_barcode_segments: Some(vec![FastqSegment {
                seq: b"GATTACA".to_vec(),
                quals: b";;;;;;;".to_vec(),
            }]),
        };

        assert_eq!(mol_set.molecular_barcode_header_component(), b":GATTACA");

        let mol_set = ReadSet {
            before_mol_barcode_header_chunk: b"doesnt matter".to_vec(),
            after_record_number_before_sample_barcode_header_chunk: b"doesnt matter".to_vec(),
            sample_barcode_sequence: b"doesnt matter".to_vec(),
            template_segments: None,
            sample_barcode_segments: None,
            molecular_barcode_segments: Some(vec![
                FastqSegment { seq: b"GATTACA".to_vec(), quals: b";;;;;;;".to_vec() },
                FastqSegment { seq: b"TAT".to_vec(), quals: b";;;".to_vec() },
            ]),
        };

        assert_eq!(mol_set.molecular_barcode_header_component(), b":GATTACA+TAT");

        let mol_set = ReadSet {
            before_mol_barcode_header_chunk: b"doesnt matter".to_vec(),
            after_record_number_before_sample_barcode_header_chunk: b"doesnt matter".to_vec(),
            sample_barcode_sequence: b"doesnt matter".to_vec(),
            template_segments: None,
            sample_barcode_segments: None,
            molecular_barcode_segments: Some(vec![
                FastqSegment { seq: b"GATTACA".to_vec(), quals: b";;;;;;;".to_vec() },
                FastqSegment { seq: b"TAT".to_vec(), quals: b";;;".to_vec() },
                FastqSegment { seq: b"T".to_vec(), quals: b";".to_vec() },
            ]),
        };

        assert_eq!(mol_set.molecular_barcode_header_component(), b":GATTACA+TAT+T");
    }

    #[test]
    fn test_barcode_insertion() {
        let tmp = TempDir::new().unwrap();
        let read_structures = vec![ReadStructure::from_str("8M17B100T").unwrap()];
        let mol_barcode = "TAGGAAAA";
        let s1_barcode = "AAAAAAAAGATTACAGA";
        let sample_metadata = metadata_file(
            &tmp,
            &[s1_barcode, "CCCCCCCCGATTACAGA", "GGGGGGGGGATTACAGA", "GGGGGGTTGATTACAGA"],
        );
        let input_files = vec![fastq_file(
            &tmp,
            "ex",
            "ex",
            &[&(mol_barcode.to_owned() + s1_barcode + &"A".repeat(100))],
        )];

        let output_dir = tmp.path().to_path_buf().join("output");

        // Inserting both
        let demux_inputs = Demux {
            inputs: input_files.clone(),
            read_structures: read_structures.clone(),
            sample_metadata: sample_metadata.clone(),
            output_types: vec!['T'],
            output: output_dir.clone(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 3,
            compression_level: 5,
            no_cache: false,
            no_metrics: true,
            metrics_filename: "fqtk-demux-metrics.tsv".to_string(),
            no_molecular_barcode_in_header: false,
            no_sample_barcode_in_header: false,
        };
        demux_inputs.execute().unwrap();

        let expected_output = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads = read_fastq(&expected_output);

        assert_eq!(fq_reads.len(), 1);
        assert_eq!(
            fq_reads[0],
            OwnedRecord {
                head: b"ex_0:::::::TAGGAAAA 1:N:0:AAAAAAAAGATTACAGA".to_vec(),
                seq: "A".repeat(100).as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            }
        );

        // Inserting only mol barcode
        let demux_inputs = Demux {
            inputs: input_files.clone(),
            read_structures: read_structures.clone(),
            sample_metadata: sample_metadata.clone(),
            output_types: vec!['T'],
            output: output_dir.clone(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 3,
            compression_level: 5,
            no_cache: false,
            no_metrics: true,
            metrics_filename: "fqtk-demux-metrics.tsv".to_string(),
            no_molecular_barcode_in_header: false,
            no_sample_barcode_in_header: true,
        };
        demux_inputs.execute().unwrap();

        let expected_output = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads = read_fastq(&expected_output);

        assert_eq!(fq_reads.len(), 1);
        print_bytes(fq_reads[0].head());
        assert_eq!(
            fq_reads[0],
            OwnedRecord {
                head: b"ex_0:::::::TAGGAAAA".to_vec(),
                seq: "A".repeat(100).as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            }
        );

        // Inserting only sample barcode
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
            compression_level: 5,
            no_cache: false,
            no_metrics: true,
            metrics_filename: "fqtk-demux-metrics.tsv".to_string(),
            no_molecular_barcode_in_header: true,
            no_sample_barcode_in_header: false,
        };
        demux_inputs.execute().unwrap();

        let expected_output = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads = read_fastq(&expected_output);

        assert_eq!(fq_reads.len(), 1);
        assert_eq!(
            fq_reads[0],
            OwnedRecord {
                head: b"ex_0 1:N:0:AAAAAAAAGATTACAGA".to_vec(),
                seq: "A".repeat(100).as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            }
        );
    }

    #[test]
    fn test_determine_boundaries_in_existing_header() {
        //                                    will trim here if setting new mol barcode, otherwise at end
        //                                    |
        //                                    |       extra colon
        //                                    |       |
        //                                    v       v
        let existing: &[u8; 19] = b"ex_0:::::::GATTACA:";
        let (end, num_colons, opt_after_space) =
            ReadSet::determine_boundaries_in_existing_header(existing, false);
        assert_eq!(existing.len(), end);
        assert_eq!(bytecount::count(existing, b':'), num_colons);
        assert!(opt_after_space.is_none());

        let (end, num_colons, opt_after_space) =
            ReadSet::determine_boundaries_in_existing_header(existing, true);
        assert_eq!(10, end);
        assert_eq!(bytecount::count(existing, b':'), num_colons);
        assert!(opt_after_space.is_none());

        let existing: &[u8; 29] = b"ex_0:::::::GATTACA: 1:T:1:AAA";
        if let (end, num_colons_before_space, Some((after_space_start, after_space_end))) =
            ReadSet::determine_boundaries_in_existing_header(existing, false)
        {
            assert_eq!("ex_0:::::::GATTACA:".len(), end);
            assert_eq!(num_colons_before_space, 8);
            assert_eq!(&existing[after_space_start..after_space_end], b":T:1:");
        } else {
            panic!("error reading after space")
        }
    }

    #[test]
    fn test_process_existing_read_header() {
        let existing: &[u8; 19] = b"ex_0:::::::GATTACA:";
        let (start, end) = ReadSet::process_existing_read_header(existing, true, true, true);
        assert_eq!(&start, &existing[..10]);
        assert_eq!(&end, b":N:0:");

        // padding bool doesnt affect it if theres no need for new padding.
        let (start, end) = ReadSet::process_existing_read_header(existing, false, true, true);
        assert_eq!(&start, &existing[..10]);
        assert_eq!(&end, b":N:0:");

        // not adding mol barcode to header reads entire before space string.
        let (start, end) = ReadSet::process_existing_read_header(existing, true, false, true);
        assert_eq!(&start, &existing);
        assert_eq!(&end, b":N:0:");

        // padding bool doesnt affect it if theres no need for new padding
        let (start, end) = ReadSet::process_existing_read_header(existing, false, false, true);
        assert_eq!(&start, &existing);
        assert_eq!(&end, b":N:0:");

        // Not adding sample barcode will produce empty slice if nothing after slice
        let (_, end) = ReadSet::process_existing_read_header(existing, true, true, false);
        assert_eq!(&end, b"");

        let existing: &[u8; 33] = b"ex_0::::::test:GATTACA: 1:T:1:AAA";
        // Not adding sample barcode will trim read num but otherwise leave after space untouched
        // if something after slice
        let (_, end) = ReadSet::process_existing_read_header(existing, true, true, false);
        assert_eq!(&end, &existing[25..]);

        let short_existing: &[u8; 8] = b"ex_0::::";

        // if not padding will not add additional colons.
        let (start, _) = ReadSet::process_existing_read_header(short_existing, false, true, true);
        assert_eq!(&start, &short_existing);

        let (start, _) = ReadSet::process_existing_read_header(short_existing, false, false, true);
        print_bytes(&start);
        assert_eq!(&start, &short_existing);

        // if not adding new molecular barcode also will not add additional colons.
        let (start, _) = ReadSet::process_existing_read_header(short_existing, true, false, true);
        assert_eq!(&start, &short_existing);

        // if not padding and theres a barcode to add will add additional colons.
        let (start, end) = ReadSet::process_existing_read_header(short_existing, true, true, true);
        assert_eq!(&start, &existing[..10]);
        assert_eq!(&end, b":N:0:");
    }

    #[test]
    #[should_panic(expected = "Invalid FASTQ contents after space")]
    fn test_process_existing_read_header_panics_if_improper_after_space_segment() {
        let existing: &[u8; 42] = b"ex_0::::::test:GATTACA: not_valid_illumina";
        // Not adding sample barcode will trim read num but otherwise leave after space untouched
        // if something after slice
        let (_, _) = ReadSet::process_existing_read_header(existing, true, true, false);
    }
}
