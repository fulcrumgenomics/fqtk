use crate::commands::command::Command;
use anyhow::{anyhow, ensure, Result};
use bstr::ByteSlice;
use clap::Parser;
use fgoxide::io::{DelimFile, Io};
use fgoxide::iter::IntoChunkedReadAheadIterator;
use fqtk_lib::barcode_matching::BarcodeMatcher;
use fqtk_lib::samples::SampleGroup;
use itertools::Itertools;
use log::info;
use pooled_writer::{bgzf::BgzfCompressor, Pool, PoolBuilder, PooledWriter};
use proglog::{CountFormatterKind, ProgLogBuilder};
use read_structure::ReadStructure;
use read_structure::ReadStructureError;
use read_structure::SegmentType;
use seq_io::fastq::Reader as FastqReader;
use seq_io::fastq::Record;
use serde::{Deserialize, Serialize};
use std::collections::HashSet;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::iter::Filter;
use std::slice::Iter;
use std::{
    fs,
    path::{Path, PathBuf},
};

/// Type alias to prevent clippy complaining about type complexity
type VecOfReaders = Vec<Box<dyn BufRead + Send>>;

/// Type alias for segment type iter functions, which iterate over the segments of a ``ReadSet``
/// filtering for a specific type.
type SegmentTypeIter<'a> = Filter<Iter<'a, FastqSegment>, fn(&&FastqSegment) -> bool>;

const BUFFER_SIZE: usize = 1024 * 1024;

/// The bases and qualities associated with a segment of a FASTQ record.
#[derive(Debug, Clone)]
struct FastqSegment {
    /// bases of the FASTQ subsection
    seq: Vec<u8>,
    /// qualities of the FASTQ subsection
    quals: Vec<u8>,
    /// the type of segment being stored
    segment_type: SegmentType,
}

////////////////////////////////////////////////////////////////////////////////
// ReadSet and it's impls
////////////////////////////////////////////////////////////////////////////////

/// One unit of FASTQ records separated into their component read segments.
#[derive(Debug, Clone)]
struct ReadSet {
    /// Header of the FASTQ record
    header: Vec<u8>,
    /// Segments of reads
    segments: Vec<FastqSegment>,
}

impl ReadSet {
    const PREFIX: u8 = b'@';
    const SPACE: u8 = b' ';
    const COLON: u8 = b':';
    const PLUS: u8 = b'+';

    /// Produces an iterator over references to the template segments stored in this ``ReadSet``.
    fn template_segments(&self) -> SegmentTypeIter {
        self.segments.iter().filter(|s| s.segment_type == SegmentType::Template)
    }

    /// Produces an iterator over references to the sample barcode segments stored in this
    /// ``ReadSet``.
    fn sample_barcode_segments(&self) -> SegmentTypeIter {
        self.segments.iter().filter(|s| s.segment_type == SegmentType::SampleBarcode)
    }

    /// Produces an iterator over references to the molecular barcode segments stored in this
    /// ``ReadSet``.
    fn molecular_barcode_segments(&self) -> SegmentTypeIter {
        self.segments.iter().filter(|s| s.segment_type == SegmentType::MolecularBarcode)
    }

    /// Generates the sample barcode sequence for this read set and returns it as a Vec of bytes.
    fn sample_barcode_sequence(&self) -> Vec<u8> {
        self.sample_barcode_segments().flat_map(|s| &s.seq).copied().collect()
    }

    /// Combines ``ReadSet`` structs together into a single ``ReadSet``
    fn combine_readsets(readsets: Vec<Self>) -> Self {
        assert!(!readsets.is_empty(), "Cannot call combine readsets on an empty vec!");
        let mut readset_iter = readsets.into_iter();
        let mut first = readset_iter.next().expect("Cannot call combine readsets on an empty vec!");
        for next_readset in readset_iter {
            first.segments.extend(next_readset.segments);
        }
        first
    }

    /// Writes the FASTQ header to the given writer.  Substitutes in the given read number into
    /// the comment section.  Also adds in the UMI(s) from the UMI segments and the sample barcodes
    /// into the appropriate places in the header.
    ///
    /// UMI and sample barcode segments and (separately) concatenated with `+`s in between. If
    /// there is an existing UMI or sample barcode in the header, the segments are appended after
    /// adding a `+`, else the segments are inserted.
    ///
    /// Supports headers that are just the `name` segment, or `name comment`.  The name must have
    /// at most 8 colon-separated parts.  If there are seven or fewer parts in the name, the
    /// UMI is appended as the last part.  If there are eight, the eighth is assumed to be an
    /// existing UMI, and any UMI is appended.
    ///
    /// If the comment is present is must have exactly four colon-separated parts.
    ///
    /// Format of the header is:
    ///   @name comment
    /// Where
    ///   name = @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI>
    ///   comment = <read>:<is filtered>:<control number>:<index>
    fn write_header<W: Write>(&self, writer: &mut W, read_num: usize) -> Result<()> {
        Self::write_header_internal(
            writer,
            read_num,
            self.header.as_slice(),
            self.sample_barcode_segments(),
            self.molecular_barcode_segments(),
        )
    }

    fn write_header_internal<W: Write>(
        writer: &mut W,
        read_num: usize,
        header: &[u8],
        sample_barcode_segments: SegmentTypeIter,
        mut molecular_barcode_segments: SegmentTypeIter,
    ) -> Result<()> {
        // Extract the name and optionally the comment
        let (name, comment) = match header.find_byte(Self::SPACE) {
            Some(x) => (&header[0..x], Some(&header[x + 1..])),
            None => (header, None),
        };

        writer.write_all(&[Self::PREFIX])?;

        // Handle the 'name' component of the header.  If we don't have any UMI segments
        // we can emit the name part as is.  Otherwise we need to append the UMIs to the name.
        if let Some(first_seg) = molecular_barcode_segments.next() {
            let sep_count = name.iter().filter(|c| **c == Self::COLON).count();
            ensure!(
                sep_count <= 7,
                "Can't handle read name with more than 8 segments: {}",
                String::from_utf8(header.to_vec())?
            );

            writer.write_all(name)?;
            if sep_count == 7 {
                // UMI already present, append to it with a UMI separator first
                writer.write_all(&[Self::PLUS])?;
            } else {
                // UMI not present yet, insert a minor separator
                writer.write_all(&[Self::COLON])?;
            }

            writer.write_all(first_seg.seq.as_slice())?;
            // Append all the UMI segments with pluses in between.
            for seg in molecular_barcode_segments {
                writer.write_all(&[Self::PLUS])?;
                writer.write_all(seg.seq.as_slice())?;
            }
        } else {
            writer.write_all(name)?;
        }

        writer.write_all(&[Self::SPACE])?;

        // Then the 'comment' section
        match comment {
            None => {
                // If no pre-existing comment, assume the read is a passing filter, non-control
                // read and generate a comment for it (sample barcode is added below).
                write!(writer, "{}:N:0:", read_num)?;
            }
            Some(chars) => {
                // Else check it's a 4-part name... fix the read number at the front and
                // check to see if there's a real sample barcode on the back
                let sep_count = chars.iter().filter(|c| **c == Self::COLON).count();
                ensure!(
                    sep_count == 3,
                    "Comment in did not have 4 segments: {}",
                    String::from_utf8(header.to_vec())?
                );
                let first_colon_idx = chars.iter().position(|ch| *ch == Self::COLON).unwrap();

                // Illumina, in the unmatched FASTQs, can place a "0" in the index position, sigh
                let remainder = if chars.last().unwrap().is_ascii_digit() {
                    &chars[first_colon_idx + 1..chars.len() - 1]
                } else {
                    &chars[first_colon_idx + 1..chars.len()]
                };

                write!(writer, "{}:", read_num)?;
                writer.write_all(remainder)?;

                if *remainder.last().unwrap() != Self::COLON {
                    writer.write_all(&[Self::PLUS])?;
                }
            }
        }

        // Append all the sample barcode segments to the new comment
        for (idx, seg) in sample_barcode_segments.enumerate() {
            if idx > 0 {
                writer.write_all(&[Self::PLUS])?;
            }
            writer.write_all(seg.seq.as_slice())?;
        }

        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////
// ReadSetIterator and it's impls
////////////////////////////////////////////////////////////////////////////////

/// A struct for iterating over the records in multiple FASTQ files simultaneously, destructuring
/// them according to the provided read structures, yielding ``ReadSet`` objects on each iteration.
struct ReadSetIterator {
    /// Read structure object describing the structure of the reads in this file.
    read_structure: ReadStructure,
    /// The file containing FASTQ records.
    source: FastqReader<Box<dyn BufRead + Send>>,
}

impl Iterator for ReadSetIterator {
    type Item = ReadSet;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(rec) = self.source.next() {
            let mut segments = Vec::with_capacity(self.read_structure.number_of_segments());
            let next_fq_rec = rec.expect("Unexpected error parsing FASTQs");
            let read_name = next_fq_rec.head().to_vec();

            for read_segment in self.read_structure.iter() {
                let (seq, quals) = read_segment
                    .extract_bases_and_quals(next_fq_rec.seq(), next_fq_rec.qual())
                    .unwrap_or_else(|e| {
                        panic!(
                            "Error extracting read segment bases or quals from FASTQ record {}; {}",
                            String::from_utf8(next_fq_rec.head().to_vec()).unwrap(),
                            e
                        )
                    });
                segments.push(FastqSegment {
                    seq: seq.to_vec(),
                    quals: quals.to_vec(),
                    segment_type: read_segment.kind,
                });
            }
            Some(ReadSet { header: read_name, segments })
        } else {
            None
        }
    }
}

impl ReadSetIterator {
    /// Instantiates a new iterator over the read sets for a set of FASTQs with defined read
    /// structures
    pub fn new(
        read_structure: ReadStructure,
        source: FastqReader<Box<dyn BufRead + Send>>,
    ) -> Self {
        Self { read_structure, source }
    }
}

////////////////////////////////////////////////////////////////////////////////
// SampleWriters and it's impls
////////////////////////////////////////////////////////////////////////////////

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
    fn write(&mut self, read_set: &ReadSet) -> Result<()> {
        for (writers_opt, segments) in [
            (&mut self.template_writers, &mut read_set.template_segments()),
            (&mut self.sample_barcode_writers, &mut read_set.template_segments()),
            (&mut self.molecular_barcode_writers, &mut read_set.template_segments()),
        ] {
            if let Some(writers) = writers_opt {
                for (read_idx, (writer, segment)) in writers.iter_mut().zip(segments).enumerate() {
                    read_set.write_header(writer, read_idx + 1)?;
                    writer.write_all(b"\n")?;
                    writer.write_all(segment.seq.as_slice())?;
                    writer.write_all(b"\n+\n")?;
                    writer.write_all(segment.quals.as_slice())?;
                    writer.write_all(b"\n")?;
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

////////////////////////////////////////////////////////////////////////////////
// DemuxMetric and it's impls
////////////////////////////////////////////////////////////////////////////////

/// A set of metrics for a single sample from demultiplexing.  "Template" in this context
/// refers to the set of reads that share a read name - i.e. a set of one read each from all
/// of the input FASTQ files.
///
/// The `ratio_*` fields are calculated using the mean and max template counts across all samples
/// *excluding* the unmatched pseudo-sample.
#[allow(clippy::module_name_repetitions)]
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct DemuxMetric {
    /// The ID of the sample being reported on.
    sample_id: String,
    /// The expected barcode sequence associated with the sample.
    barcode: String,
    /// The number of templates (querynames, inserts) assigned to the sample.
    templates: usize,
    /// The fraction of all templates in the input that were assigned to the sample.
    frac_templates: f64,
    /// The ratio of this sample's `templates` to the mean across all samples.
    ratio_to_mean: f64,
    /// The ratio of this sample's `templates` to the best (max) across all samples.
    ratio_to_best: f64,
}

impl DemuxMetric {
    /// Create a new ``DemuxMetric`` with the given sample id and barcode.
    fn new(sample: &str, barcode: &str) -> DemuxMetric {
        DemuxMetric {
            sample_id: sample.to_string(),
            barcode: barcode.to_string(),
            templates: 0,
            frac_templates: 0.0,
            ratio_to_mean: 0.0,
            ratio_to_best: 0.0,
        }
    }

    /// Update all the derived fields in all the provided metrics objects.
    fn update(samples: &mut [DemuxMetric], unmatched: &mut DemuxMetric) {
        let sample_total: f64 = samples.iter().map(|s| s.templates as f64).sum();
        let total = sample_total + unmatched.templates as f64;
        let mean = sample_total / samples.len() as f64;
        let best = samples.iter().map(|s| s.templates).max().unwrap_or(0) as f64;

        for sample in samples {
            sample.frac_templates = sample.templates as f64 / total;
            sample.ratio_to_mean = sample.templates as f64 / mean;
            sample.ratio_to_best = sample.templates as f64 / best;
        }

        unmatched.frac_templates = unmatched.templates as f64 / total;
        unmatched.ratio_to_mean = unmatched.templates as f64 / mean;
        unmatched.ratio_to_best = unmatched.templates as f64 / best;
    }
}

////////////////////////////////////////////////////////////////////////////////
// Demux (main class) and it's impls
////////////////////////////////////////////////////////////////////////////////

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
/// 1. `sample_id` - the id of the sample or library.
/// 2. `barcode` - the expected barcode sequence associated with the `sample_id`.
///
/// The read structures will be used to extract the observed sample barcode, template bases, and
/// molecular identifiers from each read.  The observed sample barcode will be matched to the
/// sample barcodes extracted from the bases in the sample metadata and associated read structures.
///
/// ## Outputs
///
/// All outputs are generated in the provided `--output` directory.  For each sample plus the
/// unmatched reads, FASTQ files are written for each read segment (specified in the read
/// structures) of one of the types supplied to `--output-types`.  FASTQ files have names
/// of the format:
///
/// ```
/// {sample_id}.{segment_type}{read_num}.fq.gz
/// ```
///
/// where `segment_type` is one of `R`, `I`, and `U` (for template, barcode/index and molecular
/// barcode/UMI reads respectively) and `read_num` is a number starting at 1 for each segment
/// type.
///
/// In addition a `demux-metrics.txt` file is written that is a tab-delimited file with counts
/// of how many reads were assigned to each sample and derived metrics.
///
/// ## Example Command Line
///
/// As an example, if the sequencing run was 2x100bp (paired end) with two 8bp index reads both
/// reading a sample barcode, as well as an in-line 8bp sample barcode in read one, the command
/// line would be:
///
/// ```
/// fqtk demux \
///     --inputs r1.fq.gz i1.fq.gz i2.fq.gz r2.fq.gz \
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
    #[clap(long, short = 't', default_value = "5")]
    threads: usize,

    /// The level of compression to use to compress outputs.
    #[clap(long, short = 'c', default_value = "5")]
    compression_level: usize,

    /// The chunk size to use when reading FASTQ files
    // TODO remove me when we've figured out good defaults
    #[clap(long, default_value = "1000")]
    chunk_size: usize,

    /// The buffer size to use when reading ahead on FASTQ files
    // TODO remove me when we've figured out good defaults
    #[clap(long, default_value = "1000")]
    buffer_size: usize,

    /// If true, caching of barcodes will be disabled
    #[clap(long, hide = true)]
    no_cache: bool,
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
                &sample.sample_id,
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
    fn build_writer_pool(
        sample_writers: Vec<SampleWriters<BufWriter<File>>>,
        compression_level: usize,
        threads: usize,
    ) -> Result<(Pool, Vec<SampleWriters<PooledWriter>>)> {
        let mut new_sample_writers = Vec::with_capacity(sample_writers.len());
        let mut pool_builder = PoolBuilder::<_, BgzfCompressor>::new()
            .threads(threads)
            .queue_size(threads * 50)
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
        Ok((pool, new_sample_writers))
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
        let fgio = Io::new(5, BUFFER_SIZE);
        let fq_readers_result = self
            .inputs
            .iter()
            .map(|p| fgio.new_reader(p))
            .collect::<Result<VecOfReaders, fgoxide::FgError>>();
        if let Err(e) = &fq_readers_result {
            constraint_errors.push(format!("Error opening input files for reading: {}", e));
        }

        if self.threads < 5 {
            constraint_errors
                .push(format!("Threads provided {} was too low! Must be 5 or more.", self.threads));
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
        let fq_sources =
            fq_readers.into_iter().map(|fq| FastqReader::with_capacity(fq, BUFFER_SIZE));

        let (mut pool, mut sample_writers) = Self::build_writer_pool(
            Self::create_writers(
                &self.read_structures,
                &sample_group,
                &output_segment_types,
                &self.unmatched_prefix,
                &self.output,
            )?,
            self.compression_level,
            // reserve 1 for main thread and 2 for read ahead, remaining on writing.
            self.threads - 3,
        )?;
        let unmatched_index = sample_writers.len() - 1;
        info!("Created sample and {} writers.", self.unmatched_prefix);

        // Create the metrics as a Vec that is index sync'd with the sample group
        let mut sample_metrics = sample_group
            .samples
            .iter()
            .map(|s| DemuxMetric::new(s.sample_id.as_str(), s.barcode.as_str()))
            .collect_vec();
        let mut unmatched_metric = DemuxMetric::new(self.unmatched_prefix.as_str(), ".");

        // Setup the barcode matcher - the primary return from here is the index of the samp
        let mut barcode_matcher = BarcodeMatcher::new(
            &sample_group.samples.iter().map(|s| s.barcode.as_str()).collect::<Vec<_>>(),
            u8::try_from(self.max_mismatches)?,
            u8::try_from(self.min_mismatch_delta)?,
            !self.no_cache,
        );

        let mut fq_iterators = fq_sources
            .zip(self.read_structures.clone().into_iter())
            .map(|(source, read_structure)| {
                ReadSetIterator::new(read_structure, source)
                    .read_ahead(self.chunk_size, self.buffer_size)
            })
            .collect::<Vec<_>>();

        let logger = ProgLogBuilder::new()
            .name("fqtk")
            .noun("records")
            .verb("demultiplexed")
            .unit(1_000_000)
            .count_formatter(CountFormatterKind::Comma)
            .level(log::Level::Info)
            .build();

        loop {
            let mut next_read_sets = Vec::with_capacity(fq_iterators.len());
            for iter in &mut fq_iterators {
                if let Some(rec) = iter.next() {
                    next_read_sets.push(rec);
                }
            }
            if next_read_sets.is_empty() {
                break;
            }
            assert_eq!(
                next_read_sets.len(),
                fq_iterators.len(),
                "FASTQ sources out of sync at records: {:?}",
                next_read_sets
            );
            let read_set = ReadSet::combine_readsets(next_read_sets);
            if let Some(barcode_match) = barcode_matcher.assign(&read_set.sample_barcode_sequence())
            {
                sample_writers[barcode_match.best_match].write(&read_set)?;
                sample_metrics[barcode_match.best_match].templates += 1;
            } else {
                sample_writers[unmatched_index].write(&read_set)?;
                unmatched_metric.templates += 1;
            }
            logger.record();
        }

        // Shut down the pool
        info!("Finished reading input FASTQs.");
        sample_writers.into_iter().map(SampleWriters::close).collect::<Result<Vec<_>>>()?;
        pool.stop_pool()?;
        info!("Output FASTQ writing complete.");

        // Write out the metrics
        let metrics_path = self.output.join("demux-metrics.txt");
        DemuxMetric::update(sample_metrics.as_mut_slice(), &mut unmatched_metric);
        sample_metrics.push(unmatched_metric);
        DelimFile::default().write_tsv(&metrics_path, sample_metrics)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fqtk_lib::samples::Sample;
    use rstest::rstest;
    use seq_io::fastq::OwnedRecord;
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

    fn assert_equal(actual: &impl seq_io::fastq::Record, expected: &impl seq_io::fastq::Record) {
        assert_eq!(
            String::from_utf8(actual.head().to_vec()).unwrap(),
            String::from_utf8(expected.head().to_vec()).unwrap()
        );

        assert_eq!(
            String::from_utf8(actual.seq().to_vec()).unwrap(),
            String::from_utf8(expected.seq().to_vec()).unwrap()
        );

        assert_eq!(
            String::from_utf8(actual.qual().to_vec()).unwrap(),
            String::from_utf8(expected.qual().to_vec()).unwrap()
        );
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
            threads: 5,
            compression_level: 5,
            chunk_size: 1000,
            buffer_size: 1000,
            no_cache: false,
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
            threads: 5,
            compression_level: 5,
            chunk_size: 1000,
            buffer_size: 1000,
            no_cache: false,
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
            threads: 5,
            compression_level: 5,
            chunk_size: 1000,
            buffer_size: 1000,
            no_cache: false,
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
            chunk_size: 1000,
            buffer_size: 1000,
            no_cache: false,
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
            chunk_size: 1000,
            buffer_size: 1000,
            no_cache: false,
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
            threads: 5,
            compression_level: 5,
            chunk_size: 1000,
            buffer_size: 1000,
            no_cache,
        };
        demux_inputs.execute().unwrap();

        let output_path = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads = read_fastq(&output_path);

        assert_eq!(fq_reads.len(), 1);
        assert_equal(
            &fq_reads[0],
            &OwnedRecord {
                head: b"ex_0 1:N:0:AAAAAAAAGATTACAGA".to_vec(),
                seq: "A".repeat(100).as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            },
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
            threads: 5,
            compression_level: 5,
            chunk_size: 1000,
            buffer_size: 1000,
            no_cache,
        };
        demux_inputs.execute().unwrap();

        let r1_path = output_dir.join("Sample0000.R1.fq.gz");
        let r1_reads = read_fastq(&r1_path);

        assert_eq!(r1_reads.len(), 1);
        assert_equal(
            &r1_reads[0],
            &OwnedRecord {
                head: b"ex_0 1:N:0:AAAAAAAA+GATTACAGA".to_vec(),
                seq: "A".repeat(100).as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            },
        );

        let r2_path = output_dir.join("Sample0000.R2.fq.gz");
        let r2_reads = read_fastq(&r2_path);

        assert_eq!(r2_reads.len(), 1);
        assert_equal(
            &r2_reads[0],
            &OwnedRecord {
                head: b"ex_0 2:N:0:AAAAAAAA+GATTACAGA".to_vec(),
                seq: "T".repeat(100).as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            },
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

        let demux = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output_types: vec!['T'],
            output: output_dir.clone(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 5,
            compression_level: 5,
            chunk_size: 1000,
            buffer_size: 1000,
            no_cache,
        };
        demux.execute().unwrap();

        let r1_path = output_dir.join("Sample0000.R1.fq.gz");
        let r1_reads = read_fastq(&r1_path);

        assert_eq!(r1_reads.len(), 1);
        assert_equal(
            &r1_reads[0],
            &OwnedRecord {
                head: b"ex_0 1:N:0:AAAAAAAA+GATTACAGA".to_vec(),
                seq: "A".repeat(100).as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            },
        );
        let r2_path = output_dir.join("Sample0000.R2.fq.gz");
        let r2_reads = read_fastq(&r2_path);

        assert_eq!(r2_reads.len(), 1);
        assert_equal(
            &r2_reads[0],
            &OwnedRecord {
                head: b"ex_0 2:N:0:AAAAAAAA+GATTACAGA".to_vec(),
                seq: "T".repeat(100).as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            },
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

        let demux = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output_types: vec!['T'],
            output: output_dir.clone(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 5,
            compression_level: 5,
            chunk_size: 1000,
            buffer_size: 1000,
            no_cache,
        };
        demux.execute().unwrap();

        let r1_path = output_dir.join("Sample0000.R1.fq.gz");
        let r1_reads = read_fastq(&r1_path);

        assert_eq!(r1_reads.len(), 1);
        assert_equal(
            &r1_reads[0],
            &OwnedRecord {
                head: b"ex_0:CCCC+A 1:N:0:AAAA+AAAA+GAT+TACAGA".to_vec(),
                seq: "A".repeat(100).as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            },
        );
        let r2_path = output_dir.join("Sample0000.R2.fq.gz");
        let r2_reads = read_fastq(&r2_path);

        assert_eq!(r2_reads.len(), 1);
        assert_equal(
            &r2_reads[0],
            &OwnedRecord {
                head: b"ex_0:CCCC+A 2:N:0:AAAA+AAAA+GAT+TACAGA".to_vec(),
                seq: "T".as_bytes().to_vec(),
                qual: ";".as_bytes().to_vec(),
            },
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
            threads: 5,
            compression_level: 5,
            chunk_size: 1000,
            buffer_size: 1000,
            no_cache,
        };
        demux_inputs.execute().unwrap();

        let r1_path = output_dir.join("Sample0000.R1.fq.gz");
        let r1_reads = read_fastq(&r1_path);
        assert_eq!(r1_reads.len(), 1);
        assert_equal(
            &r1_reads[0],
            &OwnedRecord {
                head: b"ex_0 1:N:0:AAAAAAAAGATTACAGA".to_vec(),
                seq: "A".repeat(20).as_bytes().to_vec(),
                qual: ";".repeat(20).as_bytes().to_vec(),
            },
        );

        let r2_path = output_dir.join("Sample0000.R2.fq.gz");
        let r2_reads = read_fastq(&r2_path);
        assert_eq!(r2_reads.len(), 1);
        assert_equal(
            &r2_reads[0],
            &OwnedRecord {
                head: b"ex_0 2:N:0:AAAAAAAAGATTACAGA".to_vec(),
                seq: "T".repeat(20).as_bytes().to_vec(),
                qual: ";".repeat(20).as_bytes().to_vec(),
            },
        );

        let r3_path = output_dir.join("Sample0000.R3.fq.gz");
        let r3_reads = read_fastq(&r3_path);

        assert_eq!(r3_reads.len(), 1);
        assert_equal(
            &r3_reads[0],
            &OwnedRecord {
                head: b"ex_0 3:N:0:AAAAAAAAGATTACAGA".to_vec(),
                seq: "G".repeat(20).as_bytes().to_vec(),
                qual: ";".repeat(20).as_bytes().to_vec(),
            },
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
            threads: 5,
            compression_level: 5,
            chunk_size: 1000,
            buffer_size: 1000,
            no_cache: false,
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
            threads: 5,
            compression_level: 5,
            chunk_size: 1000,
            buffer_size: 1000,
            no_cache: false,
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
            threads: 5,
            compression_level: 5,
            chunk_size: 1000,
            buffer_size: 1000,
            no_cache: false,
        };
        demux_inputs.execute().unwrap();
    }

    ////////////////////////////////////////////////////////////////////////////
    // Tests for ReadSet::write_header_internal
    ////////////////////////////////////////////////////////////////////////////

    fn seg(bases: &[u8], segment_type: SegmentType) -> FastqSegment {
        let quals = vec![b'#'; bases.len()];
        FastqSegment { seq: bases.to_vec(), quals, segment_type }
    }

    #[test]
    fn test_write_header_standard_no_umi() {
        let mut out = Vec::new();
        let header = b"inst:123:ABCDE:1:204:1022:2108 1:N:0:0";
        let barcode_segs =
            [seg(b"ACGT", SegmentType::SampleBarcode), seg(b"GGTT", SegmentType::SampleBarcode)];
        let umi_segs = [];
        let expected = "@inst:123:ABCDE:1:204:1022:2108 1:N:0:ACGT+GGTT".to_string();
        ReadSet::write_header_internal(
            &mut out,
            1,
            header,
            barcode_segs.iter().filter(|_| true),
            umi_segs.iter().filter(|_| true),
        )
        .unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    #[test]
    fn test_write_header_standard_with_umi() {
        let mut out = Vec::new();
        let header = b"inst:123:ABCDE:1:204:1022:2108 1:Y:0:0";
        let barcode_segs =
            [seg(b"ACGT", SegmentType::SampleBarcode), seg(b"GGTT", SegmentType::SampleBarcode)];
        let umi_segs = [seg(b"AACCGGTT", SegmentType::MolecularBarcode)];
        let expected = "@inst:123:ABCDE:1:204:1022:2108:AACCGGTT 2:Y:0:ACGT+GGTT".to_string();
        ReadSet::write_header_internal(
            &mut out,
            2,
            header,
            barcode_segs.iter().filter(|_| true),
            umi_segs.iter().filter(|_| true),
        )
        .unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    #[test]
    fn test_write_header_append_barcode_and_umi() {
        let mut out = Vec::new();
        let header = b"inst:123:ABCDE:1:204:1022:2108:AAAA 1:Y:0:TTTT";
        let barcode_segs =
            [seg(b"ACGT", SegmentType::SampleBarcode), seg(b"GGTT", SegmentType::SampleBarcode)];
        let umi_segs = [seg(b"AACCGGTT", SegmentType::MolecularBarcode)];
        let expected =
            "@inst:123:ABCDE:1:204:1022:2108:AAAA+AACCGGTT 2:Y:0:TTTT+ACGT+GGTT".to_string();
        ReadSet::write_header_internal(
            &mut out,
            2,
            header,
            barcode_segs.iter().filter(|_| true),
            umi_segs.iter().filter(|_| true),
        )
        .unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    #[test]
    fn test_write_header_short_name_no_comment() {
        let mut out = Vec::new();
        let header = b"q1";
        let barcode_segs =
            [seg(b"ACGT", SegmentType::SampleBarcode), seg(b"GGTT", SegmentType::SampleBarcode)];
        let umi_segs = [seg(b"AACCGGTT", SegmentType::MolecularBarcode)];
        let expected = "@q1:AACCGGTT 1:N:0:ACGT+GGTT".to_string();
        ReadSet::write_header_internal(
            &mut out,
            1,
            header,
            barcode_segs.iter().filter(|_| true),
            umi_segs.iter().filter(|_| true),
        )
        .unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    #[test]
    #[should_panic(expected = "8 segments")]
    fn test_write_header_name_too_many_parts() {
        let mut out = Vec::new();
        let header = b"q1:1:2:3:4:5:6:7:8:9:10";
        let barcode_segs =
            [seg(b"ACGT", SegmentType::SampleBarcode), seg(b"GGTT", SegmentType::SampleBarcode)];
        let umi_segs = [seg(b"AACCGGTT", SegmentType::MolecularBarcode)];
        ReadSet::write_header_internal(
            &mut out,
            1,
            header,
            barcode_segs.iter().filter(|_| true),
            umi_segs.iter().filter(|_| true),
        )
        .unwrap();
    }

    #[test]
    #[should_panic(expected = "4 segments")]
    fn test_write_header_comment_too_few_parts() {
        let mut out = Vec::new();
        let header = b"q1 0:0";
        let barcode_segs =
            [seg(b"ACGT", SegmentType::SampleBarcode), seg(b"GGTT", SegmentType::SampleBarcode)];
        let umi_segs = [seg(b"AACCGGTT", SegmentType::MolecularBarcode)];
        ReadSet::write_header_internal(
            &mut out,
            1,
            header,
            barcode_segs.iter().filter(|_| true),
            umi_segs.iter().filter(|_| true),
        )
        .unwrap();
    }
}
