use crate::commands::command::Command;
use anyhow::{Result, anyhow, ensure};
use bstr::ByteSlice;
use clap::{Parser, ValueEnum};
use fgoxide::io::{DelimFile, Io};
use fgoxide::iter::IntoChunkedReadAheadIterator;
use fqtk_lib::barcode_matching::BarcodeMatcher;
use fqtk_lib::samples::SampleGroup;
use itertools::Itertools;
use log::{info, warn};
use pooled_writer::{Pool, PoolBuilder, PooledWriter, bgzf::BgzfCompressor};
use proglog::{CountFormatterKind, ProgLogBuilder};
use read_structure::ReadStructure;
use read_structure::ReadStructureError;
use read_structure::SegmentType;
use seq_io::fastq::Reader as FastqReader;
use seq_io::fastq::Record;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fmt::Display;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::iter::Filter;
use std::slice::Iter;
use std::str::FromStr;
use std::{
    fs,
    path::{Path, PathBuf},
};

/// Type alias to prevent clippy complaining about type complexity
type VecOfReaders = Vec<Box<dyn BufRead + Send>>;

/// Type alias for segment type iter functions, which iterate over the segments of a ``ReadSet``
/// filtering for a specific type.
type SegmentIter<'a> = Filter<Iter<'a, FastqSegment>, fn(&&FastqSegment) -> bool>;

const BUFFER_SIZE: usize = 1024 * 1024;

/// The bases and qualities associated with a segment of a FASTQ record.
#[derive(PartialEq, Eq, Debug, Clone)]
struct FastqSegment {
    /// bases of the FASTQ subsection
    seq: Vec<u8>,
    /// qualities of the FASTQ subsection
    quals: Vec<u8>,
    /// the type of segment being stored
    segment_type: SegmentType,
    /// the index of the source FASTQ file this segment came from
    source_index: usize,
}

////////////////////////////////////////////////////////////////////////////////
// ReadSet and it's impls
////////////////////////////////////////////////////////////////////////////////

#[derive(Eq, Hash, PartialEq, Debug, Clone, Copy)]
enum SkipReason {
    /// The read had too few bases for the segment.
    TooFewBases,
}

/// The FASTQ header "shape" for each output record: whether the comment is rewritten, kept
/// verbatim, or dropped entirely.  Independent of this, `--umi-in-name` controls whether the
/// UMI(s) are folded into the read name; every format honors an explicit `--umi-in-name` override
/// that its own routing allows to take effect (see `Demux::umi_in_name`'s doc for the one
/// exception: `--template-types` can fold a UMI into the template bases instead, in which case
/// there is nothing left for `--umi-in-name` to fold into the name).
#[derive(clap::ValueEnum, Eq, Hash, PartialEq, Debug, Clone, Copy, Default)]
enum HeaderFormatKind {
    /// Casava ≥1.8 / bcl-convert-style header: the read-number field in the comment is rewritten
    /// to match the output file and the sample barcode(s) are appended to the comment.  UMI(s) are
    /// folded into the read name by default (override with `--umi-in-name`).  For example
    /// `@inst:1:FC:1:1101:5:7 1:N:0:0` becomes `@inst:1:FC:1:1101:5:7:AACCGGTT 1:N:0:ACGTACGT`.
    #[default]
    Illumina,
    /// Emit each read's original header verbatim: the read-number field is not rewritten and no
    /// sample barcode is appended.  UMI(s) are NOT folded into the name by default (override with
    /// `--umi-in-name true`).  Nothing is parsed unless a UMI is folded in, so with the default
    /// `--umi-in-name` this (and `name-only`, below) is one of the formats that passes through
    /// headers which do not follow Illumina conventions unchanged (more than eight `:`-delimited
    /// name fields, or a comment that is not four `:`-delimited fields), such as those produced by
    /// MGI, Element, and ONT instruments.
    Unmodified,
    /// Emit only the read name: any pre-existing comment is always dropped, regardless of
    /// `--umi-in-name`.  UMI(s) are NOT folded into the name by default (override with
    /// `--umi-in-name true`), in which case the same non-Illumina-header tolerance described for
    /// `unmodified` applies here too (the name is only parsed when a UMI is actually folded into
    /// it).  Intended for pipelines (e.g. `bwa mem -C`) where a Casava-style comment would break
    /// downstream parsing.  With no comment and (by default) no UMI in the name, R1/R2 headers
    /// become byte-identical (no read number anywhere) — the same tradeoff `samtools fastq -n`
    /// makes for file-paired workflows.
    NameOnly,
}

/// The full set of parameters governing how a single FASTQ header is written: whether UMIs are
/// eligible to be folded into the name, the header format, and whether the UMI is folded in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct HeaderFormat {
    /// When false, UMI (M) segments are never folded into the read name (they are written into the
    /// template bases instead, so the same bases are not emitted twice).
    include_umi: bool,
    /// The header "shape": comment rewrite/keep/drop behavior.
    kind: HeaderFormatKind,
    /// Whether to fold UMI(s) into the read name.  Already resolved from an explicit
    /// `--umi-in-name` (if given) or `kind`'s default otherwise -- by the time a `HeaderFormat` is
    /// built there is no more "unset" state.
    umi_in_name: bool,
}

impl HeaderFormatKind {
    /// Whether this format folds UMI(s) into the read name absent an explicit `--umi-in-name`.
    fn default_umi_in_name(self) -> bool {
        matches!(self, Self::Illumina)
    }

    /// True if the comment is rewritten (read-number field replaced, sample barcode appended)
    /// rather than kept verbatim or dropped.  Only `Illumina` does this.
    fn rewrites_comment(self) -> bool {
        matches!(self, Self::Illumina)
    }
}

impl Display for HeaderFormatKind {
    // Delegates to the `ValueEnum` derive's own name mapping (kebab-case of the variant name)
    // rather than hand-rolling a second copy of it, so the two can't drift out of sync -- this
    // string is also what `warn_on_discarded_segments` tells users to pass back on the CLI.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.to_possible_value().expect("no skipped variants").get_name())
    }
}

impl Display for SkipReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SkipReason::TooFewBases => write!(f, "Too few bases"),
        }
    }
}

impl FromStr for SkipReason {
    type Err = anyhow::Error;
    fn from_str(string: &str) -> Result<Self, Self::Err> {
        match string {
            "too few bases" | "too-few-bases" | "toofewbases" => Ok(SkipReason::TooFewBases),
            _ => Err(anyhow!("Invalid skip reason: {}", string)),
        }
    }
}

/// The raw bases and qualities of a single FASTQ record from one input file, retained so the
/// record can be re-segmented per-sample when per-sample read structures are in use.
#[derive(PartialEq, Eq, Debug, Clone)]
struct RawInput {
    bases: Vec<u8>,
    quals: Vec<u8>,
}

/// Error returned by [`ReadSet::resegment`].  The two variants let callers distinguish reads
/// that were too short to extract every segment (always recoverable by counting them as
/// `TooFewBases` skips) from genuinely unexpected failures (which should abort the run).
#[derive(Debug)]
enum ResegmentError {
    /// The raw bases were too short to extract one or more segments of the read structure.
    TooShort(anyhow::Error),
    /// An unexpected error other than insufficient length.
    Other(anyhow::Error),
}

/// Hint appended to a too-short-read abort, pointing the user at the flag that converts the
/// abort into a counted skip.
const TOO_FEW_BASES_HINT: &str =
    "pass `--skip-reasons too-few-bases` to skip reads that are too short instead of aborting";

/// One unit of FASTQ records separated into their component read segments.
#[derive(PartialEq, Debug, Clone)]
struct ReadSet {
    /// Headers of the FASTQ records, one per input FASTQ, in input order.  Reads from the same
    /// template share a name but differ in their comment (notably the read-number field), so each
    /// source's header must be retained to emit it verbatim for that source's outputs.
    headers: Vec<Vec<u8>>,
    /// Segments of reads (built using the global read structures).
    segments: Vec<FastqSegment>,
    /// The raw bases / qualities from each input FASTQ, in input order.  Used to rebuild
    /// segments when a matched sample has its own per-sample read structures.  Empty for the
    /// no-per-sample-read-structures case (and for skipped reads).
    raw_per_input: Vec<RawInput>,
    /// The reason the read should be skipped for demultiplexing, or None if it should be
    /// demultiplexed.
    skip_reason: Option<SkipReason>,
}

impl ReadSet {
    const PREFIX: u8 = b'@';
    const SPACE: u8 = b' ';
    const COLON: u8 = b':';
    const PLUS: u8 = b'+';
    /// Separator written between the read name and the first UMI (Illumina DRAGEN convention).
    const UMI_DELIM: u8 = b':';
    /// Separator written between consecutive UMIs, e.g. duplex strands (Illumina DRAGEN
    /// convention).
    const UMI_SEP: u8 = b'+';

    /// Produces an iterator over references to the template segments stored in this ``ReadSet``.
    fn template_segments(&self) -> SegmentIter {
        self.segments.iter().filter(|s| s.segment_type == SegmentType::Template)
    }

    /// Produces an iterator over references to the sample barcode segments stored in this
    /// ``ReadSet``.
    fn sample_barcode_segments(&self) -> SegmentIter {
        self.segments.iter().filter(|s| s.segment_type == SegmentType::SampleBarcode)
    }

    /// Produces an iterator over references to the molecular barcode segments stored in this
    /// ``ReadSet``.
    fn molecular_barcode_segments(&self) -> SegmentIter {
        self.segments.iter().filter(|s| s.segment_type == SegmentType::MolecularBarcode)
    }

    /// Produces an iterator over references to the cellular barcode segments stored in this
    /// ``ReadSet``.
    fn cellular_barcode_segments(&self) -> SegmentIter {
        self.segments.iter().filter(|s| s.segment_type == SegmentType::CellularBarcode)
    }

    /// Generates the sample barcode sequence for this read set and returns it as a Vec of bytes.
    fn sample_barcode_sequence(&self) -> Vec<u8> {
        self.sample_barcode_segments().flat_map(|s| &s.seq).copied().collect()
    }

    /// Returns a fixed-length, per-input prefix slice of the raw bases concatenated across all
    /// input FASTQs, used to perform per-sample N-padded matching when per-sample read
    /// structures are in use.  Callers are responsible for ensuring every input has at least
    /// `prefix_lens[i]` bases — `ReadSetIterator` already gates each input on this length, so
    /// any read reaching this method is guaranteed long enough.
    ///
    /// # Panics
    ///   - In debug builds, panics if any input has fewer bases than its prefix length.
    fn matching_window_bases(&self, prefix_lens: &[usize]) -> Vec<u8> {
        let total: usize = prefix_lens.iter().sum();
        let mut out = Vec::with_capacity(total);
        for (raw, &plen) in self.raw_per_input.iter().zip(prefix_lens) {
            debug_assert!(
                raw.bases.len() >= plen,
                "raw bases length {} < matching prefix length {}; ReadSetIterator should have \
                 gated this read",
                raw.bases.len(),
                plen,
            );
            out.extend_from_slice(&raw.bases[..plen]);
        }
        out
    }

    /// Consumes this read set and returns a new one with `segments` extracted from
    /// `raw_per_input` according to the supplied per-input read structures.  `raw_per_input`
    /// is dropped (no longer needed once segments are built); `header` is moved through.
    ///
    /// # Errors
    ///   - Returns [`ResegmentError::TooShort`] if any segment falls outside the available
    ///     bases for that input (i.e. the read is shorter than the read structure requires).
    ///   - Returns [`ResegmentError::Other`] for any other extraction failure.
    /// # Panics
    ///   - Will panic if `read_structures.len() != self.raw_per_input.len()`.
    fn resegment(self, read_structures: &[ReadStructure]) -> Result<Self, ResegmentError> {
        assert_eq!(
            read_structures.len(),
            self.raw_per_input.len(),
            "resegment requires one read structure per input FASTQ ({} vs {})",
            read_structures.len(),
            self.raw_per_input.len(),
        );
        let total_segments: usize = read_structures.iter().map(|rs| rs.number_of_segments()).sum();
        let mut segments = Vec::with_capacity(total_segments);
        for (source_index, (raw, rs)) in
            self.raw_per_input.iter().zip(read_structures.iter()).enumerate()
        {
            for seg in rs.iter() {
                let (s, q) = match seg.extract_bases_and_quals(&raw.bases, &raw.quals) {
                    Ok(ok) => ok,
                    Err(e) => {
                        let is_length_error = matches!(
                            e,
                            ReadStructureError::ReadEndsBeforeSegment(_)
                                | ReadStructureError::ReadEndsAfterSegment(_)
                        );
                        let wrapped = anyhow!(
                            "Error extracting bases (len: {}) for read structure ({}) from \
                             FASTQ record {}: {}",
                            raw.bases.len(),
                            rs,
                            String::from_utf8_lossy(self.first_header()),
                            e,
                        );
                        return Err(if is_length_error {
                            ResegmentError::TooShort(wrapped)
                        } else {
                            ResegmentError::Other(wrapped)
                        });
                    }
                };
                // A variable-length template (`+T`) extracts to zero bases when the read ends
                // exactly at the template's start.  The read-structure crate returns `Ok(empty)`
                // rather than a length error, but emitting a zero-length FASTQ record is invalid,
                // so treat such a read as too short (consistent with a genuinely truncated read).
                if seg.kind == SegmentType::Template && s.is_empty() {
                    return Err(ResegmentError::TooShort(anyhow!(
                        "Read structure ({}) produced a zero-length template for FASTQ record {} \
                         (len: {})",
                        rs,
                        String::from_utf8_lossy(self.first_header()),
                        raw.bases.len(),
                    )));
                }
                segments.push(FastqSegment {
                    seq: s.to_vec(),
                    quals: q.to_vec(),
                    segment_type: seg.kind,
                    source_index,
                });
            }
        }
        Ok(Self {
            headers: self.headers,
            segments,
            raw_per_input: Vec::new(),
            skip_reason: self.skip_reason,
        })
    }

    /// Combines ``ReadSet`` structs together into a single ``ReadSet``
    fn combine_readsets(readsets: Vec<Self>) -> Self {
        assert!(!readsets.is_empty(), "Cannot call combine readsets on an empty vec!");
        let total_segments: usize = readsets.iter().map(|s| s.segments.len()).sum();

        let mut readset_iter = readsets.into_iter();
        let mut first = readset_iter.next().expect("Cannot call combine readsets on an empty vec!");
        first.segments.reserve_exact(total_segments - first.segments.len());
        first.raw_per_input.reserve_exact(readset_iter.len());

        for next_readset in readset_iter {
            first.segments.extend(next_readset.segments);
            first.raw_per_input.extend(next_readset.raw_per_input);
            first.headers.extend(next_readset.headers);
        }

        first
    }

    /// Returns the header for the given input source.  A `ReadSet` retains one header per input
    /// FASTQ (see [`combine_readsets`]), so `source_index` is always in range for the segments it
    /// holds.
    fn header_for_source(&self, source_index: usize) -> &[u8] {
        self.headers[source_index].as_slice()
    }

    /// Returns the header of the first input source, used for error messages naming the record.
    fn first_header(&self) -> &[u8] {
        self.headers[0].as_slice()
    }

    /// Reconstructs a read for a given source index by concatenating segments of the specified
    /// types from that source in order. Used for --template-types output.
    ///
    /// Note: Segments are concatenated in their original order (as they appear in the read
    /// structure), preserving the physical arrangement of bases in the original read.
    fn segments_for_source(
        &self,
        source_index: usize,
        template_types: &HashSet<SegmentType>,
    ) -> (Vec<u8>, Vec<u8>) {
        let mut seq = Vec::new();
        let mut quals = Vec::new();
        for segment in &self.segments {
            if segment.source_index == source_index
                && template_types.contains(&segment.segment_type)
            {
                seq.extend_from_slice(&segment.seq);
                quals.extend_from_slice(&segment.quals);
            }
        }
        (seq, quals)
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
    ///
    /// When `include_umi` is false, UMI (M) segments are not folded into the read name; this is
    /// used when the UMI bases are instead written into the template output bases, to avoid
    /// emitting the same bases in both the header and the bases.
    ///
    /// When `format.kind` is `Unmodified` or `NameOnly`, the header of the input source that
    /// produced this output (`source_index`) is emitted rather than a rewritten one: no sample
    /// barcode is appended to the comment and the read-number field is not rewritten.  Whether the
    /// UMI is folded into the name is controlled independently by `format.umi_in_name`.
    fn write_header<W: Write>(
        &self,
        writer: &mut W,
        read_num: usize,
        source_index: usize,
        format: HeaderFormat,
    ) -> Result<()> {
        Self::write_header_internal(
            writer,
            read_num,
            self.header_for_source(source_index),
            self.sample_barcode_segments(),
            self.molecular_barcode_segments(),
            format,
        )
    }

    fn write_header_internal<W: Write>(
        writer: &mut W,
        read_num: usize,
        header: &[u8],
        sample_barcode_segments: SegmentIter,
        mut molecular_barcode_segments: SegmentIter,
        format: HeaderFormat,
    ) -> Result<()> {
        let HeaderFormat { include_umi, kind, umi_in_name } = format;
        writer.write_all(&[Self::PREFIX])?;

        let fold_umi_into_name = include_umi && umi_in_name;

        // `Unmodified` with no UMI folded into the name emits the header verbatim and stops:
        // nothing is parsed at all here.  (`NameOnly` with no UMI folded in reaches the same
        // non-Illumina-header tolerance further below, via the `else` branch of the name-writing
        // step, since the name itself is never parsed unless a UMI is being folded into it.)
        if kind == HeaderFormatKind::Unmodified && !fold_umi_into_name {
            writer.write_all(header)?;
            return Ok(());
        }

        // Extract the name and optionally the comment
        let (name, comment) = match header.find_byte(Self::SPACE) {
            Some(x) => (&header[0..x], Some(&header[x + 1..])),
            None => (header, None),
        };

        // Handle the 'name' component of the header.  If we're not folding in a UMI we can emit
        // the name part as is.  Otherwise we need to append the UMI(s) to the name.
        let first_umi_segment =
            if fold_umi_into_name { molecular_barcode_segments.next() } else { None };
        if let Some(first_seg) = first_umi_segment {
            let sep_count = name.iter().filter(|c| **c == Self::COLON).count();
            ensure!(
                sep_count <= 7,
                "Can't handle read name with more than 8 segments: {}",
                String::from_utf8(header.to_vec())?
            );

            writer.write_all(name)?;
            if sep_count == 7 {
                // UMI already present, append to it with a UMI separator first
                writer.write_all(&[Self::UMI_SEP])?;
            } else {
                // UMI not present yet, insert a minor separator
                writer.write_all(&[Self::UMI_DELIM])?;
            }

            writer.write_all(first_seg.seq.as_slice())?;
            // Append all the UMI segments with the UMI separator in between.
            for seg in molecular_barcode_segments {
                writer.write_all(&[Self::UMI_SEP])?;
                writer.write_all(seg.seq.as_slice())?;
            }
        } else {
            writer.write_all(name)?;
        }

        // Exhaustive on `kind` (rather than early-returns eliminating down to an implicit last
        // case) so a future `HeaderFormatKind` variant fails to compile here instead of silently
        // inheriting `Illumina`'s comment-rewrite behavior.
        match kind {
            // `NameOnly` never emits a comment, regardless of what the input carried.
            HeaderFormatKind::NameOnly => Ok(()),
            // `Unmodified` (reaching here only because a UMI was folded into the name) emits the
            // rest of the original comment verbatim: the read-number field is not rewritten and
            // no sample barcode is appended.  A header with no comment stays that way.
            HeaderFormatKind::Unmodified => {
                if let Some(chars) = comment {
                    writer.write_all(&[Self::SPACE])?;
                    writer.write_all(chars)?;
                }
                Ok(())
            }
            // `Illumina`: rewrite the comment (replace the read-number field, append the sample
            // barcode(s)), synthesizing one if the input had none.
            HeaderFormatKind::Illumina => {
                writer.write_all(&[Self::SPACE])?;
                match comment {
                    None => {
                        // If no pre-existing comment, assume the read is a passing filter,
                        // non-control read and generate a comment for it (sample barcode is
                        // added below).
                        write!(writer, "{}:N:0:", read_num)?;
                    }
                    Some(chars) => {
                        // Else check it's a 4-part name... fix the read number at the front and
                        // check to see if there's a real sample barcode on the back
                        let sep_count = chars.iter().filter(|c| **c == Self::COLON).count();
                        if sep_count < 3 {
                            writer.write_all(chars)?;
                            // `chars.last()` is `None` for an empty comment (e.g. a header
                            // ending in a trailing space with nothing after it); treat that the
                            // same as "didn't already end in a colon" so a separator is still
                            // written before the sample barcode below.
                            if chars.last() != Some(&Self::COLON) {
                                writer.write_all(&[Self::COLON])?;
                            }
                        } else {
                            ensure!(
                                sep_count == 3,
                                "Comment in did not have 4 segments: {}",
                                String::from_utf8(header.to_vec())?
                            );
                            let first_colon_idx =
                                chars.iter().position(|ch| *ch == Self::COLON).unwrap();

                            // Illumina, in the unmatched FASTQs, can place a "0" in the index
                            // position, sigh
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
    /// Valid reasons for skipping reads, otherwise panic!
    skip_reasons: Vec<SkipReason>,
    /// The index of this FASTQ source (0-based), used for tracking original reads.
    source_index: usize,
    /// If true, the raw bases and qualities for each record will be retained on the resulting
    /// `ReadSet` so that records can be re-segmented per-sample later.
    keep_raw: bool,
    /// The minimum number of bases required to extract every segment of `read_structure`.
    /// Pre-computed once at construction time and used to short-circuit reads that are too
    /// short.
    min_len: usize,
}

impl Iterator for ReadSetIterator {
    type Item = ReadSet;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(rec) = self.source.next() {
            let next_fq_rec = rec.expect("Unexpected error parsing FASTQs");
            let read_name = next_fq_rec.head();
            let bases = next_fq_rec.seq();
            let quals = next_fq_rec.qual();

            if bases.len() < self.min_len {
                if self.skip_reasons.contains(&SkipReason::TooFewBases) {
                    return Some(ReadSet {
                        headers: vec![read_name.to_vec()],
                        segments: vec![],
                        raw_per_input: vec![],
                        skip_reason: Some(SkipReason::TooFewBases),
                    });
                }
                panic!(
                    "Read {} had too few bases to demux {} vs. {} needed in read structure {}.",
                    String::from_utf8(read_name.to_vec()).unwrap(),
                    bases.len(),
                    self.min_len,
                    self.read_structure
                );
            }

            // In keep_raw mode the caller will re-segment per-sample (or per-globals for
            // unmatched) downstream, so skip the up-front extraction and just retain the raw
            // bases/quals.
            if self.keep_raw {
                return Some(ReadSet {
                    headers: vec![read_name.to_vec()],
                    segments: vec![],
                    raw_per_input: vec![RawInput { bases: bases.to_vec(), quals: quals.to_vec() }],
                    skip_reason: None,
                });
            }

            let mut segments = Vec::with_capacity(self.read_structure.number_of_segments());
            for (read_segment_index, read_segment) in self.read_structure.iter().enumerate() {
                let (seq, quals) =
                    read_segment.extract_bases_and_quals(bases, quals).unwrap_or_else(|e| {
                        panic!(
                            "Error extracting bases (len: {}) or quals (len: {}) for the {}th read segment ({}) in read structure ({}) from FASTQ record with name {}; {}",
                            bases.len(),
                            quals.len(),
                            read_segment_index,
                            read_segment,
                            self.read_structure,
                            String::from_utf8(read_name.to_vec()).unwrap(),
                            e
                        )
                    });
                segments.push(FastqSegment {
                    seq: seq.to_vec(),
                    quals: quals.to_vec(),
                    segment_type: read_segment.kind,
                    source_index: self.source_index,
                });
            }
            Some(ReadSet {
                headers: vec![read_name.to_vec()],
                segments,
                raw_per_input: vec![],
                skip_reason: None,
            })
        } else {
            None
        }
    }
}

impl ReadSetIterator {
    /// Instantiates a new iterator over the read sets for a set of FASTQs with defined read
    /// structures.  When `keep_raw` is true, segment extraction is deferred — the resulting
    /// `ReadSet`s carry only the raw bases/quals per input so that callers can re-segment
    /// per-sample (or per-globals) downstream.  `min_len` is the minimum number of bases a
    /// record must have for this iterator to yield it (or skip it as `TooFewBases`).
    pub fn new(
        read_structure: ReadStructure,
        source: FastqReader<Box<dyn BufRead + Send>>,
        skip_reasons: Vec<SkipReason>,
        source_index: usize,
        keep_raw: bool,
        min_len: usize,
    ) -> Self {
        Self { read_structure, source, skip_reasons, source_index, keep_raw, min_len }
    }
}

////////////////////////////////////////////////////////////////////////////////
// TemplateOutput
////////////////////////////////////////////////////////////////////////////////

/// How segments are routed into the template output bases and the read header, derived once from
/// `--template-types` and reused for every record.
struct TemplateOutput {
    /// Segment types to concatenate into the template read bases. Always includes Template (T).
    types: HashSet<SegmentType>,
    /// Whether any non-template segment types are folded into the template bases. When false, the
    /// template writer emits just the template segment (the default behavior).
    include_other_segments: bool,
    /// Whether UMI (M) segments should be written into the read header. False when UMIs are folded
    /// into the template bases, so the same bases are not written in both the header and the bases.
    umi_in_header: bool,
    /// The header "shape" for each output record.
    header_format: HeaderFormatKind,
    /// Whether to fold UMI(s) into the read name (already resolved from an explicit
    /// `--umi-in-name` or the format's default).
    umi_in_name: bool,
}

impl TemplateOutput {
    /// Builds the template output routing. `types` (the `--template-types` segment types) drives
    /// which segments are folded into the template bases and whether UMIs stay in the header;
    /// `header_format` and `umi_in_name` control how the read header is built.
    fn new(
        types: HashSet<SegmentType>,
        header_format: HeaderFormatKind,
        umi_in_name: bool,
    ) -> Self {
        let include_other_segments = types.iter().any(|t| *t != SegmentType::Template);
        let umi_in_header = !types.contains(&SegmentType::MolecularBarcode);
        Self { types, include_other_segments, umi_in_header, header_format, umi_in_name }
    }

    /// The per-record header-formatting parameters implied by this routing.
    fn header_format(&self) -> HeaderFormat {
        HeaderFormat {
            include_umi: self.umi_in_header,
            kind: self.header_format,
            umi_in_name: self.umi_in_name,
        }
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
    /// Vec of the writers for the cellular barcode read segments.
    cellular_barcode_writers: Option<Vec<W>>,
}

impl<W: Write> SampleWriters<W> {
    /// Destroys this struct and decomposes it into its component types. Used when swapping
    /// writers for pooled writers.
    #[allow(clippy::type_complexity)]
    fn into_parts(
        self,
    ) -> (String, Option<Vec<W>>, Option<Vec<W>>, Option<Vec<W>>, Option<Vec<W>>) {
        (
            self.name,
            self.template_writers,
            self.sample_barcode_writers,
            self.molecular_barcode_writers,
            self.cellular_barcode_writers,
        )
    }

    /// Writes a set of reads (defined as a ``ReadSet``) to the appropriate writers on this
    /// ``Self`` struct.
    /// Reads in the read set should be 1:1 with writers in the writer set however this is not
    /// checked at runtime as doing so substantially slows demulitplexing.
    ///
    /// The `template` parameter specifies how segments are routed into the template output. If it
    /// includes only Template, only template segments are written. If it includes additional types
    /// (e.g. MolecularBarcode), those segments are concatenated with the template, and UMIs are
    /// omitted from the read header so the same bases are not written twice.
    fn write(&mut self, read_set: &ReadSet, template: &TemplateOutput) -> Result<()> {
        // Handle template writers separately to support template output routing
        if let Some(writers) = &mut self.template_writers {
            for (read_idx, (writer, segment)) in
                writers.iter_mut().zip(read_set.template_segments()).enumerate()
            {
                read_set.write_header(
                    writer,
                    read_idx + 1,
                    segment.source_index,
                    template.header_format(),
                )?;
                writer.write_all(b"\n")?;
                if template.include_other_segments {
                    // Write segments of the requested types from this source
                    let (seq, quals) =
                        read_set.segments_for_source(segment.source_index, &template.types);
                    writer.write_all(&seq)?;
                    writer.write_all(b"\n+\n")?;
                    writer.write_all(&quals)?;
                } else {
                    // Write just the template segment
                    writer.write_all(segment.seq.as_slice())?;
                    writer.write_all(b"\n+\n")?;
                    writer.write_all(segment.quals.as_slice())?;
                }
                writer.write_all(b"\n")?;
            }
        }

        // Handle other segment types (barcodes, UMIs, etc.) - always write segments, not originals
        for (writers_opt, segments) in [
            (&mut self.sample_barcode_writers, &mut read_set.sample_barcode_segments()),
            (&mut self.molecular_barcode_writers, &mut read_set.molecular_barcode_segments()),
            (&mut self.cellular_barcode_writers, &mut read_set.cellular_barcode_segments()),
        ] {
            if let Some(writers) = writers_opt {
                for (read_idx, (writer, segment)) in writers.iter_mut().zip(segments).enumerate() {
                    read_set.write_header(
                        writer,
                        read_idx + 1,
                        segment.source_index,
                        template.header_format(),
                    )?;
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
        for writers in [
            self.template_writers,
            self.sample_barcode_writers,
            self.molecular_barcode_writers,
            self.cellular_barcode_writers,
        ]
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
/// Five kinds of operators are recognized:
///
/// 1. `T` identifies a template read
/// 2. `B` identifies a sample barcode read
/// 3. `M` identifies a unique molecular index read
/// 4. `C` identifies a unique cellular barcode read
/// 5. `S` identifies a set of bases that should be skipped or ignored
///
/// The last `<number><operator>` pair may be specified using a `+` sign instead of number to
/// denote "all remaining bases". This is useful if, e.g., fastqs have been trimmed and contain
/// reads of varying length. Both reads must have template bases.
///
/// Metadata about the samples should be given as a headered metadata TSV file with at least the
/// following two columns present:
///
/// 1. `sample_id` - the id of the sample or library.
/// 2. `barcode` - the expected barcode sequence associated with the `sample_id`.
///
/// For reads containing multiple barcodes (such as dual-indexed reads), all barcodes should be
/// concatenated together in the order they are read and stored in the `barcode` field.
///
/// IUPAC bases are supported in the (expected) `barcode` column.  An observed IUPAC base must be
/// at least as specific as the corresponding base in the expected sample barcode.  E.g. If the
/// observed base is an N, it will only match expected sample barcrods with an N.  And if the
/// observed base is an R, it will match R, V, D, and N, since the latter IUPAC codes allow both
/// A and G (R/V/D/N are a superset of the bases compare to R).
///
/// The read structures will be used to extract the observed sample barcode, template bases,
/// molecular identifiers, and cellular barcodes from each read.  The observed sample barcode will
/// be matched to the sample barcodes extracted from the bases in the sample metadata and associated
/// read structures.
///
/// An observed barcode matches an expected barcode if all the following are true:
/// 1. The number of mismatches (edits/substitutions) is less than or equal to the maximum
///    mismatches (see `--max-mismatches`).
/// 2. The difference between number of mismatches in the best and second best barcodes is greater
///    than or equal to the minimum mismatch delta (`--min-mismatch-delta`).
///
/// The expected barcode sequence may contains Ns, which are not counted as mismatches regardless
/// of the observed base (e.g. the expected barcode `AAN` will have zero mismatches relative to
/// both the observed barcodes `AAA` and `AAN`).
///
/// ## Per-sample read structures
///
/// In addition to the global `--read-structures`, the metadata TSV may include optional columns
/// `read_structure_1`, `read_structure_2`, ..., `read_structure_<N>` (where `N` is the number of
/// input FASTQs).  When any cell is non-empty for a sample, the per-sample read structures
/// replace the global `--read-structures` for that sample — both for matching and for output
/// extraction.
///
/// Per-sample structures support per-cell fall-back to the global `--read-structures`:
///
/// - A blank `read_structure_<i>` cell uses `--read-structures[i-1]` for that sample's input *i*.
/// - A row whose `read_structure_<n>` cells are all blank uses `--read-structures` entirely for
///   that sample (equivalent to omitting the columns for that sample only).
/// - The unmatched pseudo-sample always uses the global `--read-structures`.
///
/// Constraints:
///
/// 1. The number of `read_structure_<n>` columns must equal `--read-structures.len()`.
/// 2. The concatenated `B`-segment length must equal the `barcode` column length for every
///    sample (computed from each sample's *effective* read structures, i.e. with fall-backs
///    applied).
///
/// Different samples may have different per-input `(T, B, M, C)` segment counts (and hence
/// produce different sets of output files).  This supports protocols with sample-dependent
/// read structures (e.g. CODEC, where each sample may include a stagger spacer of varying
/// length so that the position of the constant ligation base shifts per sample).
///
/// During matching, each sample's expected pattern is constructed from its effective read
/// structure by filling `B`-segment positions from the `barcode` column and treating
/// `M`/`S`/`C` segment positions as `N` wildcards.  All patterns are padded with trailing `N`s
/// to a per-input matching window equal to the longest pre-template prefix across samples.
///
/// Example metadata TSV (CODEC stagger):
///
/// ```text
/// sample_id  barcode         read_structure_1   read_structure_2
/// S1         GATTACAGATTACA  3M7B1S+T           3M7B1S+T
/// S2         TTTTTTTTTTTTTT  3M1S7B1S+T         3M1S7B1S+T
/// ```
///
/// ## Outputs
///
/// All outputs are generated in the provided `--output` directory.  For each sample plus the
/// unmatched reads, FASTQ files are written for each read segment (specified in the read
/// structures) of one of the types supplied to `--output-types`.  FASTQ files have names
/// of the format:
///
/// ```bash
/// {sample_id}.{segment_type}{read_num}.fq.gz
/// ```
///
/// where `segment_type` is one of `R`, `I`, `U`, and `C` (for template, sample barcode/index,
/// molecular barcode/UMI, and cellular barcode reads, respectively) and `read_num` is a number
/// starting at 1 for each segment type.
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
/// ```bash
/// fqtk demux \
///     --inputs r1.fq.gz i1.fq.gz i2.fq.gz r2.fq.gz \
///     --read-structures 8B92T 8B 8B 100T \
///     --sample-metadata metadata.tsv \
///     --output output_folder
/// ```
///
#[derive(Parser, Debug)]
#[command(version)]
#[clap(verbatim_doc_comment)]
pub(crate) struct Demux {
    /// One or more input FASTQ files each corresponding to a sequencing read (e.g. R1, I1).
    #[clap(long, short = 'i', required = true, num_args = 1..)]
    inputs: Vec<PathBuf>,

    /// The read structures, one per input FASTQ in the same order.
    ///
    /// Per-sample read structures (see the `read_structure_<n>` metadata columns) take
    /// precedence for each matched sample, and a blank cell falls back to the corresponding
    /// `--read-structures` entry.  The unmatched pseudo-sample always uses
    /// `--read-structures` for its output extraction.  The number of `read_structure_<n>`
    /// columns must equal `--read-structures.len()`.
    #[clap(long, short = 'r' , required = true, num_args = 1..)]
    read_structures: Vec<ReadStructure>,

    /// The read structure types to write to their own files (Must be one of T, B, M, or C for
    /// template reads, sample barcode reads, molecular barcode reads, or cellular barcode reads).
    ///
    /// Multiple output types may be specified as a space-delimited list.
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
    #[clap(long, short = 't', default_value = "8")]
    threads: usize,

    /// The level of compression to use to compress outputs.
    #[clap(long, short = 'c', default_value = "5")]
    compression_level: usize,

    /// Skip demultiplexing reads for any of the following reasons, otherwise panic.
    ///
    /// 1. `too-few-bases`: there are too few bases or qualities to extract given the read
    ///    structures.  For example, if a read is 8bp long but the read structure is `10B`, or
    ///    if a read is empty and the read structure is `+T`.
    #[clap(long, short = 'S')]
    skip_reasons: Vec<SkipReason>,

    /// The read structure types to include in the template FASTQ output files.
    ///
    /// By default, only template (T) segments are included. To include additional segment types
    /// (e.g. to preserve UMIs in the output read bases), specify them here. For example,
    /// `--template-types M T` will concatenate the molecular barcode and template segments.
    ///
    /// To output the full original reads (all segments), specify all segment types present in
    /// your read structure (e.g. `--template-types B M T`).
    ///
    /// Segments are only merged *within the same physical read*: a non-`T` segment is folded into
    /// the template bases only when it is co-located with a `T` in the same read structure (e.g.
    /// `8M84T`). A segment on a separate read (e.g. a UMI on its own index read) is never merged
    /// into a template on another read; route it via `--output-types` instead, or leave it in the
    /// read header. When a UMI (M) is included here it is written into the template bases and is
    /// therefore omitted from the read header (it is not written in both places).
    ///
    /// Note: If `--template-types` includes any non-`T` type, `T` must be included in
    /// `--output-types`; each requested non-`T` type must be co-located with a `T` in every read
    /// structure where it appears; and each read structure must contain at most one `T` segment.
    #[clap(long, default_value = "T", num_args = 1..)]
    template_types: Vec<char>,

    /// The FASTQ header "shape" for each output record: whether the comment is rewritten
    /// (`illumina`), kept verbatim (`unmodified`), or dropped entirely (`name-only`).  See each
    /// value's own description below for exact per-format behavior and examples.
    ///
    /// With `unmodified` and `name-only`, sample barcode bases (and, unless `--umi-in-name true`
    /// is given, UMI bases) are retained only if routed to their own FASTQs via `--output-types`
    /// or folded into the template bases via `--template-types`; otherwise they are not present in
    /// any output.  A warning is emitted when this would silently discard bases.
    #[clap(long, value_enum, default_value_t = HeaderFormatKind::Illumina)]
    header_format: HeaderFormatKind,

    /// Whether to fold UMI(s) into the read name.
    ///
    /// If not given, this is decided by `--header-format`: `true` for `illumina`, `false` for
    /// `unmodified` and `name-only`.  An explicit value here overrides the format's default --
    /// e.g. `--header-format unmodified --umi-in-name true` keeps the original comment but still
    /// appends the UMI(s) to the name.  Has no effect when the UMI is instead folded into the
    /// template bases via `--template-types M ...`, since there is then nothing left to fold into
    /// the name.
    #[clap(long)]
    umi_in_name: Option<bool>,
}

impl Demux {
    /// The effective "fold UMI into read name" setting: an explicit `--umi-in-name` wins over the
    /// `--header-format` default, unless `--template-types` has already routed the UMI into the
    /// template bases, in which case there is nothing left to fold into the name regardless of
    /// this setting.
    fn effective_umi_in_name(&self) -> bool {
        self.umi_in_name.unwrap_or(self.header_format.default_umi_in_name())
    }

    /// Warns when the chosen `--header-format`/`--umi-in-name` combination drops segment bases
    /// that are not written to their own FASTQs (`--output-types`) or folded into the template
    /// bases (`--template-types`).
    ///
    /// The default `illumina` format carries UMIs in the read name and sample barcodes in the
    /// comment, so nothing is lost.  `unmodified` and `name-only` do not, and without this warning
    /// a run that silently discards every UMI base would look successful.
    fn warn_on_discarded_segments(
        &self,
        output_types: &HashSet<SegmentType>,
        template_types: &HashSet<SegmentType>,
    ) {
        for segment_type in Self::discarded_segment_types(
            &self.read_structures,
            self.header_format,
            self.effective_umi_in_name(),
            output_types,
            template_types,
        ) {
            match segment_type {
                // Attribute the loss correctly: an explicit `--umi-in-name false` overrides
                // `illumina` (whose own default is to write the UMI into the name), so blaming
                // `--header-format` in that case would misdescribe the cause.  `unmodified`/
                // `name-only` already default to `false`, so an explicit override there isn't the
                // cause -- fall through to the format-attributed message below, which correctly
                // tells the user to pass `--umi-in-name true`.
                SegmentType::MolecularBarcode
                    if self.umi_in_name == Some(false)
                        && self.header_format.default_umi_in_name() =>
                {
                    warn!(
                        "--umi-in-name false does not write UMIs to the read name, and M \
                         segments are not in --output-types or --template-types: the UMI bases \
                         will not appear in any output.  Add M to --output-types (or \
                         --template-types), or drop --umi-in-name false."
                    );
                }
                SegmentType::MolecularBarcode => warn!(
                    "--header-format {} does not write UMIs to the read name, and M segments are \
                     not in --output-types or --template-types: the UMI bases will not appear in \
                     any output.  Add M to --output-types (or --template-types), or pass \
                     --umi-in-name true.",
                    self.header_format
                ),
                SegmentType::SampleBarcode => warn!(
                    "--header-format {} does not write sample barcodes to the read comment, and B \
                     segments are not in --output-types or --template-types: the sample barcode \
                     bases will not appear in any output.  Add B to --output-types (or \
                     --template-types) if you need them.",
                    self.header_format
                ),
                _ => {}
            }
        }
    }

    /// Returns the segment types whose bases would be silently discarded under `header_format`/
    /// `umi_in_name` given the requested output/template routing: those present in the read
    /// structures but neither written to the read header nor routed to their own FASTQs
    /// (`output_types`) or the template bases (`template_types`).
    ///
    /// The UMI and sample-barcode axes are independent: a UMI is discarded only when it's neither
    /// folded into the name (`umi_in_name`) nor retained elsewhere; the sample barcode is
    /// discarded only when `header_format` doesn't rewrite the comment (only `illumina` does) and
    /// it isn't retained elsewhere.  `illumina` always writes both, so it never discards anything.
    fn discarded_segment_types(
        read_structures: &[ReadStructure],
        header_format: HeaderFormatKind,
        umi_in_name: bool,
        output_types: &HashSet<SegmentType>,
        template_types: &HashSet<SegmentType>,
    ) -> Vec<SegmentType> {
        // A segment type survives if it gets its own FASTQ or is folded into the template bases.
        let is_retained = |segment_type: SegmentType| {
            output_types.contains(&segment_type) || template_types.contains(&segment_type)
        };
        let present = |segment_type: SegmentType| {
            read_structures.iter().any(|rs| rs.segments_by_type(segment_type).count() > 0)
        };

        let mut discarded = Vec::new();
        if !umi_in_name
            && present(SegmentType::MolecularBarcode)
            && !is_retained(SegmentType::MolecularBarcode)
        {
            discarded.push(SegmentType::MolecularBarcode);
        }
        if !header_format.rewrites_comment()
            && present(SegmentType::SampleBarcode)
            && !is_retained(SegmentType::SampleBarcode)
        {
            discarded.push(SegmentType::SampleBarcode);
        }
        discarded
    }

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
        let mut cellular_barcode_writers = None;

        for output_type in output_types {
            let mut output_type_writers = Vec::new();

            let file_type_code = match output_type {
                SegmentType::Template => 'R',
                SegmentType::SampleBarcode => 'I',
                SegmentType::MolecularBarcode => 'U',
                SegmentType::CellularBarcode => 'C',
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
                SegmentType::CellularBarcode => {
                    cellular_barcode_writers = Some(output_type_writers);
                }
                _ => {}
            }
        }

        Ok(SampleWriters {
            name: prefix.to_string(),
            template_writers,
            sample_barcode_writers,
            molecular_barcode_writers,
            cellular_barcode_writers,
        })
    }

    /// Creates one writer per sample per read segment in the provided read structures for
    /// requested output type.  Each sample's writers are created from its per-sample read
    /// structures when present; otherwise the global `read_structures` are used.  The
    /// unmatched pseudo-sample always uses the global `read_structures`.
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
            let rs = sample.read_structures.as_deref().unwrap_or(read_structures);
            samples_writers.push(Self::create_sample_writers(
                rs,
                &sample.sample_id,
                output_types,
                output_dir,
            )?);
        }
        // Add the unmatched 'sample' using the global read structures.
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
            let (name, template_writers, barcode_writers, mol_writers, cb_writers) =
                sample.into_parts();
            let mut new_template_writers = None;
            let mut new_sample_barcode_writers = None;
            let mut new_molecular_barcode_writers = None;
            let mut new_cellular_barcode_writers = None;

            for (optional_ws, target) in [
                (template_writers, &mut new_template_writers),
                (barcode_writers, &mut new_sample_barcode_writers),
                (mol_writers, &mut new_molecular_barcode_writers),
                (cb_writers, &mut new_cellular_barcode_writers),
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
                cellular_barcode_writers: new_cellular_barcode_writers,
            });
        }
        let pool = pool_builder.build()?;
        Ok((pool, new_sample_writers))
    }

    /// Parses a list of segment type characters into a `HashSet<SegmentType>`, returning
    /// `Some(set)` on success or `None` after pushing an error message to `errors` if any
    /// character is invalid.
    fn parse_segment_types(
        chars: &[char],
        error_prefix: &str,
        errors: &mut Vec<String>,
    ) -> Option<HashSet<SegmentType>> {
        match chars.iter().map(|&c| SegmentType::try_from(c)).collect::<Result<HashSet<_>, _>>() {
            Ok(set) => Some(set),
            Err(e) => {
                errors.push(format!("{error_prefix}: {e}"));
                None
            }
        }
    }

    /// Checks that inputs to demux are valid and returns open file handles for the inputs.
    /// Checks:
    ///     - That the number of input files and number of read structs provided are the same
    ///     - That the output directory is not read-only
    ///     - That the input files exist
    ///     - That the input files have read permissions.
    ///     - That `--output-types` and `--template-types` parse to valid segment types
    ///     - That `--template-types` includes T (the template segment)
    ///     - That `--output-types` includes T when `--template-types` has non-T segments
    ///     - That no read structure has multiple T segments when `--template-types` has non-T
    fn validate_and_prepare_inputs(
        &self,
    ) -> Result<(VecOfReaders, HashSet<SegmentType>, HashSet<SegmentType>)> {
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

        let output_segment_types_result = Self::parse_segment_types(
            &self.output_types,
            "Error parsing segment types to report",
            &mut constraint_errors,
        );
        let template_segment_types_result = Self::parse_segment_types(
            &self.template_types,
            "Error parsing template segment types",
            &mut constraint_errors,
        );

        // Warn when the chosen header format would drop barcode bases that are not routed
        // anywhere else, so the loss is never silent.
        if let (Some(output_types), Some(template_types)) =
            (&output_segment_types_result, &template_segment_types_result)
        {
            self.warn_on_discarded_segments(output_types, template_types);
        }

        // Validate template_types constraints
        if let (Some(output_types), Some(template_types)) =
            (&output_segment_types_result, &template_segment_types_result)
        {
            // template_types must always contain T (the template segment)
            if !template_types.contains(&SegmentType::Template) {
                constraint_errors
                    .push("--template-types must include T (template segment)".to_owned());
            }

            // If template_types includes non-T segments, output_types must include T
            // (otherwise there's no template file to write the combined segments to)
            let has_non_template = template_types.iter().any(|t| *t != SegmentType::Template);
            if has_non_template && !output_types.contains(&SegmentType::Template) {
                constraint_errors.push(
                    "--template-types with non-T segments requires T in --output-types".to_owned(),
                );
            }

            // If template_types includes non-T segments, reject read structures with multiple
            // template segments in a single source, as segments_for_source would produce
            // duplicated output for each template writer.
            if has_non_template {
                for (idx, rs) in self.read_structures.iter().enumerate() {
                    let template_count = rs.segments_by_type(SegmentType::Template).count();
                    if template_count > 1 {
                        constraint_errors.push(format!(
                            "--template-types with non-T segments is not supported when a read structure has multiple template segments (read structure #{} has {})",
                            idx + 1, template_count
                        ));
                    }
                }

                // A non-T type is folded into the template bases only when it is co-located with
                // a T in the *same* physical read; segments are never merged across reads. Require
                // each requested non-T type to (a) appear in at least one read structure, and (b)
                // share every read structure it appears in with a T. Otherwise its bases would be
                // neither merged into a template nor written anywhere, i.e. silently dropped.
                for segment_type in template_types.iter().filter(|t| **t != SegmentType::Template) {
                    let type_char = segment_type.value();
                    let present_anywhere = self
                        .read_structures
                        .iter()
                        .any(|rs| rs.segments_by_type(*segment_type).count() > 0);
                    if !present_anywhere {
                        constraint_errors.push(format!(
                            "--template-types includes {type_char} but no read structure contains that segment type"
                        ));
                        continue;
                    }
                    for (idx, rs) in self.read_structures.iter().enumerate() {
                        let has_type = rs.segments_by_type(*segment_type).count() > 0;
                        let has_template = rs.segments_by_type(SegmentType::Template).count() > 0;
                        if has_type && !has_template {
                            constraint_errors.push(format!(
                                "--template-types includes {type_char} but read structure #{} contains {type_char} with no template (T) segment to merge it into (segments are only merged within the same read); add {type_char} to --output-types, or remove it from --template-types so it is written to the read header instead",
                                idx + 1
                            ));
                        }
                    }
                }
            }

            // Every segment type present in the read structures must have a destination so its
            // bases are not silently lost. Template (T) is forced into output, sample barcode (B)
            // and molecular barcode (M) fall back to the read header, and Skip (S) is discarded by
            // design. Cellular barcode (C) has no read-header field, so it must be routed
            // explicitly to either an output file or the template bases.
            let cellular_present = self
                .read_structures
                .iter()
                .any(|rs| rs.segments_by_type(SegmentType::CellularBarcode).count() > 0);
            if cellular_present
                && !output_types.contains(&SegmentType::CellularBarcode)
                && !template_types.contains(&SegmentType::CellularBarcode)
            {
                constraint_errors.push(
                    "Cellular barcode (C) segments are present but assigned to neither --output-types nor --template-types, and have no read-header field; add C to one of them".to_owned(),
                );
            }
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
            // Safe: parse_segment_types only returns None after pushing to constraint_errors,
            // so if constraint_errors is empty both must be Some.
            let output_segment_types =
                output_segment_types_result.expect("output segment types parsed without errors");
            let template_segment_types = template_segment_types_result
                .expect("template segment types parsed without errors");
            if output_segment_types.is_empty() {
                constraint_errors.push(
                    "No output types requested, must request at least one output segment type."
                        .to_owned(),
                );
            } else {
                return Ok((fq_readers_result?, output_segment_types, template_segment_types));
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
    #[allow(clippy::too_many_lines)]
    /// Executes the demux command
    fn execute(&self) -> Result<()> {
        let (fq_readers, output_segment_types, template_segment_types) =
            self.validate_and_prepare_inputs()?;
        let template_output = TemplateOutput::new(
            template_segment_types,
            self.header_format,
            self.effective_umi_in_name(),
        );

        let sample_group = SampleGroup::from_file(&self.sample_metadata, &self.read_structures)?;
        info!(
            "{} samples loaded from file {:?}",
            sample_group.samples.len(),
            &self.sample_metadata
        );
        let fq_sources =
            fq_readers.into_iter().map(|fq| FastqReader::with_capacity(fq, BUFFER_SIZE));

        // reserve 1 for main thread and 2 for read ahead, remaining on writing.
        let reader_threads = if self.threads <= 6 { 1 } else { 2 };
        let main_thread = 1;
        let writer_threads = self.threads - main_thread - reader_threads;

        let (mut pool, mut sample_writers) = Self::build_writer_pool(
            Self::create_writers(
                &self.read_structures,
                &sample_group,
                &output_segment_types,
                &self.unmatched_prefix,
                &self.output,
            )?,
            self.compression_level,
            writer_threads,
        )?;
        let unmatched_index = sample_writers.len() - 1;
        info!("Created sample and {} writers.", self.unmatched_prefix);

        let mut sample_metrics = sample_group
            .samples
            .iter()
            .map(|s| DemuxMetric::new(s.sample_id.as_str(), s.barcode.as_str()))
            .collect_vec();
        let mut unmatched_metric = DemuxMetric::new(self.unmatched_prefix.as_str(), ".");

        let per_sample_rs = sample_group.has_per_sample_read_structures();
        let prefix_lens = if per_sample_rs {
            sample_group.matching_prefix_lens(&self.read_structures)?
        } else {
            Vec::new()
        };
        let matching_patterns = if per_sample_rs {
            sample_group.build_matching_patterns(&self.read_structures, &prefix_lens)?
        } else {
            Vec::new()
        };

        let mut barcode_matcher = if per_sample_rs {
            BarcodeMatcher::with_patterns(
                &sample_group.samples,
                matching_patterns,
                u8::try_from(self.max_mismatches)?,
                u8::try_from(self.min_mismatch_delta)?,
                true,
            )
        } else {
            BarcodeMatcher::new(
                &sample_group.samples,
                u8::try_from(self.max_mismatches)?,
                u8::try_from(self.min_mismatch_delta)?,
                true,
            )
        };

        // In per-sample mode the iterator gates on the per-input matching-window length
        // (prefix_lens[i]); the global read structure's full min_len would over-reject reads
        // that satisfy a sample's shorter structure.  Re-segmentation downstream re-checks
        // length against the matched sample's structure.
        let mut fq_iterators = fq_sources
            .zip(self.read_structures.clone())
            .enumerate()
            .map(|(source_index, (source, read_structure))| {
                let min_len = if per_sample_rs {
                    prefix_lens[source_index]
                } else {
                    read_structure.iter().map(|s| s.length().unwrap_or(1)).sum()
                };
                ReadSetIterator::new(
                    read_structure,
                    source,
                    self.skip_reasons.clone(),
                    source_index,
                    per_sample_rs,
                    min_len,
                )
                .read_ahead(1000, 1000)
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
        let allow_too_few = self.skip_reasons.contains(&SkipReason::TooFewBases);
        let mut skip_reasons: HashMap<SkipReason, usize> = HashMap::new();
        loop {
            let mut next_read_sets = Vec::with_capacity(fq_iterators.len());
            for iter in &mut fq_iterators {
                if let Some(rec) = iter.next() {
                    next_read_sets.push(rec);
                }
            }
            if let Some(reason) = next_read_sets.iter().find_map(|read_set| read_set.skip_reason) {
                let previous = skip_reasons.get(&reason).unwrap_or(&0);
                skip_reasons.insert(reason, 1 + previous);
                continue;
            } else if next_read_sets.is_empty() {
                break;
            }
            assert_eq!(
                next_read_sets.len(),
                fq_iterators.len(),
                "FASTQ sources out of sync at records: {:?}",
                next_read_sets
            );
            let read_set = ReadSet::combine_readsets(next_read_sets);
            let observed: Vec<u8> = if per_sample_rs {
                read_set.matching_window_bases(&prefix_lens)
            } else {
                read_set.sample_barcode_sequence()
            };
            if let Some(barcode_match) = barcode_matcher.assign(&observed) {
                let matched_idx = barcode_match.best_match;
                if per_sample_rs {
                    let rs_for_sample = sample_group.samples[matched_idx]
                        .read_structures
                        .as_deref()
                        .unwrap_or(&self.read_structures);
                    match read_set.resegment(rs_for_sample) {
                        Ok(resegmented) => {
                            sample_writers[matched_idx].write(&resegmented, &template_output)?;
                            sample_metrics[matched_idx].templates += 1;
                        }
                        Err(ResegmentError::TooShort(_)) if allow_too_few => {
                            *skip_reasons.entry(SkipReason::TooFewBases).or_insert(0) += 1;
                        }
                        Err(ResegmentError::TooShort(e)) => {
                            return Err(e.context(TOO_FEW_BASES_HINT));
                        }
                        Err(ResegmentError::Other(e)) => return Err(e),
                    }
                } else {
                    sample_writers[matched_idx].write(&read_set, &template_output)?;
                    sample_metrics[matched_idx].templates += 1;
                }
            } else if per_sample_rs {
                // Iterator gating in per-sample mode is loosened to the matching prefix length,
                // so a read can clear the prefix gate while being too short for the global read
                // structure used on the unmatched path.  Honor `--skip-reasons too-few-bases`
                // exactly as the matched and legacy paths do: skip when requested, otherwise
                // abort with a remedy hint.
                match read_set.resegment(&self.read_structures) {
                    Ok(resegmented) => {
                        sample_writers[unmatched_index].write(&resegmented, &template_output)?;
                        unmatched_metric.templates += 1;
                    }
                    Err(ResegmentError::TooShort(_)) if allow_too_few => {
                        *skip_reasons.entry(SkipReason::TooFewBases).or_insert(0) += 1;
                    }
                    Err(ResegmentError::TooShort(e)) => return Err(e.context(TOO_FEW_BASES_HINT)),
                    Err(ResegmentError::Other(e)) => return Err(e),
                }
            } else {
                sample_writers[unmatched_index].write(&read_set, &template_output)?;
                unmatched_metric.templates += 1;
            }
            logger.record();
        }

        // Shut down the pool
        info!("Finished reading input FASTQs.");
        sample_writers.into_iter().map(SampleWriters::close).collect::<Result<Vec<_>>>()?;
        pool.stop_pool()?;
        info!("Output FASTQ writing complete.");

        // Write out reasons for skipping records
        if skip_reasons.is_empty() {
            info!("No records were skipped.");
        } else {
            for (reason, count) in skip_reasons.iter().sorted_by_key(|(_, c)| *c) {
                info!("{} records were skipped due to {}", count, reason);
            }
        }

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

    /// Writes a single-record FASTQ whose header line is exactly `header` (without the leading
    /// `@`) and whose bases are `bases`.  Unlike [`fastq_file`], this lets a test control the full
    /// header, including a comment, so per-source header handling can be exercised.
    fn fastq_file_with_header(
        tmpdir: &TempDir,
        filename_prefix: &str,
        header: &str,
        bases: &str,
    ) -> PathBuf {
        let io = Io::default();
        let path = tmpdir.path().join(format!("{filename_prefix}.fastq"));
        let lines =
            vec![format!("@{header}"), bases.to_owned(), "+".to_owned(), ";".repeat(bases.len())];
        io.write_lines(&path, lines).unwrap();
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

    fn read_set(segments: Vec<FastqSegment>) -> ReadSet {
        read_set_with_header(b"NOT_IMPORTANT", segments)
    }

    fn read_set_with_header(header: &[u8], segments: Vec<FastqSegment>) -> ReadSet {
        ReadSet {
            headers: vec![header.to_vec()],
            segments,
            raw_per_input: vec![],
            skip_reason: None,
        }
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
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
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
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();
    }

    #[test]
    #[should_panic(expected = "cannot be read-only")]
    #[allow(clippy::permissions_set_readonly_false)]
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
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
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
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
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
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
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
            threads: 5,
            compression_level: 5,
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
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

    #[test]
    fn test_template_types_includes_barcode() {
        let tmp = TempDir::new().unwrap();
        let read_structures = vec![ReadStructure::from_str("17B100T").unwrap()];
        let s1_barcode = "AAAAAAAAGATTACAGA";
        let template_seq = "T".repeat(100);
        let full_read = s1_barcode.to_owned() + &template_seq;
        let sample_metadata = metadata_file(&tmp, &[s1_barcode]);
        let input_files = vec![fastq_file(&tmp, "ex", "ex", &[&full_read])];

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
            skip_reasons: vec![],
            template_types: vec!['B', 'T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();

        let output_path = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads = read_fastq(&output_path);

        assert_eq!(fq_reads.len(), 1);
        // With template_types B T, the output should be the full original read (barcode + template)
        assert_equal(
            &fq_reads[0],
            &OwnedRecord {
                head: b"ex_0 1:N:0:AAAAAAAAGATTACAGA".to_vec(),
                seq: full_read.as_bytes().to_vec(),
                qual: ";".repeat(117).as_bytes().to_vec(),
            },
        );
    }

    #[test]
    fn test_template_types_with_paired_reads() {
        let tmp = TempDir::new().unwrap();
        // R1 has inline barcode, R2 is pure template
        let read_structures = vec![
            ReadStructure::from_str("8B92T").unwrap(),
            ReadStructure::from_str("100T").unwrap(),
        ];
        let s1_barcode = "AAAAAAAA";
        let r1_full = s1_barcode.to_owned() + &"A".repeat(92);
        let r2_full = "T".repeat(100);
        let sample_metadata = metadata_file(&tmp, &[s1_barcode]);
        let input_files = vec![
            fastq_file(&tmp, "ex_R1", "ex", &[&r1_full]),
            fastq_file(&tmp, "ex_R2", "ex", &[&r2_full]),
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
            skip_reasons: vec![],
            template_types: vec!['B', 'T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();

        // R1 output should be full 100bp read (barcode + template)
        let r1_path = output_dir.join("Sample0000.R1.fq.gz");
        let r1_reads = read_fastq(&r1_path);
        assert_eq!(r1_reads.len(), 1);
        assert_equal(
            &r1_reads[0],
            &OwnedRecord {
                head: b"ex_0 1:N:0:AAAAAAAA".to_vec(),
                seq: r1_full.as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            },
        );

        // R2 output should be full 100bp read (same as template since no barcode)
        let r2_path = output_dir.join("Sample0000.R2.fq.gz");
        let r2_reads = read_fastq(&r2_path);
        assert_eq!(r2_reads.len(), 1);
        assert_equal(
            &r2_reads[0],
            &OwnedRecord {
                head: b"ex_0 2:N:0:AAAAAAAA".to_vec(),
                seq: r2_full.as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            },
        );
    }

    #[test]
    fn test_template_types_with_umi() {
        let tmp = TempDir::new().unwrap();
        // Read structure with barcode, UMI, and template
        let read_structures = vec![ReadStructure::from_str("8B8M84T").unwrap()];
        let s1_barcode = "AAAAAAAA";
        let umi = "GGGGGGGG";
        let template = "T".repeat(84);
        let full_read = s1_barcode.to_owned() + umi + &template;
        let sample_metadata = metadata_file(&tmp, &[s1_barcode]);
        let input_files = vec![fastq_file(&tmp, "ex", "ex", &[&full_read])];

        let output_dir = tmp.path().to_path_buf().join("output");

        // Include only M and T (UMI + template), exclude barcode
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
            skip_reasons: vec![],
            template_types: vec!['M', 'T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();

        let output_path = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads = read_fastq(&output_path);

        assert_eq!(fq_reads.len(), 1);
        // With template_types M T, the output bases should be UMI + template (no barcode). Because
        // the UMI is written into the bases, it is *not* also appended to the read name (no
        // trailing `:GGGGGGGG` on the name), so the same bases are not emitted twice.
        let expected_seq = umi.to_owned() + &template;
        assert_equal(
            &fq_reads[0],
            &OwnedRecord {
                head: b"ex_0 1:N:0:AAAAAAAA".to_vec(),
                seq: expected_seq.as_bytes().to_vec(),
                qual: ";".repeat(92).as_bytes().to_vec(),
            },
        );
    }

    #[test]
    #[should_panic(expected = "--template-types with non-T segments requires T in --output-types")]
    fn test_template_types_requires_template_in_output_types() {
        let tmp = TempDir::new().unwrap();
        let read_structures = vec![ReadStructure::from_str("8B92T").unwrap()];
        let s1_barcode = "AAAAAAAA";
        let sample_metadata = metadata_file(&tmp, &[s1_barcode]);
        let input_files = vec![fastq_file(&tmp, "ex", "ex", &["A".repeat(100).as_str()])];

        let output_dir = tmp.path().to_path_buf().join("output");

        // This should fail: output_types has only B, but template_types has B T
        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output_types: vec!['B'],
            output: output_dir.clone(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 5,
            compression_level: 5,
            skip_reasons: vec![],
            template_types: vec!['B', 'T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();
    }

    #[test]
    #[should_panic(expected = "--template-types must include T (template segment)")]
    fn test_template_types_must_include_template() {
        let tmp = TempDir::new().unwrap();
        let read_structures = vec![ReadStructure::from_str("8B8M84T").unwrap()];
        let s1_barcode = "AAAAAAAA";
        let sample_metadata = metadata_file(&tmp, &[s1_barcode]);
        let input_files = vec![fastq_file(&tmp, "ex", "ex", &["A".repeat(100).as_str()])];

        let output_dir = tmp.path().to_path_buf().join("output");

        // This should fail: template_types doesn't include T
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
            skip_reasons: vec![],
            template_types: vec!['B', 'M'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();
    }

    #[test]
    #[should_panic(
        expected = "--template-types includes C but no read structure contains that segment type"
    )]
    fn test_template_types_absent_type_rejected() {
        let tmp = TempDir::new().unwrap();
        // Read structure has B and T but no C (cellular barcode).
        let read_structures = vec![ReadStructure::from_str("8B92T").unwrap()];
        let s1_barcode = "AAAAAAAA";
        let sample_metadata = metadata_file(&tmp, &[s1_barcode]);
        let input_files = vec![fastq_file(&tmp, "ex", "ex", &["A".repeat(100).as_str()])];

        let output_dir = tmp.path().to_path_buf().join("output");

        // This should fail: C is requested but absent from every read structure.
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
            skip_reasons: vec![],
            template_types: vec!['C', 'T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();
    }

    #[test]
    #[should_panic(expected = "Error parsing template segment types")]
    fn test_template_types_invalid_segment_type() {
        let tmp = TempDir::new().unwrap();
        let read_structures = vec![ReadStructure::from_str("8B92T").unwrap()];
        let s1_barcode = "AAAAAAAA";
        let sample_metadata = metadata_file(&tmp, &[s1_barcode]);
        let input_files = vec![fastq_file(&tmp, "ex", "ex", &["A".repeat(100).as_str()])];

        let output_dir = tmp.path().to_path_buf().join("output");

        // This should fail: 'X' is not a valid segment type
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
            skip_reasons: vec![],
            template_types: vec!['X', 'T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();
    }

    #[test]
    #[should_panic(
        expected = "--template-types with non-T segments is not supported when a read structure has multiple template segments (read structure #1 has 2)"
    )]
    fn test_template_types_non_t_rejects_multi_template_read_structure() {
        let tmp = TempDir::new().unwrap();
        // Read structure with multiple T segments in one source
        let read_structures = vec![ReadStructure::from_str("8B20T20S20T").unwrap()];
        let s1_barcode = "AAAAAAAA";
        let sample_metadata = metadata_file(&tmp, &[s1_barcode]);
        let input_files = vec![fastq_file(
            &tmp,
            "ex",
            "ex",
            &[&(s1_barcode.to_owned() + &"A".repeat(20) + &"C".repeat(20) + &"T".repeat(20))],
        )];

        let output_dir = tmp.path().to_path_buf().join("output");

        // This should fail: template_types includes B (non-T) and read structure has multiple T segments
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
            skip_reasons: vec![],
            template_types: vec!['B', 'T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();
    }

    #[test]
    #[should_panic(
        expected = "--template-types includes M but read structure #1 contains M with no template (T) segment to merge it into"
    )]
    fn test_template_types_non_colocated_segment_rejected() {
        let tmp = TempDir::new().unwrap();
        // The UMI (M) is on its own read (R1), separate from the template (R2). Since segments are
        // only merged within the same read, M cannot be folded into the template and the request
        // must be rejected rather than silently dropping the UMI bases.
        let read_structures =
            vec![ReadStructure::from_str("8M").unwrap(), ReadStructure::from_str("100T").unwrap()];
        let s1_barcode = "AAAAAAAA";
        let sample_metadata = metadata_file(&tmp, &[s1_barcode]);
        let input_files = vec![
            fastq_file(&tmp, "r1", "r1", &["GGGGGGGG"]),
            fastq_file(&tmp, "r2", "r2", &["A".repeat(100).as_str()]),
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
            skip_reasons: vec![],
            template_types: vec!['M', 'T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();
    }

    #[test]
    #[should_panic(
        expected = "Cellular barcode (C) segments are present but assigned to neither --output-types nor --template-types"
    )]
    fn test_cellular_barcode_without_destination_rejected() {
        let tmp = TempDir::new().unwrap();
        // A cellular barcode (C) is present but routed to neither --output-types nor
        // --template-types, and has no read-header field, so its bases would be silently lost.
        let read_structures = vec![ReadStructure::from_str("8B7C85T").unwrap()];
        let s1_barcode = "AAAAAAAA";
        let sample_metadata = metadata_file(&tmp, &[s1_barcode]);
        let input_files = vec![fastq_file(&tmp, "ex", "ex", &["A".repeat(100).as_str()])];

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
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();
    }

    #[test]
    fn test_template_types_all_segments_with_skip() {
        let tmp = TempDir::new().unwrap();
        // Read structure with barcode, skip, UMI, and template
        let read_structures = vec![ReadStructure::from_str("8B4S8M80T").unwrap()];
        let s1_barcode = "AAAAAAAA";
        let skip = "NNNN";
        let umi = "GGGGGGGG";
        let template = "T".repeat(80);
        let full_read = s1_barcode.to_owned() + skip + umi + &template;
        let sample_metadata = metadata_file(&tmp, &[s1_barcode]);
        let input_files = vec![fastq_file(&tmp, "ex", "ex", &[&full_read])];

        let output_dir = tmp.path().to_path_buf().join("output");

        // Include all segment types including skip
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
            skip_reasons: vec![],
            template_types: vec!['B', 'S', 'M', 'T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();

        let output_path = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads = read_fastq(&output_path);

        assert_eq!(fq_reads.len(), 1);
        // With all template_types, the output should be the full original read. The UMI (M) is
        // folded into the bases, so it is not also appended to the read name.
        assert_equal(
            &fq_reads[0],
            &OwnedRecord {
                head: b"ex_0 1:N:0:AAAAAAAA".to_vec(),
                seq: full_read.as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            },
        );
    }

    #[test]
    fn test_template_types_with_cellular_barcode() {
        let tmp = TempDir::new().unwrap();
        let read_structures = vec![ReadStructure::from_str("10M8B7C75T").unwrap()];
        let s1_barcode = "AAAAAAAA";
        let s1_umi = "ATCGATCGAT";
        let s1_cellular_barcode = "GATTACA";
        let template = "T".repeat(75);
        let sample_metadata = metadata_file(&tmp, &[s1_barcode]);
        let input_files = vec![fastq_file(
            &tmp,
            "ex",
            "ex",
            &[&format!("{s1_umi}{s1_barcode}{s1_cellular_barcode}{template}")],
        )];

        let output_dir = tmp.path().to_path_buf().join("output");

        // template_types includes C (cellular barcode) and T; output should be C+T (82 bp).
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
            skip_reasons: vec![],
            template_types: vec!['C', 'T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();

        let output_path = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads = read_fastq(&output_path);

        assert_eq!(fq_reads.len(), 1);
        let expected_seq = s1_cellular_barcode.to_owned() + &template;
        assert_equal(
            &fq_reads[0],
            &OwnedRecord {
                head: b"ex_0:ATCGATCGAT 1:N:0:AAAAAAAA".to_vec(),
                seq: expected_seq.as_bytes().to_vec(),
                qual: ";".repeat(82).as_bytes().to_vec(),
            },
        );
    }

    #[test]
    fn test_header_format_unmodified_leaves_header_unchanged() {
        // Identical setup to test_template_types_with_cellular_barcode, but with
        // --header-format unmodified. The output read name must be the input header verbatim (no
        // UMI appended to the name, no sample barcode in the comment) while the read bases are
        // unaffected by the format.
        let tmp = TempDir::new().unwrap();
        let read_structures = vec![ReadStructure::from_str("10M8B7C75T").unwrap()];
        let s1_barcode = "AAAAAAAA";
        let s1_umi = "ATCGATCGAT";
        let s1_cellular_barcode = "GATTACA";
        let template = "T".repeat(75);
        let sample_metadata = metadata_file(&tmp, &[s1_barcode]);
        let input_files = vec![fastq_file(
            &tmp,
            "ex",
            "ex",
            &[&format!("{s1_umi}{s1_barcode}{s1_cellular_barcode}{template}")],
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
            skip_reasons: vec![],
            template_types: vec!['C', 'T'],
            header_format: HeaderFormatKind::Unmodified,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();

        let output_path = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads = read_fastq(&output_path);

        assert_eq!(fq_reads.len(), 1);
        let expected_seq = s1_cellular_barcode.to_owned() + &template;
        assert_equal(
            &fq_reads[0],
            &OwnedRecord {
                head: b"ex_0".to_vec(),
                seq: expected_seq.as_bytes().to_vec(),
                qual: ";".repeat(82).as_bytes().to_vec(),
            },
        );
    }

    /// Builds a paired-end, single-sample `Demux` over two reads whose headers differ only in the
    /// read-number field, with an 8bp UMI on each end (duplex) and an inline sample barcode on R1.
    /// Returns the configured `Demux` and its output directory; the caller sets `header_format`
    /// and `umi_in_name`.
    fn paired_demux_with_header_format(
        tmp: &TempDir,
        header_format: HeaderFormatKind,
        umi_in_name: Option<bool>,
    ) -> (Demux, PathBuf) {
        let read_structures = vec![
            ReadStructure::from_str("8M8B20T").unwrap(),
            ReadStructure::from_str("8M20T").unwrap(),
        ];
        let s1_barcode = "ACGTACGT";
        let (r1_umi, r2_umi) = ("AAAAAAAA", "CCCCCCCC");
        let template = "T".repeat(20);
        let sample_metadata = metadata_file(tmp, &[s1_barcode]);
        let input_files = vec![
            fastq_file_with_header(
                tmp,
                "r1",
                "inst:1:FC:1:1101:5:7 1:N:0:0",
                &format!("{r1_umi}{s1_barcode}{template}"),
            ),
            fastq_file_with_header(
                tmp,
                "r2",
                "inst:1:FC:1:1101:5:7 2:N:0:0",
                &format!("{r2_umi}{template}"),
            ),
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
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format,
            umi_in_name,
        };
        (demux, output_dir)
    }

    /// Covers every `(--header-format, --umi-in-name)` combination whose expected output is a
    /// single pair of R1/R2 header assertions (the `name-only` cases don't fit this shape --
    /// see `test_header_format_name_only_drops_comment_even_with_umi_in_name` below -- since one
    /// of them needs two sub-scenarios in the same test to show the comment stays dropped even
    /// once a UMI starts appearing in the name).
    #[rstest]
    // unmodified, no override: regression test that each output carries the header of ITS OWN
    // source, so R2 keeps its "2:N" read-number field rather than inheriting R1's "1:N".
    // Verbatim: no UMI in the name, no barcode in the comment, read-number field unchanged.
    #[case(
        HeaderFormatKind::Unmodified,
        None,
        b"inst:1:FC:1:1101:5:7 1:N:0:0",
        b"inst:1:FC:1:1101:5:7 2:N:0:0"
    )]
    // unmodified + --umi-in-name true: appends the (duplex) UMI to the name, but leaves the
    // comment verbatim -- including R2's "2:N" -- and does not add the sample barcode.
    #[case(
        HeaderFormatKind::Unmodified,
        Some(true),
        b"inst:1:FC:1:1101:5:7:AAAAAAAA+CCCCCCCC 1:N:0:0",
        b"inst:1:FC:1:1101:5:7:AAAAAAAA+CCCCCCCC 2:N:0:0"
    )]
    // The default `illumina` format (with `--umi-in-name` unset) appends the UMI to the name AND
    // the barcode to the comment, and rewrites the read-number field per output file.  This is
    // the byte-identical-to-0.4.x invariant's expected output; see
    // `test_bare_cli_invocation_defaults_to_illumina_with_umi_in_name_unset` for the companion
    // test that pins `HeaderFormatKind::Illumina`/`None` as the actual clap defaults a bare
    // `fqtk demux` invocation gets.
    #[case(
        HeaderFormatKind::default(),
        None,
        b"inst:1:FC:1:1101:5:7:AAAAAAAA+CCCCCCCC 1:N:0:ACGTACGT",
        b"inst:1:FC:1:1101:5:7:AAAAAAAA+CCCCCCCC 2:N:0:ACGTACGT"
    )]
    // illumina + explicit --umi-in-name false: only reachable now that the UMI-in-name axis is
    // independent of the header format.  The comment is still rewritten and the sample barcode
    // still appended, but the UMI is omitted from the name entirely.
    #[case(
        HeaderFormatKind::Illumina,
        Some(false),
        b"inst:1:FC:1:1101:5:7 1:N:0:ACGTACGT",
        b"inst:1:FC:1:1101:5:7 2:N:0:ACGTACGT"
    )]
    fn test_header_format_paired_end_headers(
        #[case] header_format: HeaderFormatKind,
        #[case] umi_in_name: Option<bool>,
        #[case] expected_r1: &[u8],
        #[case] expected_r2: &[u8],
    ) {
        let tmp = TempDir::new().unwrap();
        let (demux, output_dir) = paired_demux_with_header_format(&tmp, header_format, umi_in_name);
        demux.execute().unwrap();

        let r1 = read_fastq(&output_dir.join("Sample0000.R1.fq.gz"));
        let r2 = read_fastq(&output_dir.join("Sample0000.R2.fq.gz"));
        assert_eq!(r1.len(), 1);
        assert_eq!(r2.len(), 1);
        assert_eq!(r1[0].head, expected_r1);
        assert_eq!(r2[0].head, expected_r2);
    }

    #[test]
    fn test_header_format_name_only_drops_comment_even_with_umi_in_name() {
        // `name-only` always drops the comment; with no UMI in the name either (the default),
        // R1/R2 headers become byte-identical.
        let tmp = TempDir::new().unwrap();
        let (demux, output_dir) =
            paired_demux_with_header_format(&tmp, HeaderFormatKind::NameOnly, None);
        demux.execute().unwrap();

        let r1 = read_fastq(&output_dir.join("Sample0000.R1.fq.gz"));
        let r2 = read_fastq(&output_dir.join("Sample0000.R2.fq.gz"));
        assert_eq!(r1[0].head, b"inst:1:FC:1:1101:5:7");
        assert_eq!(r2[0].head, b"inst:1:FC:1:1101:5:7");

        // With `--umi-in-name true`, the UMI is folded into the (now-identical) name, but the
        // comment is still always dropped.
        let (demux2, output_dir2) =
            paired_demux_with_header_format(&tmp, HeaderFormatKind::NameOnly, Some(true));
        demux2.execute().unwrap();
        let r1_with_umi = read_fastq(&output_dir2.join("Sample0000.R1.fq.gz"));
        assert_eq!(r1_with_umi[0].head, b"inst:1:FC:1:1101:5:7:AAAAAAAA+CCCCCCCC");
    }

    #[test]
    fn test_bare_cli_invocation_defaults_to_illumina_with_umi_in_name_unset() {
        // Pins the actual `#[clap(...)]` defaults a bare `fqtk demux` invocation gets, which the
        // `HeaderFormatKind::default()`/`None` case of `test_header_format_paired_end_headers`
        // does not: that case builds a `Demux` struct literal directly and would not catch a
        // regression in the clap attributes themselves (e.g. an accidentally changed
        // `default_value_t`).
        let demux = Demux::try_parse_from([
            "demux", "-i", "r1.fq", "-r", "8M8B20T", "-s", "meta.tsv", "-o", "outdir",
        ])
        .unwrap();
        assert_eq!(demux.header_format, HeaderFormatKind::Illumina);
        assert_eq!(demux.umi_in_name, None);
    }

    #[test]
    fn test_header_format_unmodified_accepts_non_illumina_header() {
        // A header with more than eight `:`-delimited name fields (e.g. MGI/ONT) aborts the
        // default illumina path but is emitted verbatim under `unmodified`.
        let tmp = TempDir::new().unwrap();
        let read_structures = vec![ReadStructure::from_str("8B20T").unwrap()];
        let s1_barcode = "ACGTACGT";
        let sample_metadata = metadata_file(&tmp, &[s1_barcode]);
        let weird_header = "a:b:c:d:e:f:g:h:i";
        let input_files = vec![fastq_file_with_header(
            &tmp,
            "n1",
            weird_header,
            &format!("{s1_barcode}{}", "T".repeat(20)),
        )];
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
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Unmodified,
            umi_in_name: None,
        };
        demux.execute().unwrap();

        let r1 = read_fastq(&output_dir.join("Sample0000.R1.fq.gz"));
        assert_eq!(r1[0].head, weird_header.as_bytes());
    }

    #[test]
    fn test_header_format_kind_display_round_trips_through_value_enum_parsing() {
        // `Display` feeds the user-facing warning text that recommends a `--header-format` value
        // (see `warn_on_discarded_segments`); pin that its output is always something
        // `clap::ValueEnum` will parse back to the same variant, so a future rename of a variant
        // or the `Display` impl can't silently recommend a value clap would reject.
        for kind in HeaderFormatKind::value_variants() {
            let rendered = kind.to_string();
            assert_eq!(HeaderFormatKind::from_str(&rendered, false).unwrap(), *kind);
        }
    }

    #[test]
    fn test_discarded_segment_types() {
        // Read structure with a UMI (M), a sample barcode (B), and a template (T).
        let rs = vec![ReadStructure::from_str("8M8B84T").unwrap()];
        let none: HashSet<SegmentType> = HashSet::new();
        let discarded = |header_format, umi_in_name, output: &[char], template: &[char]| {
            let to_set =
                |cs: &[char]| cs.iter().map(|c| SegmentType::try_from(*c).unwrap()).collect();
            let mut d = Demux::discarded_segment_types(
                &rs,
                header_format,
                umi_in_name,
                &to_set(output),
                &to_set(template),
            );
            d.sort_by_key(|s| *s as u8);
            d
        };

        // illumina discards nothing: UMI goes in the name, barcode in the comment.
        assert_eq!(discarded(HeaderFormatKind::Illumina, true, &['T'], &['T']), vec![]);

        // illumina with `--umi-in-name false` still writes the barcode (the comment is always
        // rewritten), but now drops the UMI -- unreachable under the old fused enum, where
        // `rewrite` always implied UMI-in-name.
        assert_eq!(
            discarded(HeaderFormatKind::Illumina, false, &['T'], &['T']),
            vec![SegmentType::MolecularBarcode]
        );

        // unmodified with no UMI in the name and only T output drops both the UMI and the barcode.
        // Sorted ascending by discriminant byte: SampleBarcode (b'B') before MolecularBarcode (b'M').
        assert_eq!(
            discarded(HeaderFormatKind::Unmodified, false, &['T'], &['T']),
            vec![SegmentType::SampleBarcode, SegmentType::MolecularBarcode]
        );

        // unmodified + umi-in-name keeps the UMI in the name, so only the sample barcode is
        // dropped.
        assert_eq!(
            discarded(HeaderFormatKind::Unmodified, true, &['T'], &['T']),
            vec![SegmentType::SampleBarcode]
        );

        // name-only takes the same `!rewrites_comment()` branch as unmodified: no UMI in the
        // name and only T output drops both the UMI and the barcode.
        assert_eq!(
            discarded(HeaderFormatKind::NameOnly, false, &['T'], &['T']),
            vec![SegmentType::SampleBarcode, SegmentType::MolecularBarcode]
        );

        // Routing M and B to their own FASTQs rescues them under unmodified.
        assert_eq!(
            discarded(HeaderFormatKind::Unmodified, false, &['T', 'M', 'B'], &['T']),
            vec![]
        );

        // Folding M and B into the template bases also rescues them.
        assert_eq!(
            discarded(HeaderFormatKind::Unmodified, false, &['T'], &['T', 'M', 'B']),
            vec![]
        );

        // A read structure without a UMI or barcode discards nothing regardless of format.
        let template_only = vec![ReadStructure::from_str("100T").unwrap()];
        assert_eq!(
            Demux::discarded_segment_types(
                &template_only,
                HeaderFormatKind::Unmodified,
                false,
                &none,
                &none
            ),
            vec![]
        );
    }

    #[test]
    fn test_output_type_reads() {
        let tmp = TempDir::new().unwrap();
        let read_structures = vec![ReadStructure::from_str("10M8B7C100T").unwrap()];
        let s1_barcode = "AAAAAAAA";
        let s1_umi = "ATCGATCGAT";
        let s1_cellular_barcode = "GATTACA";
        let sample_metadata =
            metadata_file(&tmp, &[s1_barcode, "CCCCCCCC", "GGGGGGGG", "TTTTTTTT"]);
        let template_sequence = "A".repeat(100);

        let input_files = vec![fastq_file(
            &tmp,
            "ex",
            "ex",
            &[&format!("{s1_umi}{s1_barcode}{s1_cellular_barcode}{template_sequence}")],
        )];

        let output_dir = tmp.path().to_path_buf().join("output");

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output_types: vec!['T', 'B', 'M', 'C'],
            output: output_dir.clone(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 5,
            compression_level: 5,
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();

        let output_path = output_dir.join("Sample0000.R1.fq.gz");
        let barcode_output_path = output_dir.join("Sample0000.I1.fq.gz");
        let umi_output_path = output_dir.join("Sample0000.U1.fq.gz");
        let cellular_barcode_output_path = output_dir.join("Sample0000.C1.fq.gz");
        let fq_reads = read_fastq(&output_path);
        let barcode_fq_reads = read_fastq(&barcode_output_path);
        let umi_fq_reads = read_fastq(&umi_output_path);
        let cellular_barcode_fq_reads = read_fastq(&cellular_barcode_output_path);

        assert_eq!(fq_reads.len(), 1);
        assert_equal(
            &fq_reads[0],
            &OwnedRecord {
                head: b"ex_0:ATCGATCGAT 1:N:0:AAAAAAAA".to_vec(),
                seq: "A".repeat(100).as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            },
        );

        assert_eq!(barcode_fq_reads.len(), 1);
        assert_equal(
            &barcode_fq_reads[0],
            &OwnedRecord {
                head: b"ex_0:ATCGATCGAT 1:N:0:AAAAAAAA".to_vec(),
                seq: b"AAAAAAAA".to_vec(),
                qual: ";".repeat(8).as_bytes().to_vec(),
            },
        );

        assert_eq!(umi_fq_reads.len(), 1);
        assert_equal(
            &umi_fq_reads[0],
            &OwnedRecord {
                head: b"ex_0:ATCGATCGAT 1:N:0:AAAAAAAA".to_vec(),
                seq: b"ATCGATCGAT".to_vec(),
                qual: ";".repeat(10).as_bytes().to_vec(),
            },
        );

        assert_eq!(cellular_barcode_fq_reads.len(), 1);
        assert_equal(
            &cellular_barcode_fq_reads[0],
            &OwnedRecord {
                head: b"ex_0:ATCGATCGAT 1:N:0:AAAAAAAA".to_vec(),
                seq: b"GATTACA".to_vec(),
                qual: ";".repeat(7).as_bytes().to_vec(),
            },
        );
    }

    #[test]
    fn test_demux_with_catchall_barcode() {
        let tmp = TempDir::new().unwrap();
        let read_structures = vec![ReadStructure::from_str("7B+T").unwrap()];
        let s1_barcode = "NNNNNNN";
        let sample_metadata = metadata_file(&tmp, &[s1_barcode]);
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
            max_mismatches: 0,
            min_mismatch_delta: 2,
            threads: 5,
            compression_level: 5,
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();

        let unmatched_path = output_dir.join("unmatched.R1.fq.gz");
        let unmatched_reads = read_fastq(&unmatched_path);
        assert_eq!(unmatched_reads.len(), 0);

        let output_path = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads = read_fastq(&output_path);

        assert_eq!(fq_reads.len(), 1);
        assert_equal(
            &fq_reads[0],
            &OwnedRecord {
                head: b"ex_0 1:N:0:NNNNNNN".to_vec(),
                seq: "A".repeat(100).as_bytes().to_vec(),
                qual: ";".repeat(100).as_bytes().to_vec(),
            },
        );
    }

    #[test]
    fn test_demux_with_iupac_bases_in_barcode() {
        let tmp = TempDir::new().unwrap();
        let read_structures = vec![ReadStructure::from_str("7B+T").unwrap()];
        let s1_barcode = "MMMMMMM";
        let s2_barcode = "KKKKKKK";
        let sample_metadata = metadata_file(&tmp, &[s1_barcode, s2_barcode]);
        let input_files = vec![fastq_file(
            &tmp,
            "ex",
            "ex",
            &[
                &("AAAAAAA".to_owned() + &"A".repeat(5)), // barcode s1
                &("CCCCCCC".to_owned() + &"A".repeat(5)), // barcode s1
                &("ACACACA".to_owned() + &"A".repeat(5)), // barcode s1
                &("GTGTGTG".to_owned() + &"C".repeat(5)), // barcode s2
                &("TGTGTGT".to_owned() + &"C".repeat(5)), // barcode s2
                &("CGCGCGC".to_owned() + &"T".repeat(5)), // unmatched
            ],
        )];

        let output_dir = tmp.path().to_path_buf().join("output");

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output_types: vec!['T'],
            output: output_dir.clone(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 0,
            min_mismatch_delta: 0,
            threads: 5,
            compression_level: 5,
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();

        let output_path = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads = read_fastq(&output_path);
        assert_eq!(fq_reads.len(), 3);
        assert_equal(
            &fq_reads[0],
            &OwnedRecord {
                head: b"ex_0 1:N:0:AAAAAAA".to_vec(),
                seq: "A".repeat(5).as_bytes().to_vec(),
                qual: ";".repeat(5).as_bytes().to_vec(),
            },
        );

        let output_path = output_dir.join("Sample0001.R1.fq.gz");
        let fq_reads = read_fastq(&output_path);
        assert_eq!(fq_reads.len(), 2);
        assert_equal(
            &fq_reads[0],
            &OwnedRecord {
                head: b"ex_3 1:N:0:GTGTGTG".to_vec(),
                seq: "C".repeat(5).as_bytes().to_vec(),
                qual: ";".repeat(5).as_bytes().to_vec(),
            },
        );

        // Should not match since it has 3 no calls, and barcodes have at maximum 1 no-call
        let unmatched_path = output_dir.join("unmatched.R1.fq.gz");
        let unmatched_reads = read_fastq(&unmatched_path);
        assert_eq!(unmatched_reads.len(), 1);
        assert_equal(
            &unmatched_reads[0],
            &OwnedRecord {
                head: b"ex_5 1:N:0:CGCGCGC".to_vec(),
                seq: "T".repeat(5).as_bytes().to_vec(),
                qual: ";".repeat(5).as_bytes().to_vec(),
            },
        );
    }

    #[test]
    fn test_demux_with_ns_in_barcode() {
        let tmp = TempDir::new().unwrap();
        let read_structures = vec![ReadStructure::from_str("7B+T").unwrap()];
        let s1_barcode = "NNAAAAA";
        let s2_barcode = "NNCCCCC";
        let sample_metadata = metadata_file(&tmp, &[s1_barcode, s2_barcode]);
        let input_files = vec![fastq_file(
            &tmp,
            "ex",
            "ex",
            &[
                &("ANAAAAA".to_owned() + &"A".repeat(5)),
                &("ANCCCCC".to_owned() + &"C".repeat(5)),
                &("NNNAAAA".to_owned() + &"T".repeat(5)),
            ],
        )];

        let output_dir = tmp.path().to_path_buf().join("output");

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output_types: vec!['T'],
            output: output_dir.clone(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 0,
            min_mismatch_delta: 0,
            threads: 5,
            compression_level: 5,
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();

        let output_path = output_dir.join("Sample0000.R1.fq.gz");
        let fq_reads = read_fastq(&output_path);
        assert_eq!(fq_reads.len(), 1);
        assert_equal(
            &fq_reads[0],
            &OwnedRecord {
                head: b"ex_0 1:N:0:ANAAAAA".to_vec(),
                seq: "A".repeat(5).as_bytes().to_vec(),
                qual: ";".repeat(5).as_bytes().to_vec(),
            },
        );

        let output_path = output_dir.join("Sample0001.R1.fq.gz");
        let fq_reads = read_fastq(&output_path);
        assert_eq!(fq_reads.len(), 1);
        assert_equal(
            &fq_reads[0],
            &OwnedRecord {
                head: b"ex_1 1:N:0:ANCCCCC".to_vec(),
                seq: "C".repeat(5).as_bytes().to_vec(),
                qual: ";".repeat(5).as_bytes().to_vec(),
            },
        );

        // Should not match since it has 3 no calls, and barcodes have at maximum 1 no-call
        let unmatched_path = output_dir.join("unmatched.R1.fq.gz");
        let unmatched_reads = read_fastq(&unmatched_path);
        assert_eq!(unmatched_reads.len(), 1);
        assert_equal(
            &unmatched_reads[0],
            &OwnedRecord {
                head: b"ex_2 1:N:0:NNNAAAA".to_vec(),
                seq: "T".repeat(5).as_bytes().to_vec(),
                qual: ";".repeat(5).as_bytes().to_vec(),
            },
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
            threads: 5,
            compression_level: 5,
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
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

        let output_dir: PathBuf = tmp.path().to_path_buf().join("output");

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
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
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
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
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
            threads: 5,
            compression_level: 5,
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
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
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
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
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
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
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();
    }

    #[test]
    #[should_panic(
        expected = "Read ex_2 had too few bases to demux 0 vs. 1 needed in read structure +T."
    )]
    fn test_fails_if_reads_too_short() {
        let tmp = TempDir::new().unwrap();
        let read_structures =
            vec![ReadStructure::from_str("+T").unwrap(), ReadStructure::from_str("7B").unwrap()];

        let records = [
            vec!["AAAAAAA", &SAMPLE1_BARCODE[0..7]], // barcode too short
            vec!["CCCCCCC", SAMPLE1_BARCODE],        // barcode the correct length
            vec!["", SAMPLE1_BARCODE],               // template basese too short
            vec!["G", SAMPLE1_BARCODE],              // barcode the correct length
        ];

        let input_files = vec![
            fastq_file(&tmp, "read1", "ex", &[records[0][0], records[1][0], records[2][0]]),
            fastq_file(&tmp, "index1", "ex", &[records[0][1], records[1][1], records[2][1]]),
        ];
        let sample_metadata = metadata(&tmp);
        let output_dir: PathBuf = tmp.path().to_path_buf().join("output");

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output_types: vec!['T', 'B'],
            output: output_dir,
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 5,
            compression_level: 5,
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();
    }

    #[test]
    fn test_skip_reads_too_short() {
        let tmp = TempDir::new().unwrap();
        let read_structures =
            vec![ReadStructure::from_str("+T").unwrap(), ReadStructure::from_str("7B").unwrap()];

        let records = [
            vec!["AAAAAAA", &SAMPLE1_BARCODE[0..7]], // barcode too short
            vec!["CCCCCCC", SAMPLE1_BARCODE],        // barcode the correct length
            vec!["", SAMPLE1_BARCODE],               // template basese too short
            vec!["G", SAMPLE1_BARCODE],              // barcode the correct length
        ];

        let input_files = vec![
            fastq_file(&tmp, "read1", "ex", &[records[0][0], records[1][0], records[2][0]]),
            fastq_file(&tmp, "index1", "ex", &[records[0][1], records[1][1], records[2][1]]),
        ];
        let sample_metadata = metadata(&tmp);
        let output_dir: PathBuf = tmp.path().to_path_buf().join("output");

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            sample_metadata,
            output_types: vec!['T', 'B'],
            output: output_dir.clone(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 5,
            compression_level: 5,
            skip_reasons: vec![SkipReason::TooFewBases],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();

        // Write out the metrics
        let metrics_path = output_dir.join("demux-metrics.txt");
        let metrics: Vec<DemuxMetric> = DelimFile::default().read_tsv(&metrics_path).unwrap();
        let demuxed_reads: usize = metrics.iter().map(|m| m.templates).sum();
        let sample0000_reads =
            metrics.iter().find(|m| m.sample_id == "Sample0000").unwrap().templates;
        assert_eq!(demuxed_reads, 2);
        assert_eq!(sample0000_reads, 2);

        let r1_path = output_dir.join("Sample0000.R1.fq.gz");
        let r1_reads = read_fastq(&r1_path);
        assert_eq!(r1_reads.len(), 2);

        let r2_path = output_dir.join("Sample0000.I1.fq.gz");
        let r2_reads = read_fastq(&r2_path);
        assert_eq!(r2_reads.len(), 2);
    }

    ////////////////////////////////////////////////////////////////////////////
    // Tests for ReadSet::write_header_internal
    ////////////////////////////////////////////////////////////////////////////

    fn seg(bases: &[u8], segment_type: SegmentType) -> FastqSegment {
        seg_with_source(bases, segment_type, 0)
    }

    fn seg_with_source(
        bases: &[u8],
        segment_type: SegmentType,
        source_index: usize,
    ) -> FastqSegment {
        let quals = vec![b'#'; bases.len()];
        FastqSegment { seq: bases.to_vec(), quals, segment_type, source_index }
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
            HeaderFormat { include_umi: true, kind: HeaderFormatKind::Illumina, umi_in_name: true },
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
            HeaderFormat { include_umi: true, kind: HeaderFormatKind::Illumina, umi_in_name: true },
        )
        .unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    #[test]
    fn test_write_header_omits_umi_when_in_template_bases() {
        // When include_umi is false (UMI bases are written into the template output instead), the
        // UMI must not be appended to the read name, but sample barcodes are still appended.
        let mut out = Vec::new();
        let header = b"inst:123:ABCDE:1:204:1022:2108 1:Y:0:0";
        let barcode_segs =
            [seg(b"ACGT", SegmentType::SampleBarcode), seg(b"GGTT", SegmentType::SampleBarcode)];
        let umi_segs = [seg(b"AACCGGTT", SegmentType::MolecularBarcode)];
        let expected = "@inst:123:ABCDE:1:204:1022:2108 2:Y:0:ACGT+GGTT".to_string();
        ReadSet::write_header_internal(
            &mut out,
            2,
            header,
            barcode_segs.iter().filter(|_| true),
            umi_segs.iter().filter(|_| true),
            HeaderFormat {
                include_umi: false,
                kind: HeaderFormatKind::Illumina,
                umi_in_name: true,
            },
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
            HeaderFormat { include_umi: true, kind: HeaderFormatKind::Illumina, umi_in_name: true },
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
            HeaderFormat { include_umi: true, kind: HeaderFormatKind::Illumina, umi_in_name: true },
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
            HeaderFormat { include_umi: true, kind: HeaderFormatKind::Illumina, umi_in_name: true },
        )
        .unwrap();
    }

    #[test]
    fn test_write_header_comment_too_few_parts() {
        let mut out = Vec::new();
        let header = b"q1 0:0";
        let barcode_segs =
            [seg(b"ACGT", SegmentType::SampleBarcode), seg(b"GGTT", SegmentType::SampleBarcode)];
        let umi_segs = [seg(b"AACCGGTT", SegmentType::MolecularBarcode)];
        let expected = "@q1:AACCGGTT 0:0:ACGT+GGTT".to_string();
        ReadSet::write_header_internal(
            &mut out,
            1,
            header,
            barcode_segs.iter().filter(|_| true),
            umi_segs.iter().filter(|_| true),
            HeaderFormat { include_umi: true, kind: HeaderFormatKind::Illumina, umi_in_name: true },
        )
        .unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    #[test]
    fn test_write_header_illumina_with_empty_comment_after_trailing_space() {
        // A header ending in a trailing space (comment is `Some(&[])`, not `None`) used to panic
        // on `*chars.last().unwrap()`; an empty comment is now treated like "didn't already end
        // in a colon", so a separator is still written before the sample barcode.
        let mut out = Vec::new();
        let header = b"inst:123:ABCDE:1:204:1022:2108 ";
        let barcode_segs =
            [seg(b"ACGT", SegmentType::SampleBarcode), seg(b"GGTT", SegmentType::SampleBarcode)];
        let umi_segs = [];
        let expected = "@inst:123:ABCDE:1:204:1022:2108 :ACGT+GGTT".to_string();
        ReadSet::write_header_internal(
            &mut out,
            1,
            header,
            barcode_segs.iter().filter(|_| true),
            umi_segs.iter().filter(|_| true),
            HeaderFormat { include_umi: true, kind: HeaderFormatKind::Illumina, umi_in_name: true },
        )
        .unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    // Renamed from `..._keep_original_read_names` (pre-dates the `HeaderFormatKind` rename).
    #[test]
    fn test_write_header_unmodified_verbatim() {
        // `unmodified` with no UMI in the name emits the header verbatim: no UMI appended to the
        // name, no sample barcode in the comment, and no read-number rewrite in the comment.
        let mut out = Vec::new();
        let header = b"inst:123:ABCDE:1:204:1022:2108 1:N:0:0";
        let barcode_segs =
            [seg(b"ACGT", SegmentType::SampleBarcode), seg(b"GGTT", SegmentType::SampleBarcode)];
        let umi_segs = [seg(b"AACCGGTT", SegmentType::MolecularBarcode)];
        let expected = "@inst:123:ABCDE:1:204:1022:2108 1:N:0:0".to_string();
        ReadSet::write_header_internal(
            &mut out,
            2,
            header,
            barcode_segs.iter().filter(|_| true),
            umi_segs.iter().filter(|_| true),
            HeaderFormat {
                include_umi: true,
                kind: HeaderFormatKind::Unmodified,
                umi_in_name: false,
            },
        )
        .unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    // Renamed from `..._name_only` (pre-dates and is unrelated to `HeaderFormatKind::NameOnly` --
    // this exercises a header with no comment at all, under `unmodified`, not the new format).
    #[test]
    fn test_write_header_unmodified_bare_name_no_comment() {
        // A header that is just a name (no comment) is emitted verbatim with only the `@` prefix.
        let mut out = Vec::new();
        let header = b"q1";
        let barcode_segs =
            [seg(b"ACGT", SegmentType::SampleBarcode), seg(b"GGTT", SegmentType::SampleBarcode)];
        let umi_segs = [seg(b"AACCGGTT", SegmentType::MolecularBarcode)];
        let expected = "@q1".to_string();
        ReadSet::write_header_internal(
            &mut out,
            1,
            header,
            barcode_segs.iter().filter(|_| true),
            umi_segs.iter().filter(|_| true),
            HeaderFormat {
                include_umi: true,
                kind: HeaderFormatKind::Unmodified,
                umi_in_name: false,
            },
        )
        .unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    // Renamed from `..._keep_with_umi` (pre-dates the `HeaderFormatKind`/`--umi-in-name` rework).
    #[test]
    fn test_write_header_unmodified_with_umi_in_name() {
        // `unmodified` + `umi_in_name: true` appends the UMI to the name but leaves the comment
        // verbatim: the read-number field is NOT rewritten (still `1:`), and no sample barcode is
        // appended.
        let mut out = Vec::new();
        let header = b"inst:123:ABCDE:1:204:1022:2108 1:N:0:0";
        let barcode_segs =
            [seg(b"ACGT", SegmentType::SampleBarcode), seg(b"GGTT", SegmentType::SampleBarcode)];
        let umi_segs = [
            seg(b"AACC", SegmentType::MolecularBarcode),
            seg(b"GGTT", SegmentType::MolecularBarcode),
        ];
        let expected = "@inst:123:ABCDE:1:204:1022:2108:AACC+GGTT 1:N:0:0".to_string();
        ReadSet::write_header_internal(
            &mut out,
            2,
            header,
            barcode_segs.iter().filter(|_| true),
            umi_segs.iter().filter(|_| true),
            HeaderFormat {
                include_umi: true,
                kind: HeaderFormatKind::Unmodified,
                umi_in_name: true,
            },
        )
        .unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    #[test]
    fn test_write_header_name_only_drops_comment() {
        // `name-only` drops the comment unconditionally, even when the input has one.
        let mut out = Vec::new();
        let header = b"inst:123:ABCDE:1:204:1022:2108 1:N:0:0";
        let barcode_segs =
            [seg(b"ACGT", SegmentType::SampleBarcode), seg(b"GGTT", SegmentType::SampleBarcode)];
        let umi_segs = [seg(b"AACCGGTT", SegmentType::MolecularBarcode)];
        let expected = "@inst:123:ABCDE:1:204:1022:2108".to_string();
        ReadSet::write_header_internal(
            &mut out,
            1,
            header,
            barcode_segs.iter().filter(|_| true),
            umi_segs.iter().filter(|_| true),
            HeaderFormat {
                include_umi: true,
                kind: HeaderFormatKind::NameOnly,
                umi_in_name: false,
            },
        )
        .unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    #[test]
    fn test_write_header_name_only_with_umi_in_name_still_drops_comment() {
        // `name-only` + `--umi-in-name true` folds the UMI into the name but still drops the
        // comment.
        let mut out = Vec::new();
        let header = b"inst:123:ABCDE:1:204:1022:2108 1:N:0:0";
        let barcode_segs =
            [seg(b"ACGT", SegmentType::SampleBarcode), seg(b"GGTT", SegmentType::SampleBarcode)];
        let umi_segs = [seg(b"AACCGGTT", SegmentType::MolecularBarcode)];
        let expected = "@inst:123:ABCDE:1:204:1022:2108:AACCGGTT".to_string();
        ReadSet::write_header_internal(
            &mut out,
            1,
            header,
            barcode_segs.iter().filter(|_| true),
            umi_segs.iter().filter(|_| true),
            HeaderFormat { include_umi: true, kind: HeaderFormatKind::NameOnly, umi_in_name: true },
        )
        .unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    // ############################################################################################
    // Test other ``ReadSet`` functions
    // ############################################################################################

    #[test]
    fn test_sample_barcode_sequence() {
        let segments = vec![
            seg("AGCT".as_bytes(), SegmentType::Template),
            seg("GATA".as_bytes(), SegmentType::SampleBarcode),
            seg("CAC".as_bytes(), SegmentType::SampleBarcode),
            seg("GACCCC".as_bytes(), SegmentType::MolecularBarcode),
        ];

        let read_set = read_set(segments);

        assert_eq!(read_set.sample_barcode_sequence(), "GATACAC".as_bytes().to_owned());
    }

    #[test]
    fn test_cellular_barcode_segments() {
        let segments = vec![
            seg("AGCT".as_bytes(), SegmentType::Template),
            seg("GATA".as_bytes(), SegmentType::SampleBarcode),
            seg("CAC".as_bytes(), SegmentType::SampleBarcode),
            seg("GACCCC".as_bytes(), SegmentType::MolecularBarcode),
            seg("ACGTACGT".as_bytes(), SegmentType::CellularBarcode),
        ];

        let read_set = read_set(segments);

        let expected = vec![seg("ACGTACGT".as_bytes(), SegmentType::CellularBarcode)];

        assert_eq!(expected, read_set.cellular_barcode_segments().cloned().collect::<Vec<_>>());
    }

    #[test]
    fn test_template_segments() {
        let segments = vec![
            seg("AGCT".as_bytes(), SegmentType::SampleBarcode),
            seg("GATA".as_bytes(), SegmentType::Template),
            seg("CAC".as_bytes(), SegmentType::Template),
            seg("GACCCC".as_bytes(), SegmentType::MolecularBarcode),
        ];
        let expected = vec![
            seg("GATA".as_bytes(), SegmentType::Template),
            seg("CAC".as_bytes(), SegmentType::Template),
        ];

        let read_set = read_set(segments);

        assert_eq!(expected, read_set.template_segments().cloned().collect::<Vec<_>>());
    }
    #[test]
    fn test_sample_barcode_segments() {
        let segments = vec![
            seg("AGCT".as_bytes(), SegmentType::Template),
            seg("GATA".as_bytes(), SegmentType::SampleBarcode),
            seg("CAC".as_bytes(), SegmentType::SampleBarcode),
            seg("GACCCC".as_bytes(), SegmentType::MolecularBarcode),
        ];
        let expected = vec![
            seg("GATA".as_bytes(), SegmentType::SampleBarcode),
            seg("CAC".as_bytes(), SegmentType::SampleBarcode),
        ];

        let read_set = read_set(segments);

        assert_eq!(expected, read_set.sample_barcode_segments().cloned().collect::<Vec<_>>());
    }
    #[test]
    fn test_molecular_barcode_segments() {
        let segments = vec![
            seg("AGCT".as_bytes(), SegmentType::Template),
            seg("GATA".as_bytes(), SegmentType::MolecularBarcode),
            seg("CAC".as_bytes(), SegmentType::MolecularBarcode),
            seg("GACCCC".as_bytes(), SegmentType::SampleBarcode),
        ];
        let expected = vec![
            seg("GATA".as_bytes(), SegmentType::MolecularBarcode),
            seg("CAC".as_bytes(), SegmentType::MolecularBarcode),
        ];

        let read_set = read_set(segments);

        assert_eq!(expected, read_set.molecular_barcode_segments().cloned().collect::<Vec<_>>());
    }

    #[test]
    fn test_combine_readsets() {
        let segments1 = vec![
            seg("A".as_bytes(), SegmentType::Template),
            seg("G".as_bytes(), SegmentType::Template),
            seg("C".as_bytes(), SegmentType::MolecularBarcode),
            seg("T".as_bytes(), SegmentType::SampleBarcode),
        ];
        let read_set1 = read_set_with_header(b"read1", segments1.clone());
        let segments2 = vec![
            seg("AA".as_bytes(), SegmentType::Template),
            seg("AG".as_bytes(), SegmentType::Template),
            seg("AC".as_bytes(), SegmentType::MolecularBarcode),
            seg("AT".as_bytes(), SegmentType::SampleBarcode),
        ];
        let read_set2 = read_set_with_header(b"read2", segments2.clone());
        let segments3 = vec![
            seg("AAA".as_bytes(), SegmentType::Template),
            seg("AAG".as_bytes(), SegmentType::Template),
            seg("AAC".as_bytes(), SegmentType::MolecularBarcode),
            seg("AAT".as_bytes(), SegmentType::SampleBarcode),
        ];
        let read_set3 = read_set_with_header(b"read3", segments3.clone());

        let mut expected_segments = Vec::new();
        expected_segments.extend(segments1);
        expected_segments.extend(segments2);
        expected_segments.extend(segments3);
        // Combining retains one header per source, in input order, so a later per-source lookup
        // (e.g. keep-original-read-names for R2) can recover the right header.
        let expected = ReadSet {
            headers: vec![b"read1".to_vec(), b"read2".to_vec(), b"read3".to_vec()],
            segments: expected_segments,
            raw_per_input: vec![],
            skip_reason: None,
        };

        let result = ReadSet::combine_readsets(vec![read_set1, read_set2, read_set3]);

        assert_eq!(result, expected);
        assert_eq!(result.header_for_source(0), b"read1");
        assert_eq!(result.header_for_source(1), b"read2");
        assert_eq!(result.header_for_source(2), b"read3");
    }

    #[test]
    #[should_panic(expected = "Cannot call combine readsets on an empty vec!")]
    fn test_combine_readsets_fails_on_empty_vector() {
        let _result = ReadSet::combine_readsets(Vec::new());
    }

    #[test]
    fn test_segments_for_source() {
        // Build a ReadSet with segments from two source indices in interleaved order so we can
        // verify segments_for_source filters by source AND preserves the original ordering.
        let segments = vec![
            seg_with_source(b"AAAA", SegmentType::SampleBarcode, 0),
            seg_with_source(b"GGGGG", SegmentType::Template, 0),
            seg_with_source(b"CC", SegmentType::MolecularBarcode, 1),
            seg_with_source(b"TTT", SegmentType::Template, 1),
            seg_with_source(b"NN", SegmentType::Skip, 0),
        ];
        let rs = read_set(segments);

        let b_t: HashSet<SegmentType> =
            [SegmentType::SampleBarcode, SegmentType::Template].into_iter().collect();
        let m_t: HashSet<SegmentType> =
            [SegmentType::MolecularBarcode, SegmentType::Template].into_iter().collect();
        let t_only: HashSet<SegmentType> = [SegmentType::Template].into_iter().collect();

        // Source 0 with {B,T}: pick segments at original positions 0 and 1 in order.
        let (seq, quals) = rs.segments_for_source(0, &b_t);
        assert_eq!(seq, b"AAAAGGGGG");
        assert_eq!(quals.len(), seq.len());

        // Source 1 with {M,T}: pick segments at original positions 2 and 3 in order.
        let (seq, _) = rs.segments_for_source(1, &m_t);
        assert_eq!(seq, b"CCTTT");

        // Source 0 with {T}: only the template segment.
        let (seq, _) = rs.segments_for_source(0, &t_only);
        assert_eq!(seq, b"GGGGG");

        // Source 0 with {M,T}: source 0 has no M, so just the T.
        let (seq, _) = rs.segments_for_source(0, &m_t);
        assert_eq!(seq, b"GGGGG");

        // Unused source index returns empty.
        let (seq, quals) = rs.segments_for_source(2, &b_t);
        assert!(seq.is_empty());
        assert!(quals.is_empty());
    }

    // ############################################################################################
    // Tests for per-sample read structures (CODEC-style stagger).
    // ############################################################################################

    /// Writes a per-sample-read-structures metadata TSV with columns:
    ///     sample_id, barcode, read_structure_1, ..., read_structure_<n>
    /// `samples` is a slice of (sample_id, barcode, [read_structure_per_input...]).
    fn metadata_file_with_rs(tmpdir: &TempDir, samples: &[(&str, &str, &[&str])]) -> PathBuf {
        let n_inputs = samples[0].2.len();
        let mut header = String::from("sample_id\tbarcode");
        for i in 1..=n_inputs {
            header.push('\t');
            header.push_str(&format!("read_structure_{i}"));
        }
        let mut lines = vec![header];
        for (sid, bc, structs) in samples {
            assert_eq!(structs.len(), n_inputs);
            let mut line = format!("{sid}\t{bc}");
            for s in structs.iter() {
                line.push('\t');
                line.push_str(s);
            }
            lines.push(line);
        }
        let path = tmpdir.path().join("metadata.tsv");
        Io::default().write_lines(&path, lines).unwrap();
        path
    }

    /// End-to-end CODEC stagger test:
    ///   Sample S1 (no stagger):  3M7B1S+T   barcode = GATTACA
    ///   Sample S2 (1bp stagger): 3M1S7B1S+T barcode = TTTTTTT
    /// One paired-end read per sample, each with a 50bp template after the constant `T`
    /// stagger anchor.  Verifies that:
    ///   1. Each read is assigned to the correct sample.
    ///   2. The output template starts at the right position (no template haircut for the
    ///      shorter-stagger sample).
    ///   3. The UMI (M) and sample-barcode (B) outputs reflect the per-sample structure.
    #[test]
    fn test_demux_codec_per_sample_read_structures_paired_end() {
        let tmp = TempDir::new().unwrap();

        // The global `--read-structures` is provided but should not be used for matched
        // samples (they have their own).  Its per-input segment signature can differ from the
        // per-sample structures; the unmatched output uses globals independently.
        let read_structures =
            vec![ReadStructure::from_str("+T").unwrap(), ReadStructure::from_str("+T").unwrap()];

        let s1_template_r1 = "A".repeat(50);
        let s1_template_r2 = "C".repeat(50);
        // Sample 1, no stagger: NNN GATTACA T <template>
        let s1_r1 = format!("AAA{}T{}", "GATTACA", s1_template_r1);
        let s1_r2 = format!("CCC{}T{}", "GATTACA", s1_template_r2);

        let s2_template_r1 = "G".repeat(50);
        let s2_template_r2 = "T".repeat(50);
        // Sample 2, 1bp stagger: NNN A TTTTTTT T <template>
        let s2_r1 = format!("GGGA{}T{}", "TTTTTTT", s2_template_r1);
        let s2_r2 = format!("TTTA{}T{}", "TTTTTTT", s2_template_r2);

        let r1_file = fastq_file(&tmp, "r1", "ex", &[s1_r1.as_str(), s2_r1.as_str()]);
        let r2_file = fastq_file(&tmp, "r2", "ex", &[s1_r2.as_str(), s2_r2.as_str()]);

        let sample_metadata = metadata_file_with_rs(
            &tmp,
            &[
                ("S1", "GATTACAGATTACA", &["3M7B1S+T", "3M7B1S+T"]),
                ("S2", "TTTTTTTTTTTTTT", &["3M1S7B1S+T", "3M1S7B1S+T"]),
            ],
        );

        let output_dir = tmp.path().join("output");
        let demux_inputs = Demux {
            inputs: vec![r1_file, r2_file],
            read_structures,
            sample_metadata,
            output_types: vec!['T', 'B', 'M'],
            output: output_dir.clone(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 5,
            compression_level: 5,
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();

        // Sample 1: per-sample RS = 3M7B1S+T for both inputs.  Expect 50bp template per read,
        // 7bp barcode per read, 3bp UMI per read.
        let s1_r1_reads = read_fastq(&output_dir.join("S1.R1.fq.gz"));
        let s1_r2_reads = read_fastq(&output_dir.join("S1.R2.fq.gz"));
        let s1_b1_reads = read_fastq(&output_dir.join("S1.I1.fq.gz"));
        let s1_b2_reads = read_fastq(&output_dir.join("S1.I2.fq.gz"));
        let s1_m1_reads = read_fastq(&output_dir.join("S1.U1.fq.gz"));
        let s1_m2_reads = read_fastq(&output_dir.join("S1.U2.fq.gz"));

        assert_eq!(s1_r1_reads.len(), 1, "S1 R1 should have one read");
        assert_eq!(s1_r1_reads[0].seq, s1_template_r1.as_bytes());
        assert_eq!(s1_r2_reads[0].seq, s1_template_r2.as_bytes());
        assert_eq!(s1_b1_reads[0].seq, b"GATTACA");
        assert_eq!(s1_b2_reads[0].seq, b"GATTACA");
        assert_eq!(s1_m1_reads[0].seq, b"AAA");
        assert_eq!(s1_m2_reads[0].seq, b"CCC");

        // Sample 2: per-sample RS = 3M1S7B1S+T for both inputs.
        let s2_r1_reads = read_fastq(&output_dir.join("S2.R1.fq.gz"));
        let s2_r2_reads = read_fastq(&output_dir.join("S2.R2.fq.gz"));
        let s2_b1_reads = read_fastq(&output_dir.join("S2.I1.fq.gz"));
        let s2_b2_reads = read_fastq(&output_dir.join("S2.I2.fq.gz"));
        let s2_m1_reads = read_fastq(&output_dir.join("S2.U1.fq.gz"));
        let s2_m2_reads = read_fastq(&output_dir.join("S2.U2.fq.gz"));

        assert_eq!(s2_r1_reads.len(), 1, "S2 R1 should have one read");
        assert_eq!(s2_r1_reads[0].seq, s2_template_r1.as_bytes());
        assert_eq!(s2_r2_reads[0].seq, s2_template_r2.as_bytes());
        assert_eq!(s2_b1_reads[0].seq, b"TTTTTTT");
        // Second input's barcode/UMI must also be resegmented per the staggered structure.
        assert_eq!(s2_b2_reads[0].seq, b"TTTTTTT");
        assert_eq!(s2_m1_reads[0].seq, b"GGG");
        assert_eq!(s2_m2_reads[0].seq, b"TTT");

        // The unmatched outputs (using globals `+T +T`) only have R*.fq.gz files.  All reads
        // matched, so they should be empty.
        let unmatched_r1 = read_fastq(&output_dir.join("unmatched.R1.fq.gz"));
        let unmatched_r2 = read_fastq(&output_dir.join("unmatched.R2.fq.gz"));
        assert!(unmatched_r1.is_empty(), "unmatched R1 should be empty");
        assert!(unmatched_r2.is_empty(), "unmatched R2 should be empty");
    }

    /// End-to-end single-input CODEC-style test with two samples at different staggers (the most
    /// common inline-index layout).  Exercises heterogeneous per-sample prefix lengths within a
    /// single input and verifies that the shorter-stagger sample's template is not haircut to the
    /// longer matching window.
    #[test]
    fn test_demux_codec_per_sample_single_input() {
        let tmp = TempDir::new().unwrap();
        // Global is only used for unmatched reads.
        let read_structures = vec![ReadStructure::from_str("+T").unwrap()];

        let s1_template = "A".repeat(30);
        let s2_template = "C".repeat(30);
        // S1, no stagger:  NNN GATTACA T <template>          (pre-template = 11)
        let s1_r1 = format!("AAA{}T{}", "GATTACA", s1_template);
        // S2, 1bp stagger: NNN A TTTTTTT T <template>        (pre-template = 12)
        let s2_r1 = format!("GGGA{}T{}", "TTTTTTT", s2_template);
        let r1_file = fastq_file(&tmp, "r1", "ex", &[s1_r1.as_str(), s2_r1.as_str()]);

        let sample_metadata = metadata_file_with_rs(
            &tmp,
            &[("S1", "GATTACA", &["3M7B1S+T"]), ("S2", "TTTTTTT", &["3M1S7B1S+T"])],
        );

        let output_dir = tmp.path().join("output");
        let demux_inputs = Demux {
            inputs: vec![r1_file],
            read_structures,
            sample_metadata,
            output_types: vec!['T', 'B', 'M'],
            output: output_dir.clone(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 5,
            compression_level: 5,
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();

        // S1 (shorter stagger): template must start at offset 11, not the 12-base window.
        let s1_r1_reads = read_fastq(&output_dir.join("S1.R1.fq.gz"));
        assert_eq!(s1_r1_reads.len(), 1, "S1 should have one read");
        assert_eq!(s1_r1_reads[0].seq, s1_template.as_bytes(), "S1 template must not be haircut");
        assert_eq!(read_fastq(&output_dir.join("S1.I1.fq.gz"))[0].seq, b"GATTACA");
        assert_eq!(read_fastq(&output_dir.join("S1.U1.fq.gz"))[0].seq, b"AAA");

        // S2 (longer stagger): template starts at offset 12.
        let s2_r1_reads = read_fastq(&output_dir.join("S2.R1.fq.gz"));
        assert_eq!(s2_r1_reads.len(), 1, "S2 should have one read");
        assert_eq!(s2_r1_reads[0].seq, s2_template.as_bytes());
        assert_eq!(read_fastq(&output_dir.join("S2.I1.fq.gz"))[0].seq, b"TTTTTTT");
        assert_eq!(read_fastq(&output_dir.join("S2.U1.fq.gz"))[0].seq, b"GGG");
    }

    /// Ensures that when no per-sample read structures are present, behaviour is identical
    /// to the legacy code path (uses global `--read-structures`).
    #[test]
    fn test_demux_falls_back_to_global_read_structures() {
        let tmp = TempDir::new().unwrap();
        let read_structures = vec![ReadStructure::from_str("17B100T").unwrap()];
        let s1_barcode = "AAAAAAAAGATTACAGA";
        let sample_metadata = metadata_file(
            &tmp,
            &[s1_barcode, "CCCCCCCCGATTACAGA", "GGGGGGGGGATTACAGA", "GGGGGGTTGATTACAGA"],
        );
        let input_files =
            vec![fastq_file(&tmp, "ex", "ex", &[&(s1_barcode.to_owned() + &"A".repeat(100))])];

        let output_dir = tmp.path().join("output");
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
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();

        let fq_reads = read_fastq(&output_dir.join("Sample0000.R1.fq.gz"));
        assert_eq!(fq_reads.len(), 1);
        assert_eq!(fq_reads[0].seq, "A".repeat(100).as_bytes());
    }

    /// Errors when per-sample `read_structure_<n>` columns disagree with the number of
    /// `--read-structures` entries.
    #[test]
    fn test_demux_errors_when_per_sample_input_count_mismatches_global() {
        let tmp = TempDir::new().unwrap();
        let r1_file =
            fastq_file(&tmp, "r1", "ex", &[format!("AAAGATTACAT{}", "A".repeat(50)).as_str()]);
        // Metadata declares 1 read_structure column, but globals declare 2 inputs.
        let sample_metadata = metadata_file_with_rs(&tmp, &[("S1", "GATTACA", &["3M7B1S+T"])]);

        let demux_inputs = Demux {
            inputs: vec![r1_file.clone(), r1_file],
            read_structures: vec![
                ReadStructure::from_str("3M7B+T").unwrap(),
                ReadStructure::from_str("3M7B+T").unwrap(),
            ],
            sample_metadata,
            output_types: vec!['T'],
            output: tmp.path().join("output"),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 5,
            compression_level: 5,
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        let err = demux_inputs.execute().unwrap_err();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("`read_structure_<n>` column(s)") && msg.contains("--read-structures"),
            "got: {msg}",
        );
    }

    /// Per-sample read structures may have a different `(T, B, M, C)` segment signature than
    /// the global `--read-structures` — different samples may produce different sets of
    /// output files, and the unmatched output uses the global structures independently.
    #[test]
    fn test_demux_per_sample_signature_may_differ_from_global() {
        let tmp = TempDir::new().unwrap();
        // Per-sample has (T=1, B=1, M=1, C=0); global has (T=1, B=1, M=0, C=0) — M counts differ.
        let r1 = format!("AAAGATTACAT{}", "A".repeat(50));
        let r1_file = fastq_file(&tmp, "r1", "ex", &[r1.as_str()]);
        let sample_metadata = metadata_file_with_rs(&tmp, &[("S1", "GATTACA", &["3M7B1S+T"])]);

        let demux_inputs = Demux {
            inputs: vec![r1_file],
            read_structures: vec![ReadStructure::from_str("7B+T").unwrap()],
            sample_metadata,
            output_types: vec!['T', 'M'],
            output: tmp.path().join("output").clone(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 5,
            compression_level: 5,
            skip_reasons: vec![],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        let output_dir = demux_inputs.output.clone();
        demux_inputs.execute().unwrap();
        // S1 produces R1 + U1 (matches per-sample 3M7B1S+T → 1T, 1M).
        let s1_r1 = read_fastq(&output_dir.join("S1.R1.fq.gz"));
        let s1_u1 = read_fastq(&output_dir.join("S1.U1.fq.gz"));
        assert_eq!(s1_r1.len(), 1);
        assert_eq!(s1_u1.len(), 1);
        assert_eq!(s1_u1[0].seq, b"AAA");
        // Globals (7B+T) have no M segment, so unmatched only writes R1; the file may not even
        // exist if no unmatched reads occurred.  In this test, the read matches S1, so there
        // are no unmatched records.
    }

    /// Builds a per-sample demux scenario in which a read clears the (loosened) per-sample
    /// prefix gate, matches no sample, and is too short for the global `--read-structures` when
    /// resegmented on the unmatched path.  `skip_reasons` is supplied by the caller.
    fn unmatched_short_read_demux(tmp: &TempDir, skip_reasons: Vec<SkipReason>) -> Demux {
        // Global is 8B+T per input; sample S1 has 3B+T per input (prefix 3).
        let read_structures = vec![
            ReadStructure::from_str("8B+T").unwrap(),
            ReadStructure::from_str("8B+T").unwrap(),
        ];
        let sample_metadata = metadata_file_with_rs(tmp, &[("S1", "AAAAAA", &["3B+T", "3B+T"])]);
        // 3-byte reads per input that don't match S1's "AAA"+"AAA" matching pattern.
        let r1_file = fastq_file(tmp, "r1", "ex", &["TTT"]);
        let r2_file = fastq_file(tmp, "r2", "ex", &["TTT"]);

        Demux {
            inputs: vec![r1_file, r2_file],
            read_structures,
            sample_metadata,
            output_types: vec!['T', 'B'],
            output: tmp.path().join("output"),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 5,
            compression_level: 5,
            skip_reasons,
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        }
    }

    /// In per-sample mode the iterator gates on the per-sample matching prefix, which can be
    /// shorter than the global `--read-structures` minimum.  A read that clears the prefix gate
    /// but matches no sample reaches the unmatched branch, where resegmentation against the
    /// global structure fails on length.  Without `--skip-reasons too-few-bases` this is fatal
    /// (consistent with the matched and legacy paths) and the error names the remedy.
    #[test]
    fn test_demux_per_sample_unmatched_short_read_aborts_without_flag() {
        let tmp = TempDir::new().unwrap();
        let demux_inputs = unmatched_short_read_demux(&tmp, vec![]);
        let err = demux_inputs.execute().unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("--skip-reasons too-few-bases"), "got: {msg}");
    }

    /// With `--skip-reasons too-few-bases`, the same too-short unmatched read is counted as a
    /// skip rather than aborting, and the run completes with no reads demultiplexed.
    #[test]
    fn test_demux_per_sample_unmatched_short_read_skipped_with_flag() {
        let tmp = TempDir::new().unwrap();
        let demux_inputs = unmatched_short_read_demux(&tmp, vec![SkipReason::TooFewBases]);
        let output_dir = demux_inputs.output.clone();
        demux_inputs.execute().unwrap();

        let metrics_path = output_dir.join("demux-metrics.txt");
        let metrics: Vec<DemuxMetric> = DelimFile::default().read_tsv(&metrics_path).unwrap();
        let demuxed: usize = metrics.iter().map(|m| m.templates).sum();
        assert_eq!(demuxed, 0, "no reads should be demuxed; the unmatched read is too short");
    }

    /// A read whose bases end exactly at the start of a variable-length template (`+T`) has zero
    /// template bases.  The underlying read-structure crate returns an empty slice (not an error)
    /// for the trailing `+T`, so without an explicit guard `resegment` would emit a zero-length
    /// FASTQ record.  Such a read must instead be treated as `TooShort`.
    #[test]
    fn test_resegment_zero_length_template_is_too_short() {
        let read_set = ReadSet {
            headers: vec![b"ex".to_vec()],
            segments: vec![],
            raw_per_input: vec![RawInput {
                bases: b"ACGTACGT".to_vec(),
                quals: b"FFFFFFFF".to_vec(),
            }],
            skip_reason: None,
        };
        // 8B consumes all 8 bases, leaving the `+T` template with zero bases.
        let read_structures = vec![ReadStructure::from_str("8B+T").unwrap()];
        let result = read_set.resegment(&read_structures);
        assert!(
            matches!(result, Err(ResegmentError::TooShort(_))),
            "expected TooShort for a read with no template bases, got {result:?}",
        );
    }

    /// A variable-length template with at least one base resegments normally (guards against the
    /// zero-length fix over-rejecting non-empty templates).
    #[test]
    fn test_resegment_one_base_template_is_ok() {
        let read_set = ReadSet {
            headers: vec![b"ex".to_vec()],
            segments: vec![],
            raw_per_input: vec![RawInput {
                bases: b"ACGTACGTT".to_vec(),
                quals: b"FFFFFFFFF".to_vec(),
            }],
            skip_reason: None,
        };
        // 8B consumes 8 bases, leaving a single `T` base for the template.
        let read_structures = vec![ReadStructure::from_str("8B+T").unwrap()];
        let resegmented = read_set.resegment(&read_structures).expect("should resegment");
        let template: Vec<u8> =
            resegmented.template_segments().flat_map(|s| s.seq.clone()).collect();
        assert_eq!(template, b"T");
    }

    /// End-to-end: a matched per-sample read whose bases end exactly at the template start must be
    /// counted as a `TooFewBases` skip, not written as a zero-length template record.
    #[test]
    fn test_demux_per_sample_zero_length_template_is_skipped() {
        let tmp = TempDir::new().unwrap();
        // Per-sample RS 8B+T (prefix 8); the read is exactly 8 bases (barcode only, no template).
        let sample_metadata = metadata_file_with_rs(&tmp, &[("S1", "ACGTACGT", &["8B+T"])]);
        let r1_file = fastq_file(&tmp, "r1", "ex", &["ACGTACGT"]);

        let output_dir = tmp.path().join("output");
        let demux_inputs = Demux {
            inputs: vec![r1_file],
            read_structures: vec![ReadStructure::from_str("8B+T").unwrap()],
            sample_metadata,
            output_types: vec!['T', 'B'],
            output: output_dir.clone(),
            unmatched_prefix: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: 5,
            compression_level: 5,
            // Tolerate the too-short read so the matched path skips rather than aborts.
            skip_reasons: vec![SkipReason::TooFewBases],
            template_types: vec!['T'],
            header_format: HeaderFormatKind::Illumina,
            umi_in_name: None,
        };
        demux_inputs.execute().unwrap();

        // No template record should have been written for S1 (no zero-length record).
        let s1_r1 = read_fastq(&output_dir.join("S1.R1.fq.gz"));
        assert!(s1_r1.is_empty(), "no zero-length template record should be written");

        let metrics_path = output_dir.join("demux-metrics.txt");
        let metrics: Vec<DemuxMetric> = DelimFile::default().read_tsv(&metrics_path).unwrap();
        let demuxed: usize = metrics.iter().map(|m| m.templates).sum();
        assert_eq!(demuxed, 0, "the boundary read must be skipped, not demuxed");
    }
}
