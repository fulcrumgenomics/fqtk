use crate::is_valid_iupac;

use anyhow::{Context, Result, anyhow, bail, ensure};
use fgoxide::io::Io;
use itertools::Itertools;
use read_structure::{ReadStructure, SegmentType};
use std::fmt::{self, Display};
use std::ops::Range;
use std::path::Path;
use std::str::FromStr;

const DEFAULT_FILE_DELIMETER: u8 = b'\t';
const SAMPLE_ID_HEADER: &str = "sample_id";
const BARCODE_HEADER: &str = "barcode";
const READ_STRUCTURE_PREFIX: &str = "read_structure_";

/// Struct for describing a single sample and metadata associated with that sample.
#[derive(Clone, Debug, PartialEq)]
pub struct Sample {
    /// ID of the sample or library
    pub sample_id: String,
    /// DNA barcode associated with the sample
    pub barcode: String,
    /// Optional per-sample read structures (one per input FASTQ).  When present, these
    /// override the global `--read-structures` for this sample, both for matching pattern
    /// construction and for output extraction.
    pub read_structures: Option<Vec<ReadStructure>>,
    /// Index of the sample in the [`SampleGroup`] object, used for syncing indices across
    /// different structs.
    pub(crate) ordinal: usize,
}

impl Display for Sample {
    /// Implements a nice format display for the [`Sample`] struct.
    /// E.g. A sample with ordinal 2, name test-sample, and barcode GATTACA would look like:
    /// Sample(0002) - { name: test-sample    barcode: GATTACA }
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Sample({:04}) - {{ name: {}\tbarcode: {} }}",
            self.ordinal, self.sample_id, self.barcode
        )
    }
}

impl Sample {
    /// Validates inputs to generate a [`Self`] struct and instantiates the struct if they are
    /// valid.
    /// # Panics
    ///   - Panics if sample name is empty string.
    ///   - Panics if barcode is empty string.
    ///   - Panics if barcode has bases other than A, C, G, T, U, R, Y, S, W, K, M, D, V, H, B, or
    ///     N/n/.
    #[must_use]
    pub fn new(ordinal: usize, name: String, barcode: String) -> Self {
        Self::with_read_structures(ordinal, name, barcode, None)
    }

    /// Like [`Sample::new`] but allows attaching per-sample read structures.
    #[must_use]
    pub fn with_read_structures(
        ordinal: usize,
        name: String,
        barcode: String,
        read_structures: Option<Vec<ReadStructure>>,
    ) -> Self {
        assert!(!name.is_empty(), "Sample name cannot be empty");
        assert!(!barcode.is_empty(), "Sample barcode cannot be empty");
        assert!(
            barcode.as_bytes().iter().all(|&b| is_valid_iupac(b)),
            "All sample barcode bases must be one of A, C, G, T, U, R, Y, S, W, K, M, D, V, H, B, N"
        );
        Self { sample_id: name, barcode, read_structures, ordinal }
    }

    /// Returns the header line expected by the metadata file deserializer.  Only the two required
    /// columns are reported; per-sample `read_structure_<n>` columns are optional.
    #[must_use]
    pub fn deserialize_header_line() -> String {
        format!("{SAMPLE_ID_HEADER}\t{BARCODE_HEADER}")
    }
}

/// Struct for storing information about multiple samples and for defining functions associated
/// with groups of [`Sample`]s, rather than individual structs.
#[derive(Clone, Debug, PartialEq)]
pub struct SampleGroup {
    /// A group of samples
    pub samples: Vec<Sample>,
}

impl Display for SampleGroup {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "SampleGroup {{")?;
        for sample in &self.samples {
            writeln!(f, "    {sample}")?;
        }
        writeln!(f, "}}")
    }
}

impl SampleGroup {
    /// Validates a group of [`Sample`]s and instantiates a [`Self`] struct if they are valid. Will
    /// clone the [`Sample`] structs and change the number on the `ordinal` field on those cloned
    /// to match the order in which they are stored in this [`Self`].
    ///
    /// Per-sample read structures may differ in their per-input `(T, B, M, C)` segment counts
    /// across samples (which lets each sample produce a different set of output files), but
    /// the total length of all sample-barcode (`B`) segments per sample must equal the
    /// `barcode` column for that sample (and barcodes are required to be the same length).
    ///
    /// # Errors
    ///   - Will error if no samples are provided.
    ///   - Will error if there are duplicate sample names provided.
    ///   - Will error if there are duplicate barcodes provided.
    ///   - Will error if barcodes don't all have the same length.
    ///   - Will error if any sample's per-sample sample-barcode (`B`) segment is not fixed length.
    ///   - Will error if any sample's per-sample sample-barcode (`B`) segment lengths don't
    ///     sum to that sample's `barcode` column length.
    pub fn from_samples(samples: &[Sample]) -> Result<Self> {
        ensure!(!samples.is_empty(), "Must provide one or more sample");

        ensure!(
            samples.iter().map(|s| &s.sample_id).all_unique(),
            "Each sample name must be unique, duplicate identified"
        );

        ensure!(
            samples.iter().map(|s| &s.barcode).all_unique(),
            "Each sample barcode must be unique, duplicate identified",
        );

        let first_barcode_length = samples[0].barcode.len();
        ensure!(
            samples.iter().map(|s| &s.barcode).all(|b| b.len() == first_barcode_length),
            "All barcodes must have the same length",
        );

        // Per-sample read structures (when present) must have B-segments whose lengths sum to
        // the sample's `barcode` column length.  Sample-barcode segments must be fixed-length.
        for sample in samples {
            let Some(rs) = sample.read_structures.as_ref() else { continue };
            let mut b_len: usize = 0;
            for seg in rs.iter().flat_map(|r| r.segments_by_type(SegmentType::SampleBarcode)) {
                let len = seg.length.ok_or_else(|| {
                    anyhow!(
                        "Sample {}: sample-barcode (B) segments in per-sample read structures \
                         must be fixed length",
                        sample.sample_id,
                    )
                })?;
                b_len += len;
            }
            ensure!(
                b_len == sample.barcode.len(),
                "Sample {}: total sample-barcode (B) length across per-sample read structures \
                 is {} but barcode column has {} bases",
                sample.sample_id,
                b_len,
                sample.barcode.len(),
            );
        }

        Ok(Self {
            samples: samples
                .iter()
                .enumerate()
                .map(|(ordinal, sample)| {
                    Sample::with_read_structures(
                        ordinal,
                        sample.sample_id.clone(),
                        sample.barcode.clone(),
                        sample.read_structures.clone(),
                    )
                })
                .collect(),
        })
    }

    /// Attempts to load a [`Self`] object from a tab-delimited file.  The file must have a header
    /// row containing at least the columns `sample_id` and `barcode`.  Optional additional columns
    /// `read_structure_1`, `read_structure_2`, ..., `read_structure_<N>` attach per-sample read
    /// structures to each sample (one per input FASTQ in the order given to `--inputs`).  When
    /// these columns are present:
    ///
    /// - `globals` (the `--read-structures` argument) must have exactly `N` entries, one per
    ///   input FASTQ.
    /// - A blank cell in `read_structure_<i>` falls back to `globals[i-1]` for that sample.
    /// - A row whose `read_structure_<n>` cells are all blank uses `globals` entirely (i.e. is
    ///   equivalent to omitting the per-sample columns for that sample only).
    ///
    /// When no `read_structure_<n>` columns are present, `globals` is unused and every sample
    /// uses the global structures.
    ///
    /// # Errors
    ///   - Will error if the file cannot be read.
    ///   - Will error if the header row is missing required columns.
    ///   - Will error if `read_structure_<n>` columns are non-contiguous or don't match
    ///     `globals.len()` (when present).
    ///   - Will error if any `read_structure_<n>` cell has a value that fails to parse.
    ///   - Will error if [`SampleGroup::from_samples`] rejects the parsed records (e.g. duplicate
    ///     names/barcodes, inconsistent barcode lengths, or barcode/`B`-segment length mismatch).
    pub fn from_file<P: AsRef<Path>>(path: P, globals: &[ReadStructure]) -> Result<SampleGroup> {
        let path = path.as_ref();
        let io = Io::default();
        let lines = io
            .read_lines(path)
            .with_context(|| format!("failed to read sample metadata file {path:?}"))?;
        let mut iter = lines.into_iter().filter(|l| !l.trim().is_empty());

        let header = iter.next().ok_or_else(|| {
            anyhow!(
                "sample metadata file {path:?} is empty (expected header line {})",
                Sample::deserialize_header_line(),
            )
        })?;
        // Strip a UTF-8 BOM and any trailing CR so files saved on Windows or with a BOM
        // (e.g. by Excel/Notepad) parse correctly.
        let header = strip_bom_and_cr(&header);
        let header_fields: Vec<&str> = header.split(DEFAULT_FILE_DELIMETER as char).collect();

        let sample_id_idx =
            header_fields.iter().position(|c| *c == SAMPLE_ID_HEADER).ok_or_else(|| {
                anyhow!("sample metadata header is missing column `{SAMPLE_ID_HEADER}`")
            })?;
        let barcode_idx =
            header_fields.iter().position(|c| *c == BARCODE_HEADER).ok_or_else(|| {
                anyhow!("sample metadata header is missing column `{BARCODE_HEADER}`")
            })?;

        // Discover read_structure_<n> columns in header order; require contiguous indexing
        // starting at 1 and a total count matching `globals.len()` (when present).
        let mut rs_columns: Vec<(usize, usize)> = Vec::new(); // (n, header_idx)
        for (idx, name) in header_fields.iter().enumerate() {
            if let Some(suffix) = name.strip_prefix(READ_STRUCTURE_PREFIX) {
                let n: usize = suffix
                    .parse()
                    .with_context(|| format!("metadata column `{name}` has non-integer suffix"))?;
                ensure!(n >= 1, "metadata column `{name}` must use 1-based indexing");
                rs_columns.push((n, idx));
            }
        }
        rs_columns.sort_by_key(|(n, _)| *n);
        for (i, (n, _)) in rs_columns.iter().enumerate() {
            ensure!(
                *n == i + 1,
                "per-sample read structure columns must be contiguous starting at \
                 `{READ_STRUCTURE_PREFIX}1` (found `{READ_STRUCTURE_PREFIX}{n}` at position {})",
                i + 1,
            );
        }
        if !rs_columns.is_empty() {
            ensure!(
                rs_columns.len() == globals.len(),
                "metadata file has {} `{READ_STRUCTURE_PREFIX}<n>` column(s) but \
                 `--read-structures` has {} entry/entries",
                rs_columns.len(),
                globals.len(),
            );
        }

        let mut samples: Vec<Sample> = Vec::new();
        for (line_no, line) in iter.enumerate() {
            let row_no = line_no + 2; // header is line 1
            let line = line.trim_end_matches('\r');
            let cols: Vec<&str> = line.split(DEFAULT_FILE_DELIMETER as char).collect();
            ensure!(
                cols.len() == header_fields.len(),
                "sample metadata row {row_no} has {} columns but header has {}",
                cols.len(),
                header_fields.len(),
            );
            let sample_id = cols[sample_id_idx].to_owned();
            let barcode = cols[barcode_idx].to_owned();
            let read_structures =
                parse_per_sample_read_structures(row_no, &cols, &rs_columns, globals)?;
            samples.push(Sample::with_read_structures(
                samples.len(),
                sample_id,
                barcode,
                read_structures,
            ));
        }

        if samples.is_empty() {
            bail!("sample metadata file {path:?} contained no sample rows");
        }
        Self::from_samples(&samples)
    }

    /// Returns true if this group has per-sample read structures (i.e. at least one sample
    /// carries custom read structures that override `--read-structures`).
    #[must_use]
    pub fn has_per_sample_read_structures(&self) -> bool {
        self.samples.iter().any(|s| s.read_structures.is_some())
    }

    /// Returns the per-input matching prefix lengths needed to demultiplex this sample group.
    /// For each input FASTQ index, this is the maximum across samples of the number of bases
    /// before the first Template segment or the end of the last sample-barcode (`B`) segment,
    /// whichever is larger, in that sample's read structure.  Falls back to the corresponding
    /// entry in `default_structures` for samples without per-sample read structures.
    ///
    /// # Errors
    ///   - Returns an error if any non-Template segment before the first Template segment, or
    ///     any sample-barcode segment, has a non-fixed length.
    ///   - Returns an error if a sample's per-sample read-structure count differs from
    ///     `default_structures.len()`.
    pub fn matching_prefix_lens(&self, default_structures: &[ReadStructure]) -> Result<Vec<usize>> {
        let n = default_structures.len();
        let mut maxes = vec![0usize; n];
        for sample in &self.samples {
            let rs_for_sample = sample.read_structures.as_deref().unwrap_or(default_structures);
            ensure!(
                rs_for_sample.len() == n,
                "sample {}: number of read structures ({}) does not match number of inputs ({})",
                sample.sample_id,
                rs_for_sample.len(),
                n,
            );
            for (i, rs) in rs_for_sample.iter().enumerate() {
                let plen = pre_template_fixed_len(rs).with_context(|| {
                    format!(
                        "per-sample demultiplexing requires a fixed-length matching window, so \
                         every segment before the template in {} must have a fixed length",
                        read_structure_source(sample, i)
                    )
                })?;
                let barcode_end = sample_barcode_spans(rs)
                    .with_context(|| {
                        format!("invalid read structure in {}", read_structure_source(sample, i))
                    })?
                    .last()
                    .map_or(0, |span| span.end);
                maxes[i] = maxes[i].max(plen).max(barcode_end);
            }
        }
        Ok(maxes)
    }

    /// Builds the matching pattern bytes for each sample.  The pattern is a fixed-length
    /// byte sequence (concatenated across all input FASTQs) of total length
    /// `prefix_lens.iter().sum()`.  For each input, B (sample-barcode) segment positions are
    /// filled from `sample.barcode` (taken in read order across inputs), and all other
    /// positions up to the prefix length (M/S/C/T segments and trailing bases) are filled with
    /// `N` (treated as a wildcard by the matcher).  `B` segments may appear anywhere in the
    /// read, including after a fixed-length template, since only the last segment of a read
    /// structure may be variable length and every `B` therefore has a fixed offset.
    ///
    /// # Errors
    ///   - Returns an error if any sample-barcode (`B`) segment has a non-fixed length.
    ///   - Returns an error if a sample-barcode (`B`) segment crosses the matching window
    ///     boundary, since silently truncating a barcode would corrupt matching.
    pub fn build_matching_patterns(
        &self,
        default_structures: &[ReadStructure],
        prefix_lens: &[usize],
    ) -> Result<Vec<Vec<u8>>> {
        ensure!(
            default_structures.len() == prefix_lens.len(),
            "expected one prefix length per input FASTQ"
        );
        let total_len: usize = prefix_lens.iter().sum();
        let mut patterns = Vec::with_capacity(self.samples.len());
        for sample in &self.samples {
            let rs_for_sample = sample.read_structures.as_deref().unwrap_or(default_structures);
            ensure!(
                rs_for_sample.len() == prefix_lens.len(),
                "sample {}: number of read structures ({}) does not match number of inputs ({})",
                sample.sample_id,
                rs_for_sample.len(),
                prefix_lens.len(),
            );

            // The effective read structures' total sample-barcode (B) length must equal the
            // barcode column length.  This also validates global-only samples in mixed mode, whose
            // effective structure is the `--read-structures` fallback (a case `from_samples` cannot
            // check because it has no access to the global structures).
            let expected_barcode_len: usize = rs_for_sample
                .iter()
                .flat_map(|rs| rs.segments_by_type(SegmentType::SampleBarcode))
                .map(|seg| seg.length.unwrap_or(0))
                .sum();
            ensure!(
                expected_barcode_len == sample.barcode.len(),
                "sample {}: {}read structure(s) declare {} sample-barcode (B) base(s) but the \
                 barcode column has {} base(s)",
                sample.sample_id,
                if sample.read_structures.is_none() {
                    "(using --read-structures fallback) "
                } else {
                    ""
                },
                expected_barcode_len,
                sample.barcode.len(),
            );

            let mut pattern = Vec::with_capacity(total_len);
            let mut barcode_cursor = 0usize;
            let barcode_bytes = sample.barcode.as_bytes();
            for (i, (rs, &prefix_len)) in rs_for_sample.iter().zip(prefix_lens).enumerate() {
                let mut input_pattern = vec![b'N'; prefix_len];
                let spans = sample_barcode_spans(rs).with_context(|| {
                    format!("invalid read structure in {}", read_structure_source(sample, i))
                })?;
                for span in spans {
                    let end = span.end;
                    ensure!(
                        end <= prefix_len,
                        "sample {}: a sample-barcode segment in read structure {rs} crosses the \
                         matching window boundary (segment ends at base {end}, but the window is \
                         {prefix_len} bases)",
                        sample.sample_id,
                    );
                    let len = span.len();
                    ensure!(
                        barcode_cursor + len <= barcode_bytes.len(),
                        "sample {}: barcode ({} bases) is shorter than the total \
                         sample-barcode length required by its read structure",
                        sample.sample_id,
                        barcode_bytes.len(),
                    );
                    input_pattern[span]
                        .copy_from_slice(&barcode_bytes[barcode_cursor..barcode_cursor + len]);
                    barcode_cursor += len;
                }
                pattern.extend_from_slice(&input_pattern);
            }
            ensure!(
                barcode_cursor == sample.barcode.len(),
                "sample {}: only consumed {} of {} barcode bases when building matching pattern",
                sample.sample_id,
                barcode_cursor,
                sample.barcode.len(),
            );
            patterns.push(pattern);
        }
        Ok(patterns)
    }

    /// Builds the [`MatchingWindow`] used to assign reads to samples in per-sample read
    /// structure mode, keeping only the positions where some sample's pattern is not a no-call.
    /// The other positions can never change a match, and dropping them keeps template bases out
    /// of the barcode matcher's cache keys.  Errors as [`Self::build_matching_patterns`] does.
    pub fn matching_window(&self, default_structures: &[ReadStructure]) -> Result<MatchingWindow> {
        let prefix_lens = self.matching_prefix_lens(default_structures)?;
        let full_patterns = self.build_matching_patterns(default_structures, &prefix_lens)?;
        let coordinates: Vec<(usize, usize)> = prefix_lens
            .iter()
            .enumerate()
            .flat_map(|(input, &len)| (0..len).map(move |offset| (input, offset)))
            .collect();
        let mut kept: Vec<usize> = (0..coordinates.len())
            .filter(|&j| full_patterns.iter().any(|p| !crate::byte_is_nocall(p[j])))
            .collect();
        if kept.is_empty() {
            // No position is constrained; keep them all so the patterns are never empty.
            kept = (0..coordinates.len()).collect();
        }
        Ok(MatchingWindow {
            positions: kept.iter().map(|&j| coordinates[j]).collect(),
            patterns: full_patterns
                .iter()
                .map(|pattern| kept.iter().map(|&j| pattern[j]).collect())
                .collect(),
        })
    }
}

/// The read positions and per-sample expected bases used to assign reads to samples in
/// per-sample read structure mode.  Built by [`SampleGroup::matching_window`].
#[derive(Clone, Debug, PartialEq)]
pub struct MatchingWindow {
    /// The `(input index, offset in the read)` of each compared base, in read order.
    positions: Vec<(usize, usize)>,
    /// For each sample (in sample order), its expected base at each of `positions`.  A no-call
    /// (`N`) is a wildcard.
    patterns: Vec<Vec<u8>>,
}

impl MatchingWindow {
    /// Returns each sample's expected bases at the window positions, in sample order.
    #[must_use]
    pub fn patterns(&self) -> &[Vec<u8>] {
        &self.patterns
    }

    /// Returns the bases of `reads` (one per input FASTQ, in input order) at the window
    /// positions.  A position past the end of its read is returned as `N`, so a short read still
    /// matches the samples whose barcode bases it covers.
    #[must_use]
    pub fn observed_bases<R: AsRef<[u8]>>(&self, reads: &[R]) -> Vec<u8> {
        self.positions
            .iter()
            .map(|&(input, offset)| reads[input].as_ref().get(offset).copied().unwrap_or(b'N'))
            .collect()
    }
}

/// Resolves the per-sample read structure cells for one row of the metadata file.  Returns
/// `None` if `rs_columns` is empty or every cell is blank (sample uses globals entirely);
/// otherwise returns `Some(vec)` with each blank cell replaced by `globals[i]`.
///
/// Sample-barcode (`B`) segments in any parsed cell must be fixed length: a variable-length
/// `+B` is rejected here so the error surfaces as a normal `Result` rather than panicking
/// later during `from_samples` validation.
fn parse_per_sample_read_structures(
    row_no: usize,
    cols: &[&str],
    rs_columns: &[(usize, usize)],
    globals: &[ReadStructure],
) -> Result<Option<Vec<ReadStructure>>> {
    if rs_columns.is_empty() {
        return Ok(None);
    }
    let mut entries: Vec<Option<ReadStructure>> = Vec::with_capacity(rs_columns.len());
    for (n, idx) in rs_columns {
        let raw = cols[*idx].trim();
        if raw.is_empty() {
            entries.push(None);
        } else {
            let rs = ReadStructure::from_str(raw).with_context(|| {
                format!(
                    "sample metadata row {row_no} column `{READ_STRUCTURE_PREFIX}{n}` has \
                     invalid read structure `{raw}`",
                )
            })?;
            for seg in rs.segments_by_type(SegmentType::SampleBarcode) {
                ensure!(
                    seg.length.is_some(),
                    "sample metadata row {row_no} column `{READ_STRUCTURE_PREFIX}{n}`: \
                     sample-barcode segment {seg} must be fixed length (variable-length `+B` \
                     is not supported in per-sample read structures)",
                );
            }
            entries.push(Some(rs));
        }
    }
    if entries.iter().all(Option::is_none) {
        return Ok(None);
    }
    let resolved: Vec<ReadStructure> = entries
        .into_iter()
        .enumerate()
        .map(|(i, e)| e.unwrap_or_else(|| globals[i].clone()))
        .collect();
    Ok(Some(resolved))
}

/// Strips a leading UTF-8 byte-order mark and any trailing carriage return from a line so that
/// files saved with a BOM or CRLF line endings parse correctly.
fn strip_bom_and_cr(s: &str) -> &str {
    s.strip_prefix('\u{FEFF}').unwrap_or(s).trim_end_matches('\r')
}

/// Returns the read positions covered by each sample-barcode (`B`) segment in `rs`, in order.
/// Only the last segment of a read structure may be variable length, so every fixed-length `B`
/// segment has a fixed offset.  Errors if a `B` segment is variable length.
fn sample_barcode_spans(rs: &ReadStructure) -> Result<Vec<Range<usize>>> {
    let mut spans = Vec::new();
    let mut offset = 0;
    for seg in rs.iter() {
        if seg.kind == SegmentType::SampleBarcode {
            let len = seg.length.ok_or_else(|| {
                anyhow!(
                    "sample-barcode segment {seg} in read structure {rs} must have a fixed length"
                )
            })?;
            spans.push(offset..offset + len);
        }
        match seg.length {
            Some(len) => offset += len,
            None => break,
        }
    }
    Ok(spans)
}

/// Describes where a sample's read structure for input `input_idx` (0-based) came from, for
/// error messages: its own `read_structure_<n>` cell or the `--read-structures` fallback.
fn read_structure_source(sample: &Sample, input_idx: usize) -> String {
    if sample.read_structures.is_some() {
        format!("sample {}'s `{READ_STRUCTURE_PREFIX}{}`", sample.sample_id, input_idx + 1)
    } else {
        format!(
            "the `--read-structures` fallback for sample {} (input {})",
            sample.sample_id,
            input_idx + 1,
        )
    }
}

/// Returns the offset of the first Template segment in a read structure.  If the read structure
/// has no Template segments, returns the sum of all (fixed-length) segment lengths.
///
/// # Errors
///   - Returns an error if a non-Template segment lacks a fixed length.
fn pre_template_fixed_len(rs: &ReadStructure) -> Result<usize> {
    let mut len = 0;
    for seg in rs.iter() {
        if seg.kind == SegmentType::Template {
            return Ok(len);
        }
        len += seg.length.ok_or_else(|| {
            anyhow!("non-template segment {seg} in read structure {rs} must have a fixed length")
        })?;
    }
    Ok(len)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fgoxide::io::Io;
    use std::str::FromStr;
    use tempfile::TempDir;

    // ############################################################################################
    // Test [`SampleGroup::from_file`] - Expected to pass
    // ############################################################################################
    #[test]
    fn test_reading_from_tsv_file() {
        let lines = vec![
            Sample::deserialize_header_line(),
            "sample1\tGATTACA".to_owned(),
            "sample2\tCATGCTA".to_owned(),
        ];
        let tempdir = TempDir::new().unwrap();
        let f1 = tempdir.path().join("sample_metadata.tsv");

        let io = Io::default();
        io.write_lines(&f1, &lines).unwrap();
        let samples_metadata = SampleGroup::from_file(&f1, &[]).unwrap();

        assert!(samples_metadata.samples[0].sample_id == "sample1");
        assert!(samples_metadata.samples[1].sample_id == "sample2");
        assert!(samples_metadata.samples[0].barcode == "GATTACA");
        assert!(samples_metadata.samples[1].barcode == "CATGCTA");
        assert!(!samples_metadata.has_per_sample_read_structures());
    }

    #[test]
    fn test_reading_from_file_with_empty_lines_at_end() {
        let lines = vec![
            Sample::deserialize_header_line(),
            "sample1\tGATTACA".to_owned(),
            "sample2\tCATGCTA".to_owned(),
            String::new(),
            String::new(),
        ];
        let tempdir = TempDir::new().unwrap();
        let f1 = tempdir.path().join("sample_metadata.tsv");

        let io = Io::default();
        io.write_lines(&f1, &lines).unwrap();
        let samples_metadata = SampleGroup::from_file(&f1, &[]).unwrap();

        assert!(samples_metadata.samples[0].sample_id == "sample1");
        assert!(samples_metadata.samples[1].sample_id == "sample2");
        assert!(samples_metadata.samples[0].barcode == "GATTACA");
        assert!(samples_metadata.samples[1].barcode == "CATGCTA");
    }

    #[test]
    fn test_new_sample_non_agct_bases_in_barcode_allowed() {
        let name = "s_1_example_name".to_owned();
        let barcode = "GATTANN".to_owned();
        let ordinal = 0;
        let _sample = Sample::new(ordinal, name, barcode);
    }

    #[test]
    fn test_tsv_file_delim_error() {
        let lines: Vec<String> = ["sample_id,barcode", "sample1,GATTACA", "sample2,CATGCTA"]
            .iter()
            .map(|&s| s.into())
            .collect();
        let tempdir = TempDir::new().unwrap();
        let f1 = tempdir.path().join("sample_metadata.tsv");

        let io = Io::default();
        io.write_lines(&f1, &lines).unwrap();
        let err = SampleGroup::from_file(&f1, &[]).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("missing column `sample_id`"), "got: {msg}");
    }

    // ############################################################################################
    // Test [`SampleGroup::from_file`] - Expected to error or panic
    // ############################################################################################
    #[test]
    fn test_reading_from_file_with_no_header() {
        let lines = vec!["sample1\tGATTACA", "sample2\tCATGCTA"];
        let tempdir = TempDir::new().unwrap();
        let f1 = tempdir.path().join("sample_metadata.tsv");

        let io = Io::default();
        io.write_lines(&f1, &lines).unwrap();
        let err = SampleGroup::from_file(&f1, &[]).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("missing column `sample_id`"), "got: {msg}");
    }

    #[test]
    fn test_reading_header_only_file() {
        let lines = vec![Sample::deserialize_header_line()];
        let tempdir = TempDir::new().unwrap();
        let f1 = tempdir.path().join("sample_metadata.tsv");

        let io = Io::default();
        io.write_lines(&f1, &lines).unwrap();
        let err = SampleGroup::from_file(&f1, &[]).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("contained no sample rows"), "got: {msg}");
    }

    #[test]
    fn test_reading_empty_file() {
        let lines = vec![""];
        let tempdir = TempDir::new().unwrap();
        let f1 = tempdir.path().join("sample_metadata.tsv");

        let io = Io::default();
        io.write_lines(&f1, &lines).unwrap();
        let err = SampleGroup::from_file(&f1, &[]).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("is empty") || msg.contains("missing column"), "got: {msg}");
    }

    #[test]
    fn test_reading_from_file_with_duplicate_barcodes_errors() {
        // A malformed user TSV (duplicate barcodes) must return a clean error, not panic.
        let lines = vec!["sample_id\tbarcode", "sample1\tGATTACA", "sample2\tGATTACA"];
        let tempdir = TempDir::new().unwrap();
        let f1 = tempdir.path().join("sample_metadata.tsv");

        let io = Io::default();
        io.write_lines(&f1, &lines).unwrap();
        let err = SampleGroup::from_file(&f1, &[]).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("Each sample barcode must be unique"), "got: {msg}");
    }

    #[test]
    fn test_reading_non_existent_file() {
        let tempdir = TempDir::new().unwrap();
        let f1 = tempdir.path().join("sample_metadata.tsv");
        let err = SampleGroup::from_file(&f1, &[]).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("failed to read sample metadata file"), "got: {msg}");
    }

    // ############################################################################################
    // Test [`Sample::new`] - Expected to pass
    // ############################################################################################
    #[test]
    fn test_new_sample_success() {
        let name = "s_1_example_name".to_owned();
        let barcode = "GATTACA".to_owned();
        let ordinal = 0;
        let sample = Sample::new(ordinal, name.clone(), barcode.clone());
        assert_eq!(
            Sample { sample_id: name, barcode, read_structures: None, ordinal },
            sample,
            "Sample differed from expectation"
        );
    }

    // ############################################################################################
    // Test [`Sample::new`] - Expected to panic
    // ############################################################################################
    #[test]
    #[should_panic(expected = "Sample name cannot be empty")]
    fn test_new_sample_fail1_empty_sample_name() {
        let name = String::new();
        let barcode = "GATTACA".to_owned();
        let ordinal = 0;
        let _sample = Sample::new(ordinal, name, barcode);
    }

    #[test]
    #[should_panic(expected = "Sample barcode cannot be empty")]
    fn test_new_sample_fail2_empty_barcode() {
        let name = "s_1_example_name".to_owned();
        let barcode = String::new();
        let ordinal = 0;
        let _sample = Sample::new(ordinal, name, barcode);
    }

    // ############################################################################################
    // Test [`SampleGroup::from_samples`] - expected to pass
    // ############################################################################################
    #[test]
    fn test_from_samples_sample_group_pass1_single_sample() {
        let sample1 = Sample::new(0, "sample_1".to_owned(), "GATTACA".to_owned());
        let samples_vec = vec![sample1.clone()];
        let sample_group = SampleGroup::from_samples(&samples_vec).unwrap();

        assert_eq!(sample_group, SampleGroup { samples: vec![sample1] });
    }

    #[test]
    fn test_from_samples_sample_group_pass2_multi_unique_samples() {
        let sample1 = Sample::new(0, "sample_1".to_owned(), "GATTACA".to_owned());
        let sample2 = Sample::new(1, "sample_2".to_owned(), "CATGGAT".to_owned());
        let samples_vec = vec![sample1.clone(), sample2.clone()];
        let sample_group = SampleGroup::from_samples(&samples_vec).unwrap();

        assert_eq!(sample_group, SampleGroup { samples: vec![sample1, sample2] });
    }

    #[test]
    fn test_from_samples_sample_group_pass3_ordinal_values_will_be_changed_by_new() {
        let sample1 = Sample::new(0, "sample_1".to_owned(), "GATTACA".to_owned());
        let sample2_before = Sample::new(2, "sample_2".to_owned(), "CATGGAT".to_owned());
        let sample2_after = Sample::new(1, "sample_2".to_owned(), "CATGGAT".to_owned());
        let samples_vec = vec![sample1.clone(), sample2_before];
        let sample_group = SampleGroup::from_samples(&samples_vec).unwrap();

        assert_eq!(sample_group, SampleGroup { samples: vec![sample1, sample2_after] });
    }

    // ############################################################################################
    // Test [`SampleGroup::from_samples`] - expected to error
    // ############################################################################################
    #[test]
    fn test_from_samples_sample_group_fail1_no_samples() {
        let samples = vec![];
        let err = SampleGroup::from_samples(&samples).unwrap_err();
        assert!(err.to_string().contains("Must provide one or more sample"), "got: {err:#}");
    }

    #[test]
    fn test_from_samples_sample_group_fail2_duplicate_sample_names() {
        let samples = vec![
            Sample::new(0, "sample_1".to_owned(), "GATTACA".to_owned()),
            Sample::new(0, "sample_1".to_owned(), "CATGGAT".to_owned()),
        ];
        let err = SampleGroup::from_samples(&samples).unwrap_err();
        assert!(
            err.to_string().contains("Each sample name must be unique, duplicate identified"),
            "got: {err:#}",
        );
    }

    #[test]
    fn test_from_samples_sample_group_fail3_duplicate_barcodes() {
        let samples = vec![
            Sample::new(0, "sample_1".to_owned(), "GATTACA".to_owned()),
            Sample::new(0, "sample_2".to_owned(), "GATTACA".to_owned()),
        ];
        let err = SampleGroup::from_samples(&samples).unwrap_err();
        assert!(
            err.to_string().contains("Each sample barcode must be unique, duplicate identified"),
            "got: {err:#}",
        );
    }

    #[test]
    fn test_from_samples_sample_group_fail4_barcodes_of_different_lengths() {
        let samples = vec![
            Sample::new(0, "sample_1".to_owned(), "GATTACA".to_owned()),
            Sample::new(0, "sample_2".to_owned(), "CATGGA".to_owned()),
        ];
        let err = SampleGroup::from_samples(&samples).unwrap_err();
        assert!(err.to_string().contains("All barcodes must have the same length"), "got: {err:#}");
    }

    // ############################################################################################
    // Tests for per-sample read structures.
    // ############################################################################################
    fn make_rs(s: &str) -> ReadStructure {
        ReadStructure::from_str(s).unwrap()
    }

    fn sample_with_rs(name: &str, barcode: &str, structures: &[&str]) -> Sample {
        let rs = structures.iter().map(|s| make_rs(s)).collect();
        Sample::with_read_structures(0, name.to_owned(), barcode.to_owned(), Some(rs))
    }

    fn write_metadata(tempdir: &TempDir, lines: &[String]) -> std::path::PathBuf {
        let f1 = tempdir.path().join("metadata.tsv");
        Io::default().write_lines(&f1, lines).unwrap();
        f1
    }

    #[test]
    fn test_per_sample_read_structures_round_trip_via_metadata() {
        // Two-input case: each input contributes 7B, so the barcode column carries 14 bases.
        let lines = vec![
            "sample_id\tbarcode\tread_structure_1\tread_structure_2".to_owned(),
            "S1\tGATTACAACGTACG\t3M7B1S+T\t3M7B1S+T".to_owned(),
            "S2\tGGGGGGGTTTTTTT\t3M1S7B1S+T\t3M1S7B1S+T".to_owned(),
        ];
        let tempdir = TempDir::new().unwrap();
        let f1 = write_metadata(&tempdir, &lines);
        let globals = vec![make_rs("3M9B+T"), make_rs("9B+T")];
        let group = SampleGroup::from_file(&f1, &globals).unwrap();
        assert!(group.has_per_sample_read_structures());
        let s1_rs = group.samples[0].read_structures.as_ref().unwrap();
        assert_eq!(s1_rs.len(), 2);
        let s2_rs = group.samples[1].read_structures.as_ref().unwrap();
        assert_eq!(s2_rs.len(), 2);
    }

    /// Per-sample read structures may have differing per-input `(T, B, M, C)` segment counts
    /// across samples, as long as each sample's B-segment lengths sum to its `barcode` column.
    #[test]
    fn test_per_sample_read_structures_signatures_may_differ_across_samples() {
        // S1 has one B-segment of length 14; S2 has two B-segments (7+7) of total length 14.
        let s1 = sample_with_rs("S1", "GATTACAACGTACG", &["3M14B+T"]);
        let s2 = sample_with_rs("S2", "TTTTTTTGGGGGGG", &["3M7B7B+T"]);
        let group = SampleGroup::from_samples(&[s1, s2]).unwrap();
        assert!(group.has_per_sample_read_structures());
    }

    /// Mixed mode: some samples have per-sample read structures, others fall back to globals.
    #[test]
    fn test_per_sample_read_structures_mixed_with_global_only_samples() {
        let s1 = sample_with_rs("S1", "GATTACA", &["3M7B1S+T"]);
        let s2 = Sample::new(0, "S2".to_owned(), "CCCCCCC".to_owned());
        let group = SampleGroup::from_samples(&[s1, s2]).unwrap();
        assert!(group.has_per_sample_read_structures());
        assert!(group.samples[0].read_structures.is_some());
        assert!(group.samples[1].read_structures.is_none());
    }

    #[test]
    fn test_per_sample_read_structures_barcode_length_mismatch() {
        // S1 declares 7B + 7B = 14 bases of barcode but provides only 7 in the column.
        let s1 = sample_with_rs("S1", "GATTACA", &["3M7B1S+T", "3M7B1S+T"]);
        let err = SampleGroup::from_samples(&[s1]).unwrap_err();
        assert!(err.to_string().contains("barcode column has"), "got: {err:#}");
    }

    #[test]
    fn test_matching_prefix_lens_uses_max_across_samples() {
        // S1 prefix per input: 3+7+1 = 11; S2 prefix: 3+1+7+1 = 12
        let s1 = sample_with_rs("S1", "GATTACA", &["3M7B1S+T"]);
        let s2 = sample_with_rs("S2", "GGGGGGG", &["3M1S7B1S+T"]);
        let group = SampleGroup::from_samples(&[s1, s2]).unwrap();
        let defaults = vec![make_rs("3M9B+T")];
        let lens = group.matching_prefix_lens(&defaults).unwrap();
        assert_eq!(lens, vec![12]);
    }

    #[test]
    fn test_build_matching_patterns_codec_two_samples() {
        let s1 = sample_with_rs("S1", "GATTACA", &["3M7B1S+T"]);
        let s2 = sample_with_rs("S2", "GGGGGGG", &["3M1S7B1S+T"]);
        let group = SampleGroup::from_samples(&[s1, s2]).unwrap();
        let defaults = vec![make_rs("3M9B+T")];
        let lens = group.matching_prefix_lens(&defaults).unwrap();
        let patterns = group.build_matching_patterns(&defaults, &lens).unwrap();
        // Total length is 12 for each pattern.
        // S1: NNN + GATTACA + N (S=skip wildcarded) + N (padding to 12)
        assert_eq!(patterns[0], b"NNNGATTACANN");
        // S2: NNN + N (S=stagger wildcarded) + GGGGGGG + N (S=trailing wildcarded)
        assert_eq!(patterns[1], b"NNNNGGGGGGGN");
    }

    #[test]
    fn test_build_matching_patterns_falls_back_to_defaults_when_no_per_sample() {
        let s1 = Sample::new(0, "S1".to_owned(), "GATTACA".to_owned());
        let group = SampleGroup::from_samples(&[s1]).unwrap();
        let defaults = vec![make_rs("3M7B1S+T")];
        let lens = group.matching_prefix_lens(&defaults).unwrap();
        assert_eq!(lens, vec![11]);
        let patterns = group.build_matching_patterns(&defaults, &lens).unwrap();
        assert_eq!(patterns[0], b"NNNGATTACAN");
    }

    #[test]
    fn test_build_matching_patterns_dual_input_concatenated() {
        // Two inputs, each with its own staggered structure; barcode column concatenates the
        // two B-segments left-to-right across inputs.
        let s1 = sample_with_rs("S1", "GATTACAACGTACG", &["3M7B1S+T", "7B+T"]);
        let s2 = sample_with_rs("S2", "GGGGGGGTTTTTTT", &["3M1S7B1S+T", "1S7B+T"]);
        let group = SampleGroup::from_samples(&[s1, s2]).unwrap();
        let defaults = vec![make_rs("3M9B+T"), make_rs("9B+T")];
        let lens = group.matching_prefix_lens(&defaults).unwrap();
        // Input 1: max(11, 12) = 12; Input 2: max(7, 8) = 8.
        assert_eq!(lens, vec![12, 8]);
        let patterns = group.build_matching_patterns(&defaults, &lens).unwrap();
        // Total pattern length: 12 + 8 = 20.
        let mut expected_s1 = b"NNNGATTACANN".to_vec();
        expected_s1.extend_from_slice(b"ACGTACGN");
        assert_eq!(patterns[0], expected_s1);
        let mut expected_s2 = b"NNNNGGGGGGGN".to_vec();
        expected_s2.extend_from_slice(b"NTTTTTTT");
        assert_eq!(patterns[1], expected_s2);
    }

    /// A sample-barcode segment crossing the matching window is rejected with a clear error
    /// rather than silently truncated.
    #[test]
    fn test_build_matching_patterns_b_segment_crossing_window_errors() {
        let s1 = sample_with_rs("S1", "GATTACA", &["3M7B1S+T"]);
        let group = SampleGroup::from_samples(&[s1]).unwrap();
        let defaults = vec![make_rs("3M7B1S+T")];
        let lens = vec![5usize]; // forced narrow window for the test
        let err = group.build_matching_patterns(&defaults, &lens).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("crosses the matching window boundary"), "got: {msg}");
    }

    /// A sample-barcode segment after the template that crosses a (forced) narrow window is
    /// rejected too, naming where the segment ends.
    #[test]
    fn test_build_matching_patterns_post_template_b_crossing_window_errors() {
        let s1 = sample_with_rs("S1", "ACGTCCCCCCCC", &["4B10T8B"]);
        let group = SampleGroup::from_samples(&[s1]).unwrap();
        let defaults = vec![make_rs("4B10T8B")];
        let err = group.build_matching_patterns(&defaults, &[20]).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("segment ends at base 22, but the window is 20"), "{msg}");
    }

    /// A `SampleGroup` constructed directly (bypassing `from_samples`' barcode-length
    /// validation) whose barcode is shorter than its read structure's sample-barcode
    /// segments surfaces a clean error naming the declared vs. provided barcode lengths.
    #[test]
    fn test_build_matching_patterns_barcode_shorter_than_b_segments_errors() {
        // `3M7B1S+T` requires 7 barcode bases, but the barcode column only has 3.
        let s1 = sample_with_rs("S1", "GAT", &["3M7B1S+T"]);
        let group = SampleGroup { samples: vec![s1] };
        let defaults = vec![make_rs("3M7B1S+T")];
        let lens = vec![11usize];
        let err = group.build_matching_patterns(&defaults, &lens).unwrap_err();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("declare 7 sample-barcode (B) base(s) but the barcode column has 3"),
            "got: {msg}",
        );
    }

    /// Sample-barcode (B) segments after a fixed-length template extend the matching window to
    /// the end of the last B, and their bases are matched in read order (#107).
    #[test]
    fn test_build_matching_patterns_barcode_after_template() {
        // 4B10T8B: barcode segments straddle the template.
        let s1 = sample_with_rs("S1", "ACGTCCCCCCCC", &["4B10T8B"]);
        let group = SampleGroup::from_samples(&[s1]).unwrap();
        let defaults = vec![make_rs("4B10T8B")];
        let lens = group.matching_prefix_lens(&defaults).unwrap();
        assert_eq!(lens, vec![22]);
        let patterns = group.build_matching_patterns(&defaults, &lens).unwrap();
        assert_eq!(patterns[0], b"ACGTNNNNNNNNNNCCCCCCCC".to_vec());
    }

    /// Barcode bases are taken in read order across B segments separated by a UMI, and the
    /// window stops at the end of the last B even when other segments follow it.
    #[test]
    fn test_build_matching_patterns_barcodes_split_by_umi_after_template() {
        // R2 layout from #106 (scaled down): template, index, UMI, index, then trailing skip.
        let s1 = sample_with_rs("S1", "AAAACCCC", &["6T", "6T4B3M4B2S"]);
        let group = SampleGroup::from_samples(&[s1]).unwrap();
        let defaults = vec![make_rs("6T"), make_rs("6T4B3M4B2S")];
        let lens = group.matching_prefix_lens(&defaults).unwrap();
        assert_eq!(lens, vec![0, 17]);
        let patterns = group.build_matching_patterns(&defaults, &lens).unwrap();
        assert_eq!(patterns[0], b"NNNNNNAAAANNNCCCC".to_vec());
    }

    /// A variable-length sample-barcode (B) segment has no fixed end, so a global fallback
    /// structure ending in one is rejected in per-sample mode, naming the sample that uses it.
    #[test]
    fn test_matching_prefix_lens_variable_length_barcode_errors() {
        let s1 = sample_with_rs("S1", "ACGT", &["4B+T"]);
        let s2 = Sample::new(1, "S2".to_owned(), "TTTT".to_owned());
        let group = SampleGroup { samples: vec![s1, s2] };
        let defaults = vec![make_rs("4T+B")];
        let err = group.matching_prefix_lens(&defaults).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("the `--read-structures` fallback for sample S2 (input 1)"), "{msg}");
        assert!(msg.contains("sample-barcode segment +B"), "{msg}");
    }

    /// The matching window keeps only positions some sample's barcode occupies, across inputs
    /// and on both sides of a template, and each sample's pattern has wildcards at the others'.
    #[test]
    fn test_matching_window_keeps_only_barcode_positions() {
        let s1 = sample_with_rs("S1", "AAAACCCC", &["6T", "6T4B3M4B"]);
        let s2 = sample_with_rs("S2", "GGGGTTTT", &["4B+T", "4B+T"]);
        let group = SampleGroup::from_samples(&[s1, s2]).unwrap();
        let defaults = vec![make_rs("6T"), make_rs("6T4B3M4B")];
        let window = group.matching_window(&defaults).unwrap();
        assert_eq!(
            window.patterns(),
            [b"NNNNNNNNAAAACCCC".to_vec(), b"GGGGTTTTNNNNNNNN".to_vec()].as_slice()
        );
        assert_eq!(
            window.observed_bases(&[b"ACGTAC".as_slice(), b"TTTTTTAAAAGTACCCC".as_slice()]),
            b"ACGTTTTTAAAACCCC".to_vec()
        );
    }

    /// Positions past the end of a short read are observed as `N`, so the read can still be
    /// compared against every sample instead of failing a length check.
    #[test]
    fn test_matching_window_pads_short_reads_with_no_calls() {
        let s1 = sample_with_rs("S1", "AAAACCCC", &["6T", "6T4B3M4B"]);
        let s2 = sample_with_rs("S2", "GGGGTTTT", &["4B+T", "4B+T"]);
        let group = SampleGroup::from_samples(&[s1, s2]).unwrap();
        let defaults = vec![make_rs("6T"), make_rs("6T4B3M4B")];
        let window = group.matching_window(&defaults).unwrap();
        assert_eq!(
            window.observed_bases(&[b"GGGGA".as_slice(), b"TTTTTT".as_slice()]),
            b"GGGGTTTTNNNNNNNN".to_vec()
        );
    }

    /// If no position is constrained (every barcode base is a no-call), every position is kept
    /// so the patterns are never empty.
    #[test]
    fn test_matching_window_all_no_call_barcodes_keeps_every_position() {
        let s1 = sample_with_rs("S1", "NNnN", &["2M4B+T"]);
        let group = SampleGroup::from_samples(&[s1]).unwrap();
        let window = group.matching_window(&[make_rs("2M4B+T")]).unwrap();
        assert_eq!(window.patterns(), [b"NNNNnN".to_vec()].as_slice());
    }

    /// A global-only sample (mixed mode) whose barcode length doesn't match the global
    /// `--read-structures` B-length is rejected with a message pointing at the fallback.
    #[test]
    fn test_build_matching_patterns_global_only_barcode_mismatch_errors() {
        // S1 (per-sample 6B) and S2 (global-only) both have 6-base barcodes, but the global
        // structure declares 7 B bases, so S2's effective structure doesn't match its barcode.
        let s1 = sample_with_rs("S1", "GATTAC", &["3M6B1S+T"]);
        let s2 = Sample::new(0, "S2".to_owned(), "CCCCCC".to_owned());
        let group = SampleGroup::from_samples(&[s1, s2]).unwrap();
        let defaults = vec![make_rs("3M7B1S+T")];
        let lens = group.matching_prefix_lens(&defaults).unwrap();
        let err = group.build_matching_patterns(&defaults, &lens).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("--read-structures fallback") && msg.contains("S2"), "got: {msg}",);
    }

    /// In per-sample mode every input's pre-template segments must be fixed length; a variable
    /// pre-template segment yields a message explaining the fixed-window requirement.
    #[test]
    fn test_matching_prefix_lens_variable_pre_template_errors() {
        let s1 = sample_with_rs("S1", "ACGTACGT", &["8B+M"]);
        let group = SampleGroup::from_samples(&[s1]).unwrap();
        let defaults = vec![make_rs("8B+M")];
        let err = group.matching_prefix_lens(&defaults).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("fixed-length matching window"), "got: {msg}");
    }

    // ############################################################################################
    // Per-cell fallback tests for `from_file`.
    // ############################################################################################
    /// A blank `read_structure_<n>` cell falls back to `globals[n-1]` for that sample.
    #[test]
    fn test_per_cell_fallback_uses_globals_for_blank_cells() {
        let lines = vec![
            "sample_id\tbarcode\tread_structure_1\tread_structure_2".to_owned(),
            // S1 overrides only input 1; input 2 is blank → uses globals[1].
            "S1\tGATTACAGGGGGGG\t3M7B1S+T\t".to_owned(),
            // S2 overrides only input 2; input 1 is blank → uses globals[0].
            "S2\tCCCCCCCAAAAAAA\t\t1S7B+T".to_owned(),
        ];
        let tempdir = TempDir::new().unwrap();
        let f1 = write_metadata(&tempdir, &lines);
        let globals = vec![make_rs("3M7B+T"), make_rs("7B+T")];
        let group = SampleGroup::from_file(&f1, &globals).unwrap();
        let s1_rs = group.samples[0].read_structures.as_ref().unwrap();
        assert_eq!(s1_rs.len(), 2);
        assert_eq!(s1_rs[0].to_string(), "3M7B1S+T");
        assert_eq!(s1_rs[1].to_string(), "7B+T");
        let s2_rs = group.samples[1].read_structures.as_ref().unwrap();
        assert_eq!(s2_rs[0].to_string(), "3M7B+T");
        assert_eq!(s2_rs[1].to_string(), "1S7B+T");
    }

    /// A row whose `read_structure_<n>` cells are all blank uses globals entirely (the
    /// sample's stored `read_structures` is `None`).
    #[test]
    fn test_per_cell_all_blank_row_falls_back_to_globals_entirely() {
        let lines = vec![
            "sample_id\tbarcode\tread_structure_1\tread_structure_2".to_owned(),
            "S1\tGATTACAGGGGGGG\t3M7B1S+T\t3M7B1S+T".to_owned(),
            "S2\tCCCCCCCAAAAAAA\t\t".to_owned(),
        ];
        let tempdir = TempDir::new().unwrap();
        let f1 = write_metadata(&tempdir, &lines);
        let globals = vec![make_rs("3M7B+T"), make_rs("3M7B+T")];
        let group = SampleGroup::from_file(&f1, &globals).unwrap();
        assert!(group.samples[0].read_structures.is_some());
        assert!(group.samples[1].read_structures.is_none());
    }

    /// `read_structure_<n>` column count must match `globals.len()` when columns are present.
    #[test]
    fn test_per_sample_column_count_must_match_globals() {
        let lines = vec![
            "sample_id\tbarcode\tread_structure_1".to_owned(),
            "S1\tGATTACA\t3M7B1S+T".to_owned(),
        ];
        let tempdir = TempDir::new().unwrap();
        let f1 = write_metadata(&tempdir, &lines);
        let globals = vec![make_rs("3M7B+T"), make_rs("100T")];
        let err = SampleGroup::from_file(&f1, &globals).unwrap_err();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("`read_structure_<n>` column(s)") && msg.contains("--read-structures"),
            "got: {msg}",
        );
    }

    /// A variable-length sample-barcode (`+B`) in a per-sample column must surface as a
    /// normal `Result` error rather than panic during downstream validation.
    #[test]
    fn test_per_sample_variable_length_b_segment_errors() {
        let lines =
            vec!["sample_id\tbarcode\tread_structure_1".to_owned(), "S1\tGATTACA\t3M+B".to_owned()];
        let tempdir = TempDir::new().unwrap();
        let f1 = write_metadata(&tempdir, &lines);
        let globals = vec![make_rs("3M7B+T")];
        let err = SampleGroup::from_file(&f1, &globals).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("must be fixed length"), "got: {msg}");
    }

    /// A UTF-8 BOM on the first header field should not cause a confusing
    /// "missing column `sample_id`" error.
    #[test]
    fn test_header_with_utf8_bom_is_handled() {
        let lines = vec![
            format!("\u{FEFF}{}", Sample::deserialize_header_line()),
            "sample1\tGATTACA".to_owned(),
        ];
        let tempdir = TempDir::new().unwrap();
        let f1 = write_metadata(&tempdir, &lines);
        let group = SampleGroup::from_file(&f1, &[]).unwrap();
        assert_eq!(group.samples[0].sample_id, "sample1");
        assert_eq!(group.samples[0].barcode, "GATTACA");
    }

    /// CRLF line endings on a row should not cause the trailing `\r` to bleed into the
    /// barcode column and break downstream validation.
    #[test]
    fn test_rows_with_crlf_endings_are_handled() {
        let header = format!("{}\r", Sample::deserialize_header_line());
        let lines = vec![header, "sample1\tGATTACA\r".to_owned(), "sample2\tCATGCTA\r".to_owned()];
        let tempdir = TempDir::new().unwrap();
        let f1 = write_metadata(&tempdir, &lines);
        let group = SampleGroup::from_file(&f1, &[]).unwrap();
        assert_eq!(group.samples[0].barcode, "GATTACA");
        assert_eq!(group.samples[1].barcode, "CATGCTA");
    }

    /// `matching_prefix_lens` must error (not panic) when a sample's per-sample read
    /// structure count differs from `default_structures.len()`.
    #[test]
    fn test_matching_prefix_lens_errors_on_rs_count_mismatch() {
        // Sample declares two read structures; defaults provide only one.
        let s1 = sample_with_rs("S1", "GATTACAGGGGGGG", &["3M7B1S+T", "7B+T"]);
        let group = SampleGroup::from_samples(&[s1]).unwrap();
        let defaults = vec![make_rs("3M7B+T")];
        let err = group.matching_prefix_lens(&defaults).unwrap_err();
        let msg = format!("{err:#}");
        assert!(
            msg.contains("number of read structures") && msg.contains("number of inputs"),
            "got: {msg}",
        );
    }

    /// Non-contiguous `read_structure_<n>` columns must error.
    #[test]
    fn test_per_sample_columns_must_be_contiguous() {
        let lines = vec![
            "sample_id\tbarcode\tread_structure_1\tread_structure_3".to_owned(),
            "S1\tGATTACA\t3M7B1S+T\t1S7B+T".to_owned(),
        ];
        let tempdir = TempDir::new().unwrap();
        let f1 = write_metadata(&tempdir, &lines);
        let globals = vec![make_rs("3M7B+T"), make_rs("7B+T")];
        let err = SampleGroup::from_file(&f1, &globals).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("contiguous"), "got: {msg}");
    }

    /// A non-integer suffix on a `read_structure_<n>` column must error.
    #[test]
    fn test_per_sample_columns_must_have_integer_suffix() {
        let lines = vec![
            "sample_id\tbarcode\tread_structure_abc".to_owned(),
            "S1\tGATTACA\t3M7B1S+T".to_owned(),
        ];
        let tempdir = TempDir::new().unwrap();
        let f1 = write_metadata(&tempdir, &lines);
        let globals = vec![make_rs("3M7B+T")];
        let err = SampleGroup::from_file(&f1, &globals).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("non-integer suffix"), "got: {msg}");
    }

    /// A zero-indexed `read_structure_0` column must error since indexing is 1-based.
    #[test]
    fn test_per_sample_columns_must_be_one_indexed() {
        let lines = vec![
            "sample_id\tbarcode\tread_structure_0".to_owned(),
            "S1\tGATTACA\t3M7B1S+T".to_owned(),
        ];
        let tempdir = TempDir::new().unwrap();
        let f1 = write_metadata(&tempdir, &lines);
        let globals = vec![make_rs("3M7B+T")];
        let err = SampleGroup::from_file(&f1, &globals).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("1-based indexing"), "got: {msg}");
    }
}
