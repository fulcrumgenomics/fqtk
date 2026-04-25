use crate::is_valid_iupac;

use anyhow::{Context, Result, anyhow, bail, ensure};
use fgoxide::io::Io;
use itertools::Itertools;
use read_structure::{ReadStructure, SegmentType};
use std::fmt::{self, Display};
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
    /// # Panics
    ///   - Will panic if sample metadata sheet is improperly formatted.
    ///   - Will panic if there are duplicate sample names provided.
    ///   - Will panic if there are duplicate barcodes provided.
    ///   - Will panic if barcodes don't all have the same length.
    ///   - Will panic if any sample's per-sample sample-barcode (`B`) segment lengths don't
    ///     sum to that sample's `barcode` column length.
    #[must_use]
    pub fn from_samples(samples: &[Sample]) -> Self {
        assert!(!samples.is_empty(), "Must provide one or more sample");

        assert!(
            samples.iter().map(|s| &s.sample_id).all_unique(),
            "Each sample name must be unique, duplicate identified"
        );

        assert!(
            samples.iter().map(|s| &s.barcode).all_unique(),
            "Each sample barcode must be unique, duplicate identified",
        );

        let first_barcode_length = samples[0].barcode.len();
        assert!(
            samples.iter().map(|s| &s.barcode).all(|b| b.len() == first_barcode_length),
            "All barcodes must have the same length",
        );

        // Per-sample read structures (when present) must have B-segments whose lengths sum to
        // the sample's `barcode` column length.  Sample-barcode segments must be fixed-length.
        for sample in samples {
            let Some(rs) = sample.read_structures.as_ref() else { continue };
            let b_len: usize = rs
                .iter()
                .flat_map(|r| r.segments_by_type(SegmentType::SampleBarcode))
                .map(|seg| {
                    seg.length.expect(
                        "sample barcode segments in per-sample read structures must be fixed \
                         length",
                    )
                })
                .sum();
            assert!(
                b_len == sample.barcode.len(),
                "Sample {}: total sample-barcode (B) length across per-sample read structures \
                 is {} but barcode column has {} bases",
                sample.sample_id,
                b_len,
                sample.barcode.len(),
            );
        }

        Self {
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
        }
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
    /// # Panics
    ///   - Will panic if [`SampleGroup::from_samples`] panics on the parsed records.
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
        Ok(Self::from_samples(&samples))
    }

    /// Returns true if this group has per-sample read structures (i.e. at least one sample
    /// carries custom read structures that override `--read-structures`).
    #[must_use]
    pub fn has_per_sample_read_structures(&self) -> bool {
        self.samples.iter().any(|s| s.read_structures.is_some())
    }

    /// Returns the per-input matching prefix lengths needed to demultiplex this sample group.
    /// For each input FASTQ index, this is the maximum across samples of the number of bases
    /// before the first Template segment in that sample's read structure.  Falls back to the
    /// corresponding entry in `default_structures` for samples without per-sample read
    /// structures.
    ///
    /// # Errors
    ///   - Returns an error if any non-Template segment in the matching window has a
    ///     non-fixed length.
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
                let plen = pre_template_fixed_len(rs)
                    .with_context(|| format!("sample {} input {}", sample.sample_id, i + 1))?;
                if plen > maxes[i] {
                    maxes[i] = plen;
                }
            }
        }
        Ok(maxes)
    }

    /// Builds the matching pattern bytes for each sample.  The pattern is a fixed-length
    /// byte sequence (concatenated across all input FASTQs) of total length
    /// `prefix_lens.iter().sum()`.  For each input, the sample's read structure is walked
    /// position-by-position; B (sample-barcode) segment positions are filled from
    /// `sample.barcode` (taken in order across inputs), and M/S/C/T segment positions and any
    /// trailing positions up to the prefix length are filled with `N` (treated as a wildcard
    /// by the matcher).
    ///
    /// # Errors
    ///   - Returns an error if any sample's read structure has a non-fixed-length non-template
    ///     segment in the matching window.
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
            let mut pattern = Vec::with_capacity(total_len);
            let mut barcode_cursor = 0usize;
            let barcode_bytes = sample.barcode.as_bytes();
            for (rs, &prefix_len) in rs_for_sample.iter().zip(prefix_lens) {
                let mut filled = 0usize;
                for seg in rs.iter() {
                    if seg.kind == SegmentType::Template {
                        break;
                    }
                    let len = seg.length.ok_or_else(|| {
                        anyhow!(
                            "sample {}: non-template segment {seg} in read structure must have a \
                             fixed length",
                            sample.sample_id,
                        )
                    })?;
                    let remaining = prefix_len - filled;
                    let take = len.min(remaining);
                    if seg.kind == SegmentType::SampleBarcode {
                        ensure!(
                            take == len,
                            "sample {}: sample-barcode segment {seg} crosses the matching \
                             window boundary (segment length {}, but only {} bases remain in \
                             the {}-base window)",
                            sample.sample_id,
                            len,
                            remaining,
                            prefix_len,
                        );
                        pattern.extend_from_slice(
                            &barcode_bytes[barcode_cursor..barcode_cursor + len],
                        );
                        barcode_cursor += len;
                    } else {
                        pattern.extend(std::iter::repeat_n(b'N', take));
                    }
                    filled += take;
                    if take < len {
                        break;
                    }
                }
                if filled < prefix_len {
                    pattern.extend(std::iter::repeat_n(b'N', prefix_len - filled));
                }
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
        let sample_group = SampleGroup::from_samples(&samples_vec);

        assert_eq!(sample_group, SampleGroup { samples: vec![sample1] });
    }

    #[test]
    fn test_from_samples_sample_group_pass2_multi_unique_samples() {
        let sample1 = Sample::new(0, "sample_1".to_owned(), "GATTACA".to_owned());
        let sample2 = Sample::new(1, "sample_2".to_owned(), "CATGGAT".to_owned());
        let samples_vec = vec![sample1.clone(), sample2.clone()];
        let sample_group = SampleGroup::from_samples(&samples_vec);

        assert_eq!(sample_group, SampleGroup { samples: vec![sample1, sample2] });
    }

    #[test]
    fn test_from_samples_sample_group_pass3_ordinal_values_will_be_changed_by_new() {
        let sample1 = Sample::new(0, "sample_1".to_owned(), "GATTACA".to_owned());
        let sample2_before = Sample::new(2, "sample_2".to_owned(), "CATGGAT".to_owned());
        let sample2_after = Sample::new(1, "sample_2".to_owned(), "CATGGAT".to_owned());
        let samples_vec = vec![sample1.clone(), sample2_before];
        let sample_group = SampleGroup::from_samples(&samples_vec);

        assert_eq!(sample_group, SampleGroup { samples: vec![sample1, sample2_after] });
    }

    // ############################################################################################
    // Test [`SampleGroup::from_samples`] - expected to panic
    // ############################################################################################
    #[test]
    #[should_panic(expected = "Must provide one or more sample")]
    fn test_from_samples_sample_group_fail1_no_samples() {
        let samples = vec![];
        let _sample_group = SampleGroup::from_samples(&samples);
    }

    #[test]
    #[should_panic(expected = "Each sample name must be unique, duplicate identified")]
    fn test_from_samples_sample_group_fail2_duplicate_barcodes() {
        let samples = vec![
            Sample::new(0, "sample_1".to_owned(), "GATTACA".to_owned()),
            Sample::new(0, "sample_1".to_owned(), "CATGGAT".to_owned()),
        ];
        let _sample_group = SampleGroup::from_samples(&samples);
    }

    #[test]
    #[should_panic(expected = "Each sample barcode must be unique, duplicate identified")]
    fn test_from_samples_sample_group_fail3_duplicate_sample_names() {
        let samples = vec![
            Sample::new(0, "sample_1".to_owned(), "GATTACA".to_owned()),
            Sample::new(0, "sample_2".to_owned(), "GATTACA".to_owned()),
        ];
        let _sample_group = SampleGroup::from_samples(&samples);
    }

    #[test]
    #[should_panic(expected = "All barcodes must have the same length")]
    fn test_from_samples_sample_group_fail4_barcodes_of_different_lengths() {
        let samples = vec![
            Sample::new(0, "sample_1".to_owned(), "GATTACA".to_owned()),
            Sample::new(0, "sample_2".to_owned(), "CATGGA".to_owned()),
        ];
        let _sample_group = SampleGroup::from_samples(&samples);
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
        let group = SampleGroup::from_samples(&[s1, s2]);
        assert!(group.has_per_sample_read_structures());
    }

    /// Mixed mode: some samples have per-sample read structures, others fall back to globals.
    #[test]
    fn test_per_sample_read_structures_mixed_with_global_only_samples() {
        let s1 = sample_with_rs("S1", "GATTACA", &["3M7B1S+T"]);
        let s2 = Sample::new(0, "S2".to_owned(), "CCCCCCC".to_owned());
        let group = SampleGroup::from_samples(&[s1, s2]);
        assert!(group.has_per_sample_read_structures());
        assert!(group.samples[0].read_structures.is_some());
        assert!(group.samples[1].read_structures.is_none());
    }

    #[test]
    #[should_panic(expected = "barcode column has")]
    fn test_per_sample_read_structures_barcode_length_mismatch() {
        // S1 declares 7B + 7B = 14 bases of barcode but provides only 7 in the column.
        let s1 = sample_with_rs("S1", "GATTACA", &["3M7B1S+T", "3M7B1S+T"]);
        let _ = SampleGroup::from_samples(&[s1]);
    }

    #[test]
    fn test_matching_prefix_lens_uses_max_across_samples() {
        // S1 prefix per input: 3+7+1 = 11; S2 prefix: 3+1+7+1 = 12
        let s1 = sample_with_rs("S1", "GATTACA", &["3M7B1S+T"]);
        let s2 = sample_with_rs("S2", "GGGGGGG", &["3M1S7B1S+T"]);
        let group = SampleGroup::from_samples(&[s1, s2]);
        let defaults = vec![make_rs("3M9B+T")];
        let lens = group.matching_prefix_lens(&defaults).unwrap();
        assert_eq!(lens, vec![12]);
    }

    #[test]
    fn test_build_matching_patterns_codec_two_samples() {
        let s1 = sample_with_rs("S1", "GATTACA", &["3M7B1S+T"]);
        let s2 = sample_with_rs("S2", "GGGGGGG", &["3M1S7B1S+T"]);
        let group = SampleGroup::from_samples(&[s1, s2]);
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
        let group = SampleGroup::from_samples(&[s1]);
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
        let group = SampleGroup::from_samples(&[s1, s2]);
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
        let group = SampleGroup::from_samples(&[s1]);
        let defaults = vec![make_rs("3M7B1S+T")];
        let lens = vec![5usize]; // forced narrow window for the test
        let err = group.build_matching_patterns(&defaults, &lens).unwrap_err();
        let msg = format!("{err:#}");
        assert!(msg.contains("crosses the matching window boundary"), "got: {msg}");
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
        let group = SampleGroup::from_samples(&[s1]);
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
