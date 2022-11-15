use anyhow::Result;
use bstr::BString;
use fgoxide::io::DelimFile;
use itertools::Itertools;
use serde::Deserialize;
use std::fmt::{self, Display};
use std::path::Path;

/// Checks whether a provided byte is an A, G, C, or T.
fn is_valid_base(byte: u8) -> bool {
    byte == b'A' || byte == b'C' || byte == b'G' || byte == b'T'
}

/// Struct for describing a single sample and metadata associated with that sample.
#[derive(Clone, Deserialize, Debug)]
pub struct Sample {
    /// name of the sample
    pub name: String,
    /// DNA barcode associated with the sample
    pub barcode: String,
    /// index of the sample in the [`SampleGroup`] object, used for syncing indices across
    /// different structs
    #[serde(skip_deserializing)]
    ordinal: usize,
}

impl Display for Sample {
    /// Implements a nice format display for the [`Sample`] struct.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Sample({:04}) - {{ name: {}\tbarcode: {} }}",
            self.ordinal, self.name, self.barcode
        )
    }
}

impl Sample {
    /// Validates inputs to generate a [`Self`] struct and instantiates the struct if they are
    /// valid.
    /// # Panics
    ///   - Panics if sample name is empty string.
    ///   - Panics if barcode is empty string.
    ///   - Panics if barcode has bases other than A, C, G, or T.
    #[must_use]
    pub fn new(ordinal: usize, name: String, barcode: String) -> Self {
        assert!(!name.is_empty(), "Sample name cannot be empty");
        assert!(!barcode.is_empty(), "Sample barcode cannot be empty");
        assert!(
            barcode.as_bytes().iter().all(|&b| is_valid_base(b)),
            "All sample barcode bases must be one of A, C, G, or T"
        );
        Self { name, barcode, ordinal }
    }
}

/// Struct for storing information about multiple samples and for defining functions associated
/// with groups of [`Sample`]s, rather than individual structs.
#[derive(Clone, Debug)]
#[allow(clippy::module_name_repetitions)]
pub struct SampleGroup {
    /// A group of samples
    pub samples: Vec<Sample>,
}

impl Display for SampleGroup {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "SampleGroup {{")?;
        for sample in &self.samples {
            writeln!(f, "    {}", sample)?;
        }
        writeln!(f, "}}")
    }
}

impl SampleGroup {
    /// Validates a group of [`Sample`]s and instantiates a [`Self`] struct if they are
    /// valid. Will clone the [`Sample`] structs and change the number on the `ordinal` field on
    /// those cloneto match the order in which they are stored in this [`Self`]
    /// # Panics
    ///   - Will panic if sample metadata sheet is improperly formatted
    ///   - Will panic if a different number of names and barcodes are provided
    ///   - Will panic if each
    #[must_use]
    pub fn from_samples(samples: &[Sample]) -> Self {
        // Validate that we have at least one name
        assert!(!samples.is_empty(), "Must provide one or more sample");

        // Validate that all the sample names are unique
        assert!(
            samples.iter().map(|s| &s.name).all_unique(),
            "Each sample name must be unique, duplicate identified"
        );

        // Convert barcodes to BString
        let bstr_barcodes: Vec<BString> =
            samples.iter().map(|b| BString::from(b.barcode.as_bytes())).collect();

        // Validate that the barcodes are all unique
        assert!(
            bstr_barcodes.iter().all_unique(),
            "Each sample barcode must be unique, duplicate identified",
        );

        let first_barcode_length = bstr_barcodes[0].len();
        assert!(
            bstr_barcodes.iter().all(|b| b.len() == first_barcode_length),
            "All barcodes must have the same length",
        );

        Self {
            samples: samples
                .iter()
                .enumerate()
                .map(|(ordinal, sample)| {
                    Sample::new(ordinal, sample.name.clone(), sample.barcode.clone())
                })
                .collect(),
        }
    }

    /// Attempts to load a [`Self`] object from a file. File should be delimeted with
    /// `delimiter`, should have a header with `name` and `barcode` fields present.
    /// # Errors
    ///   - Will error if file cannot be read, either due to not the file not existing or due to
    ///     the format being different from the format expected.
    /// # Panics
    ///   - Will panic if sample metadata sheet is improperly formatted
    ///   - Will panic if a different number of names and barcodes are provided
    ///   - Will panic if each
    pub fn from_file<P: AsRef<Path>>(
        path: &P,
        delimiter: u8,
    ) -> Result<SampleGroup, fgoxide::FgError> {
        let reader = DelimFile::default();
        Ok(Self::from_samples(&reader.read(path, delimiter, false)?))
    }
}

#[cfg(test)]
mod tests {
    use core::panic;

    use super::*;
    use fgoxide::{self, io::Io};
    use tempfile::TempDir;

    // ############################################################################################
    // Test [`SampleGroup::from_file`] - Expected to pass
    // ############################################################################################
    #[test]
    fn test_reading_from_tsv_file() {
        let lines = vec!["name\tbarcode", "sample1\tGATTACA", "sample2\tCATGCTA"];
        let tempdir = TempDir::new().unwrap();
        let f1 = tempdir.path().join("sample_metadata.tsv");

        let io = Io::default();
        io.write_lines(&f1, &lines).unwrap();
        let samples_metadata = SampleGroup::from_file(&f1, b'\t').unwrap();

        assert!(samples_metadata.samples[0].name == *"sample1");
        assert!(samples_metadata.samples[1].name == *"sample2");
        assert!(samples_metadata.samples[0].barcode == "GATTACA");
        assert!(samples_metadata.samples[1].barcode == "CATGCTA");
    }
    #[test]
    fn test_reading_from_csv_file() {
        let lines = vec!["name,barcode", "sample1,GATTACA", "sample2,CATGCTA"];
        let tempdir = TempDir::new().unwrap();
        let f1 = tempdir.path().join("sample_metadata.tsv");

        let io = Io::default();
        io.write_lines(&f1, &lines).unwrap();
        let samples_metadata = SampleGroup::from_file(&f1, b',').unwrap();

        assert!(samples_metadata.samples[0].name == *"sample1");
        assert!(samples_metadata.samples[1].name == *"sample2");
        assert!(samples_metadata.samples[0].barcode == "GATTACA");
        assert!(samples_metadata.samples[1].barcode == "CATGCTA");
    }

    #[test]
    fn test_reading_from_file_with_empty_lines_at_end() {
        let lines = vec!["name\tbarcode", "sample1\tGATTACA", "sample2\tCATGCTA", "", ""];
        let tempdir = TempDir::new().unwrap();
        let f1 = tempdir.path().join("sample_metadata.tsv");

        let io = Io::default();
        io.write_lines(&f1, &lines).unwrap();
        let samples_metadata = SampleGroup::from_file(&f1, b'\t').unwrap();

        assert!(samples_metadata.samples[0].name == *"sample1");
        assert!(samples_metadata.samples[1].name == *"sample2");
        assert!(samples_metadata.samples[0].barcode == "GATTACA");
        assert!(samples_metadata.samples[1].barcode == "CATGCTA");
    }

    // ############################################################################################
    // Test [`SampleGroup::from_file`] - Expected to panic
    // ############################################################################################
    #[test]
    fn test_reading_from_file_with_no_header() {
        let expected_error_message =
            "CSV deserialize error: record 1 (line: 2, byte: 16): missing field `name`";
        let lines = vec!["sample1\tGATTACA", "sample2\tCATGCTA"];
        let tempdir = TempDir::new().unwrap();
        let f1 = tempdir.path().join("sample_metadata.tsv");

        let io = Io::default();
        io.write_lines(&f1, &lines).unwrap();
        println!("{}", SampleGroup::from_file(&f1, b'\t').unwrap_err());
        if let fgoxide::FgError::ConversionError(e) =
            SampleGroup::from_file(&f1, b'\t').unwrap_err()
        {
            assert_eq!(e.to_string(), expected_error_message);
        } else {
            panic!("Different error type than expected reading from headerless file.")
        }
    }

    #[test]
    #[should_panic(expected = "Must provide one or more sample")]
    fn test_reading_header_only_file() {
        let lines = vec!["name\tbarcode"];
        let tempdir = TempDir::new().unwrap();
        let f1 = tempdir.path().join("sample_metadata.tsv");

        let io = Io::default();
        io.write_lines(&f1, &lines).unwrap();
        let _sm = SampleGroup::from_file(&f1, b'\t').unwrap();
    }

    #[test]
    #[should_panic(expected = "Must provide one or more sample")]
    fn test_reading_empty_file() {
        let lines = vec![""];
        let tempdir = TempDir::new().unwrap();
        let f1 = tempdir.path().join("sample_metadata.tsv");

        let io = Io::default();
        io.write_lines(&f1, &lines).unwrap();
        let _sm = SampleGroup::from_file(&f1, b'\t').unwrap();
    }

    #[test]
    fn test_reading_non_existent_file() {
        let tempdir = TempDir::new().unwrap();
        let f1 = tempdir.path().join("sample_metadata.tsv");
        if let fgoxide::FgError::IoError(e) = SampleGroup::from_file(&f1, b'\t').unwrap_err() {
            assert_eq!(e.to_string(), "No such file or directory (os error 2)");
        } else {
            panic!("Different error than expected reading non-existent file")
        }
    }

    // ############################################################################################
    // Test is_valid_base
    // ############################################################################################
    #[test]
    fn test_is_valid_base() {
        assert!(!is_valid_base(b'N'));
        assert!(!is_valid_base(b'n'));
        assert!(!is_valid_base(b'.'));
        assert!(!is_valid_base(b'a'));
        assert!(!is_valid_base(b'c'));
        assert!(!is_valid_base(b'g'));
        assert!(!is_valid_base(b't'));
        assert!(is_valid_base(b'A'));
        assert!(is_valid_base(b'C'));
        assert!(is_valid_base(b'G'));
        assert!(is_valid_base(b'T'));
    }

    // ############################################################################################
    // Test [`Sample::new`] - Expected to pass
    // ############################################################################################
    #[test]
    fn test_new_sample_success() {
        let name = "s_1_example_name".to_owned();
        let barcode = "GATTACA".to_owned();
        let ordinal = 0;
        let sample = Sample::new(ordinal, name, barcode);
        assert_eq!(
            format!("{}", sample),
            "Sample(0000) - { name: s_1_example_name\tbarcode: GATTACA }"
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

    #[test]
    #[should_panic(expected = "All sample barcode bases must be one of A, C, G, or T")]
    fn test_new_sample_fail3_non_agct_bases_in_barcode() {
        let name = "s_1_example_name".to_owned();
        let barcode = "GATTANN".to_owned();
        let ordinal = 0;
        let sample = Sample::new(ordinal, name, barcode);
        assert_eq!(
            format!("{}", sample),
            "Sample(0000) - { name: s_1_example_name\tbarcode: GATTACA }"
        );
    }

    // ############################################################################################
    // Test [`SampleGroup::from_samples`] - expected to pass
    // ############################################################################################
    #[test]
    fn test_from_samples_sample_group_pass1_single_sample() {
        let sample1 = Sample::new(0, "sample_1".to_owned(), "GATTACA".to_owned());
        let samples_vec = vec![sample1.clone()];
        let samples_metadata = SampleGroup::from_samples(&samples_vec);

        let expected_formatted_string = format!("SampleGroup {{\n    {sample1}\n}}\n");
        assert_eq!(format!("{samples_metadata}"), expected_formatted_string);
    }

    #[test]
    fn test_from_samples_sample_group_pass2_multi_unique_samples() {
        let sample1 = Sample::new(0, "sample_1".to_owned(), "GATTACA".to_owned());
        let sample2 = Sample::new(1, "sample_2".to_owned(), "CATGGAT".to_owned());
        let samples_vec = vec![sample1.clone(), sample2.clone()];
        let samples_metadata = SampleGroup::from_samples(&samples_vec);

        let expected_formatted_string =
            format!("SampleGroup {{\n    {sample1}\n    {sample2}\n}}\n");
        assert_eq!(format!("{samples_metadata}"), expected_formatted_string);
    }

    // ############################################################################################
    // Test [`SampleGroup::from_samples`] - expected to panic
    // ############################################################################################
    #[test]
    #[should_panic(expected = "Must provide one or more sample")]
    fn test_from_samples_sample_group_fail1_no_samples() {
        let samples = vec![];
        let _samples_metadata = SampleGroup::from_samples(&samples);
    }

    #[test]
    #[should_panic(expected = "Each sample name must be unique, duplicate identified")]
    fn test_from_samples_sample_group_fail2_duplicate_barcodes() {
        let samples = vec![
            Sample::new(0, "sample_1".to_owned(), "GATTACA".to_owned()),
            Sample::new(0, "sample_1".to_owned(), "CATGGAT".to_owned()),
        ];
        let _samples_metadata = SampleGroup::from_samples(&samples);
    }

    #[test]
    #[should_panic(expected = "Each sample barcode must be unique, duplicate identified")]
    fn test_from_samples_sample_group_fail3_duplicate_sample_names() {
        let samples = vec![
            Sample::new(0, "sample_1".to_owned(), "GATTACA".to_owned()),
            Sample::new(0, "sample_2".to_owned(), "GATTACA".to_owned()),
        ];
        let _samples_metadata = SampleGroup::from_samples(&samples);
    }

    #[test]
    #[should_panic(expected = "All barcodes must have the same length")]
    fn test_from_samples_sample_group_fail4_barcodes_of_different_lengths() {
        let samples = vec![
            Sample::new(0, "sample_1".to_owned(), "GATTACA".to_owned()),
            Sample::new(0, "sample_2".to_owned(), "CATGGA".to_owned()),
        ];
        let _samples_metadata = SampleGroup::from_samples(&samples);
    }
}
