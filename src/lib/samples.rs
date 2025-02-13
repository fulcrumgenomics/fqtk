use crate::is_valid_iupac;

use anyhow::Result;
use fgoxide::io::DelimFile;
use itertools::Itertools;
use serde::Deserialize;
use serde_aux::prelude::*;
use std::collections::HashSet;
use std::collections::hash_map::RandomState;
use std::fmt::{self, Display};
use std::path::Path;

const DEFAULT_FILE_DELIMETER: u8 = b'\t';

/// Struct for describing a single sample and metadata associated with that sample.
#[derive(Clone, Deserialize, Debug, PartialEq, Eq)]
pub struct Sample {
    /// ID of the sample or library
    pub sample_id: String,
    /// DNA barcode associated with the sample
    pub barcode: String,
    /// index of the sample in the [`SampleGroup`] object, used for syncing indices across
    /// different structs
    #[serde(skip_deserializing)]
    ordinal: usize,
}

impl Display for Sample {
    /// Implements a nice format display for the [`Sample`] struct.
    /// E.g. A sample with ordinal 2, name test-sample, and barcode GATTACA would look like:
    /// Sample(0002) - { name: test-sample    barcode: GATTACA}
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
    ///   - Panics if barcode has bases other than A, C, G, T, or N/n/.
    #[must_use]
    pub fn new(ordinal: usize, name: String, barcode: String) -> Self {
        assert!(!name.is_empty(), "Sample name cannot be empty");
        assert!(!barcode.is_empty(), "Sample barcode cannot be empty");
        assert!(
            barcode.as_bytes().iter().all(|&b| is_valid_iupac(b)),
            "All sample barcode bases must be one of A, C, G, T, U, R, Y, S, W, K, M, D, V, H, B, N"
        );
        Self { sample_id: name, barcode, ordinal }
    }

    /// Returns the header line expected by serde when deserializing
    #[must_use]
    pub fn deserialize_header_line() -> String {
        let field_names = serde_introspect::<Self>();
        let skip_deserialize_fields: HashSet<&str, RandomState> = HashSet::from_iter(["ordinal"]);
        let final_field_names: Vec<String> = field_names
            .iter()
            .filter(|&&f| !skip_deserialize_fields.contains(f))
            .map(|&f| f.to_owned())
            .collect();
        final_field_names.join("\t")
    }
}

/// Struct for storing information about multiple samples and for defining functions associated
/// with groups of [`Sample`]s, rather than individual structs.
#[derive(Clone, Debug, PartialEq, Eq)]
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
    ///   - Will panic if there are duplicate sample names provided
    ///   - Will panic if there are duplicate barcodes provided
    ///   - Will panic if barcodes don't all have the same length
    #[must_use]
    pub fn from_samples(samples: &[Sample]) -> Self {
        // Validate that we have at least one name
        assert!(!samples.is_empty(), "Must provide one or more sample");

        // Validate that all the sample names are unique
        assert!(
            samples.iter().map(|s| &s.sample_id).all_unique(),
            "Each sample name must be unique, duplicate identified"
        );

        // Validate that the barcodes are all unique
        assert!(
            samples.iter().map(|s| &s.barcode).all_unique(),
            "Each sample barcode must be unique, duplicate identified",
        );

        // Validate that the barcodes are all of the same length
        let first_barcode_length = samples[0].barcode.len();
        assert!(
            samples.iter().map(|s| &s.barcode).all(|b| b.len() == first_barcode_length),
            "All barcodes must have the same length",
        );

        Self {
            samples: samples
                .iter()
                .enumerate()
                .map(|(ordinal, sample)| {
                    Sample::new(ordinal, sample.sample_id.clone(), sample.barcode.clone())
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
    pub fn from_file<P: AsRef<Path>>(path: &P) -> Result<SampleGroup, fgoxide::FgError> {
        let reader = DelimFile::default();
        Ok(Self::from_samples(&reader.read(path, DEFAULT_FILE_DELIMETER, false)?))
    }
}

#[cfg(test)]
mod tests {
    use core::panic;

    use super::*;
    use csv::DeserializeErrorKind as CsvDeserializeErrorEnum;
    use csv::ErrorKind as CsvErrorEnum;
    use fgoxide::{self, io::Io};
    use serde::de::Error;
    use serde::de::value::Error as SerdeError;
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
        let samples_metadata = SampleGroup::from_file(&f1).unwrap();

        assert!(samples_metadata.samples[0].sample_id == "sample1");
        assert!(samples_metadata.samples[1].sample_id == "sample2");
        assert!(samples_metadata.samples[0].barcode == "GATTACA");
        assert!(samples_metadata.samples[1].barcode == "CATGCTA");
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
        let samples_metadata = SampleGroup::from_file(&f1).unwrap();

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

    // ############################################################################################
    // Test [`SampleGroup::from_file`] - Expected to panic
    // ############################################################################################
    #[test]
    fn test_reading_from_file_with_no_header() {
        let lines = vec!["sample1\tGATTACA", "sample2\tCATGCTA"];
        let tempdir = TempDir::new().unwrap();
        let f1 = tempdir.path().join("sample_metadata.tsv");

        let io = Io::default();
        io.write_lines(&f1, &lines).unwrap();
        let mut to_panic = true;
        if let fgoxide::FgError::ConversionError(csv_e) = SampleGroup::from_file(&f1).unwrap_err() {
            if let CsvErrorEnum::Deserialize { pos: _, err: csv_de_err } = csv_e.into_kind() {
                if let CsvDeserializeErrorEnum::Message(s) = csv_de_err.kind() {
                    to_panic = false;
                    assert_eq!(s, &SerdeError::missing_field("sample_id").to_string());
                }
            }
        }
        assert!(!to_panic, "Different error type than expected reading from headerless file.");
    }

    #[test]
    #[should_panic(expected = "Must provide one or more sample")]
    fn test_reading_header_only_file() {
        let lines = vec![Sample::deserialize_header_line()];
        let tempdir = TempDir::new().unwrap();
        let f1 = tempdir.path().join("sample_metadata.tsv");

        let io = Io::default();
        io.write_lines(&f1, &lines).unwrap();
        let _sm = SampleGroup::from_file(&f1).unwrap();
    }

    #[test]
    #[should_panic(expected = "Must provide one or more sample")]
    fn test_reading_empty_file() {
        let lines = vec![""];
        let tempdir = TempDir::new().unwrap();
        let f1 = tempdir.path().join("sample_metadata.tsv");

        let io = Io::default();
        io.write_lines(&f1, &lines).unwrap();
        let _sm = SampleGroup::from_file(&f1).unwrap();
    }

    #[test]
    fn test_reading_non_existent_file() {
        let tempdir = TempDir::new().unwrap();
        let f1 = tempdir.path().join("sample_metadata.tsv");
        if let fgoxide::FgError::IoError(e) = SampleGroup::from_file(&f1).unwrap_err() {
            assert_eq!(e.to_string(), "No such file or directory (os error 2)");
        } else {
            panic!("Different error than expected reading non-existent file")
        }
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
            Sample { sample_id: name, barcode, ordinal },
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
}
