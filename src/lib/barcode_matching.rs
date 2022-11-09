use bstr::{BString, ByteSlice};

/// Checks whether a given u8 byte is a "No-call"-ed base, signified by the bytes 'N', 'n' and '.'
fn byte_is_nocall(byte: u8) -> bool {
    byte == b'N' || byte == b'n' || byte == b'.'
}

/// The struct that contains the info related to the best and next best sample barcode match.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct BarcodeMatch {
    /// Index of the best bardcode match in the corresponding BarcodeMatcher struct that generated
    /// this match.
    best_match: usize,
    /// The number of mismatches to the best matching barcode for the read described by this match.
    best_mismatches: u8,
    /// The number of mismatches to the second best matching barcode for the read described by this
    /// match
    next_best_mismatches: u8,
}

/// The struct responsible for matching barcodes to a ``Vec`` of sample barcodes.
#[derive(Clone, Debug)]
pub struct BarcodeMatcher {
    /// Vec of the barcodes for each sample
    /// Note - this is to be replaced by a sample struct in task 3. For now we're keeping things
    /// very simple.
    sample_barcodes: Vec<BString>,
    /// The maximum number of no calls (Ns) in the sample barcode bases allowed for matching.
    max_no_calls: usize,
    /// The maximum mismatches to match a sample barcode.
    max_mismatches: u8,
    /// The minimum difference between number of mismatches in the best and second best barcodes
    /// for a barcode to be considered a match.
    min_mismatch_delta: u8,
}

impl BarcodeMatcher {
    /// Instantiates a new ``BarcodeMatcher`` struct. Checks that the sample barcodes vector is not
    /// empty and that none of the barcodes provided are the empty string.
    ///
    /// # Panics
    /// - Will panic if provided an empty vec of sample barcodes.
    /// - Will panic if any provided barcode is length zero.
    #[must_use]
    pub fn new(
        sample_barcodes: &[&str],
        max_no_calls: Option<usize>,
        max_mismatches: Option<u8>,
        min_mismatch_delta: Option<u8>,
    ) -> Self {
        assert!(!sample_barcodes.is_empty(), "Must provide at least one sample barcode");

        let modified_sample_barcodes = sample_barcodes.iter().map(
            |barcode|
            if barcode.is_empty() {
                None
            } else {
                Some(BString::from(barcode.to_ascii_uppercase()))
            }
        ).collect::<Option<Vec<_>>>().expect("Sample barcode cannot be empty string");

        Self {
            sample_barcodes: modified_sample_barcodes,
            max_no_calls: max_no_calls.unwrap_or(2),
            max_mismatches: max_mismatches.unwrap_or(2),
            min_mismatch_delta: min_mismatch_delta.unwrap_or(1),
        }
    }

    /// Counts the number of bases that differ between two byte arrays.
    fn count_mismatches(observed_bases: &[u8], expected_bases: &[u8]) -> u8 {
        assert_eq!(
            observed_bases.len(),
            expected_bases.len(),
            "observed_bases: {}, expected_bases: {}",
            observed_bases.len(),
            expected_bases.len()
        );
        let mut count: usize = 0;
        for (&expected_base, &observed_base) in expected_bases.iter().zip(observed_bases.iter()) {
            if !byte_is_nocall(expected_base) && expected_base != observed_base {
                count += 1;
            }
        }
        u8::try_from(count).expect("Overflow on number of mismatch bases")
    }

    /// Assigns the barcode that best matches the provided ``read_bases``.
    #[must_use]
    pub fn assign(&self, read_bases: &[u8]) -> Option<BarcodeMatch> {
        let num_no_call = read_bases
            .iter()
            .filter(|&&byte| byte_is_nocall(byte))
            .count();

        if num_no_call > self.max_no_calls {
            return None;
        }

        let mut best_barcode_index = self.sample_barcodes.len();
        let mut best_mismatches = 255u8;
        let mut next_best_mismatches = 255u8;
        for (index, sample_barcode) in self.sample_barcodes.iter().enumerate() {
            let mismatches = Self::count_mismatches(read_bases, sample_barcode.as_bstr());

            if mismatches < best_mismatches {
                next_best_mismatches = best_mismatches;
                best_mismatches = mismatches;
                best_barcode_index = index;
            } else if  mismatches < next_best_mismatches {
                next_best_mismatches = mismatches;
            }
        }
        if best_mismatches > self.max_mismatches || (next_best_mismatches - best_mismatches) < self.min_mismatch_delta
        {
            None
        } else {
           Some(
                BarcodeMatch {
                    best_match: best_barcode_index,
                    best_mismatches,
                    next_best_mismatches,
                }
           )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use bstr::BStr;

    // ############################################################################################
    // Test ``BarcodeMatcher`` instantiation panics.
    // ############################################################################################
    #[test]
    fn test_barcode_matcher_instantiation_can_succeed() {
        let sample_barcodes = vec!["AGCT"];
        let _matcher = BarcodeMatcher::new(
            &sample_barcodes,
            None,
            None,
            None,
        );
    }

    #[test]
    #[should_panic(expected = "Must provide at least one sample barcode")]
    fn test_barcode_matcher_fails_if_no_sample_barcodes_provided() {
        let sample_barcodes: Vec<&str> = vec![];
        let _matcher = BarcodeMatcher::new(
            &sample_barcodes,
            None,
            None,
            None,
        );
    }

    #[test]
    #[should_panic(expected = "Sample barcode cannot be empty string")]
    fn test_barcode_matcher_fails_if_empty_sample_barcode_provided() {
        let sample_barcodes = vec!["AGCT", ""];
        let _matcher = BarcodeMatcher::new(
            &sample_barcodes,
            None,
            None,
            None,
        );
    }

    // ############################################################################################
    // Test byte_is_no_call
    // ############################################################################################
    #[test]
    fn test_byte_is_no_call() {
        assert!(byte_is_nocall(b'N'));
        assert!(byte_is_nocall(b'n'));
        assert!(byte_is_nocall(b'.'));
        assert!(!byte_is_nocall(b'A'));
        assert!(!byte_is_nocall(b'C'));
        assert!(!byte_is_nocall(b'G'));
        assert!(!byte_is_nocall(b'T'));
        assert!(!byte_is_nocall(b'a'));
        assert!(!byte_is_nocall(b'c'));
        assert!(!byte_is_nocall(b'g'));
        assert!(!byte_is_nocall(b't'));
    }

    // ############################################################################################
    // Test BarcodeMatcher::count_mismatches
    // ############################################################################################

    // Thought process behind not panicking on empty string is:
    //   1. sample barcodes are checked to see if they are empty and sample vs read barcodes are
    //      compared for length, so empty string will fail anyway
    //   2. the fewer operations / better optimization in this core matching function the better.
    #[test]
    fn empty_string_can_run_in_count_mismatches() {
        assert_eq!(BarcodeMatcher::count_mismatches("".as_bytes(), "".as_bytes()), 0);
    }

    #[test]
    fn find_no_mismatches() {
        assert_eq!(BarcodeMatcher::count_mismatches("GATTACA".as_bytes(), "GATTACA".as_bytes()), 0,);
    }

    #[test]
    fn ns_in_expected_barcode_dont_contribute_to_mismatch_counter() {
        assert_eq!(BarcodeMatcher::count_mismatches("GATTACA".as_bytes(), "GANNACA".as_bytes()), 0,);
    }

    #[test]
    fn find_two_mismatches() {
        assert_eq!(BarcodeMatcher::count_mismatches("GATTACA".as_bytes(), "GACCACA".as_bytes()), 2,);
    }

    #[test]
    fn not_count_no_calls() {
        assert_eq!(BarcodeMatcher::count_mismatches("GATTACA".as_bytes(), "GANNACA".as_bytes()), 0,);
    }

    #[test]
    fn find_compare_two_sequences_that_have_all_mismatches() {
        assert_eq!(BarcodeMatcher::count_mismatches("GATTACA".as_bytes(), "CTAATGT".as_bytes()), 7,);
    }

    #[test]
    #[should_panic(expected = "observed_bases: 5, expected_bases: 6")]
    fn find_compare_two_sequences_of_different_length() {
        let _mismatches = BarcodeMatcher::count_mismatches("GATTA".as_bytes(), "CTATGT".as_bytes());
    }

    // ############################################################################################
    // Test BarcodeMatcher::assign
    // ############################################################################################
    // Some returns
    #[test]
    fn test_assign_exact_match() {
        const EXPECTED_BARCODE_INDEX: usize = 0;
        let sample_barcodes = vec!["ACGT", "AAAG", "CACA"];
        let matcher = BarcodeMatcher::new(&sample_barcodes, Some(2), Some(2), Some(2));
        assert_eq!(
            matcher.assign(sample_barcodes[EXPECTED_BARCODE_INDEX].as_bytes()),
            Some(
                BarcodeMatch{
                    best_match: EXPECTED_BARCODE_INDEX,
                    best_mismatches: 0,
                    next_best_mismatches: 3,
                }
            ),
        );
    }

    #[test]
    fn test_assign_imprecise_match() {
        let sample_barcodes = vec!["AAAT", "AGAG", "CACA"];
        let matcher = BarcodeMatcher::new(
            &sample_barcodes, Some(2), Some(2), Some(2)
        );
        //                                   1 different base
        //                                   |
        //                                   v
        let test_barcode: &BStr = BStr::new("GAAT".as_bytes());
        let expected = BarcodeMatch{best_match: 0, best_mismatches: 1, next_best_mismatches: 3};
        assert_eq!(matcher.assign(test_barcode), Some(expected));
    }

    #[test]
    fn test_assign_precise_match_with_no_call() {
        let sample_barcodes = vec!["AAAT", "AGAG", "CACA"];
        let matcher = BarcodeMatcher::new(
            &sample_barcodes, Some(2), Some(2), Some(2)
        );
        //                                   1 no-call
        //                                   |
        //                                   v
        let test_barcode: &BStr = BStr::new("NAAT".as_bytes());
        let expected = BarcodeMatch{best_match: 0, best_mismatches: 1, next_best_mismatches: 3};
        assert_eq!(matcher.assign(test_barcode), Some(expected));
    }

    #[test]
    fn test_assign_imprecise_match_with_no_call() {
        let sample_barcodes = vec!["AAATTT", "AGAGGG", "CACAGG"];
        let matcher = BarcodeMatcher::new(
            &sample_barcodes, Some(2), Some(2), Some(2)
        );
        //                             1 no-call
        //                             |
        //                             | 1 different base
        //                             | |
        //                             v v
        let test_barcode: &[u8; 6] = b"NAGTTT";
        let expected = BarcodeMatch{best_match: 0, best_mismatches: 2, next_best_mismatches: 5};
        assert_eq!(matcher.assign(test_barcode), Some(expected));
    }

    #[test]
    fn test_sample_no_call_doesnt_contribute_to_mismatch_number() {
        let sample_barcodes = vec!["NAGTTT", "AGAGGG", "CACAGG"];
        let matcher = BarcodeMatcher::new(
            &sample_barcodes, Some(2), Some(1), Some(2)
        );
        //                             1 no-call
        //                             |
        //                             | 1 different base
        //                             | |
        //                             v v
        let test_barcode: &[u8; 6] = b"AAATTT";
        let expected = BarcodeMatch{best_match: 0, best_mismatches: 1, next_best_mismatches: 4};
        assert_eq!(matcher.assign(test_barcode), Some(expected));
    }

    // None returns
    #[test]
    fn test_read_no_call_contributes_to_mismatch_number() {
        let sample_barcodes = vec!["AAATTT", "AGAGGG", "CACAGG"];
        let matcher = BarcodeMatcher::new(
            &sample_barcodes, Some(2), Some(1), Some(2)
        );
        //                             1 no-call
        //                             |
        //                             | 1 different base
        //                             | |
        //                             v v
        let test_barcode: &[u8; 6] = b"NAGTTT";
        assert_eq!(matcher.assign(test_barcode), None);
    }

    #[test]
    fn test_produce_no_match_if_too_many_mismatches() {
        let sample_barcodes = vec![
            "AAGCTAG",
            "CAGCTAG",
            "GAGCTAG",
            "TAGCTAG",
        ];
        let assignment_barcode: &[u8] = b"ATCGATC";
        let matcher = BarcodeMatcher::new(&sample_barcodes, Some(100), Some(0), Some(100));
        assert_eq!(matcher.assign(assignment_barcode), None);
    }

    #[test]
    fn test_produce_no_match_if_within_mismatch_delta() {
        let sample_barcodes = vec![
            "AAAAAAAA",
            "CCCCCCCC",
            "GGGGGGGG",
            "GGGGGGTT",
        ];
        let assignment_barcode: &[u8] = sample_barcodes[3].as_bytes();
        let matcher = BarcodeMatcher::new(&sample_barcodes, Some(100), Some(100), Some(3));
        assert_eq!(matcher.assign(assignment_barcode), None);
    }

    #[test]
    fn test_produce_no_match_if_too_many_no_calls() {
        let sample_barcodes = vec![
            "AAAAAAAA",
            "CCCCCCCC",
            "GGGGGGGG",
            "GGGGGGTT",
        ];
        let assignment_barcode: &[u8] = b"GGGGGGTN";
        let matcher = BarcodeMatcher::new(&sample_barcodes, Some(0), Some(100), Some(100));
        assert_eq!(matcher.assign(assignment_barcode), None);
    }
}
