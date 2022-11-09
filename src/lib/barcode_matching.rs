use bstr::{BStr, BString, ByteSlice};

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
    pub fn assign(&self, read_bases: &BStr) -> Option<BarcodeMatch> {
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
    use bstr::{BStr, BString};

    const SAMPLE_BARCODE_1: &str = "AAAAAAAAGATTACAGA";
    const SAMPLE_BARCODE_2: &str = "CCCCCCCCGATTACAGA";
    const SAMPLE_BARCODE_3: &str = "GGGGGGGGGATTACAGA";
    const SAMPLE_BARCODE_4: &str = "GGGGGGTTGATTACAGA";

    fn default_barcodes() -> Vec<BString> {
        vec![
            BString::from(SAMPLE_BARCODE_1.as_bytes()),
            BString::from(SAMPLE_BARCODE_2.as_bytes()),
            BString::from(SAMPLE_BARCODE_3.as_bytes()),
            BString::from(SAMPLE_BARCODE_4.as_bytes()),
        ]
    }

    fn expected_barcode_match(
        expected_barcode_index: usize,
        best_mismatches: u8,
        next_best_mismatches: u8,
    ) -> BarcodeMatch {
        BarcodeMatch {
            best_match: expected_barcode_index,
            best_mismatches,
            next_best_mismatches,
        }
    }

    fn get_matcher(
        sample_barcodes: Vec<BString>,
        max_no_calls: Option<usize>,
        max_mismatches: Option<u8>,
        min_mismatch_delta: Option<u8>,
    ) -> BarcodeMatcher {
        BarcodeMatcher {
            sample_barcodes,
            max_no_calls: max_no_calls.unwrap_or(2),
            max_mismatches: max_mismatches.unwrap_or(2),
            min_mismatch_delta: min_mismatch_delta.unwrap_or(1),
        }
    }

    fn default_matcher() -> BarcodeMatcher {
        get_matcher(default_barcodes(), None, None, None)
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
    }

    // ############################################################################################
    // Test BarcodeMatcher::count_mismatches
    // ############################################################################################
    #[test]
    fn find_no_mismatches() {
        assert_eq!(BarcodeMatcher::count_mismatches("GATTACA".as_bytes(), "GATTACA".as_bytes()), 0,);
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

    // ############################################################################################
    // Test BarcodeMatcher::assign
    // ############################################################################################
    // Some returns
    #[test]
    fn test_assign_fragment_reads_demuxes_exact_match() {
        let matcher = default_matcher();
        assert_eq!(
            matcher.assign(BStr::new(SAMPLE_BARCODE_1.as_bytes())),
            Some(expected_barcode_match(0, 0, 8)),
        );
    }

    #[test]
    fn test_assign_fragment_reads_demuxes_imprecise_match() {
        let matcher = default_matcher();
        //                                   1 different base
        //                                   |
        //                                   v
        let test_barcode: &BStr = BStr::new("GAAAAAAAGATTACAGA".as_bytes());
        let expected = expected_barcode_match(0, 1, 7);
        assert_eq!(matcher.assign(test_barcode), Some(expected));
    }

    #[test]
    fn test_assign_fragment_reads_demuxes_precise_match_with_no_call() {
        let matcher = default_matcher();
        //                                   1 no-call
        //                                   |
        //                                   v
        let test_barcode: &BStr = BStr::new("NAAAAAAAGATTACAGA".as_bytes());
        let expected = expected_barcode_match(0, 1, 8);
        assert_eq!(matcher.assign(test_barcode), Some(expected));
    }

    #[test]
    fn test_assign_fragment_reads_demuxes_imprecise_match_with_no_call() {
        let matcher = default_matcher();
        //                                   1 no-call
        //                                   |
        //                                   | 1 different base
        //                                   | |
        //                                   v v
        let test_barcode: &BStr = BStr::new("NAGAAAAAGATTACAGA".as_bytes());
        let expected = expected_barcode_match(0, 2, 7);
        assert_eq!(matcher.assign(test_barcode), Some(expected));
    }

    // None returns
    #[test]
    fn test_produce_no_match_if_too_many_mismatches() {
        let matcher = get_matcher(default_barcodes(), None, Some(0), None);
        assert_eq!(matcher.assign(BStr::new(b"AAAAAAAAGATTACAGT")), None,);
    }

    #[test]
    fn test_produce_no_match_if_within_mismatch_delta() {
        let matcher = get_matcher(default_barcodes(), None, Some(10), Some(3));
        assert_eq!(matcher.assign(BStr::new(SAMPLE_BARCODE_4.as_bytes())), None,);
    }

    #[test]
    fn test_produce_no_match_if_too_many_no_calls() {
        let matcher = get_matcher(default_barcodes(), Some(0), Some(10), None);
        assert_eq!(matcher.assign(BStr::new("GGGGGGTTGATTACAGN".as_bytes())), None,);
    }
}
