use crate::encode_expected;
use crate::encode_observed;

use super::byte_is_nocall;
use ahash::HashMap as AHashMap;
use ahash::HashMapExt;
use bio_seq::prelude::*;
use itertools::Itertools;

const STARTING_CACHE_SIZE: usize = 1_000_000;

/// The struct that contains the info related to the best and next best sample barcode match.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct BarcodeMatch {
    /// Index of the best bardcode match in the corresponding ``BarcodeMatcher`` struct that generated
    /// this match.
    pub best_match: usize,
    /// The number of mismatches to the best matching barcode for the read described by this match.
    pub best_mismatches: u8,
    /// The number of mismatches to the second best matching barcode for the read described by this
    /// match
    pub next_best_mismatches: u8,
}

/// The struct responsible for matching barcodes to a ``Vec`` of sample barcodes.
#[derive(Clone, Debug)]
pub struct BarcodeMatcher {
    /// Vec of the barcodes for each sample
    sample_barcodes: Vec<Seq<Iupac>>,
    /// The maxium number of Ns in any barcode in set of sample barcodes
    max_ns_in_barcodes: usize,
    /// The maximum mismatches to match a sample barcode.
    max_mismatches: u8,
    /// The minimum difference between number of mismatches in the best and second best barcodes
    /// for a barcode to be considered a match.
    min_mismatch_delta: u8,
    /// If true will attempt to use the cache when matching.
    use_cache: bool,
    /// Caching struct for storing results of previous matches
    cache: AHashMap<Vec<u8>, BarcodeMatch>,
    /// Vec of booleans for each sample indicating whether the sample barcode has any Ns.
    has_ns: Vec<bool>,
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
        max_mismatches: u8,
        min_mismatch_delta: u8,
        use_cache: bool,
    ) -> Self {
        assert!(!sample_barcodes.is_empty(), "Must provide at least one sample barcode");
        assert!(
            sample_barcodes.iter().all(|b| !b.is_empty()),
            "Sample barcode cannot be empty string"
        );

        let mut max_ns_in_barcodes = 0;
        let modified_sample_barcodes = sample_barcodes
            .iter()
            .map(|barcode| {
                let num_ns: usize = barcode.chars().filter(|&b| byte_is_nocall(b as u8)).count();
                let barcode: Seq<Iupac> = encode_expected(barcode.to_ascii_uppercase().as_bytes());
                max_ns_in_barcodes = max_ns_in_barcodes.max(num_ns);
                barcode
            })
            .collect::<Vec<_>>();
        let has_ns =
            modified_sample_barcodes.iter().map(|b| (*b).iter().contains(&Iupac::N)).collect();

        Self {
            sample_barcodes: modified_sample_barcodes,
            max_ns_in_barcodes,
            max_mismatches,
            min_mismatch_delta,
            use_cache,
            cache: AHashMap::with_capacity(STARTING_CACHE_SIZE),
            has_ns,
        }
    }

    /// Counts the number of bases that differ between two byte arrays.
    fn count_mismatches(
        observed_bases: &Seq<Iupac>,
        expected_bases: &Seq<Iupac>,
        expected_has_ns: bool,
    ) -> u8 {
        assert_eq!(
            observed_bases.len(),
            expected_bases.len(),
            "observed_bases: {}, expected_bases: {}",
            observed_bases.len(),
            expected_bases.len()
        );

        let observed_bases: &SeqSlice<Iupac> = observed_bases;
        let expected_bases: &SeqSlice<Iupac> = expected_bases;
        let intersection: Seq<Iupac> = observed_bases & expected_bases;
        let mut count: usize = intersection.iter().filter(|b| b.to_char() as u8 == b'-').count();

        // If observed bases contains '-' and the expected bases contains 'N', treat
        // '-'/'N' as a match.
        // TODO: would just need to complement observed, then take the intersection, and count the
        // # of Ns.  Needs bio-seq to implement the complement of Iupac bases:
        // https://github.com/jeff-k/bio-seq/blame/c5c657603ecb845ba7398f5d78cb31bdf19351f5/bio-seq/src/codec/iupac.rs#L153
        if expected_has_ns && observed_bases.iter().contains(&Iupac::X) {
            for (obs, exp) in observed_bases.iter().zip(expected_bases.iter()) {
                if obs == Iupac::X && exp == Iupac::N {
                    count -= 1;
                }
            }
        }
        u8::try_from(count).expect("Overflow on number of mismatch bases")
    }

    /// Returns the expected barcode length, assuming a fixed length for all samples.
    fn expected_barcode_length(&self) -> usize {
        self.sample_barcodes[0].len()
    }

    /// Assigns the barcode that best matches the provided ``read_bases``.
    #[must_use]
    fn assign_internal(&self, read_bases: &[u8]) -> Option<BarcodeMatch> {
        let mut best_barcode_index = self.sample_barcodes.len();
        let mut best_mismatches = 255u8;
        let mut next_best_mismatches = 255u8;
        let read_bases = encode_observed(&read_bases.to_ascii_uppercase());
        for (index, sample_barcode) in self.sample_barcodes.iter().enumerate() {
            let mismatches =
                Self::count_mismatches(&read_bases, sample_barcode, self.has_ns[index]);
            if mismatches < best_mismatches {
                next_best_mismatches = best_mismatches;
                best_mismatches = mismatches;
                best_barcode_index = index;
            } else if mismatches < next_best_mismatches {
                next_best_mismatches = mismatches;
            }
        }

        if best_mismatches > self.max_mismatches
            || (next_best_mismatches - best_mismatches) < self.min_mismatch_delta
        {
            None
        } else {
            Some(BarcodeMatch {
                best_match: best_barcode_index,
                best_mismatches,
                next_best_mismatches,
            })
        }
    }

    /// Assigns the barcode that best matches the provided ``read_bases``, using internal caching
    /// if configured to do so and skipping calculation for reads that cannot match any barcode (
    /// due to having too many no-called bases).
    pub fn assign(&mut self, read_bases: &[u8]) -> Option<BarcodeMatch> {
        // do not try matching if there are not enough bases
        if read_bases.len() < self.expected_barcode_length() {
            return None;
        }
        let num_no_calls = read_bases.iter().filter(|&&b| byte_is_nocall(b)).count();
        if num_no_calls > (self.max_mismatches as usize) + self.max_ns_in_barcodes {
            None
        } else if self.use_cache {
            if let Some(cached_match) = self.cache.get(read_bases) {
                Some(*cached_match)
            } else {
                let maybe_match = self.assign_internal(read_bases);
                if let Some(internal_val) = maybe_match {
                    self.cache.insert(read_bases.to_vec(), internal_val);
                };
                maybe_match
            }
        } else {
            self.assign_internal(read_bases)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    // ############################################################################################
    // Test ``BarcodeMatcher`` instantiation panics.
    // ############################################################################################
    #[rstest]
    #[case(true)]
    #[case(false)]
    fn test_barcode_matcher_instantiation_can_succeed(#[case] use_cache: bool) {
        let sample_barcodes = vec!["AGCT"];
        let _matcher = BarcodeMatcher::new(&sample_barcodes, 2, 1, use_cache);
    }

    #[rstest]
    #[case(true)]
    #[case(false)]
    #[should_panic(expected = "Must provide at least one sample barcode")]
    fn test_barcode_matcher_fails_if_no_sample_barcodes_provided(#[case] use_cache: bool) {
        let sample_barcodes: Vec<&str> = vec![];
        let _matcher = BarcodeMatcher::new(&sample_barcodes, 2, 1, use_cache);
    }

    #[rstest]
    #[case(true)]
    #[case(false)]
    #[should_panic(expected = "Sample barcode cannot be empty string")]
    fn test_barcode_matcher_fails_if_empty_sample_barcode_provided(#[case] use_cache: bool) {
        let sample_barcodes = vec!["AGCT", ""];
        let _matcher = BarcodeMatcher::new(&sample_barcodes, 2, 1, use_cache);
    }

    fn count_mismatches(observed_bases: &str, expected_bases: &str) -> u8 {
        let expected_bases = encode_expected(expected_bases.as_bytes());
        let has_ns = expected_bases.iter().contains(&Iupac::N);

        BarcodeMatcher::count_mismatches(
            &encode_observed(observed_bases.as_bytes()),
            &expected_bases,
            has_ns,
        )
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
        assert_eq!(count_mismatches("", ""), 0);
    }

    #[test]

    fn find_no_mismatches() {
        assert_eq!(count_mismatches("GATTACA", "GATTACA"), 0,);
    }

    #[test]
    fn ns_in_expected_barcode_dont_contribute_to_mismatch_counter() {
        assert_eq!(count_mismatches("GATTACA", "GANNACA"), 0,);
    }

    #[test]
    fn all_ns_barcode_have_no_mismatches() {
        assert_eq!(count_mismatches("GANNACA", "NNNNNNN"), 0,);
    }

    #[test]
    fn find_two_mismatches() {
        assert_eq!(count_mismatches("GATTACA", "GACCACA"), 2,);
    }

    #[test]
    fn not_count_no_calls() {
        assert_eq!(count_mismatches("GATTACA", "GANNACA"), 0,);
    }

    #[test]
    fn find_compare_two_sequences_that_have_all_mismatches() {
        assert_eq!(count_mismatches("GATTACA", "CTAATGT"), 7,);
    }

    #[test]
    fn find_compare_iupac_barcode() {
        assert_eq!(count_mismatches("ACGTAAACCGAAACA", "ACGTMRWSYKVHDBN"), 0,);
        // IUPAC bases are mismatches in the observed barcodes
        assert_eq!(count_mismatches("ACGTMRWSYKVHDBN", "ACGTAAACCGAAACA"), 11,);
    }

    #[test]
    fn count_mismatches_iupac_bases_assymetry() {
        // if the observed base is an N, it will not match anything but an N
        assert_eq!(count_mismatches("N", "R"), 1,);
        assert_eq!(count_mismatches("N", "N"), 0,);
        // if the observed base is not A, C, G, T, U, or N, it will not match anything except an N
        assert_eq!(count_mismatches("R", "R"), 1,);
        assert_eq!(count_mismatches("R", "V"), 1,);
        assert_eq!(count_mismatches("R", "D"), 1,);
        assert_eq!(count_mismatches("R", "B"), 1,);
        assert_eq!(count_mismatches("R", "N"), 0,);
    }

    #[test]
    #[should_panic(expected = "observed_bases: 5, expected_bases: 6")]
    fn find_compare_two_sequences_of_different_length() {
        let _mismatches = count_mismatches("GATTA", "CTATGT");
    }

    // ############################################################################################
    // Test BarcodeMatcher::assign
    // ############################################################################################
    // Some returns
    #[rstest]
    #[case(true)]
    #[case(false)]
    fn test_assign_exact_match(#[case] use_cache: bool) {
        const EXPECTED_BARCODE_INDEX: usize = 0;
        let sample_barcodes = vec!["ACGT", "AAAG", "CACA"];
        let mut matcher = BarcodeMatcher::new(&sample_barcodes, 2, 2, use_cache);
        assert_eq!(
            matcher.assign(sample_barcodes[EXPECTED_BARCODE_INDEX].as_bytes()),
            Some(BarcodeMatch {
                best_match: EXPECTED_BARCODE_INDEX,
                best_mismatches: 0,
                next_best_mismatches: 3,
            }),
        );
    }

    #[rstest]
    #[case(true)]
    #[case(false)]
    fn test_assign_imprecise_match(#[case] use_cache: bool) {
        let sample_barcodes = vec!["AAAT", "AGAG", "CACA"];
        let mut matcher = BarcodeMatcher::new(&sample_barcodes, 2, 2, use_cache);
        //                          1 different base
        //                          |
        //                          v
        let test_barcode: &[u8] = b"GAAT";
        let expected = BarcodeMatch { best_match: 0, best_mismatches: 1, next_best_mismatches: 3 };
        assert_eq!(matcher.assign(test_barcode), Some(expected));
    }

    #[rstest]
    #[case(true)]
    #[case(false)]
    fn test_assign_precise_match_with_no_call(#[case] use_cache: bool) {
        let sample_barcodes = vec!["AAAT", "AGAG", "CACA"];
        let mut matcher = BarcodeMatcher::new(&sample_barcodes, 2, 2, use_cache);
        //                             1 no-call
        //                             |
        //                             v
        let test_barcode: &[u8; 4] = b"NAAT";
        let expected = BarcodeMatch { best_match: 0, best_mismatches: 1, next_best_mismatches: 3 };
        assert_eq!(matcher.assign(test_barcode), Some(expected));
    }

    #[rstest]
    #[case(true)]
    #[case(false)]
    fn test_assign_imprecise_match_with_no_call(#[case] use_cache: bool) {
        let sample_barcodes = vec!["AAATTT", "AGAGGG", "CACAGG"];
        let mut matcher = BarcodeMatcher::new(&sample_barcodes, 2, 2, use_cache);
        //                             1 no-call
        //                             |
        //                             | 1 different base
        //                             | |
        //                             v v
        let test_barcode: &[u8; 6] = b"NAGTTT";
        let expected = BarcodeMatch { best_match: 0, best_mismatches: 2, next_best_mismatches: 5 };
        assert_eq!(matcher.assign(test_barcode), Some(expected));
    }

    #[rstest]
    #[case(true)]
    #[case(false)]
    fn test_sample_no_call_doesnt_contribute_to_mismatch_number(#[case] use_cache: bool) {
        let sample_barcodes = vec!["NAGTTT", "AGAGGG", "CACAGG"];
        let mut matcher = BarcodeMatcher::new(&sample_barcodes, 1, 2, use_cache);
        //                             1 no-call
        //                             |
        //                             | 1 different base
        //                             | |
        //                             v v
        let test_barcode: &[u8; 6] = b"AAATTT";
        let expected = BarcodeMatch { best_match: 0, best_mismatches: 1, next_best_mismatches: 4 };
        assert_eq!(matcher.assign(test_barcode), Some(expected));
    }

    // None returns
    #[rstest]
    #[case(true)]
    #[case(false)]
    fn test_read_no_call_contributes_to_mismatch_number(#[case] use_cache: bool) {
        let sample_barcodes = vec!["AAATTT", "AGAGGG", "CACAGG"];
        let mut matcher = BarcodeMatcher::new(&sample_barcodes, 1, 2, use_cache);
        //                             1 no-call
        //                             |
        //                             | 1 different base
        //                             | |
        //                             v v
        let test_barcode: &[u8; 6] = b"NAGTTT";
        assert_eq!(matcher.assign(test_barcode), None);
    }

    #[rstest]
    #[case(true)]
    #[case(false)]
    fn test_produce_no_match_if_too_many_mismatches(#[case] use_cache: bool) {
        let sample_barcodes = vec!["AAGCTAG", "CAGCTAG", "GAGCTAG", "TAGCTAG"];
        let assignment_barcode: &[u8] = b"ATCGATC";
        let mut matcher = BarcodeMatcher::new(&sample_barcodes, 0, 100, use_cache);
        assert_eq!(matcher.assign(assignment_barcode), None);
    }

    #[rstest]
    #[case(true)]
    #[case(false)]
    fn test_produce_no_match_if_within_mismatch_delta(#[case] use_cache: bool) {
        let sample_barcodes = vec!["AAAAAAAA", "CCCCCCCC", "GGGGGGGG", "GGGGGGTT"];
        let assignment_barcode: &[u8] = sample_barcodes[3].as_bytes();
        let mut matcher = BarcodeMatcher::new(&sample_barcodes, 100, 3, use_cache);
        assert_eq!(matcher.assign(assignment_barcode), None);
    }

    #[rstest]
    #[case(true)]
    #[case(false)]
    fn test_produce_no_match_if_too_many_mismatches_via_nocalls(#[case] use_cache: bool) {
        let sample_barcodes = vec!["AAAAAAAA", "CCCCCCCC", "GGGGGGGG", "GGGGGGTT"];
        let assignment_barcode: &[u8] = b"GGGGGGTN";
        let mut matcher = BarcodeMatcher::new(&sample_barcodes, 0, 100, use_cache);
        assert_eq!(matcher.assign(assignment_barcode), None);
    }
}
