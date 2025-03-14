use std::cmp::min;

use crate::decode;
use crate::encode;

use super::byte_is_nocall;
use super::samples::Sample;
use crate::bitenc::BitEnc;
use rustc_hash::FxHashMap;

const STARTING_CACHE_SIZE: usize = 1_000_000;

/// The struct that contains the info related to the best and next best sample barcode match.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct BarcodeMatch {
    /// Index of the best barcode match in the corresponding ``BarcodeMatcher`` struct that
    /// generated this match.
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
    /// Vec storing each sample
    samples: Vec<Sample>,
    /// Vec of the barcodes for each sample
    sample_barcodes: Vec<BitEnc>,
    /// The maxium number of Ns in any barcode in set of sample barcodes
    max_ns_in_barcodes: usize,
    /// The maximum mismatches to match a sample barcode.
    max_mismatches: u8,
    /// The minimum difference between number of mismatches in the best and second best barcodes
    /// for a barcode to be considered a match.
    min_mismatch_delta: u8,
    /// Cache of previous matches; `None` disables caching, `Some(map)` enables it.
    cache: Option<FxHashMap<Vec<u8>, BarcodeMatch>>,
}

impl BarcodeMatcher {
    /// Instantiates a new ``BarcodeMatcher`` struct. Checks that the samples vector is not
    /// empty and that none of the barcodes provided are the empty string.
    ///
    /// # Panics
    /// - Will panic if provided an empty vec of samples.
    /// - Will panic if any provided barcode is length zero.
    #[must_use]
    pub fn new(
        samples: &[Sample],
        max_mismatches: u8,
        min_mismatch_delta: u8,
        use_cache: bool,
    ) -> Self {
        assert!(!samples.is_empty(), "Must provide at least one sample");
        assert!(
            samples.iter().all(|b| !b.barcode.is_empty()),
            "Sample barcode cannot be empty string"
        );

        let mut max_ns_in_barcodes = 0;
        let mut modified_samples = samples.to_vec();
        let mut sample_barcodes = Vec::with_capacity(samples.len());
        for sample in &mut modified_samples {
            sample.barcode = sample.barcode.to_ascii_uppercase();
            let bytes = sample.barcode.as_bytes();
            let num_ns: usize = bytes.iter().filter(|&&b| byte_is_nocall(b)).count();
            max_ns_in_barcodes = max_ns_in_barcodes.max(num_ns);
            sample_barcodes.push(encode(bytes));
        }
        let cache = use_cache
            .then(|| FxHashMap::with_capacity_and_hasher(STARTING_CACHE_SIZE, Default::default()));
        Self {
            samples: modified_samples,
            sample_barcodes,
            max_ns_in_barcodes,
            max_mismatches,
            min_mismatch_delta,
            cache,
        }
    }

    /// Counts the number of bases that differ between two equal-length encoded byte arrays.
    /// Callers must ensure lengths match (validated once in [`Self::assign_internal`]).
    fn count_mismatches(
        observed_bases: &BitEnc,
        expected_bases: &BitEnc,
        max_mismatches: u8,
    ) -> u8 {
        let count = observed_bases.hamming(expected_bases, u32::from(max_mismatches));
        u8::try_from(count).expect("Overflow on number of mismatch bases")
    }

    /// Returns the expected barcode length, assuming a fixed length for all samples.
    fn expected_barcode_length(&self) -> usize {
        self.samples[0].barcode.len()
    }

    /// Assigns the barcode that best matches the provided ``read_bases``.
    #[must_use]
    fn assign_internal(&self, read_bases: &[u8]) -> Option<BarcodeMatch> {
        let mut best_barcode_index = self.samples.len();
        let mut best_mismatches = 255u8;
        let mut next_best_mismatches = 255u8;
        let mut max_mismatches = 255u8;
        let read_bases = encode(read_bases); // NB: this encodes IUPAC bases in the read, but count_mismatches will treat them as no-calls.
        let expected_len = self.expected_barcode_length();
        assert!(
            read_bases.nr_symbols() == expected_len,
            "Read barcode ({}) length ({}) differs from expected barcode length ({})",
            decode(&read_bases),
            read_bases.nr_symbols(),
            expected_len,
        );
        for (index, sample_barcode) in self.sample_barcodes.iter().enumerate() {
            let mismatches = Self::count_mismatches(&read_bases, sample_barcode, max_mismatches);
            if mismatches < best_mismatches {
                next_best_mismatches = best_mismatches;
                best_mismatches = mismatches;
                best_barcode_index = index;
                if next_best_mismatches < 255u8 - self.min_mismatch_delta {
                    max_mismatches =
                        min(max_mismatches, next_best_mismatches + self.min_mismatch_delta);
                }
            } else if mismatches < next_best_mismatches {
                next_best_mismatches = mismatches;
                if next_best_mismatches < 255u8 - self.min_mismatch_delta {
                    max_mismatches =
                        min(max_mismatches, next_best_mismatches + self.min_mismatch_delta);
                }
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
            return None;
        }
        if let Some(cache) = &self.cache {
            if let Some(cached_match) = cache.get(read_bases) {
                return Some(*cached_match);
            }
        }
        let maybe_match = self.assign_internal(read_bases);
        if let (Some(cache), Some(m)) = (self.cache.as_mut(), maybe_match) {
            cache.insert(read_bases.to_vec(), m);
        }
        maybe_match
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rstest::rstest;

    /// Given a barcode and integer identifier, generate a sample instance. This helper function
    /// allows creating more succinct testing code.
    fn barcode_to_sample(barcode: &str, idx: usize) -> Sample {
        Sample::new(idx, format!("sample_{idx}"), barcode.to_string())
    }

    /// Create a vector of samples from a list of barcodes, for more succinct testing code.
    fn barcodes_to_samples(barcodes: &[&str]) -> Vec<Sample> {
        barcodes
            .iter()
            .enumerate()
            .map(|(idx, barcode)| barcode_to_sample(barcode, idx))
            .collect::<Vec<_>>()
    }

    /// Helper for running ``BarcodeMatcher::count_mismatches``.  Encodes the bases.
    fn count_mismatches(observed_bases: &str, expected_bases: &str) -> u8 {
        BarcodeMatcher::count_mismatches(
            &encode(observed_bases.as_bytes()),
            &encode(expected_bases.as_bytes()),
            255,
        )
    }

    // ############################################################################################
    // Test ``BarcodeMatcher`` instantiation
    // ############################################################################################
    #[rstest]
    #[case(true)]
    #[case(false)]
    fn test_barcode_matcher_instantiation_can_succeed(#[case] use_cache: bool) {
        let samples = barcodes_to_samples(&["ACGT"]);
        let _matcher = BarcodeMatcher::new(&samples, 2, 1, use_cache);
    }

    #[rstest]
    #[case(true)]
    #[case(false)]
    #[should_panic(expected = "Must provide at least one sample")]
    fn test_barcode_matcher_fails_if_no_samples_provided(#[case] use_cache: bool) {
        let samples = barcodes_to_samples(&[]);
        let _matcher = BarcodeMatcher::new(&samples, 2, 1, use_cache);
    }

    // ############################################################################################
    // Test BarcodeMatcher::count_mismatches
    // ############################################################################################

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
        assert_eq!(count_mismatches("ACGTTAAACCGAAACA", "ACGTUMRWSYKVHDBN"), 0,);
        // IUPAC bases are mismatches in the observed barcodes
        assert_eq!(count_mismatches("ACGTUMRWSYKVHDBN", "ACGTTAAACCGAAACA"), 11,);
    }

    #[test]
    fn count_mismatches_iupac_bases_assymetry() {
        // if the observed base is an N, it will not match anything but an N
        assert_eq!(count_mismatches("N", "R"), 1,);
        assert_eq!(count_mismatches("N", "N"), 0,);
        // if the observed base is an R, it will match R, V, D, and N
        assert_eq!(count_mismatches("R", "R"), 0,);
        assert_eq!(count_mismatches("R", "V"), 0,);
        assert_eq!(count_mismatches("R", "D"), 0,);
        assert_eq!(count_mismatches("R", "N"), 0,);
        assert_eq!(count_mismatches("R", "B"), 1,);
    }

    #[test]
    #[should_panic(
        expected = "Read barcode (GATTAGA) length (7) differs from expected barcode length (6)"
    )]
    fn assign_panics_when_read_is_longer_than_expected_barcode() {
        let samples = barcodes_to_samples(&["CTATGT"]);
        let mut matcher = BarcodeMatcher::new(&samples, 1, 2, false);
        let _ = matcher.assign(b"GATTAGA");
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
        let samples = barcodes_to_samples(&["ACGT", "AAAG", "CACA"]);
        let mut matcher = BarcodeMatcher::new(&samples, 2, 2, use_cache);
        assert_eq!(
            matcher.assign(samples[EXPECTED_BARCODE_INDEX].barcode.as_bytes()),
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
        let samples = barcodes_to_samples(&["AAAT", "AGAG", "CACA"]);
        let mut matcher = BarcodeMatcher::new(&samples, 2, 2, use_cache);
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
        let samples = barcodes_to_samples(&["AAAT", "AGAG", "CACA"]);
        let mut matcher = BarcodeMatcher::new(&samples, 2, 2, use_cache);
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
        let samples = barcodes_to_samples(&["AAATTT", "AGAGGG", "CACAGG"]);
        let mut matcher = BarcodeMatcher::new(&samples, 2, 2, use_cache);
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
        let samples = barcodes_to_samples(&["NAGTTT", "AGAGGG", "CACAGG"]);
        let mut matcher = BarcodeMatcher::new(&samples, 1, 2, use_cache);
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
        let samples = barcodes_to_samples(&["AAATTT", "AGAGGG", "CACAGG"]);
        let mut matcher = BarcodeMatcher::new(&samples, 1, 2, use_cache);
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
        let samples = barcodes_to_samples(&["AAGCTAG", "CAGCTAG", "GAGCTAG", "TAGCTAG"]);
        let assignment_barcode: &[u8] = b"ATCGATC";
        let mut matcher = BarcodeMatcher::new(&samples, 0, 100, use_cache);
        assert_eq!(matcher.assign(assignment_barcode), None);
    }

    #[rstest]
    #[case(true)]
    #[case(false)]
    fn test_produce_no_match_if_within_mismatch_delta(#[case] use_cache: bool) {
        let samples = barcodes_to_samples(&["AAAAAAAA", "CCCCCCCC", "GGGGGGGG", "GGGGGGTT"]);
        let assignment_barcode: &[u8] = samples[3].barcode.as_bytes();
        let mut matcher = BarcodeMatcher::new(&samples, 100, 3, use_cache);
        assert_eq!(matcher.assign(assignment_barcode), None);
    }

    #[rstest]
    #[case(true)]
    #[case(false)]
    fn test_produce_no_match_if_too_many_mismatches_via_nocalls(#[case] use_cache: bool) {
        let samples = barcodes_to_samples(&["AAAAAAAA", "CCCCCCCC", "GGGGGGGG", "GGGGGGTT"]);
        let assignment_barcode: &[u8] = b"GGGGGGTN";
        let mut matcher = BarcodeMatcher::new(&samples, 0, 100, use_cache);
        assert_eq!(matcher.assign(assignment_barcode), None);
    }
}
