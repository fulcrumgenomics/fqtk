use super::moving_average::MovingAverage;
use seq_io::fastq::OwnedRecord;

use std::ops::Range;

/// Returns an optional point at which to trim the read due to quality oscillations.
///
/// # Arguments
///
/// * `quals` - The quality scores for the read.
/// * `osc_delta` - The minimum difference in quality scores to be considered an oscillation.
/// * `window_size` - The number of bases to consider in a window.
/// * `max_oscillations` - The maximum number of oscillations allowed in a window.
pub fn identify_trim_point(
    quals: &[u8],
    osc_delta: i32,
    window_size: usize,
    max_oscillations: usize,
) -> Option<usize> {
    // Compute a vector of booleans to represent whether each position is immediately
    // after an oscillation
    let mut is_osc = vec![false; quals.len()];

    for i in 1..quals.len() {
        is_osc[i] = (quals[i] as i32 - quals[i - 1] as i32).abs() >= osc_delta;
    }

    // TODO [brent]: this is ~O(n^2/2) and could be ~O(n)
    // instead, convert is_osc into a vector of indices of oscillations. then can use a
    // window of `max_oscillations` and check if there are at least `window_size` oscillations in the window.
    for i in 1..is_osc.len() {
        if !is_osc[i] {
            continue;
        }
        let mut n = 1;
        #[allow(clippy::needless_range_loop)] // will refactor this later anyway
        for j in (i + 1)..std::cmp::min(i + window_size, quals.len()) {
            if is_osc[j] {
                n += 1;
                if n > max_oscillations {
                    return Some(i);
                }
            }
        }
    }

    None
}

// Indicates which tail(s) to clip.
#[derive(clap::ValueEnum, Debug, Copy, Clone)]
pub enum Tail {
    Left,
    Right,
    Both,
}

/// Uses a moving average to return a range of high quality bases.
/// If all bases are high-quality, the range is the full read.
pub fn find_high_quality_bases(
    bqs: &[u8],
    min_quality: u8,
    window: u8,
    tail: Tail,
) -> Range<usize> {
    let mut left = 0;
    let mut right = bqs.len();
    if matches!(tail, Tail::Left | Tail::Both) {
        let mut ma = MovingAverage::<u8>::new(window as usize);
        for &bq in bqs {
            let mean = ma.push(bq);
            if mean >= min_quality as f64 {
                break;
            }
            left += 1;
        }
    }
    if matches!(tail, Tail::Right | Tail::Both) {
        let mut ma = MovingAverage::<u8>::new(window as usize);
        for &bq in bqs.iter().rev() {
            let mean = ma.push(bq);
            if mean >= min_quality as f64 {
                break;
            }
            right -= 1;
        }
    }
    left..right
}

/// Hard clip the Record to the given range.
pub fn clip_read(record: &mut OwnedRecord, range: Range<usize>) {
    if range.start == 0 && range.end == record.seq.len() {
        return;
    }
    record.seq = record.seq[range.clone()].to_vec();
    record.qual = record.qual[range].to_vec();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_hq_all() {
        let bqs = b"IIIIIIII";
        let range = find_high_quality_bases(bqs, 'I' as u8, 3, Tail::Both);
        assert_eq!(range, 0..bqs.len());
    }

    #[test]
    fn test_find_hq_ends() {
        let bqs = b"EIIIIIIE";
        let range = find_high_quality_bases(bqs, 'I' as u8, 1, Tail::Both);
        assert_eq!(range, 1..bqs.len() - 1);

        let bqs = b"EIIIIIIE";
        let range = find_high_quality_bases(bqs, 'I' as u8, 1, Tail::Left);
        assert_eq!(range, 1..bqs.len());

        let bqs = b"EIIIIIIE";
        let range = find_high_quality_bases(bqs, 'I' as u8, 1, Tail::Right);
        assert_eq!(range, 0..bqs.len() - 1);
    }

    #[test]
    fn test_identify_trim_point() {
        let quals = vec![33; 100];
        let oscillation = 10;
        let window_size = 20;
        let max_oscillations = 3;

        let trim_point = identify_trim_point(&quals, oscillation, window_size, max_oscillations);

        assert_eq!(trim_point, None);
    }

    #[test]
    fn test_identify_trim_point_not_quite_big_enough_changes() {
        let base = [9u8, 18, 27, 36, 27, 18];
        let quals: Vec<u8> = base.into_iter().cycle().take(base.len() * 20).collect();

        let oscillation = 10;
        let window_size = 20;
        let max_oscillations = 3;

        let trim_point = identify_trim_point(&quals, oscillation, window_size, max_oscillations);
        assert_eq!(trim_point, None);
    }

    #[test]
    fn test_identify_trim_point_finds_point() {
        let left = vec![30u8; 50];
        let right = [15u8, 22, 35, 20, 32];
        let quals =
            [&left[..], &right.into_iter().cycle().take(right.len() * 5).collect::<Vec<u8>>()]
                .concat();

        let trim_point = identify_trim_point(&quals, 10, 20, 3);
        assert_eq!(trim_point, Some(50));
    }
}
