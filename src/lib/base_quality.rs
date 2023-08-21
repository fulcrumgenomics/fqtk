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
    // collect indices of oscillations where quals[i] - quals[i-1] >= osc_delta
    let osc = quals
        .windows(2)
        .enumerate()
        .filter(|(_, w)| (w[1] as i32 - w[0] as i32).abs() >= osc_delta)
        .map(|(i, _)| i + 1)
        .collect::<Vec<_>>();
    // NOTE [performance]: we can do this without allocating but the .windows() method used below is only
    // available on a slice and this only allocates for as many oscillations as are found.

    // here we have indices of oscillations. e.g.: [50, 52, 53, 54, 55, 57, 58, 59, 60, 62, 63, 64, 65, 67, 68, 69, 70, 72, 73, 74]

    // use a window of `max_oscilations` and check if there are at least `max_oscillations` oscillations in the window.
    osc.windows(max_oscillations).find(|w| 
        // given a e.g. [50, 52, 55] means we found 3 oscillations in 5 bases (55 - 50)
        // and we had window_size to find that many. so if the last - first < window_size, we found an osc window.
        w[w.len() - 1] - w[0] < window_size
    ).map(|w| w[0])
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
    #[test]
    fn test_identify_trim_point_finds_point_within() {
        let left = vec![20u8; 50];
        let right = [15u8, 22, 35, 20, 32];
        let quals =
            [&left[..], &right.into_iter().cycle().take(right.len() * 5).collect::<Vec<u8>>()]
                .concat();

        let trim_point = identify_trim_point(&quals, 10, 20, 3);
        assert_eq!(trim_point, Some(52));
    }
}
