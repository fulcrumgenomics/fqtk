use super::moving_average::MovingAverage;

use std::ops::Range;

pub(crate) fn find_oscillating_quals(bqs: &[u8]) -> Range<usize> {
    return 0..0;
}

pub(crate) enum Tail {
    Left,
    Right,
    Both,
}

/// Uses a moving average to return a range of high quality bases.
/// If all bases are high-quality, the range is the full read.
pub(crate) fn find_high_quality_bases(
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
}
