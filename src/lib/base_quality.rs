use super::moving_average::MovingAverage;
use seq_io::fastq::OwnedRecord;

use std::ops::Range;

pub fn find_oscillating_quals(_bqs: &[u8]) -> Range<usize> {
    unimplemented!("find_oscillating_quals");
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
}
