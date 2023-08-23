use seq_io::fastq::OwnedRecord;

fn complement(base: u8) -> u8 {
    match base {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        _ => b'N',
    }
}

/// Overlap describes the overlap between two reads.
/// if `adapter` is then `shift` is the amount of adapter sequence at ends of the reads.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub struct PairOverlap {
    // shift required to get read1 to overlap with read2
    pub shift: usize,
    // length of the overlap
    pub overlap: usize,
    // number of single-base differences within the overlap
    pub diffs: usize,
    // true if the end of read1 is "aligned" right of end of read2.
    pub adapter: bool,
}

impl PairOverlap {
    /// Corrects bases and adjusts base-qualities in the pair given the
    /// previously calculated overlap.
    /// Only bases where the base-quality difference between the reads is higher than `bq_delta` are corrected.
    /// (currently) Base-qualities are summed for bases that agree.
    /// If mask_quality is None, bases are hard-clipped,
    ///  otherwise, they are set to the given quality.
    /// returns the number of corrections performed in read1 and read2 respectively.
    pub fn correct(
        &self,
        r1: &mut OwnedRecord,
        r2: &mut OwnedRecord,
        bq_delta: u8,
        mask_quality: Option<u8>,
    ) -> (usize, usize) {
        // NOTE: copying everything for now. This is not ideal but simplifies code.
        // Only applies to subset of reads that actually overlap.
        let mut r1_seq = r1.seq.clone();
        let mut r2_seq = reverse_complement(&r2.seq);
        let mut r1_quals = r1.qual.clone();
        let mut r2_quals = r2.qual.clone();
        r2_quals.reverse();

        if self.adapter {
            std::mem::swap(&mut r1_seq, &mut r2_seq);
            std::mem::swap(&mut r1_quals, &mut r2_quals);
        }
        let mut correction_counts = (0usize, 0);

        for i in 0..self.overlap {
            let agree = r1_seq[i + self.shift] == r2_seq[i];

            // when bases agree, set both to the sum of the qualities.
            if agree {
                // avoid overflow and cap bq to 127 (last ascii character)
                r2_quals[i] = (r1_quals[i + self.shift] as u16 + r2_quals[i] as u16).min(126) as u8;
                r1_quals[i + self.shift] = r2_quals[i];
            } else {
                // when bases disagree, set the lower quality base to the higher quality base.
                if r1_quals[i + self.shift] >= r2_quals[i] + bq_delta {
                    r2_seq[i] = r1_seq[i + self.shift];
                    r2_quals[i] = r1_quals[i + self.shift];
                    correction_counts.1 += 1;
                } else if r2_quals[i] >= r1_quals[i + self.shift] + bq_delta {
                    r1_seq[i + self.shift] = r2_seq[i];
                    r1_quals[i + self.shift] = r2_quals[i];
                    correction_counts.0 += 1;
                }
            }
        }

        if self.adapter {
            std::mem::swap(&mut r1_seq, &mut r2_seq);
            std::mem::swap(&mut r1_quals, &mut r2_quals);
            correction_counts = (correction_counts.1, correction_counts.0);
            if let Some(maskq) = mask_quality {
                // mask qualities
                let l = r1_quals.len();
                for q in &mut r1_quals[l - self.shift..] {
                    *q = maskq;
                }
                for q in &mut r2_quals[..self.shift] {
                    *q = maskq;
                }
            } else {
                // clip right end of r1
                r1_seq = r1_seq[..r1_seq.len() - self.shift].to_vec();
                r1_quals = r1_quals[..r1_quals.len() - self.shift].to_vec();
                // we're still reversed here, so we also clip of left end of r2
                r2_seq = r2_seq[self.shift..].to_vec();
                r2_quals = r2_quals[self.shift..].to_vec();
            }
        }

        r1.seq = r1_seq;
        r1.qual = r1_quals;

        r2.seq = reverse_complement(&r2_seq);
        r2.qual = r2_quals.iter().rev().copied().collect::<Vec<u8>>();
        correction_counts
    }
}

fn reverse_complement(s: &[u8]) -> Vec<u8> {
    s.iter().rev().map(|x| complement(*x)).collect::<Vec<u8>>()
}

/// find overlap given the sequence of a pair of reads.
/// max_overlap_error_rate of e.g. 0.1 means that we allow 1 error per 10 bases of overlap.
pub fn find_overlap(
    s1: &[u8],
    s2: &[u8],
    min_overlap: usize,
    max_overlap_error_rate: f64,
) -> Option<PairOverlap> {
    // note that it's possible to do this without allocating, but then we must (re)complement many times.
    let s2 = s2.iter().rev().map(|x| complement(*x)).collect::<Vec<u8>>();

    'shift_loop: for shift in 0..s1.len().max(min_overlap) - min_overlap {
        let mut diffs = 0;
        let overlap_len = (s1.len().max(shift) - shift).min(s2.len());
        let max_diffs = (overlap_len as f64 * max_overlap_error_rate).floor() as usize;
        for i in 0..overlap_len {
            if s1[i + shift] != s2[i] {
                diffs += 1;
                if diffs > max_diffs {
                    continue 'shift_loop;
                }
            }
        }
        // since were here, we have a large enough overlap
        return Some(PairOverlap { shift, overlap: overlap_len, diffs, adapter: false });
    }

    // now handle case where there is adapter sequenced.
    // shift is in s2 this time.
    'shift_loop: for shift in 0..s2.len() - min_overlap {
        let mut diffs = 0;
        let overlap_len = (s2.len() - shift).min(s1.len());
        let max_diffs = (overlap_len as f64 * max_overlap_error_rate).floor() as usize;
        for i in 0..overlap_len {
            if s1[i] != s2[i + shift] {
                diffs += 1;
                if diffs > max_diffs {
                    continue 'shift_loop;
                }
            }
        }
        // since were here, we have a large enough overlap
        return Some(PairOverlap { shift, overlap: overlap_len, diffs, adapter: true });
    }

    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_overlap() {
        // shift: 0
        //  i:  0123456789
        // s1:  ACGAAAAA
        // s2:  AAAAAGGC
        // shift: 3
        // i+shift    34567
        // s1:     ACGAAAAA
        // i:         01234567
        // s2:        AAAAAGGC
        let r1 = b"ACGAAAAA";
        let r2 = b"AAAAAGGC".iter().rev().map(|x| complement(*x)).collect::<Vec<u8>>();
        let overlap = find_overlap(r1, &r2, 4, 0f64);
        assert_eq!(overlap, Some(PairOverlap { shift: 3, overlap: 5, diffs: 0, adapter: false }));
    }

    #[test]
    fn test_overlap_with_adapters() {
        // sequenced adapter
        // we know this is the case when end of read1 is beyond end of read2 (here beyond == to the right of)
        // shift: 0
        //  i:  0123456789
        // s1:  AAAAAGGC
        // s2:  ACGAAAAA
        // shift: 3
        // s1:     AAAAAGGC
        // s2:  ACGAAAAA
        let r1 = b"AAAAAGGC";
        let r2 = b"ACGAAAAA".iter().rev().map(|x| complement(*x)).collect::<Vec<u8>>();
        let overlap = find_overlap(r1, &r2, 4, 0f64);
        assert_eq!(overlap, Some(PairOverlap { shift: 3, overlap: 5, diffs: 0, adapter: true }));
    }

    #[test]
    fn test_overlap_correction() {
        let overlap = PairOverlap { shift: 3, overlap: 5, diffs: 0, adapter: false };
        let mut r1 = OwnedRecord {
            head: b"r1".to_vec(),
            seq: b"ACGAAATA".to_vec(),
            qual: b"IIIIIIII".to_vec(),
        };
        let mut r2 = OwnedRecord {
            head: b"r2".to_vec(),
            seq: reverse_complement(b"AAAAAGGC"),
            qual: b"IIIIIIII".to_vec(),
        };

        overlap.correct(&mut r1, &mut r2, 5, None);

        // the sequences are different, but not corrected becase BQs are the same.
        assert_eq!(r1.seq, b"ACGAAATA".to_vec());
        assert_eq!(r2.seq, reverse_complement(b"AAAAAGGC"));

        // now adjust r2 to have higher BQs than r1 for that base.
        // ACGAAATA
        //    AAAAAGGC
        //    IIIQIIII

        r2.qual = b"IIIQIIII".iter().rev().map(|b| *b).collect::<Vec<u8>>();
        overlap.correct(&mut r1, &mut r2, 5, None);
        assert_eq!(r1.seq, b"ACGAAAAA".to_vec());
    }

    #[test]
    fn test_overlap_correction_with_adapter() {
        let overlap = PairOverlap { shift: 3, overlap: 5, diffs: 0, adapter: true };
        //       IIIQIIII
        // r1:   AAAAAGGC
        // r2: ACGAAATA

        let mut r1 = OwnedRecord {
            head: b"r1".to_vec(),
            seq: b"AAAAAGGC".to_vec(),
            qual: b"IIIQIIII".to_vec(),
        };
        let mut r2 = OwnedRecord {
            head: b"r2".to_vec(),
            seq: reverse_complement(b"ACGAAATA").to_vec(),
            qual: b"IIIIIIII".to_vec(),
        };

        overlap.correct(&mut r1, &mut r2, 5, None);

        // sequence 2 is corrected.
        assert_eq!(r2.seq, reverse_complement(b"AAAAA"));
        assert_eq!(r1.seq, b"AAAAA");

        assert_eq!(r2.qual, b"~Q~~~");
        assert_eq!(r1.qual, b"~~~Q~");
    }
}
