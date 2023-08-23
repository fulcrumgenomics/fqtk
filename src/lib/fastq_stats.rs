use std::fmt;

use crate::pair_overlap::PairOverlap;
pub struct Stat<T> {
    pub(crate) pre: T,
    pub(crate) post: T,
}

pub enum When {
    Pre,
    Post,
}

pub enum ReadI {
    R1 = 0,
    R2 = 1,
}

impl From<usize> for ReadI {
    fn from(i: usize) -> Self {
        match i {
            0 => Self::R1,
            1 => Self::R2,
            _ => panic!("invalid read index"),
        }
    }
}

/// Stats related to pair overlaps.
pub struct OverlapStats {
    overlap_hist: Vec<usize>,
    // total bases of overlap
    bases: usize,
    // mismatches in the overlap
    differences: usize,
    // corrections for ead read
    corrections: [usize; 2],
}

impl Default for OverlapStats {
    fn default() -> Self {
        Self::new()
    }
}

impl OverlapStats {
    /// create a new OverlapStats struct.
    pub fn new() -> Self {
        Self { overlap_hist: vec![0; DEFAULT_LEN], bases: 0, differences: 0, corrections: [0, 0] }
    }

    /// update stats given the PairOverlap.
    pub fn update(&mut self, po: PairOverlap) {
        self.update_overlap(po.overlap);
        self.differences += po.diffs;
    }

    /// update with the amount of overlap.
    fn update_overlap(&mut self, overlap: usize) {
        if overlap >= self.overlap_hist.len() {
            self.overlap_hist.resize(overlap + 1, 0);
        }
        self.overlap_hist[overlap] += 1;
        self.bases += overlap;
    }

    /// update with corrections between the overlapping reads.
    pub fn update_corrections(&mut self, n: usize, read: ReadI) {
        self.corrections[read as usize] += n;
    }
}

/// Stats for each pair before and after processing
pub struct Stats {
    // r1, r2
    pub(crate) length_hist: [Stat<Vec<usize>>; 2],
    pub overlap_stats: OverlapStats,
    pub(crate) errors_corrected: [usize; 2],
    // number of oscillations in each read.
    pub(crate) oscillations: [usize; 2],

    // number of reads filtered by length
    pub(crate) length_filtered: usize,
    max_read_length: usize,
}

const DEFAULT_LEN: usize = 400;

impl Default for Stats {
    fn default() -> Self {
        Self::new()
    }
}

impl Stats {
    /// create a new Stats struct.
    pub fn new() -> Self {
        Self {
            length_hist: [
                Stat { pre: vec![0; DEFAULT_LEN], post: vec![0; DEFAULT_LEN] },
                Stat { pre: vec![0; DEFAULT_LEN], post: vec![0; DEFAULT_LEN] },
            ],
            overlap_stats: OverlapStats::default(),
            errors_corrected: [0, 0],
            oscillations: [0, 0],
            length_filtered: 0,
            max_read_length: 0,
        }
    }
    pub fn update_oscillations(&mut self, n: usize, read: ReadI) {
        self.oscillations[read as usize] += n;
    }

    pub fn increment_length_filter(&mut self) {
        self.length_filtered += 1;
    }

    /// update length histogram for a read.
    pub fn update_length(&mut self, length: usize, p: When, read: ReadI) {
        let index = read as usize;
        self.max_read_length = self.max_read_length.max(length);
        match p {
            When::Pre => {
                if length >= self.length_hist[index].pre.len() {
                    self.length_hist[index].pre.resize(length + 1, 0);
                }
                self.length_hist[index].pre[length] += 1;
            }
            When::Post => {
                if length >= self.length_hist[index].post.len() {
                    self.length_hist[index].post.resize(length + 1, 0);
                }
                self.length_hist[index].post[length] += 1;
            }
        }
    }

    pub fn update_errors_corrected(&mut self, n: usize, read: ReadI) {
        self.errors_corrected[read as usize] += n;
    }
}

impl fmt::Display for OverlapStats {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "Overlap length histogram:\noverlap\tcount")?;
        for (i, v) in self.overlap_hist.iter().enumerate() {
            if *v > 0 {
                writeln!(f, "{}\t{}", i, v)?;
            }
        }
        writeln!(f, "Differences in overlapped bases: {}", self.differences)?;

        writeln!(f, "Errors corrected:")?;
        writeln!(f, "Read 1: {}", self.corrections[0])?;
        writeln!(f, "Read 2: {}", self.corrections[1])
    }
}

impl fmt::Display for Stats {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "Oscillations: R1={} R2={}", self.oscillations[0], self.oscillations[1])?;
        writeln!(f, "Reads Filtered based on length: {}", self.length_filtered)?;
        writeln!(f, "Read-length histogram:")?;
        writeln!(f, "length\tr1_pre\tr1_post\tr2_pre\tr2_post")?;
        for i in 0..self.max_read_length + 1 {
            let pre1 = self.length_hist[0].pre.get(i).unwrap_or(&0);
            let post1 = self.length_hist[0].post.get(i).unwrap_or(&0);
            let pre2 = self.length_hist[1].pre.get(i).unwrap_or(&0);
            let post2 = self.length_hist[1].post.get(i).unwrap_or(&0);
            writeln!(f, "{}\t{}\t{}\t{}\t{}", i, pre1, post1, pre2, post2)?;
        }
        self.overlap_stats.fmt(f)
    }
}
