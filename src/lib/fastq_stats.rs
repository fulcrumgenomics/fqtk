use std::fmt;
use std::ops::Index;
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

/// Stats for each pair before and after processing
pub struct Stats {
    // r1, r2
    pub(crate) length_hist: [Stat<Vec<usize>>; 2],
    pub(crate) overlap_hist: Stat<Vec<usize>>,
    pub(crate) errors_corrected: [usize; 2],
}

impl Stats {
    /// create a new Stats struct.
    pub fn new() -> Self {
        Self {
            length_hist: [
                Stat { pre: vec![0; 1000], post: vec![0; 1000] },
                Stat { pre: vec![0; 1000], post: vec![0; 1000] },
            ],
            overlap_hist: Stat { pre: vec![0; 1000], post: vec![0; 1000] },
            errors_corrected: [0, 0],
        }
    }

    /// update length histogram for a read.
    pub fn update_length(&mut self, length: usize, p: When, read: ReadI) {
        let index = read as usize;
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

    pub fn update_bases_overlap(&mut self, overlap: usize, p: When) {
        match p {
            When::Pre => {
                if overlap >= self.overlap_hist.pre.len() {
                    self.overlap_hist.pre.resize(overlap + 1, 0);
                }
                self.overlap_hist.pre[overlap] += 1;
            }
            When::Post => {
                if overlap >= self.overlap_hist.post.len() {
                    self.overlap_hist.post.resize(overlap + 1, 0);
                }
                self.overlap_hist.post[overlap] += 1;
            }
        }
    }
}

impl fmt::Display for Stats {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Read 1 length histogram:\n")?;
        for (i, (pre, post)) in
            self.length_hist[0].pre.iter().zip(self.length_hist[0].post.iter()).enumerate()
        {
            if *pre > 0 || *post > 0 {
                write!(f, "{}\t{}\t{}\n", i, pre, post)?;
            }
        }
        write!(f, "Read 2 length histogram:\n")?;
        for (i, (pre, post)) in
            self.length_hist[1].pre.iter().zip(self.length_hist[1].post.iter()).enumerate()
        {
            if *pre > 0 || *post > 0 {
                write!(f, "{}\t{}\t{}\n", i, pre, post)?;
            }
        }
        write!(f, "Overlap length histogram:\n")?;
        for (i, (pre, post)) in
            self.overlap_hist.pre.iter().zip(self.overlap_hist.post.iter()).enumerate()
        {
            if *pre > 0 || *post > 0 {
                write!(f, "{}\t{}\t{}\n", i, pre, post)?;
            }
        }
        write!(f, "Errors corrected:\n")?;
        write!(f, "Read 1: {}\n", self.errors_corrected[0])?;
        write!(f, "Read 2: {}\n", self.errors_corrected[1])?;
        Ok(())
    }
}
