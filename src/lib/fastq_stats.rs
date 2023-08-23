use std::fmt;
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
    pub(crate) overlap_hist: Vec<usize>,
    pub(crate) errors_corrected: [usize; 2],
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
            overlap_hist: vec![0; DEFAULT_LEN],
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

    pub fn update_bases_overlap(&mut self, overlap: usize) {
        if overlap >= self.overlap_hist.len() {
            self.overlap_hist.resize(overlap + 1, 0);
        }
        self.overlap_hist[overlap] += 1;
    }
}

impl fmt::Display for Stats {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "Read-length histogram:")?;
        writeln!(f, "length\tr1_pre\tr1_post\tr2_pre\tr2_post")?;
        let max_len = self.length_hist[0]
            .pre
            .len()
            .max(self.length_hist[0].post.len())
            .max(self.length_hist[1].pre.len().max(self.length_hist[1].post.len()));

        for i in 0..max_len {
            let pre1 = self.length_hist[0].pre.get(i).unwrap_or(&0);
            let post1 = self.length_hist[0].post.get(i).unwrap_or(&0);
            let pre2 = self.length_hist[1].pre.get(i).unwrap_or(&0);
            let post2 = self.length_hist[1].post.get(i).unwrap_or(&0);
            writeln!(f, "{}\t{}\t{}\t{}\t{}", i, pre1, post1, pre2, post2)?;
        }

        writeln!(f, "Overlap length histogram:")?;
        for (i, v) in self.overlap_hist.iter().enumerate() {
            if *v > 0 {
                writeln!(f, "{}\t{}", i, v)?;
            }
        }
        writeln!(f, "Errors corrected:")?;
        writeln!(f, "Read 1: {}", self.errors_corrected[0])?;
        writeln!(f, "Read 2: {}", self.errors_corrected[1])?;
        Ok(())
    }
}
