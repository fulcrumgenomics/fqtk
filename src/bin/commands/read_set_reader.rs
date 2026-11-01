//! Lockstep reading of synchronized FASTQ inputs (e.g. R1 and R2), one record from each input at
//! a time, shared by the subcommands that process paired or multi-read data.

use crate::commands::fastq_readahead::{ChunkReader, ReadAheadBuilder, open_readers};
use anyhow::{Result, bail};
use memchr::memchr2;
use seq_io::fastq::{Error as FastqError, Reader as FastqReader, Record, RefRecord};
use std::io;
use std::path::PathBuf;

/// Reads synchronized FASTQ inputs one read set at a time: the next record from each input.
///
/// Errors if the inputs end at different points and, when name checking is on, if the records in
/// a read set don't share a [`base_read_name`].  Errors name the inputs responsible.
///
/// Not an [`Iterator`]: each read set borrows the readers' buffers, so it must be dropped before
/// the next one is read.
pub(crate) struct ReadSetReader {
    /// One reader per input, in the order the inputs were given.
    readers: Vec<FastqReader<ChunkReader>>,
    /// The path of each input, for error messages.
    paths: Vec<PathBuf>,
    /// Whether to check that the records in each read set share a name.
    check_names: bool,
    /// Number of read sets returned so far.
    sets_read: u64,
}

impl ReadSetReader {
    /// Opens every input on its own read-ahead thread, failing if any input can't be opened or
    /// decoded (see [`open_readers`]).  `check_names` turns on the read-name check.
    pub(crate) fn open(inputs: Vec<ReadAheadBuilder>, check_names: bool) -> io::Result<Self> {
        let paths = inputs.iter().map(|input| input.path.clone()).collect();
        let readers = open_readers(inputs)?;
        Ok(Self { readers, paths, check_names, sets_read: 0 })
    }

    /// Returns the next read set, one record per input in input order, or `None` once every
    /// input has ended.  Errors if an input can't be read or parsed, if only some inputs have
    /// ended, or, when checking names, if the records' names differ.
    pub(crate) fn next(&mut self) -> Result<Option<Vec<RefRecord<'_>>>> {
        let set_number = self.sets_read + 1;
        let mut records = Vec::with_capacity(self.readers.len());
        let mut ended = Vec::new();
        for (reader, path) in self.readers.iter_mut().zip(&self.paths) {
            match reader.next() {
                Some(Ok(record)) => records.push(record),
                // I/O errors from the read-ahead reader already name the input.
                Some(Err(FastqError::Io(e))) => bail!("Error reading FASTQ input: {e}"),
                Some(Err(e)) => bail!("Error reading FASTQ input {path:?}: {e}"),
                None => ended.push(path),
            }
        }
        if records.is_empty() {
            return Ok(None);
        }
        if !ended.is_empty() {
            bail!(
                "FASTQ inputs are out of sync: {ended:?} ended before read {set_number}, but the \
                 other inputs did not."
            );
        }
        if self.check_names {
            check_read_names(&records, &self.paths, set_number)?;
        }
        self.sets_read = set_number;
        Ok(Some(records))
    }
}

/// Returns the read name portion of a FASTQ header, stripping any trailing `/1` or `/2`
/// paired-end suffix and any comment after whitespace.
pub(crate) fn base_read_name(head: &[u8]) -> &[u8] {
    let name_end = memchr2(b' ', b'\t', head).unwrap_or(head.len());
    let name = &head[..name_end];
    if name.len() >= 2
        && name[name.len() - 2] == b'/'
        && matches!(name[name.len() - 1], b'1' | b'2')
    {
        &name[..name.len() - 2]
    } else {
        name
    }
}

/// Checks that every record in a read set has the same [`base_read_name`] as the first.
/// `paths` and `set_number` (1-based) are used in errors.
fn check_read_names(records: &[RefRecord], paths: &[PathBuf], set_number: u64) -> Result<()> {
    let expected = base_read_name(records[0].head());
    for (record, path) in records.iter().zip(paths).skip(1) {
        let name = base_read_name(record.head());
        if name != expected {
            bail!(
                "Read name mismatch at read {set_number}: {:?} in {:?} but {:?} in {path:?}.",
                String::from_utf8_lossy(expected),
                paths[0],
                String::from_utf8_lossy(name),
            );
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;
    use tempfile::TempDir;

    /// Writes a FASTQ file at `path` with one record per header in `heads`.
    fn write_fastq(path: &Path, heads: &[&str]) -> PathBuf {
        let text = heads.iter().fold(String::new(), |mut text, head| {
            text.push_str(&format!("@{head}\nACGT\n+\nIIII\n"));
            text
        });
        std::fs::write(path, text).unwrap();
        path.to_path_buf()
    }

    /// Opens `paths` with the default chunking.
    fn open(paths: &[PathBuf], check_names: bool) -> ReadSetReader {
        let inputs =
            paths.iter().map(|path| ReadAheadBuilder { path: path.clone(), ..Default::default() });
        ReadSetReader::open(inputs.collect(), check_names).unwrap()
    }

    /// Reads every read set, returning each one's record names.
    fn read_all(reader: &mut ReadSetReader) -> Result<Vec<Vec<String>>> {
        let mut sets = Vec::new();
        while let Some(records) = reader.next()? {
            sets.push(
                records.iter().map(|r| String::from_utf8(r.head().to_vec()).unwrap()).collect(),
            );
        }
        Ok(sets)
    }

    #[test]
    fn each_read_set_has_the_next_record_from_every_input_in_input_order() {
        let dir = TempDir::new().unwrap();
        let r1 = write_fastq(&dir.path().join("r1.fq"), &["a/1", "b/1"]);
        let r2 = write_fastq(&dir.path().join("r2.fq"), &["a/2", "b/2"]);
        let i1 = write_fastq(&dir.path().join("i1.fq"), &["a", "b"]);

        let sets = read_all(&mut open(&[r1, r2, i1], true)).unwrap();
        assert_eq!(sets, [["a/1", "a/2", "a"], ["b/1", "b/2", "b"]]);
    }

    #[test]
    fn empty_inputs_have_no_read_sets() {
        let dir = TempDir::new().unwrap();
        let r1 = write_fastq(&dir.path().join("r1.fq"), &[]);
        let r2 = write_fastq(&dir.path().join("r2.fq"), &[]);

        assert!(read_all(&mut open(&[r1, r2], true)).unwrap().is_empty());
    }

    #[test]
    fn an_input_that_ends_first_is_out_of_sync_and_named() {
        let dir = TempDir::new().unwrap();
        let r1 = write_fastq(&dir.path().join("r1.fq"), &["a", "b", "c"]);
        let r2 = write_fastq(&dir.path().join("r2.fq"), &["a", "b"]);

        let err = read_all(&mut open(&[r1, r2.clone()], true)).unwrap_err().to_string();
        assert!(err.contains("out of sync"), "{err}");
        assert!(err.contains(&format!("{r2:?}")), "error does not name r2: {err}");
        assert!(err.contains("before read 3"), "error does not give the read number: {err}");
    }

    #[test]
    fn the_first_input_ending_first_is_out_of_sync_and_named() {
        let dir = TempDir::new().unwrap();
        let r1 = write_fastq(&dir.path().join("r1.fq"), &["a"]);
        let r2 = write_fastq(&dir.path().join("r2.fq"), &["a", "b"]);

        let err = read_all(&mut open(&[r1.clone(), r2], true)).unwrap_err().to_string();
        assert!(err.contains("out of sync"), "{err}");
        assert!(err.contains(&format!("{r1:?}")), "error does not name r1: {err}");
    }

    #[test]
    fn mismatched_names_are_an_error_naming_the_read_and_inputs() {
        let dir = TempDir::new().unwrap();
        let r1 = write_fastq(&dir.path().join("r1.fq"), &["a/1", "b/1"]);
        let r2 = write_fastq(&dir.path().join("r2.fq"), &["a/2", "c/2"]);

        let err = read_all(&mut open(&[r1.clone(), r2.clone()], true)).unwrap_err().to_string();
        assert!(err.contains("Read name mismatch at read 2"), "{err}");
        assert!(err.contains(&format!("{r1:?}")) && err.contains(&format!("{r2:?}")), "{err}");
    }

    #[test]
    fn mismatched_names_are_read_when_name_checking_is_off() {
        let dir = TempDir::new().unwrap();
        let r1 = write_fastq(&dir.path().join("r1.fq"), &["a"]);
        let r2 = write_fastq(&dir.path().join("r2.fq"), &["b"]);

        assert_eq!(read_all(&mut open(&[r1, r2], false)).unwrap(), [["a", "b"]]);
    }

    #[test]
    fn pair_suffixes_and_comments_are_not_name_differences() {
        let dir = TempDir::new().unwrap();
        let r1 = write_fastq(&dir.path().join("r1.fq"), &["a/1 1:N:0:ACGT"]);
        let r2 = write_fastq(&dir.path().join("r2.fq"), &["a/2\t2:N:0:ACGT"]);

        assert_eq!(read_all(&mut open(&[r1, r2], true)).unwrap().len(), 1);
    }

    #[test]
    fn a_malformed_record_is_an_error_naming_the_input() {
        let dir = TempDir::new().unwrap();
        let r1 = dir.path().join("r1.fq");
        std::fs::write(&r1, "@a\nACGT\n-\nIIII\n").unwrap();

        let err = read_all(&mut open(&[r1.clone()], true)).unwrap_err().to_string();
        assert!(err.contains(&format!("{r1:?}")), "error does not name the input: {err}");
    }

    #[test]
    fn base_read_name_keeps_a_plain_name() {
        assert_eq!(base_read_name(b"read1"), b"read1");
    }

    #[test]
    fn base_read_name_strips_a_slash_1_suffix() {
        assert_eq!(base_read_name(b"read1/1"), b"read1");
    }

    #[test]
    fn base_read_name_strips_a_slash_2_suffix() {
        assert_eq!(base_read_name(b"read1/2"), b"read1");
    }

    #[test]
    fn base_read_name_strips_a_comment() {
        assert_eq!(base_read_name(b"read1 1:N:0:ATCACG"), b"read1");
    }

    #[test]
    fn base_read_name_strips_a_suffix_and_a_comment() {
        assert_eq!(base_read_name(b"read1/1 1:N:0:ATCACG"), b"read1");
    }

    #[test]
    fn base_read_name_strips_a_tab_separated_comment() {
        assert_eq!(base_read_name(b"read1\tBC:Z:ACGT"), b"read1");
    }

    #[test]
    fn base_read_name_keeps_a_slash_3_suffix() {
        assert_eq!(base_read_name(b"readname/3"), b"readname/3");
    }
}
