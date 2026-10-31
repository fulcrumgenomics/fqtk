//! Background-decompressing FASTQ input readers shared across subcommands.
//!
//! Splits gzip inflation off the consuming thread: for each input, a background read-ahead thread
//! opens and decompresses the file into fixed-size byte chunks ([`ByteChunker`]), and the returned
//! readers parse records from those bytes on the consuming (main) thread ([`ChunkReader`]).
//! Keeping inflation off the main thread leaves it free to parse records and route them
//! downstream (e.g. to a compression pool).

use fgoxide::FgError;
use fgoxide::io::Io;
use fgoxide::iter::{ChunkedReadAheadIterator, IntoChunkedReadAheadIterator};
use seq_io::fastq::Reader as FastqReader;
use std::cmp::min;
use std::io::{self, BufRead, Read};
use std::path::{Path, PathBuf};

/// Buffer size used when opening BufReads for input files and sizing the FASTQ parser.
const BUFFER_SIZE: usize = 1024 * 1024;

/// Default size, in bytes, of each decompressed input chunk handed from the decompression thread
/// to the consuming thread.
const DEFAULT_CHUNK_SIZE: usize = 128 * 1024;

/// Default number of decompressed input chunks kept in flight on the read-ahead thread.
const DEFAULT_CHUNK_COUNT: usize = 32;

/// A decompressed chunk of input bytes shuttled from a decompression thread to the consuming
/// (parsing) thread, or an I/O error encountered while decompressing.
type ByteChunk = io::Result<Vec<u8>>;

/// Read-ahead stream of decompressed input byte chunks produced on a background thread.
type ByteChunkStream = ChunkedReadAheadIterator<ByteChunk>;

/// Builds a single background-decompressing FASTQ reader for one input file.
///
/// Construct with the input `path` and, optionally, override `chunk_size`/`chunk_count`; the
/// `Default` impl supplies sensible values for the latter two, so the common case is:
///
/// ```ignore
/// let reader = ReadAheadBuilder { path, ..Default::default() }.build();
/// ```
pub(crate) struct ReadAheadBuilder {
    /// Path to the input FASTQ (may be uncompressed, gzipped, or block-gzipped).
    pub(crate) path: PathBuf,
    /// Target size, in bytes, of each decompressed chunk handed to the consuming thread.
    pub(crate) chunk_size: usize,
    /// Number of decompressed chunks kept in flight on the read-ahead thread.
    pub(crate) chunk_count: usize,
}

impl ReadAheadBuilder {
    /// Returns a FASTQ reader that parses records on the consuming thread while a background
    /// read-ahead thread opens the input and decompresses it into fixed-size byte chunks.
    ///
    /// Never blocks: the input is opened on the read-ahead thread, so a failure to open or decode
    /// it surfaces as an error, naming the path, on the reader's first read.
    pub(crate) fn build(self) -> FastqReader<ChunkReader> {
        let chunker =
            ByteChunker { path: self.path, reader: None, chunk_size: self.chunk_size, done: false };
        let chunks = chunker.read_ahead(1, self.chunk_count);
        let chunk_reader = ChunkReader { chunks, current: Vec::new(), pos: 0 };
        FastqReader::with_capacity(chunk_reader, BUFFER_SIZE)
    }
}

impl Default for ReadAheadBuilder {
    fn default() -> Self {
        Self {
            path: PathBuf::new(),
            chunk_size: DEFAULT_CHUNK_SIZE,
            chunk_count: DEFAULT_CHUNK_COUNT,
        }
    }
}

/// Iterator that reads fixed-size chunks of decompressed bytes from an input file.
///
/// Run on a background (read-ahead) thread, this keeps gzip inflation off the consuming thread.
/// Each chunk ends on an arbitrary byte boundary; records that span two chunks are stitched back
/// together downstream by [`ChunkReader`] before `seq_io` parses them.
///
/// The input is opened on the first call to `next()` rather than on construction, so that it is
/// opened on the read-ahead thread.  Opening a FIFO blocks until its writer opens it, and a gzip
/// decoder reads the header as soon as it is built; doing either on the caller's thread would
/// open inputs one at a time, deadlocking a producer that opens all its FIFOs before writing.
struct ByteChunker {
    /// Path to the input file.
    path: PathBuf,
    /// The decompressing input reader, once opened.
    reader: Option<Box<dyn BufRead + Send>>,
    /// Target size, in bytes, of each emitted chunk.
    chunk_size: usize,
    /// Set once the underlying reader is exhausted or has errored.
    done: bool,
}

impl ByteChunker {
    /// Reads the next chunk, opening the input first if this is the first call.  Returns an empty
    /// chunk once the input is exhausted.
    fn read_chunk(&mut self) -> io::Result<Vec<u8>> {
        let reader = match &mut self.reader {
            Some(reader) => reader,
            None => self.reader.insert(open_input(&self.path)?),
        };
        // Fill a fresh buffer as fully as possible; a single `read` on a decompressing reader
        // often returns less than requested, so loop until full or the input is exhausted.
        let mut buf = vec![0u8; self.chunk_size];
        let mut filled = 0;
        while filled < self.chunk_size {
            match reader.read(&mut buf[filled..])? {
                0 => break,
                n => filled += n,
            }
        }
        buf.truncate(filled);
        Ok(buf)
    }
}

impl Iterator for ByteChunker {
    type Item = ByteChunk;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        match self.read_chunk() {
            Ok(chunk) if chunk.is_empty() => {
                self.done = true;
                None
            }
            Ok(chunk) => Some(Ok(chunk)),
            Err(e) => {
                self.done = true;
                Some(Err(io::Error::new(e.kind(), format!("{}: {e}", self.path.display()))))
            }
        }
    }
}

/// A [`Read`] adapter over a [`ByteChunkStream`] that serves the decompressed bytes received
/// from the decompression thread.  This is the seam that lets FASTQ parsing run on the consuming
/// thread while decompression runs upstream.
pub(crate) struct ChunkReader {
    /// The stream of decompressed byte chunks from the decompression thread.
    chunks: ByteChunkStream,
    /// The chunk currently being served from.
    current: Vec<u8>,
    /// Offset of the next unserved byte within `current`.
    pos: usize,
}

impl Read for ChunkReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
        // Per the `Read` contract a zero-length buffer must return `Ok(0)` without consuming
        // input; without this guard the loop below could block fetching the next chunk.
        if buf.is_empty() {
            return Ok(0);
        }
        loop {
            if self.pos < self.current.len() {
                let n = min(buf.len(), self.current.len() - self.pos);
                buf[..n].copy_from_slice(&self.current[self.pos..self.pos + n]);
                self.pos += n;
                return Ok(n);
            }
            match self.chunks.next() {
                Some(Ok(chunk)) => {
                    self.current = chunk;
                    self.pos = 0;
                }
                Some(Err(e)) => return Err(e),
                None => return Ok(0),
            }
        }
    }
}

/// Opens `path` for reading, decompressing it if its extension marks it as gzipped.
fn open_input(path: &Path) -> io::Result<Box<dyn BufRead + Send>> {
    // Unwrap the underlying `io::Error`: `FgError`'s message hides the cause (e.g. "No such
    // file or directory").
    Io::new(5, BUFFER_SIZE).new_reader(path).map_err(|e| match e {
        FgError::IoError(e) => e,
        e => io::Error::other(e),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    /// Returns `count` FASTQ records named `{prefix}1` through `{prefix}{count}`.
    fn fastq_text(prefix: &str, count: usize) -> String {
        let mut text = String::new();
        for i in 1..=count {
            text.push_str(&format!("@{prefix}{i}\nACGT\n+\nIIII\n"));
        }
        text
    }

    #[test]
    fn missing_input_errors_on_first_read_naming_the_path() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("missing.fq.gz");
        let mut reader = ReadAheadBuilder { path: path.clone(), ..Default::default() }.build();

        let err = reader.next().unwrap().unwrap_err().to_string();
        assert!(err.contains(&path.display().to_string()), "error does not name the path: {err}");
    }

    #[test]
    fn invalid_gzip_input_errors_on_first_read_naming_the_path() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("plain.fq.gz");
        std::fs::write(&path, fastq_text("q", 3)).unwrap();
        let mut reader = ReadAheadBuilder { path: path.clone(), ..Default::default() }.build();

        let err = reader.next().unwrap().unwrap_err().to_string();
        assert!(err.contains(&path.display().to_string()), "error does not name the path: {err}");
    }

    /// Inputs that are FIFOs (named pipes), whose opens block until the other end is opened.
    #[cfg(unix)]
    mod fifo_inputs {
        use super::*;
        use flate2::Compression;
        use flate2::write::GzEncoder;
        use std::fs::File;
        use std::io::Write;
        use std::panic::resume_unwind;
        use std::sync::mpsc::{self, RecvTimeoutError};
        use std::thread;
        use std::time::Duration;

        /// How long a test waits for its readers; far longer than they need, so reaching it means
        /// they deadlocked.
        const DEADLINE: Duration = Duration::from_secs(10);

        /// Creates a FIFO at `path`.
        fn make_fifo(path: &Path) {
            let status = std::process::Command::new("mkfifo").arg(path).status().unwrap();
            assert!(status.success(), "mkfifo failed for {}", path.display());
        }

        /// Runs `f` on its own thread and returns its result, panicking if it takes longer than
        /// `DEADLINE` so that a deadlock fails the test instead of hanging the suite.
        fn run_with_deadline<T: Send + 'static>(f: impl FnOnce() -> T + Send + 'static) -> T {
            let (sender, receiver) = mpsc::channel();
            let handle = thread::spawn(move || {
                let _ = sender.send(f());
            });
            match receiver.recv_timeout(DEADLINE) {
                Ok(result) => result,
                Err(RecvTimeoutError::Timeout) => {
                    panic!("readers did not finish within {DEADLINE:?}")
                }
                Err(RecvTimeoutError::Disconnected) => resume_unwind(handle.join().unwrap_err()),
            }
        }

        /// Reads every record from `reader`, returning the record names.
        fn read_names(mut reader: FastqReader<ChunkReader>) -> Vec<String> {
            reader.records().map(|rec| String::from_utf8(rec.unwrap().head).unwrap()).collect()
        }

        #[test]
        fn gzipped_fifos_opened_together_before_writing_are_read() {
            let dir = TempDir::new().unwrap();
            let r1 = dir.path().join("r1.fq.gz");
            let r2 = dir.path().join("r2.fq.gz");
            make_fifo(&r1);
            make_fifo(&r2);

            // Each open blocks until the reader opens the other end, so this producer can only
            // start writing once both inputs are open.
            let (r1_out, r2_out) = (r1.clone(), r2.clone());
            let producer = thread::spawn(move || {
                let fifos = [File::create(r1_out).unwrap(), File::create(r2_out).unwrap()];
                for (fifo, prefix) in fifos.into_iter().zip(["r1_", "r2_"]) {
                    let mut gz = GzEncoder::new(fifo, Compression::fast());
                    gz.write_all(fastq_text(prefix, 3).as_bytes()).unwrap();
                    gz.finish().unwrap();
                }
            });

            let (r1_names, r2_names) = run_with_deadline(move || {
                let r1_reader = ReadAheadBuilder { path: r1, ..Default::default() }.build();
                let r2_reader = ReadAheadBuilder { path: r2, ..Default::default() }.build();
                (read_names(r1_reader), read_names(r2_reader))
            });
            producer.join().unwrap();
            assert_eq!(r1_names, ["r1_1", "r1_2", "r1_3"]);
            assert_eq!(r2_names, ["r2_1", "r2_2", "r2_3"]);
        }

        #[test]
        fn plain_fifos_opened_in_reverse_order_are_read() {
            let dir = TempDir::new().unwrap();
            let r1 = dir.path().join("r1.fq");
            let r2 = dir.path().join("r2.fq");
            make_fifo(&r1);
            make_fifo(&r2);

            // Opens R2 first, the reverse of the order in which the readers are built.
            let (r1_out, r2_out) = (r1.clone(), r2.clone());
            let producer = thread::spawn(move || {
                let mut r2_fifo = File::create(r2_out).unwrap();
                let mut r1_fifo = File::create(r1_out).unwrap();
                r1_fifo.write_all(fastq_text("r1_", 3).as_bytes()).unwrap();
                r2_fifo.write_all(fastq_text("r2_", 3).as_bytes()).unwrap();
            });

            let (r1_names, r2_names) = run_with_deadline(move || {
                let r1_reader = ReadAheadBuilder { path: r1, ..Default::default() }.build();
                let r2_reader = ReadAheadBuilder { path: r2, ..Default::default() }.build();
                (read_names(r1_reader), read_names(r2_reader))
            });
            producer.join().unwrap();
            assert_eq!(r1_names, ["r1_1", "r1_2", "r1_3"]);
            assert_eq!(r2_names, ["r2_1", "r2_2", "r2_3"]);
        }
    }
}
