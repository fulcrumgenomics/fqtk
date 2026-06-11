//! Background-decompressing FASTQ input readers shared across subcommands.
//!
//! Splits gzip inflation off the consuming thread: for each input, a background read-ahead thread
//! decompresses the file into fixed-size byte chunks ([`ByteChunker`]), and the returned readers
//! parse records from those bytes on the consuming (main) thread ([`ChunkReader`]).  Keeping
//! inflation off the main thread leaves it free to parse records and route them downstream (e.g.
//! to a compression pool).

use anyhow::Result;
use fgoxide::io::Io;
use fgoxide::iter::{ChunkedReadAheadIterator, IntoChunkedReadAheadIterator};
use seq_io::fastq::Reader as FastqReader;
use std::cmp::min;
use std::io::{self, BufRead, Read};
use std::path::PathBuf;

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
/// let reader = ReadAheadBuilder { path, ..Default::default() }.build()?;
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
    /// Opens the input and returns a FASTQ reader that parses records on the consuming thread
    /// while a background read-ahead thread decompresses the file into fixed-size byte chunks.
    /// This keeps gzip inflation off the consuming thread, leaving it free to parse and route
    /// records.
    pub(crate) fn build(self) -> Result<FastqReader<ChunkReader>> {
        let decompressed = Io::new(5, BUFFER_SIZE)
            .new_reader(&self.path)
            .map_err(|e| anyhow::anyhow!("Failed to open {:?}: {e}", self.path))?;
        let chunker =
            ByteChunker { reader: decompressed, chunk_size: self.chunk_size, done: false };
        let chunks = chunker.read_ahead(1, self.chunk_count);
        let chunk_reader = ChunkReader { chunks, current: Vec::new(), pos: 0 };
        Ok(FastqReader::with_capacity(chunk_reader, BUFFER_SIZE))
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

/// Iterator that reads fixed-size chunks of decompressed bytes from an input reader.
///
/// The wrapped reader decompresses lazily as it is read, so running this iterator on a
/// background (read-ahead) thread keeps gzip inflation off the consuming thread.  Each chunk ends
/// on an arbitrary byte boundary; records that span two chunks are stitched back together
/// downstream by [`ChunkReader`] before `seq_io` parses them.
struct ByteChunker {
    /// The decompressing input reader.
    reader: Box<dyn BufRead + Send>,
    /// Target size, in bytes, of each emitted chunk.
    chunk_size: usize,
    /// Set once the underlying reader is exhausted or has errored.
    done: bool,
}

impl Iterator for ByteChunker {
    type Item = ByteChunk;

    fn next(&mut self) -> Option<Self::Item> {
        if self.done {
            return None;
        }
        // Fill a fresh buffer as fully as possible; a single `read` on a decompressing reader
        // often returns less than requested, so loop until full or the input is exhausted.
        let mut buf = vec![0u8; self.chunk_size];
        let mut filled = 0;
        while filled < self.chunk_size {
            match self.reader.read(&mut buf[filled..]) {
                Ok(0) => break,
                Ok(n) => filled += n,
                Err(e) => {
                    self.done = true;
                    return Some(Err(e));
                }
            }
        }
        if filled == 0 {
            self.done = true;
            None
        } else {
            buf.truncate(filled);
            Some(Ok(buf))
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
