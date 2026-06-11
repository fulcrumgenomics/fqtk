use crate::commands::command::Command;
use anyhow::{Result, anyhow};
use clap::Parser;
use clap::builder::RangedU64ValueParser;
use fgoxide::io::Io;
use fgoxide::iter::{ChunkedReadAheadIterator, IntoChunkedReadAheadIterator};
use itertools::Itertools;
use log::info;
use pooled_writer::{Pool, PoolBuilder, PooledWriter, bgzf::BgzfCompressor};
use proglog::{CountFormatterKind, ProgLogBuilder};
use seq_io::fastq::Reader as FastqReader;
use std::cmp::min;
use std::fs::File;
use std::io::{self, BufRead, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

/// Buffer size used when opening BufReads for input files
const BUFFER_SIZE: usize = 1024 * 1024;

/// A decompressed chunk of input bytes shuttled from a decompression thread to the main
/// (parsing) thread, or an I/O error encountered while decompressing.
type ByteChunk = io::Result<Vec<u8>>;

/// Read-ahead stream of decompressed input byte chunks produced on a background thread.
type ByteChunkStream = ChunkedReadAheadIterator<ByteChunk>;

/// A FASTQ reader whose records are parsed on the consuming thread from a [`ByteChunkStream`].
type InputReader = FastqReader<ChunkReader>;

/// Iterator that reads fixed-size chunks of decompressed bytes from an input reader.
///
/// The wrapped reader decompresses lazily as it is read, so running this iterator on a
/// background (read-ahead) thread keeps gzip inflation off the main thread.  Each chunk ends on
/// an arbitrary byte boundary; records that span two chunks are stitched back together
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
struct ChunkReader {
    /// The stream of decompressed byte chunks from the decompression thread.
    chunks: ByteChunkStream,
    /// The chunk currently being served from.
    current: Vec<u8>,
    /// Offset of the next unserved byte within `current`.
    pos: usize,
}

impl Read for ChunkReader {
    fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
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

/// Struct to hold all the output writers for a single shard
struct ShardWriters<W: Write> {
    /// ID of the shard this set of writers is for
    shard_number: usize,
    /// The set of writers for the shard
    writers: Vec<W>,
}

impl<W: Write> ShardWriters<W> {
    /// Consumes this struct and decomposes it into its component parts. Used when swapping
    /// writers for pooled writers.
    fn into_parts(self) -> (usize, Vec<W>) {
        (self.shard_number, self.writers)
    }
}

impl ShardWriters<PooledWriter> {
    /// Attempts to gracefully shut down the writers in this struct, consuming the struct in the
    /// process. Will error if closing of the `PooledWriter`s fails for any reason.
    fn close(self) -> Result<()> {
        self.writers.into_iter().try_for_each(PooledWriter::close)?;
        Ok(())
    }
}

/// Shards a set of FASTQs into N output shards.
///
/// Shards a set of matched FASTQs (e.g. R1 and R2) into one or more set of FASTQs where each
/// input read ends up in exactly one output FASTQ. Reads are assigned to shards on a round-robin
/// basis, so e.g. if using `--shards 10` the first read in the input files will end up in the
/// first shard, the second read in the second shard ... and the tenth read in the tenth shard.
///
/// Each shard will contain one output FASTQ file per input FASTQ files.  Output files are named
/// as follows:
///
/// ```
/// {output_prefix}.{shard_prefix}{shard_num}.{read_number_prefix}{read_num}.fq.gz
/// ```
///
/// where `shard_num` is n for the nth shard (starting at 1), `read_num` corresponds to the nth
/// file in the `inputs` list (starting at 1), and all other values in `{}` are named command
/// line parameters.  The `output_prefix` may contain an absolute path, or a relative path, with
/// relative paths interpreted relative to the working directory where the command is run.
///
/// Inputs may be uncompressed, gzipped, or block-gzipped.  Output files are _always_ block gzipped.
///
#[derive(Parser, Debug)]
#[command(version)]
#[clap(verbatim_doc_comment)]
pub(crate) struct Shard {
    /// One or more input FASTQ files each corresponding to a sequencing read (e.g. R1, R2).
    #[clap(long, short = 'i', required = true, num_args = 1..)]
    inputs: Vec<PathBuf>,

    /// Output prefix for sharded FASTQ file(s).
    #[clap(long, short = 'o')]
    output_prefix: String,

    /// Prefix to place before the shard number in the generated output file names.
    #[clap(long, short = 'S', default_value = "s")]
    shard_prefix: String,

    /// Prefix to place before the read number in the generated output file names.
    #[clap(long, short = 'R', default_value = "r")]
    read_number_prefix: String,

    /// Number of shards to generate
    #[clap(long, short = 's', value_parser = RangedU64ValueParser::<usize>::new().range(1..))]
    shards: usize,

    /// The number of threads to use for compressing output files.  Minimum 2.
    #[clap(long, short = 't', default_value = "8", value_parser = RangedU64ValueParser::<usize>::new().range(2..))]
    threads: usize,

    /// The level of compression to use to compress outputs.  Defaults to 1 because sharded FASTQs
    /// are typically short-lived intermediates, where write throughput matters more than squeezing
    /// out the last few percent of file size.
    #[clap(long, short = 'c', default_value = "1")]
    compression_level: u8,

    /// (hidden) Size, in bytes, of each decompressed input chunk passed from a per-input
    /// decompression thread to the main parsing thread.
    #[clap(long, hide = true, default_value = "131072",
           value_parser = RangedU64ValueParser::<usize>::new().range(1..))]
    chunk_size: usize,

    /// (hidden) Number of decompressed input chunks to keep in flight per input file.
    #[clap(long, hide = true, default_value = "32",
           value_parser = RangedU64ValueParser::<usize>::new().range(1..))]
    chunk_count: usize,
}

impl Shard {
    /// Opens one FASTQ reader per input path, splitting decompression and parsing across threads.
    ///
    /// For each input, a background read-ahead thread decompresses the file into fixed-size byte
    /// chunks (sized by `--chunk-size`, with `--chunk-count` chunks in flight); the returned
    /// readers parse those bytes into records on the consuming thread.  This keeps gzip inflation
    /// off the main thread, leaving it free to parse and route records to the compression pool.
    fn build_record_readers(&self) -> Result<Vec<InputReader>> {
        let fgio = Io::new(5, BUFFER_SIZE);
        let mut readers = Vec::with_capacity(self.inputs.len());
        for path in &self.inputs {
            let decompressed = fgio.new_reader(path)?;
            let chunker =
                ByteChunker { reader: decompressed, chunk_size: self.chunk_size, done: false };
            let chunks = chunker.read_ahead(1, self.chunk_count);
            let chunk_reader = ChunkReader { chunks, current: Vec::new(), pos: 0 };
            readers.push(FastqReader::with_capacity(chunk_reader, BUFFER_SIZE));
        }
        Ok(readers)
    }

    /// Builds the fastq writers for each output file and shard, and then instantiates
    /// a writer pool using block-gzip compression.  Returns the writer Pool itself along
    /// with a Vec of `ShardWriters` each of which can accept and write `Vec`s of fastq records
    /// with one record from each of the input files.
    fn build_writer_pool(&self) -> Result<(Pool, Vec<ShardWriters<PooledWriter>>)> {
        // First build up the per-shard writers
        let mut shard_writers = Vec::with_capacity(self.shards);
        for shard in 1..=self.shards {
            let mut ws = Vec::with_capacity(self.inputs.len());

            for source_idx in 1..=self.inputs.len() {
                let path_str = format!(
                    "{prefix}.{shard_prefix}{shard_num}.{read_prefix}{read_num}.fq.gz",
                    prefix = self.output_prefix,
                    shard_prefix = self.shard_prefix,
                    shard_num = shard,
                    read_prefix = self.read_number_prefix,
                    read_num = source_idx
                );
                let path = Path::new(&path_str);
                let writer = BufWriter::new(File::create(path)?);
                ws.push(writer);
            }

            shard_writers.push(ShardWriters { shard_number: shard, writers: ws });
        }

        // Then construct the writer pool, reserving one thread for the main loop.
        let pool_threads = self.threads - 1;
        let mut pool_builder = PoolBuilder::<_, BgzfCompressor>::new()
            .threads(pool_threads)
            .queue_size(pool_threads * 50)
            .compression_level(self.compression_level)?;

        // Then exchange the writers
        let mut pooled_shard_writers = Vec::with_capacity(shard_writers.len());
        for shard_writer in shard_writers.into_iter() {
            let (shard, writers) = shard_writer.into_parts();
            let pooled_writers =
                writers.into_iter().map(|w| pool_builder.exchange(w)).collect_vec();
            pooled_shard_writers
                .push(ShardWriters { shard_number: shard, writers: pooled_writers });
        }

        let pool = pool_builder.build()?;
        Ok((pool, pooled_shard_writers))
    }
}

impl Command for Shard {
    fn execute(&self) -> Result<()> {
        info!("Reading {} input FASTQs and generating {} shards.", self.inputs.len(), self.shards);

        // Open one reader per input.  Decompression runs on a background read-ahead thread per
        // input and hands decompressed byte chunks to these readers, which parse records here on
        // the main thread.  Keeping inflation off the main thread lets it spend its time parsing
        // and routing records to keep the compression pool fed.
        let mut readers = self.build_record_readers()?;
        let n_inputs = readers.len();

        let (mut pool, mut shard_writers) = self.build_writer_pool()?;

        let logger = ProgLogBuilder::new()
            .name("fqtk")
            .noun("record sets")
            .verb("read")
            .unit(5_000_000)
            .count_formatter(CountFormatterKind::Comma)
            .level(log::Level::Info)
            .build();

        // Loop, consuming one record from each input file and writing the set to the current
        // shard.  Records are emitted verbatim via `write_unchanged`, which blits the original
        // record bytes rather than re-serializing field by field.  Terminate cleanly when all
        // inputs are exhausted together, or with an error if some -- but not all -- run dry.
        let mut target_shard_idx: usize = 0;
        loop {
            let target = &mut shard_writers[target_shard_idx];
            let mut n_written = 0;
            for (writer, reader) in target.writers.iter_mut().zip(readers.iter_mut()) {
                match reader.next() {
                    Some(Ok(record)) => {
                        record.write_unchanged(&mut *writer)?;
                        n_written += 1;
                    }
                    Some(Err(e)) => return Err(anyhow!("Error reading FASTQ input: {e}")),
                    None => {}
                }
            }

            if n_written == 0 {
                break;
            }
            if n_written != n_inputs {
                return Err(anyhow!(
                    "FASTQ sources out of sync; expected {} records but got {}.",
                    n_inputs,
                    n_written
                ));
            }

            target_shard_idx = (target_shard_idx + 1) % self.shards;
            logger.record();
        }

        // Shut down the pool
        info!("Finished reading input FASTQs.");
        shard_writers.into_iter().try_for_each(ShardWriters::close)?;
        pool.stop_pool()?;
        info!("Output FASTQ writing complete.");

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::commands::command::Command;
    use crate::commands::shard::Shard;
    use bstr::ByteSlice;
    use itertools::Itertools;
    use rand;
    use seq_io::fastq::{OwnedRecord, Record};
    use std::collections::HashSet;
    use std::fs::File;
    use std::io::{BufReader, BufWriter, Write};
    use std::path::{Path, PathBuf};
    use tempfile::TempDir;

    /// Writes zero or more records to the given `Write`.  Bases and quals are randomly
    /// generated.  Read names are of the format `{prefix}{suffix}{idx}` where `idx` starts
    /// from the value of the `idx` parameter, and increases by one for each record.
    fn write_fastq_records<W: Write>(
        out: &mut W,
        prefix: &str,
        suffix: &str,
        idx: usize,
        count: usize,
    ) {
        let bases = "ACGT".as_bytes();
        for i in idx..idx + count {
            let seq = (0..30).map(|_| rand::random_range(0..4)).map(|j| bases[j]).collect_vec();
            let qual = (0..30).map(|_| rand::random_range(2u8..40u8) + 33).collect_vec();
            let rec = OwnedRecord {
                head: format!("{}{}{}", prefix, i, suffix).as_bytes().to_owned(),
                seq,
                qual,
            };
            rec.write(&mut *out).unwrap();
        }
    }

    /// Writes records to a FASTQ file at `path`.  Compression is selected by file extension
    /// (`.gz` produces gzip via fgoxide; otherwise plain text).
    fn build_fastq(path: &Path, prefix: &str, suffix: &str, idx: usize, count: usize) {
        let io = fgoxide::io::Io::new(1, 8 * 1024);
        let mut out = io.new_writer(path).unwrap();
        write_fastq_records(&mut out, prefix, suffix, idx, count);
    }

    /// Writes records to a BGZF-compressed FASTQ file at `path`.
    fn build_fastq_bgzf(path: &Path, prefix: &str, suffix: &str, idx: usize, count: usize) {
        let file = BufWriter::new(File::create(path).unwrap());
        let mut writer = bgzf::Writer::new(file, bgzf::CompressionLevel::new(3).unwrap());
        write_fastq_records(&mut writer, prefix, suffix, idx, count);
        writer.finish().unwrap().flush().unwrap();
    }

    /// Runs sharding and returns the outputs.  The returned Vec is nested as follows:
    ///  - Each entry in the top level vec is a _shard_
    ///  - Each entry in the second level vec is the sharded reads from a single input fastq
    ///  - Each entry in the third level vec is an individual fastq record
    fn run_sharding(tmp: &TempDir, inputs: &[&Path], shards: usize) -> Vec<Vec<Vec<OwnedRecord>>> {
        let prefix = format!("{}/test_out", tmp.path().to_str().unwrap());
        build_sharder(inputs, &prefix, shards).execute().unwrap();
        collect_shard_outputs(&prefix, inputs.len(), shards, read_fastq)
    }

    /// Builds a `Shard` instance for testing with the given inputs, output prefix, and shard count.
    fn build_sharder(inputs: &[&Path], output_prefix: &str, shards: usize) -> Shard {
        Shard {
            inputs: inputs.iter().map(|p| p.to_path_buf()).collect_vec(),
            output_prefix: output_prefix.to_string(),
            shard_prefix: "shard".to_string(),
            read_number_prefix: "read".to_string(),
            shards,
            threads: 4,
            compression_level: 1,
            // Deliberately small so most records straddle chunk boundaries, exercising the
            // ChunkReader reassembly path across the test suite.
            chunk_size: 100,
            chunk_count: 8,
        }
    }

    /// Reads back the sharded outputs produced under `output_prefix`, using `reader` to parse each
    /// individual output file.  Returns the same nested layout as `run_sharding`.
    fn collect_shard_outputs<F>(
        output_prefix: &str,
        n_inputs: usize,
        shards: usize,
        reader: F,
    ) -> Vec<Vec<Vec<OwnedRecord>>>
    where
        F: Fn(&Path) -> Vec<OwnedRecord>,
    {
        let mut results: Vec<Vec<Vec<OwnedRecord>>> = Vec::with_capacity(shards);
        for shard in 1..=shards {
            let mut reads_vecs = Vec::with_capacity(n_inputs);
            for input_idx in 1..=n_inputs {
                let path_str = format!("{}.shard{}.read{}.fq.gz", output_prefix, shard, input_idx);
                reads_vecs.push(reader(Path::new(&path_str)));
            }
            results.push(reads_vecs);
        }
        results
    }

    /// Reads a FASTQ file into a vec of records, transparently handling gzip via fgoxide.
    fn read_fastq(path: &Path) -> Vec<OwnedRecord> {
        let io = fgoxide::io::Io::new(1, 8 * 1024);
        let mut reader = io.new_reader(path).unwrap();
        let mut fq_reader = seq_io::fastq::Reader::with_capacity(&mut reader, 8 * 1024);
        fq_reader.records().map(|r| r.unwrap()).collect_vec()
    }

    /// Reads a FASTQ file into a vec of records, decompressing with the bgzf crate's `Reader`.
    fn read_fastq_via_bgzf(path: &Path) -> Vec<OwnedRecord> {
        let file = BufReader::new(File::open(path).unwrap());
        let bgzf_reader = bgzf::Reader::new(file);
        let mut fq_reader = seq_io::fastq::Reader::with_capacity(bgzf_reader, 8 * 1024);
        fq_reader.records().map(|r| r.unwrap()).collect_vec()
    }

    /// Returns the integer index embedded in a read name by `build_fastq`.  Strips any
    /// trailing `/N` read-end marker before parsing the trailing digits.
    fn read_index(rec: &OwnedRecord) -> usize {
        let head = rec.head.to_str().unwrap();
        let trimmed = head
            .rsplit_once('/')
            .filter(|(_, tail)| tail.chars().all(|c| c.is_ascii_digit()))
            .map(|(head, _)| head)
            .unwrap_or(head);
        let digits: String = trimmed.chars().rev().take_while(|c| c.is_ascii_digit()).collect();
        digits.chars().rev().collect::<String>().parse().unwrap()
    }

    #[test]
    fn test_shard_single_file() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        build_fastq(r1.as_path(), "q", "", 1, 50);
        let outputs = run_sharding(&tmp, &[&r1], 5);

        assert_eq!(outputs.len(), 5);
        for shard in outputs.iter() {
            assert_eq!(shard.len(), 1);
            assert_eq!(shard.iter().next().unwrap().len(), 10);
        }

        let read_names: HashSet<&str> =
            outputs.iter().flatten().flatten().map(|r| r.head.to_str().unwrap()).collect();
        assert_eq!(read_names.len(), 50);
    }

    #[test]
    fn test_shard_multiple_files() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        let r2 = PathBuf::from(tmp.path()).join("r2.fq");
        build_fastq(r1.as_path(), "q", "/1", 1, 64);
        build_fastq(r2.as_path(), "q", "/2", 1, 64);
        let outputs = run_sharding(&tmp, &[&r1, &r2], 3);

        assert_eq!(outputs.len(), 3);

        for shard in outputs.iter() {
            assert_eq!(shard.len(), 2); // two reads in each shard
            assert_eq!(shard[0].len(), shard[1].len()); // both r1 and r2 have same number of reads
            assert!(shard[0].len() == 21 || shard[0].len() == 22); // 64 / 3 == 21 1/3
        }
    }

    /// Verifies that the Nth input record (1-based) lands in shard ((N-1) % shards) + 1, and
    /// that within a shard the records appear in input order.
    #[test]
    fn test_round_robin_assignment() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        let n_reads = 20;
        let n_shards = 4;
        build_fastq(r1.as_path(), "q", "", 1, n_reads);
        let outputs = run_sharding(&tmp, &[&r1], n_shards);

        for (shard_idx, shard) in outputs.iter().enumerate() {
            let expected: Vec<usize> =
                (1..=n_reads).filter(|i| (i - 1) % n_shards == shard_idx).collect();
            let actual: Vec<usize> = shard[0].iter().map(read_index).collect();
            assert_eq!(actual, expected, "shard {} contents wrong", shard_idx + 1);
        }
    }

    /// Verifies that for paired input, the j-th r1 record in any shard came from the same input
    /// record-set as the j-th r2 record (i.e., they share an index).
    #[test]
    fn test_paired_records_stay_aligned() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        let r2 = PathBuf::from(tmp.path()).join("r2.fq");
        build_fastq(r1.as_path(), "q", "/1", 1, 30);
        build_fastq(r2.as_path(), "q", "/2", 1, 30);
        let outputs = run_sharding(&tmp, &[&r1, &r2], 4);

        for shard in outputs.iter() {
            let r1_indices: Vec<usize> = shard[0].iter().map(read_index).collect();
            let r2_indices: Vec<usize> = shard[1].iter().map(read_index).collect();
            assert_eq!(r1_indices, r2_indices);
        }
    }

    /// Every input read should appear in exactly one output shard, with no losses or duplicates.
    #[test]
    fn test_no_reads_lost_or_duplicated() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        build_fastq(r1.as_path(), "q", "", 1, 100);
        let outputs = run_sharding(&tmp, &[&r1], 7);

        let all_indices: Vec<usize> =
            outputs.iter().flatten().flatten().map(read_index).sorted().collect();
        let expected: Vec<usize> = (1..=100).collect();
        assert_eq!(all_indices, expected);
    }

    #[test]
    fn test_single_shard() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        build_fastq(r1.as_path(), "q", "", 1, 12);
        let outputs = run_sharding(&tmp, &[&r1], 1);

        assert_eq!(outputs.len(), 1);
        let actual: Vec<usize> = outputs[0][0].iter().map(read_index).collect();
        assert_eq!(actual, (1..=12).collect::<Vec<_>>());
    }

    #[test]
    fn test_more_shards_than_reads() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        build_fastq(r1.as_path(), "q", "", 1, 3);
        let outputs = run_sharding(&tmp, &[&r1], 5);

        assert_eq!(outputs.len(), 5);
        assert_eq!(outputs[0][0].len(), 1);
        assert_eq!(outputs[1][0].len(), 1);
        assert_eq!(outputs[2][0].len(), 1);
        assert_eq!(outputs[3][0].len(), 0);
        assert_eq!(outputs[4][0].len(), 0);
    }

    #[test]
    fn test_empty_input() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        build_fastq(r1.as_path(), "q", "", 1, 0);
        let outputs = run_sharding(&tmp, &[&r1], 4);

        assert_eq!(outputs.len(), 4);
        for shard in outputs.iter() {
            assert_eq!(shard.len(), 1);
            assert!(shard[0].is_empty());
        }
    }

    #[test]
    fn test_mismatched_input_lengths_fails() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        let r2 = PathBuf::from(tmp.path()).join("r2.fq");
        build_fastq(r1.as_path(), "q", "/1", 1, 10);
        build_fastq(r2.as_path(), "q", "/2", 1, 8);

        let prefix = format!("{}/test_out", tmp.path().to_str().unwrap());
        let err = build_sharder(&[&r1, &r2], &prefix, 3).execute().unwrap_err();
        assert!(err.to_string().contains("out of sync"), "unexpected error message: {}", err);
    }

    #[test]
    fn test_gzip_compressed_inputs() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq.gz");
        build_fastq(r1.as_path(), "q", "", 1, 25);
        let outputs = run_sharding(&tmp, &[&r1], 5);

        let all_indices: Vec<usize> =
            outputs.iter().flatten().flatten().map(read_index).sorted().collect();
        assert_eq!(all_indices, (1..=25).collect::<Vec<_>>());
    }

    /// Verifies that BGZF inputs spanning more than one block are read in full (i.e., we don't
    /// stop after the first BGZF block).  `BGZF_BLOCK_SIZE` is ~64 KiB; we generate enough
    /// records to comfortably exceed that uncompressed.
    #[test]
    fn test_bgzf_compressed_inputs_multiple_blocks() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq.gz");
        // Each record is ~70-75 bytes; 1500 records ≈ 105 KiB uncompressed → multiple BGZF blocks.
        let n_reads = 1500;
        build_fastq_bgzf(r1.as_path(), "q", "", 1, n_reads);
        let outputs = run_sharding(&tmp, &[&r1], 3);

        let all_indices: Vec<usize> =
            outputs.iter().flatten().flatten().map(read_index).sorted().collect();
        assert_eq!(all_indices, (1..=n_reads).collect::<Vec<_>>());
    }

    /// Outputs are documented to always be BGZF — read one back with the bgzf crate's `Reader`
    /// and check that the round-tripped contents match what we sent in.
    #[test]
    fn test_output_files_are_bgzf() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        build_fastq(r1.as_path(), "q", "", 1, 30);

        let prefix = format!("{}/test_out", tmp.path().to_str().unwrap());
        build_sharder(&[&r1], &prefix, 3).execute().unwrap();

        let outputs = collect_shard_outputs(&prefix, 1, 3, read_fastq_via_bgzf);
        let all_indices: Vec<usize> =
            outputs.iter().flatten().flatten().map(read_index).sorted().collect();
        assert_eq!(all_indices, (1..=30).collect::<Vec<_>>());
    }

    #[test]
    fn test_custom_prefixes() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        build_fastq(r1.as_path(), "q", "", 1, 6);

        let prefix = format!("{}/sample", tmp.path().to_str().unwrap());
        let sharder = Shard {
            inputs: vec![r1.clone()],
            output_prefix: prefix.clone(),
            shard_prefix: "chunk_".to_string(),
            read_number_prefix: "R".to_string(),
            shards: 2,
            threads: 2,
            compression_level: 1,
            chunk_size: 64 * 1024,
            chunk_count: 4,
        };
        sharder.execute().unwrap();

        for shard in 1..=2 {
            let path = PathBuf::from(format!("{}.chunk_{}.R1.fq.gz", prefix, shard));
            assert!(path.exists(), "expected output file does not exist: {}", path.display());
            assert_eq!(read_fastq(&path).len(), 3);
        }
    }

    /// With a pathologically small chunk size, individual records (and the boundaries between
    /// them) span many decompressed chunks, so the `ChunkReader` must stitch bytes back together
    /// before the parser sees them.  Every read should still round-trip exactly once.
    #[test]
    fn test_tiny_chunk_size_reassembles_records() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        build_fastq(r1.as_path(), "q", "", 1, 40);

        let prefix = format!("{}/test_out", tmp.path().to_str().unwrap());
        let sharder = Shard {
            inputs: vec![r1.clone()],
            output_prefix: prefix.clone(),
            shard_prefix: "shard".to_string(),
            read_number_prefix: "read".to_string(),
            shards: 3,
            threads: 4,
            compression_level: 1,
            chunk_size: 3,
            chunk_count: 2,
        };
        sharder.execute().unwrap();

        let outputs = collect_shard_outputs(&prefix, 1, 3, read_fastq);
        let all_indices: Vec<usize> =
            outputs.iter().flatten().flatten().map(read_index).sorted().collect();
        assert_eq!(all_indices, (1..=40).collect::<Vec<_>>());
    }

    /// With a chunk size larger than the whole input, every record arrives within a single chunk
    /// and no cross-chunk reassembly is needed -- the opposite extreme from the tiny-chunk case.
    #[test]
    fn test_large_chunk_size_single_chunk() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        build_fastq(r1.as_path(), "q", "", 1, 20);

        let prefix = format!("{}/test_out", tmp.path().to_str().unwrap());
        let sharder = Shard {
            inputs: vec![r1.clone()],
            output_prefix: prefix.clone(),
            shard_prefix: "shard".to_string(),
            read_number_prefix: "read".to_string(),
            shards: 2,
            threads: 2,
            compression_level: 1,
            chunk_size: 8 * 1024 * 1024,
            chunk_count: 4,
        };
        sharder.execute().unwrap();

        let outputs = collect_shard_outputs(&prefix, 1, 2, read_fastq);
        let all_indices: Vec<usize> =
            outputs.iter().flatten().flatten().map(read_index).sorted().collect();
        assert_eq!(all_indices, (1..=20).collect::<Vec<_>>());
    }
}
