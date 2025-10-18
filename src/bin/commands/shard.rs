use crate::commands::command::Command;
use anyhow::{Result, anyhow};
use clap::Parser;
use fgoxide::io::Io;
use fgoxide::iter::IntoChunkedReadAheadIterator;
use itertools::Itertools;
use log::info;
use pooled_writer::{Pool, PoolBuilder, PooledWriter, bgzf::BgzfCompressor};
use proglog::{CountFormatterKind, ProgLogBuilder};
use seq_io::fastq::OwnedRecord;
use seq_io::fastq::Reader as FastqReader;
use seq_io::fastq::Record;
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::{
    io,
    path::{Path, PathBuf},
};

/// Buffer size used when opening BufReads for input files
const BUFFER_SIZE: usize = 64 * 1024;

/// Struct to hold all the output writers for a single shard
struct ShardWriters<W: Write> {
    /// Name of the sample this set of writers is for
    shard_number: usize,
    /// The set of writers for the shard
    writers: Vec<W>,
}

impl<W: Write> ShardWriters<W> {
    /// Destroys this struct and decomposes it into its component types. Used when swapping
    /// writers for pooled writers.
    fn into_parts(self) -> (usize, Vec<W>) {
        (self.shard_number, self.writers)
    }

    /// Writes a set of FASTQ records out.  The number of reads (i.e. `reads.len()`) must
    /// match the number of individual writers; if a mismatch is found an error will be raised.
    fn write(&mut self, reads: &[OwnedRecord]) -> Result<()> {
        if reads.len() != self.writers.len() {
            return Err(anyhow!("Expected {} reads, got {}", self.shard_number, reads.len()));
        }

        for (writer, read) in self.writers.iter_mut().zip(reads.iter()) {
            read.write(writer)?;
        }

        Ok(())
    }
}

impl ShardWriters<PooledWriter> {
    /// Attempts to gracefully shut down the writers in this struct, consuming the struct in the
    /// process. Will error if closing of the `PooledWriter`s fails for any reason.
    fn close(self) -> Result<()> {
        self.writers.into_iter().map(|w| w.close()).collect::<Result<Vec<_>, io::Error>>()?;
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
/// Inputs may be uncompressed gzipped, or block-gzipped.  Output files are _always_ block gzipped.
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
    #[clap(long, short = 's')]
    shards: usize,

    /// The number of threads to use.
    #[clap(long, short = 't', default_value = "8")]
    threads: usize,

    /// The level of compression to use to compress outputs.
    #[clap(long, short = 'c', default_value = "5")]
    compression_level: u8,
}

impl Shard {
    /// Opens a FASTQ reader per input path.  Handles uncompressed and gzipped FASTQ files.
    fn build_readers(paths: &[PathBuf]) -> Result<Vec<FastqReader<Box<dyn BufRead + Send>>>> {
        let fgio = Io::new(5, BUFFER_SIZE);
        let readers = paths
            .iter()
            .map(|p| fgio.new_reader(p))
            .collect::<Result<Vec<_>, fgoxide::FgError>>()?;

        let fq_readers =
            readers.into_iter().map(|r| FastqReader::with_capacity(r, 10)).collect_vec();

        Ok(fq_readers)
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

        // Then construct the writer pool
        let mut pool_builder = PoolBuilder::<_, BgzfCompressor>::new()
            .threads(self.threads)
            .queue_size(self.threads * 50)
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
    #[allow(clippy::too_many_lines)]
    fn execute(&self) -> Result<()> {
        info!("Reading {} input FASTQs and generating {} shards.", self.inputs.len(), self.shards);

        // Open the input FASTQ files and from each file generate a chunked read-ahead iterator
        // to ensure that the reading of the input files is not the bottleneck.  We don't count
        // the read-ahead threads as consuming from the thread-count as they will consume minimal
        // CPU and are mostly going to be doing/blocking on I/O.
        let fq_readers = Self::build_readers(&self.inputs)?;
        let fq_iters = fq_readers.into_iter().map(|r| r.into_records()).collect_vec();
        let mut fq_iters = fq_iters.into_iter().map(|i| i.read_ahead(128, 48)).collect_vec();

        let (mut pool, mut shard_writers) = self.build_writer_pool()?;

        let logger = ProgLogBuilder::new()
            .name("fqtk")
            .noun("record sets")
            .verb("read")
            .unit(5_000_000)
            .count_formatter(CountFormatterKind::Comma)
            .level(log::Level::Info)
            .build();

        // Loop, consuming one read from each input file, and writing to the appropriate
        // shard.  Terminate cleanly when all input iterators are exhausted, or with error
        // if at least one, but not all input iterators are exhausted.
        let mut target_shard_idx: usize = 0;
        loop {
            // Pull in the next set of reads
            let mut recs = Vec::with_capacity(fq_iters.len());
            for iter in &mut fq_iters {
                if let Some(rec) = iter.next() {
                    recs.push(rec?);
                }
            }

            if recs.is_empty() {
                break;
            }

            assert_eq!(
                recs.len(),
                fq_iters.len(),
                "FASTQ sources out of sync at records: {:?}",
                recs
            );

            shard_writers[target_shard_idx].write(&recs)?;
            target_shard_idx = (target_shard_idx + 1) % self.shards;
            logger.record();
        }

        // Shut down the pool
        info!("Finished reading input FASTQs.");
        shard_writers.into_iter().map(|w| w.close()).collect::<Result<Vec<_>>>()?;
        pool.stop_pool()?;
        info!("Output FASTQ writing complete.");

        Ok(())
    }
}
