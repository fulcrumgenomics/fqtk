use crate::commands::command::Command;
use anyhow::{Result, anyhow};
use clap::Parser;
use fgoxide::io::Io;
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

const BUFFER_SIZE: usize = 64 * 1024;

struct ShardWriters<W: Write> {
    /// Name of the sample this set of writers is for
    shard_number: usize,
    /// The set of writers for the shard
    writers: Vec<W>,
}

impl<W: Write> ShardWriters<W> {
    /// Destroys this struct and decomposes it into its component types. Used when swapping
    /// writers for pooled writers.
    #[allow(clippy::type_complexity)]
    fn into_parts(self) -> (usize, Vec<W>) {
        (self.shard_number, self.writers)
    }

    ///
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
    /// Attempts to gracefully shutdown the writers in this struct, consuming the struct in the
    /// process
    /// # Errors
    ///     - Will error if closing of the ``PooledWriter``s fails for any reason
    fn close(self) -> Result<()> {
        self.writers.into_iter().map(|w| w.close()).collect::<Result<Vec<_>, io::Error>>()?;
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////
// shard (main class) and it's impls
////////////////////////////////////////////////////////////////////////////////

/// Shards a set of FASTQs into N output shards.
///
#[derive(Parser, Debug)]
#[command(version)]
#[clap(verbatim_doc_comment)]
pub(crate) struct Shard {
    /// One or more input FASTQ files each corresponding to a sequencing read (e.g. R1, I1).
    #[clap(long, short = 'i', required = true, num_args = 1..)]
    inputs: Vec<PathBuf>,

    /// Output prefix for sharded FASTQ file(s).
    #[clap(long, short = 'o')]
    output_prefix: String,

    /// Maximum mismatches for a barcode to be considered a match.
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

    fn build_writer_pool(
        prefix: &str,
        sources: usize,
        shards: usize,
        threads: usize,
        compression: u8,
    ) -> Result<(Pool, Vec<ShardWriters<PooledWriter>>)> {
        // First build up the per-shard writers
        let mut shard_writers = Vec::with_capacity(shards);
        for shard in 1..=shards {
            let mut ws = Vec::with_capacity(sources);

            for source_idx in 1..=sources {
                let path_str = format!("{}.s{}.r{}.fq.gz", prefix, shard, source_idx);
                let path = Path::new(&path_str);
                let writer = BufWriter::new(File::create(path)?);
                ws.push(writer);
            }

            shard_writers.push(ShardWriters { shard_number: shard, writers: ws });
        }

        // Then construct the writer pool
        let mut pool_builder = PoolBuilder::<_, BgzfCompressor>::new()
            .threads(threads)
            .queue_size(threads * 50)
            .compression_level(compression)?;

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
    /// Executes the demux command
    fn execute(&self) -> Result<()> {
        info!("Reading {} input FASTQs and generating {} shards.", self.inputs.len(), self.shards);

        let mut fq_readers = Self::build_readers(&self.inputs)?;
        let mut fq_iters = fq_readers.iter_mut().map(|r| r.records()).collect_vec();

        let (mut pool, mut shard_writers) = Self::build_writer_pool(
            self.output_prefix.as_str(),
            self.inputs.len(),
            self.shards,
            self.threads,
            self.compression_level,
        )?;

        let logger = ProgLogBuilder::new()
            .name("fqtk")
            .noun("record sets")
            .verb("read")
            .unit(5_000_000)
            .count_formatter(CountFormatterKind::Comma)
            .level(log::Level::Info)
            .build();

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
