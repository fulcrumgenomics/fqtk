use crate::commands::command::Command;
use anyhow::{Error, Result};
use clap::Parser;
use fgoxide::io::Io;
use log::info;
use pooled_writer::{bgzf::BgzfCompressor, Pool, PoolBuilder, PooledWriter};
use seq_io::fastq::Reader as FastqReader;
use std::fs::{self, File};
use std::io::{BufRead, BufWriter, Write};
use std::path::Path;
use std::path::PathBuf;

/// Trimming and overlap correction of paired-end reads
#[derive(Parser, Debug)]
#[command(version)]
pub(crate) struct Trimmer {
    /// Reading/Writing threads
    #[clap(long, short = 't', default_value = "5")]
    threads: usize,

    /// Level of compression to use to compress outputs.
    #[clap(long, short = 'c', default_value = "5")]
    compression_level: usize,

    /// Minimum difference in base-quality for one read to correct an overlapping
    /// base from the other read.
    #[clap(long, short = 'd', default_value = "15")]
    overlap_min_bq_delta: usize,

    /// Minimum pair overlap length to attempt correction.
    #[clap(long, short = 'l', default_value = "50")]
    overlap_min_len: usize,

    /// Maximum error-rate allowed in the overlap.
    #[clap(long, short = 'e', default_value = "0.02")]
    overlap_max_error_rate: f64,

    /// The output directory into which to FASTQs.
    #[clap(long, short = 'o', required = true)]
    output: PathBuf,

    /// Fastq file for Read1 and Read2
    #[clap(required = true, num_args = 2)]
    fastqs: Vec<PathBuf>,
}

const BUFFER_SIZE: usize = 1024 * 1024;

/// Type alias to prevent clippy complaining about type complexity
type VecOfReaders = Vec<Box<dyn BufRead + Send>>;

fn create_writer<P: AsRef<Path>>(name: P) -> Result<BufWriter<File>, Error> {
    Ok(BufWriter::new(File::create(name)?))
}

struct TrimIO {
    pool: Pool,
    writers: Vec<PooledWriter>,
    readers: Vec<FastqReader<Box<dyn BufRead + Send>>>,
}

impl Trimmer {
    fn prepare(&self) -> Result<TrimIO, Error> {
        if !self.output.exists() {
            info!("Output directory {:#?} didn't exist, creating it.", self.output);
            fs::create_dir_all(&self.output)?;
        }

        let fgio = Io::new(5, BUFFER_SIZE);
        let fq_readers = self
            .fastqs
            .iter()
            .map(|p| fgio.new_reader(p))
            .collect::<Result<VecOfReaders, fgoxide::FgError>>()?;

        let fq_readers =
            fq_readers.into_iter().map(|fq| FastqReader::with_capacity(fq, BUFFER_SIZE)).collect();

        let writers = vec![
            create_writer(&self.output.join("reads_1.fq.gz"))?,
            create_writer(&self.output.join("reads_2.fq.gz"))?,
        ];

        let mut pool_builder = PoolBuilder::<_, BgzfCompressor>::new()
            .threads(self.threads)
            .queue_size(self.threads * 50)
            .compression_level(u8::try_from(self.compression_level)?)?;

        let pooled_writers =
            writers.into_iter().map(|w| pool_builder.exchange(w)).collect::<Vec<_>>();

        let pool = pool_builder.build()?;

        // see here for pooled writer: https://docs.rs/pooled-writer/0.3.0/pooled_writer/
        //
        // and https://github.com/fulcrumgenomics/fqtk/blob/ae91e90ecc86826fc632837bb982798a7e6b6f7a/src/bin/commands/demux.rs#L804
        // and https://github.com/fulcrumgenomics/fqtk/blob/ae91e90ecc86826fc632837bb982798a7e6b6f7a/src/bin/commands/demux.rs#L850-L851
        // for reader
        Ok(TrimIO { pool, writers: pooled_writers, readers: fq_readers })
    }
}

impl Command for Trimmer {
    fn execute(&self) -> Result<()> {
        let mut trim_io = self.prepare()?;

        writeln!(trim_io.writers[0], "Hello, world!")?;
        _ = trim_io.readers[0].next();

        trim_io.writers.into_iter().try_for_each(|w| w.close())?;
        trim_io.pool.stop_pool()?;
        Ok(())
    }
}
