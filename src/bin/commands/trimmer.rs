use crate::commands::command::Command;
use anyhow::{Error, Result};
use clap::{Parser, ValueEnum};
use fgoxide::io::Io;
use fqtk_lib::{base_quality, pair_overlap};
use log::info;
use pooled_writer::{bgzf::BgzfCompressor, Pool, PoolBuilder, PooledWriter};
use seq_io::fastq::{Reader as FastqReader, Record};
use std::fs::File;
use std::io::{BufRead, BufWriter};
use std::path::Path;
use std::path::PathBuf;

#[derive(ValueEnum, Clone, Debug)]
enum Operation {
    TrimQual,
    Overlap,
    Osc,
    LenFilter,
}

/// Trimming and overlap correction of paired-end reads
#[derive(Parser, Debug)]
#[command(version)]
pub(crate) struct TrimmerOpts {
    /// Reading/Writing threads
    #[clap(long, short = 't', default_value = "5")]
    threads: usize,

    /// Minimum base-quality to keep a base when trimming tails.
    #[clap(long, short = 'q', default_value = "20")]
    trim_tail_quality: u8,

    #[clap(value_enum)]
    trim_tail_side: base_quality::Tail,

    /// Window size for moving average when trimming tails.
    #[clap(long, short = 'w', default_value = "20")]
    trim_tail_window: u8,

    /// Level of compression to use to compress outputs.
    #[clap(long, short = 'c', default_value = "5")]
    compression_level: usize,

    /// Length requirement of shorter read.
    #[clap(long, short = 'S', default_value = "5")]
    filter_shorter: usize,

    /// Length requirement of shorter read.
    #[clap(long, short = 'L', default_value = "15")]
    filter_longer: usize,

    /// Size of window to look for oscillations.
    #[clap(long, default_value = "15")]
    osc_window: usize,

    /// Number of oscillations in a window to trigger trimming/masking.
    #[clap(long, default_value = "4")]
    osc_max_oscillations: usize,

    /// Difference between adjacent bases to be considered on oscilation.
    #[clap(long, default_value = "10")]
    osc_delta: usize,

    /// Minimum difference in base-quality for one read to correct an overlapping
    /// base from the other read.
    #[clap(long, short = 'd', default_value = "15")]
    overlap_min_bq_delta: u8,

    /// Minimum pair overlap length to attempt correction.
    #[clap(long, short = 'l', default_value = "50")]
    overlap_min_len: usize,

    /// Maximum error-rate allowed in the overlap.
    #[clap(long, short = 'e', default_value = "0.02")]
    overlap_max_error_rate: f64,

    /// The output directory into which to FASTQs.
    #[clap(long, short = 'o', required = true, num_args = 2)]
    output: Vec<PathBuf>,

    /// Fastqs file for Read1 and Read2
    #[clap(long, short = 'i', required = true, num_args = 2)]
    input: Vec<PathBuf>,

    /// Order of operations
    #[clap(value_enum, short = 'p', default_value = "overlap", num_args = 1..)]
    operations: Vec<Operation>,
}

const BUFFER_SIZE: usize = 1024 * 1024;

/// Type alias to prevent clippy complaining about type complexity
type VecOfReaders = Vec<Box<dyn BufRead + Send>>;
type VecOfFqReaders = Vec<FastqReader<Box<dyn BufRead + Send>>>;

fn create_writer<P: AsRef<Path>>(name: P) -> Result<BufWriter<File>, Error> {
    Ok(BufWriter::new(File::create(name)?))
}

impl TrimmerOpts {
    fn prepare(&self) -> Result<(Pool, Vec<PooledWriter>, VecOfFqReaders), Error> {
        let fgio = Io::new(5, BUFFER_SIZE);
        let fq_readers = self
            .input
            .iter()
            .map(|p| fgio.new_reader(p))
            .collect::<Result<VecOfReaders, fgoxide::FgError>>()?;

        let fq_readers =
            fq_readers.into_iter().map(|fq| FastqReader::with_capacity(fq, BUFFER_SIZE)).collect();

        let writers = vec![create_writer(&self.output[0])?, create_writer(&self.output[1])?];

        let mut pool_builder = PoolBuilder::<_, BgzfCompressor>::new()
            .threads(self.threads)
            .queue_size(self.threads * 50)
            .compression_level(u8::try_from(self.compression_level)?)?;

        let pooled_writers =
            writers.into_iter().map(|w| pool_builder.exchange(w)).collect::<Vec<_>>();

        let pool = pool_builder.build()?;

        Ok((pool, pooled_writers, fq_readers))
    }
}

impl Command for TrimmerOpts {
    fn execute(&self) -> Result<()> {
        let (mut pool, mut writers, mut readers) = self.prepare()?;
        let f1 = readers.remove(0);
        let f2 = readers.remove(0);

        'pair: for (r1, r2) in f1.into_records().zip(f2.into_records()) {
            let mut r1 = r1?;
            let mut r2 = r2?;

            // TODO: implement qual-masking vs clipping with enum
            for operation in &self.operations {
                match operation {
                    Operation::TrimQual => {
                        for r in [&mut r1, &mut r2].iter_mut() {
                            let hq_range = base_quality::find_high_quality_bases(
                                r.qual(),
                                self.trim_tail_quality,
                                self.trim_tail_window,
                                self.trim_tail_side,
                            );
                            base_quality::clip_read(r, hq_range);
                        }
                    }
                    Operation::Overlap => {
                        if let Some(overlap) = pair_overlap::find_overlap(
                            r1.seq(),
                            r2.seq(),
                            self.overlap_min_len,
                            self.overlap_max_error_rate,
                        ) {
                            overlap.correct(&mut r1, &mut r2, self.overlap_min_bq_delta, true);
                        }
                    }
                    Operation::Osc => {
                        if let Some(i) = base_quality::identify_trim_point(
                            r1.qual(),
                            self.osc_delta as i32,
                            self.osc_window,
                            self.osc_max_oscillations,
                        ) {
                            // TODO [Brent]: allow clip or mask with flag
                            base_quality::clip_read(&mut r1, 0usize..i);
                        }
                        if let Some(i) = base_quality::identify_trim_point(
                            r2.qual(),
                            self.osc_delta as i32,
                            self.osc_window,
                            self.osc_max_oscillations,
                        ) {
                            // TODO [Brent]: allow clip or mask with flag
                            base_quality::clip_read(&mut r2, 0usize..i);
                        }
                    }
                    Operation::LenFilter => {
                        if r1.seq().len().min(r2.seq().len()) < self.filter_shorter
                            || r1.seq().len().max(r2.seq().len()) < self.filter_longer
                        {
                            info!("Skipping pair with short read");
                            break 'pair;
                        }
                    }
                }
            }

            r1.write(&mut writers[0])?;
            r2.write(&mut writers[1])?;
        }

        writers.into_iter().try_for_each(|w| w.close())?;
        pool.stop_pool()?;

        Ok(())
    }
}
