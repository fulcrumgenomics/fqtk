use crate::commands::command::Command;
use anyhow::{Error, Result};
use clap::{Parser, ValueEnum};
use fgoxide::io::Io;
use fqtk_lib::fastq_stats as stats;
use fqtk_lib::{base_quality, pair_overlap};
use log::info;
use pooled_writer::{bgzf::BgzfCompressor, Pool, PoolBuilder, PooledWriter};
use seq_io::fastq::{Reader as FastqReader, Record};
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::path::Path;
use std::path::PathBuf;

#[derive(ValueEnum, Clone, Debug)]
enum Operation {
    Clip,
    Overlap,
    Osc,
    FilterLen,
}

/// Trimming and overlap correction of paired-end reads
#[derive(Parser, Debug)]
#[command(version)]
pub(crate) struct TrimmerOpts {
    /// Reading/Writing threads
    #[clap(long, short = 't', default_value = "5")]
    threads: usize,

    /// Clip bases with quality < this value.
    #[clap(long, short = 'q', default_value = "20")]
    clip_tail_quality: u8,

    /// Which tail(s) to clip.
    #[clap(long, value_enum, default_value = "end")]
    clip_tail_side: base_quality::Tail,

    /// Window size for moving average when clipping tails.
    #[clap(long, short = 'w', default_value = "20")]
    clip_tail_window: u8,

    /// Level of compression to use to compress outputs.
    #[clap(long, short = 'c', default_value = "5")]
    compression_level: usize,

    /// Length requirement of shorter read. Lengths below this are clipped.
    #[clap(long, short = 'S', default_value = "5")]
    filter_shorter: usize,

    /// Length requirement of longer read. Lengths below this are clipped.
    #[clap(long, short = 'L', default_value = "15")]
    filter_longer: usize,

    /// Size of window to look for oscillations.
    #[clap(long, default_value = "15")]
    osc_window: usize,

    /// Required number of oscillations in a window to trigger trimming/masking.
    #[clap(long, default_value = "4")]
    osc_max_oscillations: usize,

    /// Difference between adjacent bases to be considered on oscillation.
    #[clap(long, default_value = "10")]
    osc_delta: usize,

    /// Minimum difference in base-quality for one read to correct an overlapping
    /// base from the other read.
    #[clap(long, short = 'd', default_value = "15")]
    overlap_min_bq_delta: u8,

    /// Minimum pair overlap length to attempt correction.
    #[clap(long, short = 'l', default_value = "50")]
    overlap_min_length: usize,

    /// Maximum error-rate allowed in the overlap.
    #[clap(long, short = 'e', default_value = "0.02")]
    overlap_max_error_rate: f64,

    /// Hard clip adapter sequences from the reads detected in overlap module.
    /// If this is not specified, the adapter qualities are instead set to `mask_quality`
    #[clap(long, default_value_t = false)]
    overlap_hard_clip_adapters: bool,

    /// Quality value to use as a mask (should be lower than `clip_tail_quality`)
    #[clap(long, default_value = "0")]
    mask_quality: u8,

    /// Order of operations
    #[clap(value_enum, short = 'p', default_value = "overlap", num_args = 1..)]
    operations: Vec<Operation>,

    /// The paths for the 2 output FASTQs.
    #[clap(long, short = 'o', required = true, num_args = 2)]
    output: Vec<PathBuf>,

    /// Fastqs file for Read1 and Read2
    #[clap(long, short = 'i', required = true, num_args = 2)]
    input: Vec<PathBuf>,
}

const BUFFER_SIZE: usize = 1024 * 1024;

/// Type alias to prevent clippy complaining about type complexity
type VecOfReaders = Vec<Box<dyn BufRead + Send>>;
type VecOfFqReaders = Vec<FastqReader<Box<dyn BufRead + Send>>>;

fn create_writer<P: AsRef<Path>>(name: P) -> Result<BufWriter<File>, Error> {
    Ok(BufWriter::new(File::create(name)?))
}

fn check_extension(p: &Path) -> bool {
    let ext = p.extension().map_or("", |v| v.to_str().unwrap_or(""));
    ["bgz", "gz"].contains(&ext)
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

        let output: Vec<_> = self
            .output
            .iter()
            .map(|p| {
                if !check_extension(p) {
                    log::warn!(
                        "Output file {} does not end with .gz or .bgz. Writing bgzipped output to {} instead.",
                        p.display(),
                        p.with_extension("gz").display()
                    );
                    p.with_extension("gz")
                } else {
                    p.clone()
                }
            })
            .collect();

        let writers = vec![create_writer(&output[0])?, create_writer(&output[1])?];

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

        let mut stats = stats::Stats::new();

        'pair: for (r1, r2) in f1.into_records().zip(f2.into_records()) {
            let mut r1 = r1?;
            let mut r2 = r2?;

            stats.update_length(r1.seq.len(), stats::When::Pre, stats::ReadI::R1);
            stats.update_length(r2.seq.len(), stats::When::Pre, stats::ReadI::R2);

            for operation in &self.operations {
                match operation {
                    Operation::Clip => {
                        for r in [&mut r1, &mut r2].iter_mut() {
                            let hq_range = base_quality::find_high_quality_bases(
                                r.qual(),
                                self.clip_tail_quality,
                                self.clip_tail_window,
                                self.clip_tail_side,
                            );
                            // this is hard clip so we send None
                            base_quality::mask_read(r, hq_range, None);
                        }
                    }
                    Operation::Overlap => {
                        if let Some(overlap) = pair_overlap::find_overlap(
                            r1.seq(),
                            r2.seq(),
                            self.overlap_min_length,
                            self.overlap_max_error_rate,
                        ) {
                            log::debug!(
                                "found overlap in pair: {} shift: {}, overlap: {}, adapter: {}",
                                r1.id().unwrap_or("read"),
                                overlap.shift,
                                overlap.overlap,
                                overlap.adapter
                            );
                            stats.overlap_stats.update(overlap);
                            let corrections = overlap.correct(
                                &mut r1,
                                &mut r2,
                                self.overlap_min_bq_delta,
                                if self.overlap_hard_clip_adapters {
                                    None
                                } else {
                                    Some(self.mask_quality)
                                },
                            );
                            log::debug!("corrections: {:?}", corrections);
                            stats.overlap_stats.update_corrections(corrections.0, stats::ReadI::R1);
                            stats.overlap_stats.update_corrections(corrections.1, stats::ReadI::R2);
                        }
                    }
                    Operation::Osc => {
                        if let Some(i) = base_quality::identify_trim_point(
                            r1.qual(),
                            self.osc_delta as i32,
                            self.osc_window,
                            self.osc_max_oscillations,
                        ) {
                            base_quality::mask_read(&mut r1, 0usize..i, Some(self.mask_quality));
                            stats.update_oscillations(1, stats::ReadI::R1);
                        }
                        if let Some(i) = base_quality::identify_trim_point(
                            r2.qual(),
                            self.osc_delta as i32,
                            self.osc_window,
                            self.osc_max_oscillations,
                        ) {
                            base_quality::mask_read(&mut r2, 0usize..i, Some(self.mask_quality));
                            stats.update_oscillations(1, stats::ReadI::R2);
                        }
                    }
                    Operation::FilterLen => {
                        if r1.seq().len().min(r2.seq().len()) < self.filter_shorter
                            || r1.seq().len().max(r2.seq().len()) < self.filter_longer
                        {
                            info!("Skipping pair with short read");
                            stats.increment_length_filter();
                            continue 'pair;
                        }
                    }
                }
            }

            stats.update_length(r1.seq.len(), stats::When::Post, stats::ReadI::R1);
            stats.update_length(r2.seq.len(), stats::When::Post, stats::ReadI::R2);

            r1.write(&mut writers[0])?;
            r2.write(&mut writers[1])?;
        }

        writeln!(std::io::stderr(), "{}", stats)?;

        writers.into_iter().try_for_each(|w| w.close())?;
        pool.stop_pool()?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_check_extension_gz() {
        let path = PathBuf::from("test_file.gz");
        assert_eq!(
            check_extension(&path),
            true,
            "The function should return true for .gz extension"
        );
    }

    #[test]
    fn test_check_extension_txt() {
        let path = PathBuf::from("test_file.txt");
        assert_eq!(
            check_extension(&path),
            false,
            "The function should return false for .txt extension"
        );
    }
}
