use crate::commands::command::Command;
use anyhow::Result;
use clap::Parser;
use log::info;
use std::fs;
use std::path::PathBuf;

/// Trimming and overlap correction of paired-end reads
#[derive(Parser, Debug)]
#[command(version)]
pub(crate) struct Trimmer {
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
    fastqs: PathBuf,
}

impl Trimmer {
    fn prepare(&self) -> Result<()> {
        if !self.output.exists() {
            info!("Output directory {:#?} didn't exist, creating it.", self.output);
            fs::create_dir_all(&self.output)?;
        }
        // see here for pooled writer: https://docs.rs/pooled-writer/0.3.0/pooled_writer/
        //
        // and https://github.com/fulcrumgenomics/fqtk/blob/ae91e90ecc86826fc632837bb982798a7e6b6f7a/src/bin/commands/demux.rs#L804
        // and https://github.com/fulcrumgenomics/fqtk/blob/ae91e90ecc86826fc632837bb982798a7e6b6f7a/src/bin/commands/demux.rs#L850-L851
        // for reader
        Ok(())
    }
}

impl Command for Trimmer {
    fn execute(&self) -> Result<()> {
        self.prepare()?;
        Ok(())
    }
}
