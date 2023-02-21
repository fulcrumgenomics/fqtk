use crate::commands::command::Command;
use anyhow::{anyhow, Result};

use clap::Parser;
use fgoxide::io::Io;
use seq_io::fastq::write_to;
use seq_io::fastq::OwnedRecord;
use seq_io::fastq::Reader as FastqReader;
use seq_io::fastq::Record;

use fgoxide::iter::IntoChunkedReadAheadIterator;

use seq_io::fastq::RefRecord;
use std::io::Write;
use std::io::{BufRead, BufWriter};
use std::path::PathBuf;

/// Type aliases to prevent clippy complaining about type complexity
type VecOfReaders = Vec<Box<dyn BufRead + Send>>;
type VecOfWriters = Vec<BufWriter<Box<dyn Write + Send>>>;
/// The buffer size to use for readers and writers
const BUFFER_SIZE: usize = 1024 * 1024;

struct TrimReadIterator {
    /// The FASTQ file that is being read from by ``Self``.
    source: FastqReader<Box<dyn BufRead + Send>>,
    /// The window size over which to look for deviations in qualities.
    deviation_window_size: usize,
    /// The average deviation threshold above which to begin trimming.
    average_deviation_threshold: usize,
}
impl TrimReadIterator {
    /// Filters reads for large changes in quality scores that can occur towards the end of long
    /// Illumina reads.
    /// TODO - expand this documentation
    pub fn left_to_right_quality_deviation_filter(
        r: &RefRecord,
        window_size: usize,
        average_deviation_threshold: usize,
    ) -> Option<usize> {
        let qualities = r.qual();
        if window_size > qualities.len() {
            return None;
        }
        let mut deviations = Vec::with_capacity(qualities.len() - 1);
        for i in 0..qualities.len() - 1 {
            let j = i + 1;
            deviations.push(u8::abs_diff(qualities[i], qualities[j]));
        }
        let maximum_sum = window_size * average_deviation_threshold;

        let mut summed_deviations: usize =
            deviations[0..window_size].iter().map(|&v| v as usize).sum();

        for i in window_size..deviations.len() {
            if summed_deviations > maximum_sum {
                for (j, &d) in deviations.iter().enumerate().take(i).skip(i - window_size) {
                    if d as usize > average_deviation_threshold {
                        return Some(j + 1);
                    }
                }
            }
            if i != deviations.len() - 1 {
                summed_deviations = summed_deviations + deviations[i + 1] as usize
                    - deviations[i - window_size] as usize;
            }
        }
        None
    }
}

impl Iterator for TrimReadIterator {
    type Item = OwnedRecord;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(rec) = self.source.next() {
            let next_fq_rec = rec.expect("Unexpected error parsing FASTQs.");
            let deviation_result = Self::left_to_right_quality_deviation_filter(
                &next_fq_rec,
                self.deviation_window_size,
                self.average_deviation_threshold,
            );
            if let Some(mask_point) = deviation_result {
                Some(OwnedRecord {
                    head: next_fq_rec.head().to_owned(),
                    seq: next_fq_rec.seq()[0..mask_point].to_owned(),
                    qual: next_fq_rec.qual()[0..mask_point].to_owned(),
                })
            } else {
                Some(next_fq_rec.to_owned_record())
            }
        } else {
            None
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
// QualTrim (main class) and it's impls
////////////////////////////////////////////////////////////////////////////////

/// Performs quality trimming on FASTQs.
/// TODO - flesh out this documentation
#[derive(Parser, Debug)]
#[command(version)]
pub(crate) struct QualTrim {
    /// One or two input FASTQ file paths. Number of files provided should match the number of
    /// output files provided.
    #[clap(long, short = 'i', required = true, num_args = 1..=2)]
    inputs: Vec<PathBuf>,

    /// One or two output FASTQ file paths. Number of files provided should match the number of
    /// input files provided.
    #[clap(long, short = 'o', required = true, num_args = 1..=2)]
    outputs: Vec<PathBuf>,

    /// Minimum read length of the output FASTQ reads in order for the reads to be output.
    #[clap(long, default_value = "1", num_args = 1..=2)]
    min_read_lengths: Vec<usize>,

    /// The window size to use when filtering based on quality score oscilations that occur in
    /// longer illumina reads
    #[clap(long, default_value = "15")]
    deviation_window_size: usize,

    /// The average deviation at which a the sliding deviation window will be trimmed from the
    /// output reads.
    #[clap(long, default_value = "10")]
    average_deviation_threshold: usize,
}

impl QualTrim {
    /// Checks that inputs to demux are valid and returns open file handles for the inputs.
    /// Checks:
    ///     - That the number of input files and number of read structs provided are the same
    ///     - That the output directory is not read-only
    ///     - That the input files exist
    ///     - That the input files have read permissions.
    fn validate_and_prepare_inputs(&self) -> Result<(VecOfReaders, VecOfWriters)> {
        let mut constraint_errors = vec![];

        for input in &self.inputs {
            if !input.exists() {
                constraint_errors.push(format!("Provided input file {input:#?} doesn't exist"));
            }
        }
        // Attempt to open the files for reading.
        let fgio = Io::new(5, BUFFER_SIZE);
        let fq_readers_result = self
            .inputs
            .iter()
            .map(|p| fgio.new_reader(p))
            .collect::<Result<VecOfReaders, fgoxide::FgError>>();
        if let Err(e) = &fq_readers_result {
            constraint_errors.push(format!("Error opening input files for reading: {e}"));
        }

        let fq_writers_result = self
            .outputs
            .iter()
            .map(|p| fgio.new_writer(p))
            .collect::<Result<VecOfWriters, fgoxide::FgError>>();

        if constraint_errors.is_empty() {
            return Ok((fq_readers_result?, fq_writers_result?));
        }
        let mut details = "Inputs failed validation!\n".to_owned();
        for error_reason in constraint_errors {
            details.push_str(&format!("    - {error_reason}\n"));
        }
        Err(anyhow!("The following errors with the input(s) were detected:\n{}", details))
    }
}

impl Command for QualTrim {
    fn execute(&self) -> anyhow::Result<()> {
        assert_eq!(
            self.inputs.len(),
            self.outputs.len(),
            "Inputs and outputs should have the same length!"
        );
        assert!(self.inputs.len() <= 2, "Cannot provide more than 2 read files at this time.");

        let mut read_length_mins_sorted = if self.min_read_lengths.len() == 1 {
            vec![self.min_read_lengths[0]; self.inputs.len()]
        } else {
            self.min_read_lengths.clone()
        };
        read_length_mins_sorted.sort_unstable();

        let (inputs, mut outputs) = self.validate_and_prepare_inputs()?;

        let mut fq_sources = inputs
            .into_iter()
            .map(|fq| {
                TrimReadIterator {
                    source: FastqReader::with_capacity(fq, BUFFER_SIZE),
                    deviation_window_size: self.deviation_window_size,
                    average_deviation_threshold: self.average_deviation_threshold,
                }
                .read_ahead(1000, 1000)
            })
            .collect::<Vec<_>>();

        loop {
            let mut next_reads = Vec::with_capacity(fq_sources.len());
            for iter in &mut fq_sources {
                if let Some(rec) = iter.next() {
                    next_reads.push(rec);
                }
            }
            if next_reads.is_empty() {
                break;
            }
            assert_eq!(
                next_reads.len(),
                fq_sources.len(),
                "FASTQ sources out of sync at records {next_reads:?}"
            );

            let mut observed_read_lengths =
                next_reads.iter().map(|r| r.seq.len()).collect::<Vec<_>>();
            observed_read_lengths.sort_unstable();

            if observed_read_lengths
                .iter()
                .zip(read_length_mins_sorted.iter())
                .any(|(&observed, &minimum)| observed < minimum)
            {
                continue;
            }

            for (writer, read) in outputs.iter_mut().zip(next_reads.iter()) {
                write_to(writer, &read.head, &read.seq, &read.qual)?;
            }
        }
        Ok(())
    }
}
