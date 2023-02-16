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
use std::cmp;
use std::io::Write;
use std::io::{BufRead, BufWriter};
use std::path::PathBuf;

/// Type alias to prevent clippy complaining about type complexity
type VecOfReaders = Vec<Box<dyn BufRead + Send>>;
type VecOfWriters = Vec<BufWriter<Box<dyn Write + Send>>>;
const BUFFER_SIZE: usize = 1024 * 1024;

struct TrimReadIterator {
    source: FastqReader<Box<dyn BufRead + Send>>,
    run_eamss: bool,
    back_to_front_window_size: usize,
    front_to_back_window_size: usize,
}
impl TrimReadIterator {
    // Adjust the thresholds to the ASCII bytes that seq_io returns.
    const EAMSS_M2_GE_THRESHOLD: u8 = 30 + 33;
    const EAMSS_S1_LT_THRESHOLD: u8 = 15 + 33;

    // EAMSS is an Illumina Developed Algorithm for detecting reads whose quality has deteriorated towards
    // their end and revising the quality to the masking quality (2) if this is the case.  This algorithm
    // works as follows (with one exception):
    //
    //     Start at the end (high indices, at the right below) of the read and calculate an EAMSS
    //     tally at each location as follow:
    //     if(quality[i] < 15) tally += 1
    //     if(quality[i] >= 15 and < 30) tally = tally
    //     if(quality[i] >= 30) tally -= 2
    //
    //
    // For each location, keep track of this tally (e.g.)
    // Read Starts at <- this end
    // Cycle:       1  2  3  4  5  6  7  8  9
    // Bases:       A  C  T  G  G  G  T  C  A
    // Qualities:   32 32 16 15 8  10 32 2  2
    // Cycle Score: -2 -2 0  0  1  1  -2 1  1           //The EAMSS Score determined for this cycle alone
    // EAMSS TALLY: 0  0  2  2  2  1  0  2  1
    // X - Earliest instance of Max-Score
    // <p/>
    // You must keep track of the maximum EAMSS tally (in this case 2) and the earliest(lowest) cycle at which
    // it occurs.  If and only if, the max EAMSS tally >= 1 then from there until the end(highest cycle) of the
    // read reassign these qualities as 2 (the masking quality).  The output qualities would therefore be
    // transformed from:
    // <p/>
    // Original Qualities: 32 32 16 15 8  10 32 2  2    to
    // Final    Qualities: 32 32 2  2  2  2  2  2  2
    // X - Earliest instance of max-tally/end of masking
    // <p/>
    // IMPORTANT:
    // The one exception is: If the max EAMSS Tally is preceded by a long string of G basecalls (10 or more, with a single basecall exception per10 bases)
    // then the masking continues to the beginning of that string of G's. E.g.:
    // <p/>
    // Cycle:       1  2  3  4  5  6  7  8   9  10 11 12 13 14 15 16 17 18
    // Bases:       C  T  A  C  A  G  A  G   G  G  G  G  G  G  G  C  A  T
    // Qualities:   30 22 26 27 28 30 7  34  20 19 38 15 32 32 10 4  2  5
    // Cycle Score: -2  0  0  0  0 -2 1  -2  0  0  -2 0  -2 -2  1 1  1  1
    // EAMSS TALLY: -2 -5 -5 -5 -5 -5 -3 -4 -2 -2  -2 0   0  2  4 3  2  1
    // X- Earliest instance of Max-Tally
    // <p/>
    // Resulting Transformation:
    // Bases:                C  T  A  C  A  G  A   G   G  G  G  G  G  G  G  C  A  T
    // Original Qualities:   30 22 26 27 28 30 7  34  20 19 38 15 32 32 10  4  2  5
    // Final    Qualities:   30 22 26 27 28  2 2   2   2  2  2  2  2  2  2  2  2  2
    // X- Earliest instance of Max-Tally
    // X - Start of EAMSS masking due to G-Run
    // <p/>
    // To further clarify the exception rule here are a few examples:
    // A C G A C G G G G G G G G G G G G G G G G G G G G A C T
    // X - Earliest instance of Max-Tally
    // X - Start of EAMSS masking (with a two base call jump because we have 20 bases in the run already)
    // <p/>
    // T T G G A G G G G G G G G G G G G G G G G G G A G A C T
    // X - Earliest instance of Max-Tally
    // X - We can skip this A as well as the earlier A because we have 20 or more bases in the run already
    // X - Start of EAMSS masking (with a two base call jump because we have 20 bases in the run)
    // <p/>
    // T T G G G A A G G G G G G G G G G G G G G G G G G T T A T
    // X - Earliest instance of Max-Tally
    // X X - WE can skip these bases because the first A counts as the first skip and as far as the length of the string of G's is
    // concerned, these are both counted like G's
    // X - This A is the 20th base in the string of G's and therefore can be skipped
    // X - Note that the A's previous to the G's are only included because there are G's further on that are within the number
    // of allowable exceptions away (i.e. 2 in this instance), if there were NO G's after the A's you CANNOT count the A's
    // as part of the G strings (even if no exceptions have previously occured) In other words, the end of the string of G's
    // MUST end in a G not an "exception"
    // <p/>
    // However, if the max-tally occurs to the right of the run of Gs then this is still part of the string of G's but does count towards
    // the number of exceptions allowable
    // (e.g.)
    // T T G G G G G G G G G G A C G
    // X - Earliest instance of Max-tally
    // The first index CAN be considered as an exception, the above would be masked to
    // the following point:
    // T T G G G G G G G G G G A C G
    // X - End of EAMSS masking due to G-Run
    // <p/>
    // To sum up the final points, a string of G's CAN START with an exception but CANNOT END in an exception.
    pub fn run_eamss_on_read(r: &RefRecord) -> Option<usize> {
        let mut eamss_tally: isize = 0isize;
        let mut max_tally: isize = isize::MIN;
        let mut index_of_max: Option<usize> = None;

        for i in (0..r.seq().len()).rev() {
            eamss_tally = match r.qual()[i] {
                x if x >= Self::EAMSS_M2_GE_THRESHOLD => eamss_tally - 2,
                x if x < Self::EAMSS_S1_LT_THRESHOLD => eamss_tally + 1,
                _ => eamss_tally,
            };
            if eamss_tally >= max_tally {
                index_of_max.replace(i);
                max_tally = eamss_tally;
            }
        }

        if max_tally >= 1 {
            let mut num_gs = 0;
            let mut exceptions = 0;

            if let Some(mut i) = index_of_max {
                loop {
                    if r.seq()[i] == b'G' {
                        num_gs += 1;
                    } else if let Some(skip) = Self::skip_by(i, num_gs, exceptions, r.seq()) {
                        exceptions += skip;
                        num_gs += skip;
                        i -= skip - 1;
                    } else {
                        break;
                    }
                    if i == 0 {
                        break;
                    }
                    i -= 1;
                }
                if num_gs >= 10 {
                    index_of_max = index_of_max.map(|index| index + 1 - num_gs);
                }
            }
            index_of_max
        } else {
            None
        }
    }

    pub fn skip_by(
        index: usize,
        num_gs: usize,
        previous_exceptions: usize,
        bases: &[u8],
    ) -> Option<usize> {
        let mut skip = None;

        for backup in 1..=index {
            let exception_limit = cmp::max((num_gs + backup) / 10, 1);
            if previous_exceptions + backup > exception_limit {
                break;
            }
            if bases[index - backup] == b'G' {
                skip = Some(backup);
                break;
            }
        }
        skip
    }
}

impl Iterator for TrimReadIterator {
    type Item = (bool, OwnedRecord);

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(rec) = self.source.next() {
            let next_fq_rec = rec.expect("Unexpected error parsing FASTQs.");
            let eamss_result = Self::run_eamss_on_read(&next_fq_rec);
            let mut output_rec = next_fq_rec.to_owned_record();
            let mut eamss_filtered = false;
            if let Some(mask_point) = eamss_result {
                for i in mask_point..output_rec.seq.len() {
                    output_rec.seq[i] = u8::to_ascii_lowercase(&output_rec.seq[i]);
                }
                eamss_filtered = true;
            }
            Some((eamss_filtered, output_rec))
        } else {
            None
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
// QualTrim (main class) and it's impls
////////////////////////////////////////////////////////////////////////////////

/// Performs quality trimming on FASTQs.
#[derive(Parser, Debug)]
#[command(version)]
pub(crate) struct QualTrim {
    /// One or two input FASTQ file paths. Number of files provided should match the number of
    /// output files provided.
    #[clap(long, short = 'i', required = true, num_args = 1..2)]
    inputs: Vec<PathBuf>,

    /// One or two output FASTQ file paths. Number of files provided should match the number of
    /// input files provided.
    #[clap(long, short = 'o', required = true, num_args = 1..2)]
    outputs: Vec<PathBuf>,

    ///
    #[clap(long)]
    run_eamss: bool,

    ///
    #[clap(long, required = false)]
    back_to_front_window_size: Option<usize>,

    ///
    #[clap(long, required = false)]
    front_to_back_window_size: Option<usize>,

    /// Minimum read length of the output FASTQ reads in order for the reads to be output.
    #[clap(long, default_value = "1")]
    min_read_lengths: Vec<usize>,
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
        let (inputs, mut outputs) = self.validate_and_prepare_inputs()?;

        let mut fq_sources = inputs
            .into_iter()
            .map(|fq| {
                TrimReadIterator {
                    source: FastqReader::with_capacity(fq, BUFFER_SIZE),
                    run_eamss: true,
                    back_to_front_window_size: 5,
                    front_to_back_window_size: 5,
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
            if next_reads.iter().any(|(is_filtered, _read)| *is_filtered) {
                for (writer, (_is_filtered, read)) in outputs.iter_mut().zip(next_reads.iter()) {
                    write_to(writer, &read.head, &read.seq, &read.qual)?;
                }
            }
        }
        Ok(())
    }
}
