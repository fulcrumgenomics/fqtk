use crate::commands::command::Command;
use anyhow::anyhow;
use clap::Parser;
use fgoxide::io::Io;
use fqtk_lib::samples::SampleGroup;
use log::info;
use read_structure::ReadStructure;
use std::io::Read;
use std::str::FromStr;
use std::{fs, io::BufReader, path::PathBuf};

/// the breakdown of threads allocated to each subtask of demultiplexing.
#[allow(dead_code)]
#[derive(Copy, Clone, Debug)]
struct ThreadBreakdown {
    /// The number of threads used to demultiplex reads / assign reads to samples
    demux_threads: usize,
    /// The number of threads used to compress output before writing
    compressor_threads: usize,
    /// The number of threads used to write output to file.
    writer_threads: usize,
}

impl FromStr for ThreadBreakdown {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let threads = s.parse::<u16>().map_err(|e| anyhow!(e))?;
        if threads < 3 {
            return Err(anyhow!("Threads provided {} was too low!", threads));
        }

        // TODO - Optimize the breakdown of the thread allocation
        // TODO - figure out how to handle different number of samples to be written.
        //
        let breakdown = match threads {
            3u16 => ThreadBreakdown { demux_threads: 1, compressor_threads: 1, writer_threads: 1 },
            4u16 => ThreadBreakdown { demux_threads: 2, compressor_threads: 1, writer_threads: 1 },
            5u16 => ThreadBreakdown { demux_threads: 2, compressor_threads: 2, writer_threads: 1 },
            6u16 => ThreadBreakdown { demux_threads: 2, compressor_threads: 3, writer_threads: 1 },
            7u16 => ThreadBreakdown { demux_threads: 3, compressor_threads: 3, writer_threads: 1 },
            8u16 => ThreadBreakdown { demux_threads: 4, compressor_threads: 3, writer_threads: 1 },
            9u16 => ThreadBreakdown { demux_threads: 5, compressor_threads: 3, writer_threads: 1 },
            10u16 => ThreadBreakdown { demux_threads: 6, compressor_threads: 3, writer_threads: 1 },
            11u16 => ThreadBreakdown { demux_threads: 6, compressor_threads: 4, writer_threads: 1 },
            12u16 => ThreadBreakdown { demux_threads: 7, compressor_threads: 4, writer_threads: 2 },
            13u16 => ThreadBreakdown { demux_threads: 8, compressor_threads: 4, writer_threads: 2 },
            14u16 => ThreadBreakdown { demux_threads: 9, compressor_threads: 4, writer_threads: 2 },
            15u16 => {
                ThreadBreakdown { demux_threads: 10, compressor_threads: 4, writer_threads: 2 }
            }
            16u16 => {
                ThreadBreakdown { demux_threads: 11, compressor_threads: 4, writer_threads: 2 }
            }
            x => ThreadBreakdown {
                demux_threads: (x - 7) as usize,
                compressor_threads: 4,
                writer_threads: 3,
            },
        };
        Ok(breakdown)
    }
}

/// Demultiplexes FASTQ files.
#[derive(Parser, Debug)]
pub(crate) struct Demux {
    /// One or more input fastq files each corresponding to a sub-read (ex. index-read, read-one,
    /// read-two, fragment).
    #[clap(long, short = 'i', required = true, num_args = 1..)]
    inputs: Vec<PathBuf>,

    /// The read structure for each of the FASTQs.
    #[clap(long, short = 'r' , required = true, num_args = 1..)]
    read_structures: Vec<ReadStructure>,

    /// A file containing the metadata about the samples.
    #[clap(long, short = 'x', required = true)]
    metadata: PathBuf,

    /// The output directory in which to place sample FASTQs.
    #[clap(long, short = 'o', required = true)]
    output: PathBuf,

    /// Output FASTQ file name for the unmatched records.
    #[clap(long, short = 'u', default_value = "unmatched")]
    unmatched: String,

    /// Maximum mismatches for a barcode to be considered a match.
    #[clap(long, default_value = "1")]
    max_mismatches: usize,

    /// Minimum difference between number of mismatches in the best and second best barcodes for a
    /// barcode to be considered a match.
    #[clap(long, default_value = "2")]
    min_mismatch_delta: usize,

    /// The number of threads to use. Cannot be less than 3.
    #[clap(long, short = 't', default_value = "3")]
    threads: ThreadBreakdown,
}

impl Demux {
    /// Attempts to open input files for reading.
    ///
    /// # Errors
    ///     - Will fail if the files do not exist (which should be checked before calling this function).
    ///     - Will fail if the files do not have read permissions.
    fn open_inputs_for_reading(&self) -> anyhow::Result<Vec<BufReader<Box<dyn Read>>>> {
        let io = Io::default();
        self.inputs
            .iter()
            .map(|p| io.new_reader(p).map_err(|e| anyhow!(e)))
            .collect::<Result<Vec<_>, anyhow::Error>>()
    }

    /// Checks that inputs to demux are valid.
    /// Checks:
    ///     - That the number of input files and number of read structs provided are the same
    ///     - That the output directory is not read-only
    ///     - That the input files exist
    fn validate_inputs(&self) -> anyhow::Result<()> {
        let mut constraint_errors = vec![];

        if self.inputs.len() != self.read_structures.len() {
            let preamble = "The same number of read structures should be given as FASTQs";
            let specifics = format!(
                "(currently {} vs {} respectively)",
                self.read_structures.len(),
                self.inputs.len()
            );
            constraint_errors.push(format!("{preamble} {specifics}"));
        }

        if self.output.metadata()?.permissions().readonly() {
            constraint_errors
                .push(format!("Ouput directory {:#?} cannot be read-only", self.output));
        }

        for input in &self.inputs {
            if !input.exists() {
                constraint_errors.push(format!("Provided input file {:#?} doesn't exist", input));
            }
        }

        if constraint_errors.is_empty() {
            Ok(())
        } else {
            let mut error_string = "Inputs failed validation!\n".to_owned();
            for error_reason in constraint_errors {
                error_string.push_str(&format!("    - {}\n", error_reason));
            }
            Err(anyhow!("{}Input Validation Error!", error_string))
        }
    }
}

impl Command for Demux {
    // (For now all this does is check the inputs to the demux command, but will eventually execute
    // TODO - remove this disclaimer once working core is implemented.)
    /// Excecutes the demux command
    fn execute(&self) -> anyhow::Result<()> {
        self.validate_inputs()?;
        let _fq_readers = self.open_inputs_for_reading()?;
        let _sample_group = SampleGroup::from_file(&self.metadata)?;
        if let Some(parent_path) = self.output.parent() {
            assert!(
                parent_path.exists(),
                "Parent path ({:#?}) to output directory ({:#?}) didn't exist!\nexiting...",
                parent_path,
                self.output,
            );
        }
        if !self.output.exists() {
            info!("Output directory {:#?} didn't exist, creating it.", self.output);
            fs::create_dir(&self.output)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use rstest::rstest;
    use std::str::FromStr;

    use super::*;
    use tempfile::TempDir;

    const SAMPLE1_BARCODE: &str = "GATTGGG";

    fn fq_lines_from_bases(prefix: &str, records_bases: &[&str]) -> Vec<String> {
        let mut result = Vec::with_capacity(records_bases.len() * 4);
        for (i, &bases) in records_bases.iter().enumerate() {
            result.push(format!("Example{}_{}", prefix, i));
            result.push(bases.to_owned());
            result.push("+".to_owned());
            result.push(";".repeat(bases.len()));
        }
        result
    }

    fn fastq_file(tmpdir: &TempDir, prefix: &str, records_bases: &[&str]) -> PathBuf {
        let io = Io::default();

        let path = tmpdir.path().join(format!("{prefix}.fastq"));
        let fastq_lines = fq_lines_from_bases(prefix, records_bases);
        io.write_lines(&path, fastq_lines).unwrap();

        path
    }

    fn metadata_lines_from_barcodes(barcodes: &[&str]) -> Vec<String> {
        let mut result = Vec::with_capacity(barcodes.len() + 1);
        result.push("name\tbarcode".to_owned());
        for (i, &barcode) in barcodes.iter().enumerate() {
            result.push(format!("Sample{:04}\t{}", i, barcode));
        }
        result
    }

    fn metadata_file(tmpdir: &TempDir, barcodes: &[&str]) -> PathBuf {
        let io = Io::default();

        let path = tmpdir.path().join("metadata.tsv");
        let metadata_lines = metadata_lines_from_barcodes(barcodes);
        io.write_lines(&path, metadata_lines).unwrap();

        path
    }

    fn metadata(tmpdir: &TempDir) -> PathBuf {
        metadata_file(tmpdir, &[SAMPLE1_BARCODE])
    }

    fn read1(tmpdir: &TempDir) -> PathBuf {
        fastq_file(tmpdir, "read1", &["GATTACA"])
    }

    fn read2(tmpdir: &TempDir) -> PathBuf {
        fastq_file(tmpdir, "read2", &["TAGGATTA"])
    }

    fn index1(tmpdir: &TempDir) -> PathBuf {
        fastq_file(tmpdir, "index1", &[&SAMPLE1_BARCODE[0..3]])
    }

    fn index2(tmpdir: &TempDir) -> PathBuf {
        fastq_file(tmpdir, "index2", &[&SAMPLE1_BARCODE[3..]])
    }

    // ############################################################################################
    // Test that ``Demux:: execute`` can succeed.
    // ############################################################################################
    #[rstest]
    fn validate_inputs_can_succeed() {
        let tmpdir = TempDir::new().unwrap();
        let input_files = vec![read1(&tmpdir), read2(&tmpdir), index1(&tmpdir), index2(&tmpdir)];
        let metadata = metadata(&tmpdir);

        let read_structures = vec![
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
        ];

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            metadata,
            output: tmpdir.path().to_path_buf(),
            unmatched: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: ThreadBreakdown::from_str("3").unwrap(),
        };
        demux_inputs.execute().unwrap();
    }

    #[rstest]
    #[should_panic(expected = "The same number of read structures should be given as FASTQs")]
    #[case(vec![
        ReadStructure::from_str("+T").unwrap(),
        ReadStructure::from_str("+T").unwrap(),
        ReadStructure::from_str("+B").unwrap(),
        ])]
    #[should_panic(expected = "The same number of read structures should be given as FASTQs")]
    #[case(vec![
        ReadStructure::from_str("+T").unwrap(),
        ReadStructure::from_str("+T").unwrap(),
        ReadStructure::from_str("+B").unwrap(),
        ReadStructure::from_str("+B").unwrap(),
        ReadStructure::from_str("+B").unwrap(),
    ])]
    fn test_different_number_of_read_structs_and_inputs_fails(
        #[case] read_structures: Vec<ReadStructure>,
    ) {
        let tmpdir = TempDir::new().unwrap();
        let input_files = vec![read1(&tmpdir), read2(&tmpdir), index1(&tmpdir), index2(&tmpdir)];
        let metadata = metadata(&tmpdir);

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            metadata,
            output: tmpdir.path().to_path_buf(),
            unmatched: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: ThreadBreakdown::from_str("3").unwrap(),
        };
        demux_inputs.execute().unwrap();
    }

    #[rstest]
    #[should_panic(expected = "cannot be read-only")]
    fn test_read_only_output_dir_fails() {
        let tmpdir = TempDir::new().unwrap();
        let read_structures = vec![
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
        ];

        let input_files = vec![read1(&tmpdir), read2(&tmpdir), index1(&tmpdir), index2(&tmpdir)];
        let metadata = metadata(&tmpdir);
        let mut permissions = tmpdir.path().metadata().unwrap().permissions();
        permissions.set_readonly(true);
        fs::set_permissions(tmpdir.path(), permissions.clone()).unwrap();

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            metadata,
            output: tmpdir.path().to_path_buf(),
            unmatched: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: ThreadBreakdown::from_str("3").unwrap(),
        };
        let demux_result = demux_inputs.execute();
        permissions.set_readonly(false);
        fs::set_permissions(tmpdir.path(), permissions).unwrap();
        demux_result.unwrap();
    }

    #[rstest]
    #[should_panic(expected = "doesn't exist")]
    fn test_inputs_doesnt_exist_fails() {
        let tmpdir = TempDir::new().unwrap();
        let read_structures = vec![
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
        ];

        let input_files = vec![
            tmpdir.path().join("this_file_does_not_exist.fq"),
            read2(&tmpdir),
            index1(&tmpdir),
            index2(&tmpdir),
        ];
        let metadata = metadata(&tmpdir);

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            metadata,
            output: tmpdir.path().to_path_buf(),
            unmatched: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: ThreadBreakdown::from_str("3").unwrap(),
        };
        demux_inputs.execute().unwrap();
    }

    #[rstest]
    #[should_panic(expected = "Threads provided 2 was too low!")]
    fn test_too_few_threads_fails() {
        let tmpdir = TempDir::new().unwrap();
        let read_structures = vec![
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+T").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
            ReadStructure::from_str("+B").unwrap(),
        ];

        let input_files = vec![read1(&tmpdir), read2(&tmpdir), index1(&tmpdir), index2(&tmpdir)];
        let metadata = metadata(&tmpdir);

        let demux_inputs = Demux {
            inputs: input_files,
            read_structures,
            metadata,
            output: tmpdir.path().to_path_buf(),
            unmatched: "unmatched".to_owned(),
            max_mismatches: 1,
            min_mismatch_delta: 2,
            threads: ThreadBreakdown::from_str("2").unwrap(),
        };
        demux_inputs.execute().unwrap();
    }
}
