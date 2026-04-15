use crate::commands::command::Command;
use anyhow::{Result, ensure};
use clap::Parser;
use fgoxide::io::Io;
use log::info;
use pooled_writer::{Pool, PoolBuilder, PooledWriter, bgzf::BgzfCompressor};
use rand::Rng;
use rand::SeedableRng;
use rand_chacha::ChaCha8Rng;
use seq_io::fastq::Reader as FastqReader;
use seq_io::fastq::Record;
use std::collections::hash_map::DefaultHasher;
use std::fs::File;
use std::hash::{Hash, Hasher};
use std::io::{BufRead, BufWriter};
use std::path::PathBuf;

const BUFFER_SIZE: usize = 1024 * 1024;

/// Formats a u64 with comma separators (e.g. 1,234,567).
fn fmt_count(n: u64) -> String {
    let s = n.to_string();
    let mut result = String::with_capacity(s.len() + s.len() / 3);
    for (i, ch) in s.chars().enumerate() {
        if i > 0 && (s.len() - i) % 3 == 0 {
            result.push(',');
        }
        result.push(ch);
    }
    result
}

/// Subsamples reads from one or more synchronized FASTQ files.
///
/// Reads one or more FASTQ files (e.g. paired-end R1 and R2) and writes a
/// random subset of reads to output files. All input files must contain the
/// same number of reads in the same order; each read is either kept or
/// discarded across all files simultaneously.
///
/// Output files are named `{output}.R1.fq.gz`, `{output}.R2.fq.gz`, etc.
/// and are always BGZF compressed.
///
/// Each read is independently retained with probability equal to `--fraction`,
/// giving an approximate subsample without needing to know the total read count
/// upfront. When no explicit `--seed` is provided, a deterministic seed is
/// derived from all input parameters, so identical inputs and parameters always
/// produce identical output.
///
/// # Example
///
/// ```bash
/// fqtk subsample \
///     --input r1.fq.gz r2.fq.gz \
///     --output subsampled \
///     --fraction 0.1
/// ```
#[derive(Parser, Debug)]
#[command(version)]
#[clap(verbatim_doc_comment)]
pub(crate) struct Subsample {
    /// One or more input FASTQ files (may be gzipped). All files must have the
    /// same number of reads in the same order.
    #[clap(long, short = 'i', required = true, num_args = 1..)]
    inputs: Vec<PathBuf>,

    /// Output path prefix. Files will be named {output}.R1.fq.gz, etc.
    #[clap(long, short = 'o', required = true)]
    output: PathBuf,

    /// Fraction of reads to retain, in the range [0.0, 1.0].
    #[clap(long, short = 'f', required = true)]
    fraction: f64,

    /// Number of threads for compression. Minimum 2.
    #[clap(long, short = 't', default_value = "8")]
    threads: usize,

    /// BGZF compression level for output files.
    #[clap(long, short = 'c', default_value = "5")]
    compression_level: usize,

    /// Explicit RNG seed for reproducibility. When omitted, a deterministic
    /// seed is derived from all other parameters.
    #[clap(long, short = 's')]
    seed: Option<u64>,

    /// Disable checking that read names are in sync across input files.
    #[clap(long)]
    disable_read_name_checking: bool,
}

impl Hash for Subsample {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.inputs.hash(state);
        self.output.hash(state);
        self.fraction.to_bits().hash(state);
        self.threads.hash(state);
        self.compression_level.hash(state);
        self.seed.hash(state);
        self.disable_read_name_checking.hash(state);
    }
}

/// Returns the read name portion of a FASTQ header, stripping any trailing
/// `/1` or `/2` paired-end suffix and any comment after whitespace.
fn base_read_name(head: &[u8]) -> &[u8] {
    let name_end = head.iter().position(|&b| b == b' ' || b == b'\t').unwrap_or(head.len());
    let name = &head[..name_end];
    if name.len() >= 2
        && name[name.len() - 2] == b'/'
        && matches!(name[name.len() - 1], b'1' | b'2')
    {
        &name[..name.len() - 2]
    } else {
        name
    }
}

impl Subsample {
    /// Returns the effective RNG seed — either the explicit seed or one derived
    /// by hashing the struct.
    fn effective_seed(&self) -> u64 {
        if let Some(seed) = self.seed {
            return seed;
        }
        let mut hasher = DefaultHasher::new();
        self.hash(&mut hasher);
        hasher.finish()
    }

    /// Validates all inputs and returns a collected error if anything is wrong.
    fn validate(&self) -> Result<()> {
        let mut errors = vec![];

        if self.inputs.is_empty() {
            errors.push("At least one input file is required.".to_owned());
        }

        for input in &self.inputs {
            if !input.exists() {
                errors.push(format!("Input file {input:?} does not exist."));
            }
        }

        if !(0.0..=1.0).contains(&self.fraction) {
            errors.push(format!("Fraction must be in [0.0, 1.0], got {}.", self.fraction));
        }

        if self.threads < 2 {
            errors.push(format!("Threads must be at least 2, got {}.", self.threads));
        }

        if let Some(parent) = self.output.parent() {
            if !parent.as_os_str().is_empty() && !parent.exists() {
                errors.push(format!("Output parent directory {parent:?} does not exist."));
            }
        }

        if errors.is_empty() {
            Ok(())
        } else {
            let details = errors.iter().fold(String::new(), |mut s, e| {
                s.push_str(&format!("    - {e}\n"));
                s
            });
            Err(anyhow::anyhow!("The following errors with the input(s) were detected:\n{details}"))
        }
    }

    /// Creates BGZF-compressed pooled output writers, one per input file.
    fn create_writers(&self) -> Result<(Pool, Vec<PooledWriter>)> {
        let writer_threads = self.threads - 1;
        let mut pool_builder = PoolBuilder::<_, BgzfCompressor>::new()
            .threads(writer_threads)
            .queue_size(writer_threads * 50)
            .compression_level(u8::try_from(self.compression_level)?)?;

        let mut writers = Vec::with_capacity(self.inputs.len());
        for idx in 0..self.inputs.len() {
            let path = PathBuf::from(format!("{}.R{}.fq.gz", self.output.display(), idx + 1));
            let file_writer = BufWriter::new(File::create(&path)?);
            writers.push(pool_builder.exchange(file_writer));
        }

        let pool = pool_builder.build()?;
        Ok((pool, writers))
    }
}

impl Command for Subsample {
    fn execute(&self) -> Result<()> {
        self.validate()?;

        let seed = self.effective_seed();
        info!("Using random seed: {}", seed);
        let mut rng = ChaCha8Rng::seed_from_u64(seed);

        // Open input FASTQ readers
        let fgio = Io::new(5, BUFFER_SIZE);
        let mut sources: Vec<FastqReader<Box<dyn BufRead + Send>>> = self
            .inputs
            .iter()
            .map(|p| {
                fgio.new_reader(p)
                    .map(|r| FastqReader::with_capacity(r, BUFFER_SIZE))
                    .unwrap_or_else(|e| panic!("Failed to open {p:?}: {e}"))
            })
            .collect();

        // Create output writers
        let (mut pool, mut writers) = self.create_writers()?;

        info!(
            "Subsampling {} input file(s) at fraction {:.4} to {:?}",
            self.inputs.len(),
            self.fraction,
            self.output,
        );

        let log_unit: u64 = 5_000_000;
        let num_inputs = sources.len();
        let check_names = !self.disable_read_name_checking && num_inputs > 1;
        let mut expected_name: Vec<u8> = Vec::new();
        let mut total_read: u64 = 0;
        let mut total_kept: u64 = 0;

        loop {
            let keep = rng.random::<f64>() < self.fraction;
            let mut records_found = 0usize;

            for (i, source) in sources.iter_mut().enumerate() {
                if let Some(result) = source.next() {
                    let rec = result?;
                    records_found += 1;

                    if keep {
                        if check_names {
                            let name = base_read_name(rec.head());
                            if i == 0 {
                                expected_name.clear();
                                expected_name.extend_from_slice(name);
                            } else if name != expected_name.as_slice() {
                                anyhow::bail!(
                                    "Read name mismatch at read {}: file 0={:?}, file {i}={:?}",
                                    total_read + 1,
                                    String::from_utf8_lossy(&expected_name),
                                    String::from_utf8_lossy(name),
                                );
                            }
                        }

                        rec.write_unchanged(&mut writers[i])?;
                    }
                }
            }

            if records_found == 0 {
                break;
            }

            ensure!(
                records_found == num_inputs,
                "FASTQ files are out of sync: {} of {} files had a record at read {}",
                records_found,
                num_inputs,
                total_read + 1,
            );

            total_read += 1;
            if keep {
                total_kept += 1;
            }
            if total_read % log_unit == 0 {
                let pct = total_kept as f64 / total_read as f64 * 100.0;
                info!(
                    "[fqtk subsample] Read {} record sets and wrote {} ({:.1}%).",
                    fmt_count(total_read),
                    fmt_count(total_kept),
                    pct,
                );
            }
        }

        // Shut down writers and pool
        info!("Finished reading input FASTQs.");
        for writer in writers {
            writer.close()?;
        }
        pool.stop_pool()?;

        let pct = if total_read > 0 { total_kept as f64 / total_read as f64 * 100.0 } else { 0.0 };
        info!(
            "[fqtk subsample] Read {} record sets and wrote {} ({:.1}%).",
            fmt_count(total_read),
            fmt_count(total_kept),
            pct,
        );

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use seq_io::fastq::OwnedRecord;
    use tempfile::TempDir;

    /// Creates FASTQ content lines from the given read name prefix and base sequences.
    fn fq_lines_from_bases(prefix: &str, records_bases: &[&str]) -> Vec<String> {
        records_bases
            .iter()
            .enumerate()
            .flat_map(|(i, &bases)| {
                vec![
                    format!("@{}_{}", prefix, i),
                    bases.to_owned(),
                    "+".to_owned(),
                    ";".repeat(bases.len()),
                ]
            })
            .collect()
    }

    /// Creates FASTQ content lines with paired-end suffixes (/1 or /2) on read names.
    fn fq_lines_with_suffix(prefix: &str, suffix: &str, records_bases: &[&str]) -> Vec<String> {
        records_bases
            .iter()
            .enumerate()
            .flat_map(|(i, &bases)| {
                vec![
                    format!("@{}_{}{}", prefix, i, suffix),
                    bases.to_owned(),
                    "+".to_owned(),
                    ";".repeat(bases.len()),
                ]
            })
            .collect()
    }

    /// Writes a FASTQ file to tmpdir and returns its path.
    fn write_fastq(tmpdir: &TempDir, filename: &str, lines: &[String]) -> PathBuf {
        let io = Io::default();
        let path = tmpdir.path().join(format!("{filename}.fastq"));
        io.write_lines(&path, lines).unwrap();
        path
    }

    /// Reads all records from a (possibly gzipped) FASTQ file.
    fn read_fastq(path: &PathBuf) -> Vec<OwnedRecord> {
        let fg_io = Io::default();
        FastqReader::new(fg_io.new_reader(path).unwrap())
            .into_records()
            .collect::<Result<Vec<_>, seq_io::fastq::Error>>()
            .unwrap()
    }

    fn make_bases(n: usize) -> Vec<&'static str> {
        vec!["ACGTACGTACGT"; n]
    }

    // ---- base_read_name tests ----

    #[test]
    fn test_base_read_name_plain() {
        assert_eq!(base_read_name(b"readname"), b"readname");
    }

    #[test]
    fn test_base_read_name_with_slash_1() {
        assert_eq!(base_read_name(b"readname/1"), b"readname");
    }

    #[test]
    fn test_base_read_name_with_slash_2() {
        assert_eq!(base_read_name(b"readname/2"), b"readname");
    }

    #[test]
    fn test_base_read_name_with_comment() {
        assert_eq!(base_read_name(b"readname some comment"), b"readname");
    }

    #[test]
    fn test_base_read_name_with_suffix_and_comment() {
        assert_eq!(base_read_name(b"readname/1 some comment"), b"readname");
    }

    #[test]
    fn test_base_read_name_with_tab_comment() {
        assert_eq!(base_read_name(b"readname\tcomment"), b"readname");
    }

    #[test]
    fn test_base_read_name_slash_3_not_stripped() {
        assert_eq!(base_read_name(b"readname/3"), b"readname/3");
    }

    // ---- validation tests ----

    #[test]
    fn test_validation_missing_input() {
        let tmp = TempDir::new().unwrap();
        let cmd = Subsample {
            inputs: vec![tmp.path().join("nonexistent.fq")],
            output: tmp.path().join("out"),
            fraction: 0.5,
            threads: 2,
            compression_level: 1,
            seed: None,
            disable_read_name_checking: false,
        };
        let err = cmd.validate().unwrap_err();
        assert!(err.to_string().contains("does not exist"), "{err}");
    }

    #[test]
    fn test_validation_bad_fraction_negative() {
        let tmp = TempDir::new().unwrap();
        let lines = fq_lines_from_bases("r", &["ACGT"]);
        let input = write_fastq(&tmp, "r1", &lines);
        let cmd = Subsample {
            inputs: vec![input],
            output: tmp.path().join("out"),
            fraction: -0.1,
            threads: 2,
            compression_level: 1,
            seed: None,
            disable_read_name_checking: false,
        };
        let err = cmd.validate().unwrap_err();
        assert!(err.to_string().contains("Fraction must be"), "{err}");
    }

    #[test]
    fn test_validation_bad_fraction_over_one() {
        let tmp = TempDir::new().unwrap();
        let lines = fq_lines_from_bases("r", &["ACGT"]);
        let input = write_fastq(&tmp, "r1", &lines);
        let cmd = Subsample {
            inputs: vec![input],
            output: tmp.path().join("out"),
            fraction: 1.5,
            threads: 2,
            compression_level: 1,
            seed: None,
            disable_read_name_checking: false,
        };
        let err = cmd.validate().unwrap_err();
        assert!(err.to_string().contains("Fraction must be"), "{err}");
    }

    #[test]
    fn test_validation_insufficient_threads() {
        let tmp = TempDir::new().unwrap();
        let lines = fq_lines_from_bases("r", &["ACGT"]);
        let input = write_fastq(&tmp, "r1", &lines);
        let cmd = Subsample {
            inputs: vec![input],
            output: tmp.path().join("out"),
            fraction: 0.5,
            threads: 1,
            compression_level: 1,
            seed: None,
            disable_read_name_checking: false,
        };
        let err = cmd.validate().unwrap_err();
        assert!(err.to_string().contains("Threads must be at least 2"), "{err}");
    }

    // ---- functional tests ----

    #[test]
    fn test_single_end() {
        let tmp = TempDir::new().unwrap();
        let bases = make_bases(200);
        let lines = fq_lines_from_bases("r", &bases);
        let input = write_fastq(&tmp, "reads", &lines);
        let output = tmp.path().join("out");

        let cmd = Subsample {
            inputs: vec![input],
            output: output.clone(),
            fraction: 0.5,
            threads: 2,
            compression_level: 1,
            seed: Some(42),
            disable_read_name_checking: false,
        };
        cmd.execute().unwrap();

        let out_path = PathBuf::from(format!("{}.R1.fq.gz", output.display()));
        assert!(out_path.exists());
        let records = read_fastq(&out_path);
        // With 200 reads at 0.5 fraction, expect roughly 80-120
        assert!(records.len() > 50, "Too few records: {}", records.len());
        assert!(records.len() < 150, "Too many records: {}", records.len());
    }

    #[test]
    fn test_paired_end() {
        let tmp = TempDir::new().unwrap();
        let bases = make_bases(100);
        let lines_r1 = fq_lines_with_suffix("r", "/1", &bases);
        let lines_r2 = fq_lines_with_suffix("r", "/2", &bases);
        let input_r1 = write_fastq(&tmp, "r1", &lines_r1);
        let input_r2 = write_fastq(&tmp, "r2", &lines_r2);
        let output = tmp.path().join("out");

        let cmd = Subsample {
            inputs: vec![input_r1, input_r2],
            output: output.clone(),
            fraction: 0.3,
            threads: 2,
            compression_level: 1,
            seed: Some(42),
            disable_read_name_checking: false,
        };
        cmd.execute().unwrap();

        let out_r1 = PathBuf::from(format!("{}.R1.fq.gz", output.display()));
        let out_r2 = PathBuf::from(format!("{}.R2.fq.gz", output.display()));
        let records_r1 = read_fastq(&out_r1);
        let records_r2 = read_fastq(&out_r2);

        assert_eq!(records_r1.len(), records_r2.len(), "R1 and R2 counts differ");
        assert!(!records_r1.is_empty(), "No records kept");

        // Verify read names are in sync (same base name for each pair)
        for (r1, r2) in records_r1.iter().zip(records_r2.iter()) {
            let name1 = base_read_name(&r1.head);
            let name2 = base_read_name(&r2.head);
            assert_eq!(name1, name2, "Read names out of sync");
        }
    }

    #[test]
    fn test_three_inputs() {
        let tmp = TempDir::new().unwrap();
        let bases = make_bases(100);
        let lines1 = fq_lines_with_suffix("r", "/1", &bases);
        let lines2 = fq_lines_with_suffix("r", "/2", &bases);
        let lines3 = fq_lines_from_bases("r", &bases);
        let input1 = write_fastq(&tmp, "r1", &lines1);
        let input2 = write_fastq(&tmp, "r2", &lines2);
        let input3 = write_fastq(&tmp, "r3", &lines3);
        let output = tmp.path().join("out");

        let cmd = Subsample {
            inputs: vec![input1, input2, input3],
            output: output.clone(),
            fraction: 0.5,
            threads: 2,
            compression_level: 1,
            seed: Some(42),
            disable_read_name_checking: false,
        };
        cmd.execute().unwrap();

        let counts: Vec<usize> = (1..=3)
            .map(|i| {
                let path = PathBuf::from(format!("{}.R{i}.fq.gz", output.display()));
                assert!(path.exists(), "Missing output file R{i}");
                read_fastq(&path).len()
            })
            .collect();

        assert!(counts[0] > 0, "No records kept");
        assert_eq!(counts[0], counts[1], "R1 and R2 counts differ");
        assert_eq!(counts[0], counts[2], "R1 and R3 counts differ");
    }

    #[test]
    fn test_seed_reproducibility() {
        let tmp = TempDir::new().unwrap();
        let bases = make_bases(100);
        let lines = fq_lines_from_bases("r", &bases);
        let input = write_fastq(&tmp, "reads", &lines);

        let output1 = tmp.path().join("out1");
        let output2 = tmp.path().join("out2");

        for output in [&output1, &output2] {
            let cmd = Subsample {
                inputs: vec![input.clone()],
                output: output.clone(),
                fraction: 0.3,
                threads: 2,
                compression_level: 1,
                seed: Some(42),
                disable_read_name_checking: false,
            };
            cmd.execute().unwrap();
        }

        let records1 = read_fastq(&PathBuf::from(format!("{}.R1.fq.gz", output1.display())));
        let records2 = read_fastq(&PathBuf::from(format!("{}.R1.fq.gz", output2.display())));

        assert_eq!(records1.len(), records2.len());
        for (r1, r2) in records1.iter().zip(records2.iter()) {
            assert_eq!(r1, r2);
        }
    }

    #[test]
    fn test_deterministic_without_seed() {
        let tmp = TempDir::new().unwrap();
        let bases = make_bases(100);
        let lines = fq_lines_from_bases("r", &bases);
        let input = write_fastq(&tmp, "reads", &lines);
        let output = tmp.path().join("out");

        let cmd = Subsample {
            inputs: vec![input.clone()],
            output: output.clone(),
            fraction: 0.3,
            threads: 2,
            compression_level: 1,
            seed: None,
            disable_read_name_checking: false,
        };

        // Run once, read results, then run again and compare
        cmd.execute().unwrap();
        let out_path = PathBuf::from(format!("{}.R1.fq.gz", output.display()));
        let records1 = read_fastq(&out_path);

        cmd.execute().unwrap();
        let records2 = read_fastq(&out_path);

        assert_eq!(records1.len(), records2.len());
        for (r1, r2) in records1.iter().zip(records2.iter()) {
            assert_eq!(r1, r2);
        }
    }

    #[test]
    fn test_different_params_different_output() {
        let tmp = TempDir::new().unwrap();
        let bases = make_bases(200);
        let lines = fq_lines_from_bases("r", &bases);
        let input = write_fastq(&tmp, "reads", &lines);

        let output1 = tmp.path().join("out1");
        let output2 = tmp.path().join("out2");

        let cmd1 = Subsample {
            inputs: vec![input.clone()],
            output: output1.clone(),
            fraction: 0.5,
            threads: 2,
            compression_level: 1,
            seed: None,
            disable_read_name_checking: false,
        };
        cmd1.execute().unwrap();

        let cmd2 = Subsample {
            inputs: vec![input.clone()],
            output: output2.clone(),
            fraction: 0.51,
            threads: 2,
            compression_level: 1,
            seed: None,
            disable_read_name_checking: false,
        };
        cmd2.execute().unwrap();

        let records1 = read_fastq(&PathBuf::from(format!("{}.R1.fq.gz", output1.display())));
        let records2 = read_fastq(&PathBuf::from(format!("{}.R1.fq.gz", output2.display())));

        // Different fractions should produce different derived seeds and thus different subsets
        assert_ne!(
            records1.iter().map(|r| r.head.clone()).collect::<Vec<_>>(),
            records2.iter().map(|r| r.head.clone()).collect::<Vec<_>>(),
        );
    }

    #[test]
    fn test_fraction_zero() {
        let tmp = TempDir::new().unwrap();
        let bases = make_bases(50);
        let lines = fq_lines_from_bases("r", &bases);
        let input = write_fastq(&tmp, "reads", &lines);
        let output = tmp.path().join("out");

        let cmd = Subsample {
            inputs: vec![input],
            output: output.clone(),
            fraction: 0.0,
            threads: 2,
            compression_level: 1,
            seed: Some(42),
            disable_read_name_checking: false,
        };
        cmd.execute().unwrap();

        let records = read_fastq(&PathBuf::from(format!("{}.R1.fq.gz", output.display())));
        assert_eq!(records.len(), 0, "Expected no records with fraction=0.0");
    }

    #[test]
    fn test_fraction_one() {
        let tmp = TempDir::new().unwrap();
        let bases = make_bases(50);
        let lines = fq_lines_from_bases("r", &bases);
        let input = write_fastq(&tmp, "reads", &lines);
        let output = tmp.path().join("out");

        let cmd = Subsample {
            inputs: vec![input],
            output: output.clone(),
            fraction: 1.0,
            threads: 2,
            compression_level: 1,
            seed: Some(42),
            disable_read_name_checking: false,
        };
        cmd.execute().unwrap();

        let records = read_fastq(&PathBuf::from(format!("{}.R1.fq.gz", output.display())));
        assert_eq!(records.len(), 50, "Expected all 50 records with fraction=1.0");
    }

    #[test]
    fn test_empty_input() {
        let tmp = TempDir::new().unwrap();
        let lines: Vec<String> = vec![];
        let input = write_fastq(&tmp, "empty", &lines);
        let output = tmp.path().join("out");

        let cmd = Subsample {
            inputs: vec![input],
            output: output.clone(),
            fraction: 0.5,
            threads: 2,
            compression_level: 1,
            seed: Some(42),
            disable_read_name_checking: false,
        };
        cmd.execute().unwrap();

        let out_path = PathBuf::from(format!("{}.R1.fq.gz", output.display()));
        assert!(out_path.exists());
        let records = read_fastq(&out_path);
        assert_eq!(records.len(), 0);
    }

    #[test]
    fn test_read_name_check_passes() {
        let tmp = TempDir::new().unwrap();
        let bases = make_bases(20);
        let lines_r1 = fq_lines_with_suffix("read", "/1", &bases);
        let lines_r2 = fq_lines_with_suffix("read", "/2", &bases);
        let input_r1 = write_fastq(&tmp, "r1", &lines_r1);
        let input_r2 = write_fastq(&tmp, "r2", &lines_r2);
        let output = tmp.path().join("out");

        let cmd = Subsample {
            inputs: vec![input_r1, input_r2],
            output: output.clone(),
            fraction: 1.0,
            threads: 2,
            compression_level: 1,
            seed: Some(42),
            disable_read_name_checking: false,
        };
        // Should succeed — names match after stripping suffixes
        cmd.execute().unwrap();
    }

    #[test]
    fn test_read_name_check_fails() {
        let tmp = TempDir::new().unwrap();
        let bases_a: Vec<&str> = vec!["ACGT"; 5];
        let bases_b: Vec<&str> = vec!["ACGT"; 5];
        let lines_r1 = fq_lines_from_bases("readA", &bases_a);
        let lines_r2 = fq_lines_from_bases("readB", &bases_b);
        let input_r1 = write_fastq(&tmp, "r1", &lines_r1);
        let input_r2 = write_fastq(&tmp, "r2", &lines_r2);
        let output = tmp.path().join("out");

        let cmd = Subsample {
            inputs: vec![input_r1, input_r2],
            output: output.clone(),
            fraction: 1.0,
            threads: 2,
            compression_level: 1,
            seed: Some(42),
            disable_read_name_checking: false,
        };
        let err = cmd.execute().unwrap_err();
        assert!(err.to_string().contains("Read name mismatch"), "{err}");
    }

    #[test]
    fn test_read_name_check_disabled() {
        let tmp = TempDir::new().unwrap();
        let bases_a: Vec<&str> = vec!["ACGT"; 5];
        let bases_b: Vec<&str> = vec!["ACGT"; 5];
        let lines_r1 = fq_lines_from_bases("readA", &bases_a);
        let lines_r2 = fq_lines_from_bases("readB", &bases_b);
        let input_r1 = write_fastq(&tmp, "r1", &lines_r1);
        let input_r2 = write_fastq(&tmp, "r2", &lines_r2);
        let output = tmp.path().join("out");

        let cmd = Subsample {
            inputs: vec![input_r1, input_r2],
            output: output.clone(),
            fraction: 1.0,
            threads: 2,
            compression_level: 1,
            seed: Some(42),
            disable_read_name_checking: true,
        };
        // Should succeed with name checking disabled
        cmd.execute().unwrap();
    }
}
