use crate::commands::command::{Command, validation_result};
use crate::commands::fastq_readahead::ReadAheadBuilder;
use crate::commands::outputs::{
    check_output_paths, create_bgzf_writers, remove_incomplete_outputs,
};
use crate::commands::read_set_reader::ReadSetReader;
use anyhow::Result;
use clap::Parser;
use clap::builder::RangedU64ValueParser;
use log::info;
use pooled_writer::PooledWriter;
use proglog::{CountFormatterKind, ProgLogBuilder};
use std::path::PathBuf;

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
/// Inputs may be uncompressed, gzipped, or block-gzipped.  Output files are _always_ block gzipped.
///
/// The records read together from the inputs must share a read name, ignoring any `/1` or `/2`
/// suffix and comment, unless `--disable-read-name-checking` is given.
///
/// If sharding fails, output files are deleted so they can't be mistaken for complete ones.
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
    #[clap(long, short = 's', value_parser = RangedU64ValueParser::<usize>::new().range(1..))]
    shards: usize,

    /// The number of threads to use for compressing output files.  Minimum 2.
    #[clap(long, short = 't', default_value = "8", value_parser = RangedU64ValueParser::<usize>::new().range(2..))]
    threads: usize,

    /// The level of compression to use to compress outputs.  Defaults to 1 because sharded FASTQs
    /// are typically short-lived intermediates, where write throughput matters more than squeezing
    /// out the last few percent of file size.
    #[clap(long, short = 'c', default_value = "1",
           value_parser = RangedU64ValueParser::<u8>::new().range(1..=12))]
    compression_level: u8,

    /// Disable checking that read names are in sync across input files.
    #[clap(long)]
    disable_read_name_checking: bool,

    /// (hidden) Size, in bytes, of each decompressed input chunk passed from a per-input
    /// decompression thread to the main parsing thread.
    #[clap(long, hide = true, default_value = "131072",
           value_parser = RangedU64ValueParser::<usize>::new().range(1..))]
    chunk_size: usize,

    /// (hidden) Number of decompressed input chunks to keep in flight per input file.
    #[clap(long, hide = true, default_value = "32",
           value_parser = RangedU64ValueParser::<usize>::new().range(1..))]
    chunk_count: usize,
}

impl Shard {
    /// Checks the output paths before any input is read or output created.
    fn validate(&self) -> Result<()> {
        let mut errors = Vec::new();
        check_output_paths(&self.inputs, &self.output_paths(), &mut errors);
        validation_result(&errors)
    }

    /// Returns the output paths: one per input for the first shard, then for the second, etc.
    fn output_paths(&self) -> Vec<PathBuf> {
        let mut paths = Vec::with_capacity(self.shards * self.inputs.len());
        for shard_num in 1..=self.shards {
            for read_num in 1..=self.inputs.len() {
                paths.push(PathBuf::from(format!(
                    "{}.{}{shard_num}.{}{read_num}.fq.gz",
                    self.output_prefix, self.shard_prefix, self.read_number_prefix
                )));
            }
        }
        paths
    }
}

impl Command for Shard {
    fn execute(&self) -> Result<()> {
        self.validate()?;
        info!("Reading {} input FASTQs and generating {} shards.", self.inputs.len(), self.shards);

        // Create one reader per input.  Each input is opened and decompressed on its own
        // background read-ahead thread, and its reader parses the decompressed bytes here on the
        // main thread.  Keeping inflation off the main thread lets it spend its time parsing and
        // routing records to keep the compression pool fed.  Opening waits until every input is
        // readable, so an unreadable input fails before any shard file is created.
        let inputs = self
            .inputs
            .iter()
            .map(|path| ReadAheadBuilder {
                path: path.clone(),
                chunk_size: self.chunk_size,
                chunk_count: self.chunk_count,
            })
            .collect();
        let mut read_sets = ReadSetReader::open(inputs, !self.disable_read_name_checking)?;

        // All but this thread compress the output.
        let outputs = self.output_paths();
        let (mut pool, mut writers) =
            create_bgzf_writers(&outputs, self.threads - 1, self.compression_level)?;
        let n_inputs = self.inputs.len();

        let logger = ProgLogBuilder::new()
            .name("fqtk")
            .noun("record sets")
            .verb("read")
            .unit(5_000_000)
            .count_formatter(CountFormatterKind::Comma)
            .level(log::Level::Info)
            .build();

        // Loop, consuming one record from each input file and writing the set to the current
        // shard.  Records are emitted verbatim via `write_unchanged`, which blits the original
        // record bytes rather than re-serializing field by field.  Terminate cleanly when all
        // inputs are exhausted together, or with an error if some -- but not all -- run dry.
        //
        // The processing result is captured rather than returned directly so the writer pool is
        // always shut down afterwards (below), even on error -- flushing queued output and
        // surfacing any finalization error instead of leaving it to a panicking drop.
        let mut target_shard_idx: usize = 0;
        let process_result: Result<()> = 'process: {
            loop {
                // The whole read set is read and checked before any of it is written, so a read
                // error, an early end or a name mismatch can't leave a half-written, misaligned
                // set in the outputs.
                let records = match read_sets.next() {
                    Ok(Some(records)) => records,
                    Ok(None) => break, // all inputs ended together: a clean end of input
                    Err(e) => break 'process Err(e),
                };

                let target = &mut writers[target_shard_idx * n_inputs..][..n_inputs];
                for (record, writer) in records.iter().zip(target) {
                    if let Err(e) = record.write_unchanged(&mut *writer) {
                        break 'process Err(e.into());
                    }
                }

                target_shard_idx = (target_shard_idx + 1) % self.shards;
                logger.record();
            }
            Ok(())
        };

        // Always shut the pool down so queued output is flushed and any finalization error is
        // propagated, regardless of whether record processing succeeded.  A processing error
        // takes precedence over a cleanup error.
        info!("Finished reading input FASTQs.");
        let close_result =
            writers.into_iter().try_for_each(PooledWriter::close).map_err(anyhow::Error::from);
        let stop_result = pool.stop_pool().map_err(anyhow::Error::from);
        let result = process_result.and(close_result).and(stop_result);
        if result.is_err() {
            remove_incomplete_outputs(&outputs);
        }
        result?;
        info!("Output FASTQ writing complete.");

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use crate::commands::command::Command;
    use crate::commands::shard::Shard;
    use bstr::ByteSlice;
    use itertools::Itertools;
    use rand;
    use seq_io::fastq::{OwnedRecord, Record};
    use std::collections::HashSet;
    use std::fs::File;
    use std::io::{BufReader, BufWriter, Write};
    use std::path::{Path, PathBuf};
    use tempfile::TempDir;

    /// Writes zero or more records to the given `Write`.  Bases and quals are randomly
    /// generated.  Read names are of the format `{prefix}{idx}{suffix}` where `idx` starts
    /// from the value of the `idx` parameter, and increases by one for each record.
    fn write_fastq_records<W: Write>(
        out: &mut W,
        prefix: &str,
        suffix: &str,
        idx: usize,
        count: usize,
    ) {
        let bases = "ACGT".as_bytes();
        for i in idx..idx + count {
            let seq = (0..30).map(|_| rand::random_range(0..4)).map(|j| bases[j]).collect_vec();
            let qual = (0..30).map(|_| rand::random_range(2u8..40u8) + 33).collect_vec();
            let rec = OwnedRecord {
                head: format!("{}{}{}", prefix, i, suffix).as_bytes().to_owned(),
                seq,
                qual,
            };
            rec.write(&mut *out).unwrap();
        }
    }

    /// Writes records to a FASTQ file at `path`.  Compression is selected by file extension
    /// (`.gz` produces gzip via fgoxide; otherwise plain text).
    fn build_fastq(path: &Path, prefix: &str, suffix: &str, idx: usize, count: usize) {
        let io = fgoxide::io::Io::new(1, 8 * 1024);
        let mut out = io.new_writer(path).unwrap();
        write_fastq_records(&mut out, prefix, suffix, idx, count);
    }

    /// Writes records to a BGZF-compressed FASTQ file at `path`.
    fn build_fastq_bgzf(path: &Path, prefix: &str, suffix: &str, idx: usize, count: usize) {
        let file = BufWriter::new(File::create(path).unwrap());
        let mut writer = bgzf::Writer::new(file, bgzf::CompressionLevel::new(3).unwrap());
        write_fastq_records(&mut writer, prefix, suffix, idx, count);
        writer.finish().unwrap().flush().unwrap();
    }

    /// Runs sharding and returns the outputs.  The returned Vec is nested as follows:
    ///  - Each entry in the top level vec is a _shard_
    ///  - Each entry in the second level vec is the sharded reads from a single input fastq
    ///  - Each entry in the third level vec is an individual fastq record
    fn run_sharding(tmp: &TempDir, inputs: &[&Path], shards: usize) -> Vec<Vec<Vec<OwnedRecord>>> {
        let prefix = format!("{}/test_out", tmp.path().to_str().unwrap());
        build_sharder(inputs, &prefix, shards).execute().unwrap();
        collect_shard_outputs(&prefix, inputs.len(), shards, read_fastq)
    }

    /// Builds a `Shard` instance for testing with the given inputs, output prefix, and shard count.
    fn build_sharder(inputs: &[&Path], output_prefix: &str, shards: usize) -> Shard {
        Shard {
            inputs: inputs.iter().map(|p| p.to_path_buf()).collect_vec(),
            output_prefix: output_prefix.to_string(),
            shard_prefix: "shard".to_string(),
            read_number_prefix: "read".to_string(),
            shards,
            threads: 4,
            compression_level: 1,
            disable_read_name_checking: false,
            // Deliberately small so most records straddle chunk boundaries, exercising the
            // ChunkReader reassembly path across the test suite.
            chunk_size: 100,
            chunk_count: 8,
        }
    }

    /// Reads back the sharded outputs produced under `output_prefix`, using `reader` to parse each
    /// individual output file.  Returns the same nested layout as `run_sharding`.
    fn collect_shard_outputs<F>(
        output_prefix: &str,
        n_inputs: usize,
        shards: usize,
        reader: F,
    ) -> Vec<Vec<Vec<OwnedRecord>>>
    where
        F: Fn(&Path) -> Vec<OwnedRecord>,
    {
        let mut results: Vec<Vec<Vec<OwnedRecord>>> = Vec::with_capacity(shards);
        for shard in 1..=shards {
            let mut reads_vecs = Vec::with_capacity(n_inputs);
            for input_idx in 1..=n_inputs {
                let path_str = format!("{}.shard{}.read{}.fq.gz", output_prefix, shard, input_idx);
                reads_vecs.push(reader(Path::new(&path_str)));
            }
            results.push(reads_vecs);
        }
        results
    }

    /// Reads a FASTQ file into a vec of records, transparently handling gzip via fgoxide.
    fn read_fastq(path: &Path) -> Vec<OwnedRecord> {
        let io = fgoxide::io::Io::new(1, 8 * 1024);
        let mut reader = io.new_reader(path).unwrap();
        let mut fq_reader = seq_io::fastq::Reader::with_capacity(&mut reader, 8 * 1024);
        fq_reader.records().map(|r| r.unwrap()).collect_vec()
    }

    /// Reads a FASTQ file into a vec of records, decompressing with the bgzf crate's `Reader`.
    fn read_fastq_via_bgzf(path: &Path) -> Vec<OwnedRecord> {
        let file = BufReader::new(File::open(path).unwrap());
        let bgzf_reader = bgzf::Reader::new(file);
        let mut fq_reader = seq_io::fastq::Reader::with_capacity(bgzf_reader, 8 * 1024);
        fq_reader.records().map(|r| r.unwrap()).collect_vec()
    }

    /// Returns the integer index embedded in a read name by `build_fastq`.  Strips any
    /// trailing `/N` read-end marker before parsing the trailing digits.
    fn read_index(rec: &OwnedRecord) -> usize {
        let head = rec.head.to_str().unwrap();
        let trimmed = head
            .rsplit_once('/')
            .filter(|(_, tail)| tail.chars().all(|c| c.is_ascii_digit()))
            .map(|(head, _)| head)
            .unwrap_or(head);
        let digits: String = trimmed.chars().rev().take_while(|c| c.is_ascii_digit()).collect();
        digits.chars().rev().collect::<String>().parse().unwrap()
    }

    #[test]
    fn test_shard_single_file() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        build_fastq(r1.as_path(), "q", "", 1, 50);
        let outputs = run_sharding(&tmp, &[&r1], 5);

        assert_eq!(outputs.len(), 5);
        for shard in outputs.iter() {
            assert_eq!(shard.len(), 1);
            assert_eq!(shard.iter().next().unwrap().len(), 10);
        }

        let read_names: HashSet<&str> =
            outputs.iter().flatten().flatten().map(|r| r.head.to_str().unwrap()).collect();
        assert_eq!(read_names.len(), 50);
    }

    #[test]
    fn test_shard_multiple_files() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        let r2 = PathBuf::from(tmp.path()).join("r2.fq");
        build_fastq(r1.as_path(), "q", "/1", 1, 64);
        build_fastq(r2.as_path(), "q", "/2", 1, 64);
        let outputs = run_sharding(&tmp, &[&r1, &r2], 3);

        assert_eq!(outputs.len(), 3);

        for shard in outputs.iter() {
            assert_eq!(shard.len(), 2); // two reads in each shard
            assert_eq!(shard[0].len(), shard[1].len()); // both r1 and r2 have same number of reads
            assert!(shard[0].len() == 21 || shard[0].len() == 22); // 64 / 3 == 21 1/3
        }
    }

    /// Verifies that the Nth input record (1-based) lands in shard ((N-1) % shards) + 1, and
    /// that within a shard the records appear in input order.
    #[test]
    fn test_round_robin_assignment() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        let n_reads = 20;
        let n_shards = 4;
        build_fastq(r1.as_path(), "q", "", 1, n_reads);
        let outputs = run_sharding(&tmp, &[&r1], n_shards);

        for (shard_idx, shard) in outputs.iter().enumerate() {
            let expected: Vec<usize> =
                (1..=n_reads).filter(|i| (i - 1) % n_shards == shard_idx).collect();
            let actual: Vec<usize> = shard[0].iter().map(read_index).collect();
            assert_eq!(actual, expected, "shard {} contents wrong", shard_idx + 1);
        }
    }

    /// Verifies that for paired input, the j-th r1 record in any shard came from the same input
    /// record-set as the j-th r2 record (i.e., they share an index).
    #[test]
    fn test_paired_records_stay_aligned() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        let r2 = PathBuf::from(tmp.path()).join("r2.fq");
        build_fastq(r1.as_path(), "q", "/1", 1, 30);
        build_fastq(r2.as_path(), "q", "/2", 1, 30);
        let outputs = run_sharding(&tmp, &[&r1, &r2], 4);

        for shard in outputs.iter() {
            let r1_indices: Vec<usize> = shard[0].iter().map(read_index).collect();
            let r2_indices: Vec<usize> = shard[1].iter().map(read_index).collect();
            assert_eq!(r1_indices, r2_indices);
        }
    }

    /// Every input read should appear in exactly one output shard, with no losses or duplicates.
    #[test]
    fn test_no_reads_lost_or_duplicated() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        build_fastq(r1.as_path(), "q", "", 1, 100);
        let outputs = run_sharding(&tmp, &[&r1], 7);

        let all_indices: Vec<usize> =
            outputs.iter().flatten().flatten().map(read_index).sorted().collect();
        let expected: Vec<usize> = (1..=100).collect();
        assert_eq!(all_indices, expected);
    }

    #[test]
    fn test_single_shard() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        build_fastq(r1.as_path(), "q", "", 1, 12);
        let outputs = run_sharding(&tmp, &[&r1], 1);

        assert_eq!(outputs.len(), 1);
        let actual: Vec<usize> = outputs[0][0].iter().map(read_index).collect();
        assert_eq!(actual, (1..=12).collect::<Vec<_>>());
    }

    #[test]
    fn test_more_shards_than_reads() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        build_fastq(r1.as_path(), "q", "", 1, 3);
        let outputs = run_sharding(&tmp, &[&r1], 5);

        assert_eq!(outputs.len(), 5);
        assert_eq!(outputs[0][0].len(), 1);
        assert_eq!(outputs[1][0].len(), 1);
        assert_eq!(outputs[2][0].len(), 1);
        assert_eq!(outputs[3][0].len(), 0);
        assert_eq!(outputs[4][0].len(), 0);
    }

    #[test]
    fn test_empty_input() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        build_fastq(r1.as_path(), "q", "", 1, 0);
        let outputs = run_sharding(&tmp, &[&r1], 4);

        assert_eq!(outputs.len(), 4);
        for shard in outputs.iter() {
            assert_eq!(shard.len(), 1);
            assert!(shard[0].is_empty());
        }
    }

    #[test]
    fn test_mismatched_input_lengths_fails() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        let r2 = PathBuf::from(tmp.path()).join("r2.fq");
        build_fastq(r1.as_path(), "q", "/1", 1, 10);
        build_fastq(r2.as_path(), "q", "/2", 1, 8);

        let prefix = format!("{}/test_out", tmp.path().to_str().unwrap());
        let err = build_sharder(&[&r1, &r2], &prefix, 3).execute().unwrap_err();
        assert!(err.to_string().contains("out of sync"), "unexpected error message: {}", err);
    }

    #[test]
    fn test_out_of_sync_inputs_fail_and_remove_the_outputs() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        let r2 = PathBuf::from(tmp.path()).join("r2.fq");
        build_fastq(r1.as_path(), "q", "/1", 1, 3);
        build_fastq(r2.as_path(), "q", "/2", 1, 2);

        let prefix = format!("{}/test_out", tmp.path().to_str().unwrap());
        let sharder = build_sharder(&[&r1, &r2], &prefix, 2);
        let err = sharder.execute().unwrap_err();
        assert!(err.to_string().contains("out of sync"), "unexpected error message: {}", err);
        for output in sharder.output_paths() {
            assert!(!output.exists(), "{output:?} was left behind");
        }
    }

    #[test]
    fn test_an_output_that_is_an_input_is_rejected_before_it_is_truncated() {
        let tmp = TempDir::new().unwrap();
        let prefix = format!("{}/test_out", tmp.path().to_str().unwrap());
        let input = PathBuf::from(format!("{prefix}.shard1.read1.fq.gz"));
        build_fastq(input.as_path(), "q", "", 1, 4);
        let before = std::fs::read(&input).unwrap();

        let err = build_sharder(&[&input], &prefix, 2).execute().unwrap_err();
        assert!(err.to_string().contains("is also an input"), "unexpected error: {err}");
        assert_eq!(std::fs::read(&input).unwrap(), before);
    }

    #[test]
    fn test_a_missing_output_directory_fails_before_any_input_is_read() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        build_fastq(r1.as_path(), "q", "", 1, 4);

        let prefix = format!("{}/missing/test_out", tmp.path().to_str().unwrap());
        let err = build_sharder(&[&r1], &prefix, 2).execute().unwrap_err();
        assert!(err.to_string().contains("does not exist"), "unexpected error: {err}");
    }

    #[test]
    fn test_mismatched_read_names_fail() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        let r2 = PathBuf::from(tmp.path()).join("r2.fq");
        build_fastq(r1.as_path(), "a", "/1", 1, 5);
        build_fastq(r2.as_path(), "b", "/2", 1, 5);

        let prefix = format!("{}/test_out", tmp.path().to_str().unwrap());
        let err = build_sharder(&[&r1, &r2], &prefix, 2).execute().unwrap_err();
        assert!(err.to_string().contains("Read name mismatch"), "unexpected error: {err}");
    }

    #[test]
    fn test_mismatched_read_names_are_sharded_with_checking_disabled() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        let r2 = PathBuf::from(tmp.path()).join("r2.fq");
        build_fastq(r1.as_path(), "a", "/1", 1, 5);
        build_fastq(r2.as_path(), "b", "/2", 1, 5);

        let prefix = format!("{}/test_out", tmp.path().to_str().unwrap());
        let sharder =
            Shard { disable_read_name_checking: true, ..build_sharder(&[&r1, &r2], &prefix, 2) };
        sharder.execute().unwrap();
        let outputs = collect_shard_outputs(&prefix, 2, 2, read_fastq);
        assert_eq!(outputs.iter().flatten().flatten().count(), 10);
    }

    /// A corrupt or truncated gzip input must surface as an error rather than hanging or silently
    /// producing empty output -- this exercises the ByteChunker -> ChunkReader error propagation.
    #[test]
    fn test_decompression_error_surfaces() {
        let tmp = TempDir::new().unwrap();
        let bad = PathBuf::from(tmp.path()).join("bad.fq.gz");
        std::fs::write(&bad, b"this is definitely not valid gzip content\n").unwrap();

        let prefix = format!("{}/test_out", tmp.path().to_str().unwrap());
        let result = build_sharder(&[bad.as_path()], &prefix, 2).execute();
        assert!(result.is_err(), "expected a decompression error but the run succeeded");
    }

    #[test]
    fn test_missing_input_fails_before_any_shard_file_is_created() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        let missing = PathBuf::from(tmp.path()).join("r2.fq");
        build_fastq(r1.as_path(), "q", "/1", 1, 10);

        let prefix = format!("{}/test_out", tmp.path().to_str().unwrap());
        let err = build_sharder(&[&r1, &missing], &prefix, 2).execute().unwrap_err();
        assert!(err.to_string().contains("r2.fq"), "error does not name the input: {err}");
        let files: Vec<String> = std::fs::read_dir(tmp.path())
            .unwrap()
            .map(|entry| entry.unwrap().file_name().into_string().unwrap())
            .collect();
        assert_eq!(files, ["r1.fq"]);
    }

    #[test]
    fn test_gzip_compressed_inputs() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq.gz");
        build_fastq(r1.as_path(), "q", "", 1, 25);
        let outputs = run_sharding(&tmp, &[&r1], 5);

        let all_indices: Vec<usize> =
            outputs.iter().flatten().flatten().map(read_index).sorted().collect();
        assert_eq!(all_indices, (1..=25).collect::<Vec<_>>());
    }

    /// Verifies that BGZF inputs spanning more than one block are read in full (i.e., we don't
    /// stop after the first BGZF block).  `BGZF_BLOCK_SIZE` is ~64 KiB; we generate enough
    /// records to comfortably exceed that uncompressed.
    #[test]
    fn test_bgzf_compressed_inputs_multiple_blocks() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq.gz");
        // Each record is ~70-75 bytes; 1500 records ≈ 105 KiB uncompressed → multiple BGZF blocks.
        let n_reads = 1500;
        build_fastq_bgzf(r1.as_path(), "q", "", 1, n_reads);
        let outputs = run_sharding(&tmp, &[&r1], 3);

        let all_indices: Vec<usize> =
            outputs.iter().flatten().flatten().map(read_index).sorted().collect();
        assert_eq!(all_indices, (1..=n_reads).collect::<Vec<_>>());
    }

    /// Outputs are documented to always be BGZF — read one back with the bgzf crate's `Reader`
    /// and check that the round-tripped contents match what we sent in.
    #[test]
    fn test_output_files_are_bgzf() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        build_fastq(r1.as_path(), "q", "", 1, 30);

        let prefix = format!("{}/test_out", tmp.path().to_str().unwrap());
        build_sharder(&[&r1], &prefix, 3).execute().unwrap();

        let outputs = collect_shard_outputs(&prefix, 1, 3, read_fastq_via_bgzf);
        let all_indices: Vec<usize> =
            outputs.iter().flatten().flatten().map(read_index).sorted().collect();
        assert_eq!(all_indices, (1..=30).collect::<Vec<_>>());
    }

    #[test]
    fn test_custom_prefixes() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        build_fastq(r1.as_path(), "q", "", 1, 6);

        let prefix = format!("{}/sample", tmp.path().to_str().unwrap());
        let sharder = Shard {
            inputs: vec![r1.clone()],
            output_prefix: prefix.clone(),
            shard_prefix: "chunk_".to_string(),
            read_number_prefix: "R".to_string(),
            shards: 2,
            threads: 2,
            compression_level: 1,
            disable_read_name_checking: false,
            chunk_size: 64 * 1024,
            chunk_count: 4,
        };
        sharder.execute().unwrap();

        for shard in 1..=2 {
            let path = PathBuf::from(format!("{}.chunk_{}.R1.fq.gz", prefix, shard));
            assert!(path.exists(), "expected output file does not exist: {}", path.display());
            assert_eq!(read_fastq(&path).len(), 3);
        }
    }

    /// With a pathologically small chunk size, individual records (and the boundaries between
    /// them) span many decompressed chunks, so the `ChunkReader` must stitch bytes back together
    /// before the parser sees them.  Every read should still round-trip exactly once.
    #[test]
    fn test_tiny_chunk_size_reassembles_records() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        build_fastq(r1.as_path(), "q", "", 1, 40);

        let prefix = format!("{}/test_out", tmp.path().to_str().unwrap());
        let sharder = Shard {
            inputs: vec![r1.clone()],
            output_prefix: prefix.clone(),
            shard_prefix: "shard".to_string(),
            read_number_prefix: "read".to_string(),
            shards: 3,
            threads: 4,
            compression_level: 1,
            disable_read_name_checking: false,
            chunk_size: 3,
            chunk_count: 2,
        };
        sharder.execute().unwrap();

        let outputs = collect_shard_outputs(&prefix, 1, 3, read_fastq);
        let all_indices: Vec<usize> =
            outputs.iter().flatten().flatten().map(read_index).sorted().collect();
        assert_eq!(all_indices, (1..=40).collect::<Vec<_>>());
    }

    /// With a chunk size larger than the whole input, every record arrives within a single chunk
    /// and no cross-chunk reassembly is needed -- the opposite extreme from the tiny-chunk case.
    #[test]
    fn test_large_chunk_size_single_chunk() {
        let tmp = TempDir::new().unwrap();
        let r1 = PathBuf::from(tmp.path()).join("r1.fq");
        build_fastq(r1.as_path(), "q", "", 1, 20);

        let prefix = format!("{}/test_out", tmp.path().to_str().unwrap());
        let sharder = Shard {
            inputs: vec![r1.clone()],
            output_prefix: prefix.clone(),
            shard_prefix: "shard".to_string(),
            read_number_prefix: "read".to_string(),
            shards: 2,
            threads: 2,
            compression_level: 1,
            disable_read_name_checking: false,
            chunk_size: 8 * 1024 * 1024,
            chunk_count: 4,
        };
        sharder.execute().unwrap();

        let outputs = collect_shard_outputs(&prefix, 1, 2, read_fastq);
        let all_indices: Vec<usize> =
            outputs.iter().flatten().flatten().map(read_index).sorted().collect();
        assert_eq!(all_indices, (1..=20).collect::<Vec<_>>());
    }
}
