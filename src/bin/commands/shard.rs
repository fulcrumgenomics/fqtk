use crate::commands::command::Command;
use anyhow::{Result, anyhow};
use clap::Parser;
use clap::builder::RangedU64ValueParser;
use fgoxide::io::Io;
use fgoxide::iter::IntoChunkedReadAheadIterator;
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

/// Buffer size used when opening BufReads for input files
const BUFFER_SIZE: usize = 64 * 1024;

/// Struct to hold all the output writers for a single shard
struct ShardWriters<W: Write> {
    /// ID of the shard this set of writers is for
    shard_number: usize,
    /// The set of writers for the shard
    writers: Vec<W>,
}

impl<W: Write> ShardWriters<W> {
    /// Destroys this struct and decomposes it into its component types. Used when swapping
    /// writers for pooled writers.
    fn into_parts(self) -> (usize, Vec<W>) {
        (self.shard_number, self.writers)
    }

    /// Writes a set of FASTQ records out.  The number of reads (i.e. `reads.len()`) must
    /// match the number of individual writers; if a mismatch is found an error will be raised.
    fn write(&mut self, reads: &[OwnedRecord]) -> Result<()> {
        if reads.len() != self.writers.len() {
            return Err(anyhow!("Expected {} reads, got {}", self.writers.len(), reads.len()));
        }

        for (writer, read) in self.writers.iter_mut().zip(reads.iter()) {
            read.write(writer)?;
        }

        Ok(())
    }
}

impl ShardWriters<PooledWriter> {
    /// Attempts to gracefully shut down the writers in this struct, consuming the struct in the
    /// process. Will error if closing of the `PooledWriter`s fails for any reason.
    fn close(self) -> Result<()> {
        self.writers.into_iter().map(|w| w.close()).collect::<Result<Vec<_>, io::Error>>()?;
        Ok(())
    }
}

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

    /// The level of compression to use to compress outputs.
    #[clap(long, short = 'c', default_value = "5")]
    compression_level: u8,
}

impl Shard {
    /// Opens a FASTQ reader per input path.  Handles uncompressed and gzipped FASTQ files.
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

    /// Builds the fastq writers for each output file and shard, and then instantiates
    /// a writer pool using block-gzip compression.  Returns the writer Pool itself along
    /// with a Vec of `ShardWriters` each of which can accept and write `Vec`s of fastq records
    /// with one record from each of the input files.
    fn build_writer_pool(&self) -> Result<(Pool, Vec<ShardWriters<PooledWriter>>)> {
        // First build up the per-shard writers
        let mut shard_writers = Vec::with_capacity(self.shards);
        for shard in 1..=self.shards {
            let mut ws = Vec::with_capacity(self.inputs.len());

            for source_idx in 1..=self.inputs.len() {
                let path_str = format!(
                    "{prefix}.{shard_prefix}{shard_num}.{read_prefix}{read_num}.fq.gz",
                    prefix = self.output_prefix,
                    shard_prefix = self.shard_prefix,
                    shard_num = shard,
                    read_prefix = self.read_number_prefix,
                    read_num = source_idx
                );
                let path = Path::new(&path_str);
                let writer = BufWriter::new(File::create(path)?);
                ws.push(writer);
            }

            shard_writers.push(ShardWriters { shard_number: shard, writers: ws });
        }

        // Then construct the writer pool, reserving one thread for the main loop.
        let pool_threads = self.threads - 1;
        let mut pool_builder = PoolBuilder::<_, BgzfCompressor>::new()
            .threads(pool_threads)
            .queue_size(pool_threads * 50)
            .compression_level(self.compression_level)?;

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
    fn execute(&self) -> Result<()> {
        info!("Reading {} input FASTQs and generating {} shards.", self.inputs.len(), self.shards);

        // Open the input FASTQ files and from each file generate a chunked read-ahead iterator
        // to ensure that the reading of the input files is not the bottleneck.  We don't count
        // the read-ahead threads as consuming from the thread-count as they will consume minimal
        // CPU and are mostly going to be doing/blocking on I/O.
        let fq_readers = Self::build_readers(&self.inputs)?;
        let fq_iters = fq_readers.into_iter().map(|r| r.into_records()).collect_vec();
        let mut fq_iters = fq_iters.into_iter().map(|i| i.read_ahead(128, 48)).collect_vec();

        let (mut pool, mut shard_writers) = self.build_writer_pool()?;

        let logger = ProgLogBuilder::new()
            .name("fqtk")
            .noun("record sets")
            .verb("read")
            .unit(5_000_000)
            .count_formatter(CountFormatterKind::Comma)
            .level(log::Level::Info)
            .build();

        // Loop, consuming one read from each input file, and writing to the appropriate
        // shard.  Terminate cleanly when all input iterators are exhausted, or with error
        // if at least one, but not all input iterators are exhausted.
        let mut target_shard_idx: usize = 0;
        let mut recs = Vec::with_capacity(fq_iters.len());
        loop {
            // Pull in the next set of reads
            recs.clear();
            for iter in &mut fq_iters {
                if let Some(rec) = iter.next() {
                    recs.push(rec?);
                }
            }

            if recs.is_empty() {
                break;
            }

            if recs.len() != fq_iters.len() {
                return Err(anyhow!(
                    "FASTQ sources out of sync; expected {} records but got {}.",
                    fq_iters.len(),
                    recs.len()
                ));
            }

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
    /// generated.  Read names are of the format `{prefix}{suffix}{idx}` where `idx` starts
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
        };
        sharder.execute().unwrap();

        for shard in 1..=2 {
            let path = PathBuf::from(format!("{}.chunk_{}.R1.fq.gz", prefix, shard));
            assert!(path.exists(), "expected output file does not exist: {}", path.display());
            assert_eq!(read_fastq(&path).len(), 3);
        }
    }
}
