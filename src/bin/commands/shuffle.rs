use crate::commands::command::{Command, validation_result};
use crate::commands::fastq_readahead::ReadAheadBuilder;
use crate::commands::outputs::{
    check_output_paths, create_bgzf_writers, remove_incomplete_outputs,
};
use crate::commands::read_set_reader::{ReadSetReader, base_read_name};
use anyhow::{Context, Result, anyhow, bail};
use clap::Parser;
use clap::builder::RangedU64ValueParser;
use fgoxide::io::Io;
use flate2::bufread::MultiGzDecoder;
use log::{info, warn};
use pooled_writer::PooledWriter;
use proglog::{CountFormatterKind, ProgLog, ProgLogBuilder};
use seq_io::fastq::{Record, RefRecord};
use std::fs::{self, File, OpenOptions};
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::{Receiver, SyncSender, sync_channel};
use std::sync::{Arc, Mutex};
use std::thread::{self, JoinHandle};
use tempfile::TempDir;
use xxhash_rust::xxh3::{Xxh3, xxh3_64_with_seed};

/// Number of temporary buckets used when an input has no size to estimate from (e.g. a pipe) and
/// no `--expected-reads` is given.  Too many buckets costs little, while too few makes pass two
/// re-split them; at the default memory limit this many hold about 2 TB of FASTQ before any
/// needs re-splitting.
const STREAMING_BUCKETS: usize = 1024;

/// Frame size below which more buckets aren't worth it: the bucket count is capped at one
/// bucket's share of memory divided by this.
const MIN_FRAME_BYTES: usize = 128 * 1024;

/// Largest frame of read sets compressed as one unit in a temporary file.
const MAX_FRAME_BYTES: usize = 4 * 1024 * 1024;

/// Largest block of the in-memory buffer.  Read sets are packed into blocks so the buffer grows
/// without reallocating and copying what it already holds.
const MAX_BLOCK_BYTES: usize = 4 * 1024 * 1024;

/// Bytes of a gzipped input decompressed to estimate its compression ratio.
const RATIO_PROBE_BYTES: u64 = 4 * 1024 * 1024;

/// Magic bytes at the start of every temporary bucket file.
const BUCKET_MAGIC: [u8; 4] = *b"FQS1";

/// Bytes in an encoded read set's header: input index (u64) and length of its reads (u32).
const READ_SET_HEADER_BYTES: usize = 12;

/// Bytes of the length (u32) before each read in an encoded read set.
const READ_LEN_BYTES: usize = 4;

/// Bytes of the header before each frame in a bucket file: compressed and uncompressed lengths
/// (u32 each).
const FRAME_HEADER_BYTES: usize = 8;

/// Bytes each read set adds to the in-memory sort index.
const ENTRY_BYTES: usize = size_of::<Entry>();

/// Smallest accepted `--max-memory`.
const MIN_MAX_MEMORY: usize = 1024 * 1024;

/// A frame of encoded read sets staged for a bucket: (bucket index, bytes).
type StagedFrame = (usize, Vec<u8>);

/// A compressed frame on its way to a bucket file: (bucket index, uncompressed length, bytes).
type CompressedFrame = (usize, usize, Vec<u8>);

/// Shuffles reads in one or more synchronized FASTQ files into a random order.
///
/// Reads one or more FASTQ files (e.g. paired-end R1 and R2) and writes the same
/// records in a random order. All input files must contain the same number of
/// reads in the same order; every output file gets the same order, so paired
/// reads stay in sync. Records are written as they appear in the input, except
/// that a last record missing its newline gains one.
///
/// Writes one output file per input, in the same order, to the paths given by
/// `--outputs`. Output files are always BGZF compressed. If shuffling fails once
/// its outputs are open, those that are regular files are deleted so they can't
/// be mistaken for complete ones.
///
/// The order is set by a seeded hash of each read's name (ignoring any `/1` or
/// `/2` suffix and comment), so the same reads and seed always produce the same
/// output, whatever the thread count or memory limit. Reads that share a name
/// stay next to each other in input order, and a warning says how many repeat
/// the name before them. When no `--seed` is given, the seed is derived from the name, bases and
/// qualities of the first read set, so the same reads get the same order whether
/// they come from a plain file, a gzipped file or a pipe, and whatever their
/// comments or line endings.
///
/// Inputs that fit within `--max-memory` are shuffled in memory. Larger inputs
/// are split by hash into temporary zstd-compressed files in `--tmp-dir`, which
/// are then sorted and written out one at a time. This needs temporary disk
/// space for a compressed copy of the input, plus a little more while any
/// bucket too large for memory is split again. When reading from pipes, give
/// `--expected-reads` so the temporary files can be sized to the input. If
/// shuffle is killed, delete the `fqtk-shuffle.*` directories it leaves in the
/// temporary directory.
///
/// # Example
///
/// ```bash
/// fqtk shuffle \
///     --inputs r1.fq.gz r2.fq.gz \
///     --outputs shuffled.r1.fq.gz shuffled.r2.fq.gz
/// ```
#[derive(Parser, Debug)]
#[command(version)]
#[clap(verbatim_doc_comment)]
pub(crate) struct Shuffle {
    /// One or more input FASTQ files, gzipped if their names end in .gz or
    /// .bgz. All files must have the same number of reads in the same order.
    #[clap(long, short = 'i', required = true, num_args = 1..)]
    inputs: Vec<PathBuf>,

    /// Output FASTQ files, one per input in the same order. Always BGZF
    /// compressed.
    #[clap(long, short = 'o', required = true, num_args = 1..)]
    outputs: Vec<PathBuf>,

    /// Seed for the shuffle. When omitted, the seed is derived from the first
    /// read set.
    #[clap(long, short = 's')]
    seed: Option<u64>,

    /// Number of threads to use. Minimum 2. In addition, one thread per input
    /// for read and decompression, and up to about a quarter more run briefly
    /// while temporary files are loaded.
    #[clap(long, short = 't', default_value = "8")]
    threads: usize,

    /// BGZF compression level for output files, from 0 (BGZF blocks without
    /// compression, e.g. for a pipe to another tool) to 12.
    #[clap(long, short = 'c', default_value = "5",
           value_parser = RangedU64ValueParser::<u8>::new().range(0..=12))]
    compression_level: u8,

    /// Total memory, not per thread, for holding reads while sorting them,
    /// e.g. 512M or 4G. K, M, G and T are binary (1024-based) units.
    #[clap(long, short = 'm', default_value = "4G", value_parser = parse_memory)]
    max_memory: usize,

    /// Directory for temporary files, which need room for a compressed copy of
    /// the input. Repeat to spread them across several directories. Defaults
    /// to the system temporary directory.
    #[clap(long = "tmp-dir", short = 'T', value_name = "DIR")]
    tmp_dirs: Vec<PathBuf>,

    /// Zstd compression level for temporary files, from 1 (fastest) to 19.
    #[clap(long, default_value = "1", value_parser = clap::value_parser!(i32).range(1..=19))]
    temp_compression: i32,

    /// Disable checking that read names are in sync across input files.
    #[clap(long)]
    disable_read_name_checking: bool,

    /// Approximate number of reads in each input file (i.e. read pairs for
    /// paired-end input), used to size temporary files for large inputs. Useful
    /// when reading from pipes, whose size is unknown, and takes precedence over
    /// the estimate from input file sizes. Overestimating costs little;
    /// underestimating can add an extra pass over the data.
    #[clap(long)]
    expected_reads: Option<u64>,

    /// Number of temporary buckets to split a large input into, overriding the
    /// estimate. For testing.
    #[clap(long, hide = true,
           value_parser = RangedU64ValueParser::<usize>::new().range(2..=65_536))]
    buckets: Option<usize>,
}

impl Shuffle {
    /// Shuffles the inputs into the outputs and returns counts describing the run.  If the
    /// shuffle fails once the outputs are open, removes them (see [`remove_incomplete_outputs`]).
    fn run(&self) -> Result<ShuffleStats> {
        self.validate()?;
        for warning in self.output_warnings() {
            warn!("{warning}");
        }
        // Created now, rather than when the input first outgrows memory, so that an unusable
        // --tmp-dir fails before any input is read.
        let temp = Arc::new(TempSpace::new(&self.tmp_dirs)?);

        let inputs = self
            .inputs
            .iter()
            .map(|path| ReadAheadBuilder { path: path.clone(), ..Default::default() })
            .collect();
        // Opening waits until every input is readable, so an unreadable input fails before any
        // output is created (and truncated).
        let mut read_sets = ReadSetReader::open(inputs, !self.disable_read_name_checking)?;
        // While writing, this thread only copies records, so the rest compress; pass two's loaders
        // oversubscribe briefly, which costs less than idle cores.
        let (mut pool, mut writers) =
            create_bgzf_writers(&self.outputs, self.threads - 1, self.compression_level)?;

        let shuffle_result = self.shuffle(&mut read_sets, &mut writers, temp);
        // The writers and pool are shut down whatever happened, so queued output is flushed.
        // A failure to write output also breaks the shuffle that was writing it, so it's
        // reported in preference as the root cause.
        let close_result =
            writers.into_iter().try_for_each(PooledWriter::close).map_err(anyhow::Error::from);
        let stop_result = pool.stop_pool().map_err(anyhow::Error::from);
        let result = match stop_result.and(close_result) {
            Ok(()) => shuffle_result,
            Err(e) => Err(e.context("Failed to write the outputs")),
        };
        if result.is_err() {
            remove_incomplete_outputs(&self.outputs);
        }
        result
    }

    /// Checks the arguments, returning one error that lists everything wrong.
    fn validate(&self) -> Result<()> {
        let mut errors = vec![];

        for input in &self.inputs {
            if !input.exists() {
                errors.push(format!("Input file {input:?} does not exist."));
            } else if input.is_dir() {
                errors.push(format!("Input {input:?} is a directory."));
            }
        }

        if self.threads < 2 {
            errors.push(format!("Threads must be at least 2, got {}.", self.threads));
        }

        if self.max_memory < MIN_MAX_MEMORY {
            errors.push(format!("Max memory must be at least 1M, got {} bytes.", self.max_memory));
        }

        for dir in &self.tmp_dirs {
            if !dir.exists() {
                errors.push(format!("Temporary directory {dir:?} does not exist."));
            } else if !dir.is_dir() {
                errors.push(format!("Temporary directory {dir:?} is not a directory."));
            }
        }

        if self.outputs.len() != self.inputs.len() {
            errors.push(format!(
                "Expected one output per input: got {} inputs and {} outputs.",
                self.inputs.len(),
                self.outputs.len()
            ));
        }

        check_output_paths(&self.inputs, &self.outputs, &mut errors);
        validation_result(&errors)
    }

    /// Returns warnings about output names that suggest something other than BGZF, which is
    /// always what's written.  Outputs that already exist as something other than a regular
    /// file, e.g. /dev/stdout or a pipe, are exempt.
    fn output_warnings(&self) -> Vec<String> {
        self.outputs
            .iter()
            .filter_map(|output| {
                if output.as_os_str() == "-" {
                    Some(
                        "Output \"-\" will be written to a file named \"-\"; use /dev/stdout to \
                         write to standard output."
                            .to_owned(),
                    )
                } else if Io::is_gzip_path(output)
                    || fs::metadata(output).is_ok_and(|m| !m.is_file())
                {
                    None
                } else {
                    Some(format!(
                        "Output {output:?} doesn't end in .gz or .bgz, but will be BGZF compressed."
                    ))
                }
            })
            .collect()
    }

    /// Reads every read set from `read_sets` and writes them all to `writers` in shuffled order,
    /// using `temp` for temporary files if the input outgrows memory.
    fn shuffle(
        &self,
        read_sets: &mut ReadSetReader,
        writers: &mut [PooledWriter],
        temp: Arc<TempSpace>,
    ) -> Result<ShuffleStats> {
        match self.read_input(read_sets, &temp)? {
            PassOneResult::InMemory(mut buffer) => {
                info!("Input fits in memory; sorting and writing {} read sets.", buffer.len());
                buffer.sort();
                let logger = progress_logger("wrote");
                let shared_names = buffer.write_to(writers, &logger, &mut None)?;
                Ok(ShuffleStats {
                    read_sets: buffer.len() as u64,
                    shared_names,
                    ..Default::default()
                })
            }
            PassOneResult::Spilled { buckets, read_sets, seed } => {
                let bucket_count = buckets.len();
                let stats = self.write_buckets(buckets, temp, seed, writers)?;
                if stats.read_sets != read_sets {
                    bail!(
                        "Read {read_sets} read sets but wrote {}; the temporary files are corrupt.",
                        stats.read_sets
                    );
                }
                Ok(ShuffleStats { buckets: bucket_count, ..stats })
            }
        }
    }

    /// Pass one: reads every read set, keeping them in memory until they exceed `--max-memory`
    /// and splitting them into temporary buckets in `temp` by sort key from then on.
    fn read_input(
        &self,
        read_sets: &mut ReadSetReader,
        temp: &Arc<TempSpace>,
    ) -> Result<PassOneResult> {
        let logger = progress_logger("read");
        let mut seed = self.seed;
        if let Some(seed) = seed {
            info!("Using random seed {seed}.");
        }
        let mut buffer = SortBuffer::new(self.budget().block_bytes());
        let mut router: Option<BucketRouter> = None;
        let mut encoded = Vec::new();
        let mut text_bytes: u64 = 0;
        let mut index: u64 = 0;

        while let Some(records) = read_sets.next()? {
            encoded.clear();
            encode_read_set(&mut encoded, index, &records)?;
            let seed = *seed.get_or_insert_with(|| {
                let derived = derive_seed(&records);
                info!("Using random seed {derived}, derived from the first read set.");
                derived
            });
            let key = read_set_key(&encoded, seed)?;

            if let Some(router) = router.as_mut() {
                router.route(key, &encoded)?;
            } else {
                text_bytes +=
                    (encoded.len() - READ_SET_HEADER_BYTES - READ_LEN_BYTES * records.len()) as u64;
                buffer.push(key, &encoded)?;
                if buffer.footprint() >= self.max_memory {
                    let full = std::mem::replace(&mut buffer, SortBuffer::new(0));
                    router = Some(self.start_spilling(full, text_bytes, temp)?);
                }
            }
            index += 1;
            logger.record();
        }

        match router {
            None => Ok(PassOneResult::InMemory(buffer)),
            Some(router) => {
                let buckets = router.finish()?;
                info!("Read {index} read sets into {} temporary buckets.", buckets.len());
                let seed = seed.expect("the first read set sets the seed before any spill");
                Ok(PassOneResult::Spilled { buckets, read_sets: index, seed })
            }
        }
    }

    /// Moves the full in-memory `buffer` into new temporary buckets in `temp` and returns the
    /// router that buckets the rest of the input.  `text_bytes` is the FASTQ text the buffer was
    /// read from, used to extrapolate the size of the whole input.
    fn start_spilling(
        &self,
        mut buffer: SortBuffer,
        text_bytes: u64,
        temp: &Arc<TempSpace>,
    ) -> Result<BucketRouter> {
        let bucket_count = self.bucket_count(buffer.footprint(), buffer.len(), text_bytes)?;
        info!("Input exceeds --max-memory; splitting it into {bucket_count} temporary buckets.");
        let workers = self.split_workers();
        let mut router = BucketRouter::new(
            Arc::clone(temp),
            Partition::new(0, u128::MAX, bucket_count),
            self.budget().frame_bytes(bucket_count, workers),
            self.temp_compression,
            workers,
        );

        // Routing the read sets in sorted order fills one bucket at a time, so each bucket's
        // staging buffer can be flushed as soon as it is complete rather than every bucket's
        // filling at once on top of the full buffer.
        buffer.sort();
        let mut current = 0;
        for entry in &buffer.entries {
            let bucket = router.route(entry.key(), buffer.read_set(entry))?;
            if bucket != current {
                router.flush(current)?;
                current = bucket;
            }
        }
        router.flush(current)?;
        Ok(router)
    }

    /// Chooses how many buckets to split the input into, given the in-memory `footprint` of the
    /// first `read_sets` read sets and the `text_bytes` of FASTQ they were read from.
    fn bucket_count(&self, footprint: usize, read_sets: usize, text_bytes: u64) -> Result<usize> {
        if let Some(buckets) = self.buckets {
            return Ok(buckets);
        }
        if let Some(expected) = self.expected_reads {
            // Measuring memory per read set from the input itself keeps the estimate right
            // however long the reads are, e.g. after trimming.
            let total = expected.max(read_sets as u64) as f64;
            let total_footprint = footprint as f64 * total / read_sets.max(1) as f64;
            return Ok(self.budget().buckets_for(total_footprint));
        }
        match self.estimate_input_text_bytes()? {
            Some(total) if text_bytes > 0 => {
                let total_footprint = footprint as f64 * (total as f64 / text_bytes as f64);
                Ok(self.budget().buckets_for(total_footprint))
            }
            _ => Ok(STREAMING_BUCKETS.min(self.budget().max_buckets())),
        }
    }

    /// Estimates the total uncompressed FASTQ text across all inputs, or `None` if an input
    /// isn't a regular file (e.g. a pipe) and so has no size to go on.
    fn estimate_input_text_bytes(&self) -> Result<Option<u64>> {
        let mut total = 0.0;
        for input in &self.inputs {
            let metadata =
                fs::metadata(input).with_context(|| format!("Failed to stat {input:?}"))?;
            if !metadata.is_file() {
                return Ok(None);
            }
            let ratio = if Io::is_gzip_path(input) { gzip_ratio(input)? } else { 1.0 };
            total += metadata.len() as f64 * ratio;
        }
        Ok(Some(total as u64))
    }

    /// Pass two: loads, sorts and writes each temporary bucket in turn, its read sets hashed with
    /// `seed`.  A loader thread reads and sorts the next bucket while this thread writes the
    /// current one.
    fn write_buckets(
        &self,
        buckets: Vec<Bucket>,
        temp: Arc<TempSpace>,
        seed: u64,
        writers: &mut [PooledWriter],
    ) -> Result<ShuffleStats> {
        let loader = BucketLoader {
            temp,
            budget: self.budget(),
            seed,
            zstd_level: self.temp_compression,
            load_threads: self.load_threads(),
            resplit_workers: self.resplit_workers(),
        };
        // A rendezvous channel: the loader holds the next sorted bucket until this thread has
        // written the current one, so at most two buckets are in memory at once.
        let (sender, receiver) = sync_channel::<SortBuffer>(0);
        let handle = thread::spawn(move || loader.run(buckets, &sender));

        let logger = progress_logger("wrote");
        let mut stats = ShuffleStats::default();
        let mut previous_hash = None;
        let write_result: Result<()> = receiver.iter().try_for_each(|buffer| {
            stats.shared_names += buffer.write_to(writers, &logger, &mut previous_hash)?;
            stats.read_sets += buffer.len() as u64;
            Ok(())
        });
        // Dropping the receiver before joining lets a loader blocked in `send` return if writing
        // failed.
        drop(receiver);
        let load_result = joined(handle.join(), "The bucket loader thread");
        write_result?;
        stats.resplits = load_result?;
        Ok(stats)
    }

    /// The [`MemoryBudget`] for `--max-memory`.
    fn budget(&self) -> MemoryBudget {
        MemoryBudget { total: self.max_memory }
    }

    /// Threads compressing temporary frames during pass one: whatever `--threads` leaves after
    /// this thread and one decompression thread per gzipped input.  The thread writing
    /// temporary files only does I/O and isn't counted.
    fn split_workers(&self) -> usize {
        let gzipped = self.inputs.iter().filter(|path| Io::is_gzip_path(path)).count();
        self.threads.saturating_sub(1 + gzipped).max(1)
    }

    /// Threads decompressing and indexing frames when loading a bucket during pass two.  Loading
    /// is much cheaper than compressing the output, so one loading thread per four compressing
    /// threads keeps the compressors fed, even at compression level 1 where they are fastest.
    fn load_threads(&self) -> usize {
        (self.threads - 1).div_ceil(4)
    }

    /// Threads compressing temporary frames when re-splitting a bucket during pass two.  A
    /// re-split stalls the writing thread once it runs out of loaded buckets, idling the output
    /// compression pool, so it can use all threads but the writing and loading ones.
    fn resplit_workers(&self) -> usize {
        self.threads.saturating_sub(2).max(1)
    }
}

impl Command for Shuffle {
    fn execute(&self) -> Result<()> {
        let stats = self.run()?;
        if stats.shared_names > 0 {
            warn!(
                "{} read sets have the same read name as the one before them. Reads that share a \
                 name are kept in input order, so are not shuffled relative to each other.",
                stats.shared_names
            );
        }
        info!(
            "Shuffled {} read sets using {} temporary buckets ({} re-split).",
            stats.read_sets, stats.buckets, stats.resplits
        );
        Ok(())
    }
}

/// How the `--max-memory` budget is divided between the stages of the shuffle.
#[derive(Clone, Copy, Debug)]
struct MemoryBudget {
    /// The whole budget, in bytes.
    total: usize,
}

impl MemoryBudget {
    /// Number of buckets to split `footprint` bytes of read sets into.  Aims for buckets half
    /// the size of [`MemoryBudget::per_bucket`]: too many buckets only shrinks the frames, but a
    /// bucket that doesn't fit in memory has to be re-split, rewriting its data.
    fn buckets_for(self, footprint: f64) -> usize {
        let buckets = (2.0 * footprint / self.per_bucket() as f64).ceil() as usize;
        buckets.clamp(2, self.max_buckets())
    }

    /// Largest number of buckets worth staging at once.
    fn max_buckets(self) -> usize {
        (self.per_bucket() / MIN_FRAME_BYTES).max(2)
    }

    /// Bytes of read sets to stage per bucket before compressing them as one frame.  The staged
    /// frames of every bucket, plus up to three frames per worker in flight (queued,
    /// compressing, queued for the writer) and one being written, must fit in
    /// [`MemoryBudget::per_bucket`].
    fn frame_bytes(self, buckets: usize, workers: usize) -> usize {
        (self.per_bucket() / (buckets + 3 * workers + 1)).clamp(1, MAX_FRAME_BYTES)
    }

    /// Memory for one bucket in pass two, which holds the bucket being written alongside the
    /// next bucket being loaded.
    fn per_bucket(self) -> usize {
        self.total / 2
    }

    /// Size of each block of the in-memory buffer: small enough that the partly filled last
    /// block wastes little of the budget.
    fn block_bytes(self) -> usize {
        (self.total / 16).min(MAX_BLOCK_BYTES)
    }
}

/// The result of pass one: either every read set in memory, or temporary buckets on disk.
enum PassOneResult {
    /// Every read set, in memory.
    InMemory(SortBuffer),
    Spilled {
        /// Non-empty buckets in sort-key order.
        buckets: Vec<Bucket>,
        /// Number of read sets in all the buckets together.
        read_sets: u64,
        /// The seed the read sets were hashed with.
        seed: u64,
    },
}

/// Pass two's loading thread: loads temporary buckets in sort-key order, re-splitting any too
/// large for memory, and sends each one, sorted, to the thread writing the output.
struct BucketLoader {
    temp: Arc<TempSpace>,
    budget: MemoryBudget,
    /// The seed the read sets were hashed with in pass one.
    seed: u64,
    zstd_level: i32,
    /// Threads decompressing and indexing frames while loading a bucket.
    load_threads: usize,
    /// Threads compressing frames while re-splitting a bucket.
    resplit_workers: usize,
}

impl BucketLoader {
    /// Sends every bucket's read sets, sorted, to `sender` in sort-key order, re-splitting any
    /// bucket too large for memory, and returns the number of buckets re-split.  Stops early,
    /// without error, if the receiving thread hangs up.
    fn run(self, buckets: Vec<Bucket>, sender: &SyncSender<SortBuffer>) -> Result<usize> {
        let mut resplits = 0;
        // A stack, so that a re-split bucket's smaller buckets are all sent before the next.
        let mut pending: Vec<Bucket> = buckets.into_iter().rev().collect();
        while let Some(bucket) = pending.pop() {
            // A single read set can't be split further, however large it is.
            if bucket.footprint() <= self.budget.per_bucket() || bucket.read_sets < 2 {
                let mut buffer = bucket.load(self.load_threads, self.seed)?;
                buffer.sort();
                if sender.send(buffer).is_err() {
                    break;
                }
            } else {
                pending.extend(self.resplit(bucket)?.into_iter().rev());
                resplits += 1;
            }
        }
        Ok(resplits)
    }

    /// Splits an oversized `bucket` into smaller buckets covering its range of sort keys,
    /// returned in sort-key order, and deletes its file.  The range runs from the bucket's
    /// smallest key to its largest, so read sets sharing a name (and so a hash) split by
    /// input index.
    fn resplit(&self, bucket: Bucket) -> Result<Vec<Bucket>> {
        let count = self.budget.buckets_for(bucket.footprint() as f64);
        info!(
            "Re-splitting a bucket of {} read sets ({} bytes) that exceeds memory into {count}.",
            bucket.read_sets, bucket.raw_bytes
        );
        let corrupt = || anyhow!("Temporary file {:?} is corrupt.", bucket.path);
        let mut router = BucketRouter::new(
            Arc::clone(&self.temp),
            Partition::new(bucket.min_key, bucket.max_key, count),
            self.budget.frame_bytes(count, self.resplit_workers),
            self.zstd_level,
            self.resplit_workers,
        );
        let mut reader = FrameReader::open(&bucket.path)?;
        let mut decompressor = zstd::bulk::Decompressor::new()?;
        let mut routed: u64 = 0;
        while let Some(stored) = reader.next_stored()? {
            let frame = stored
                .decompress(&mut decompressor)
                .with_context(|| format!("Failed to read temporary file {:?}", bucket.path))?;
            let mut offset = 0;
            while offset < frame.len() {
                let len = read_set_len(&frame[offset..])?;
                let read_set = &frame[offset..offset + len];
                let key = read_set_key(read_set, self.seed)?;
                if !(bucket.min_key..=bucket.max_key).contains(&key) {
                    return Err(corrupt());
                }
                router.route(key, read_set)?;
                routed += 1;
                offset += len;
            }
        }
        if routed != bucket.read_sets {
            return Err(corrupt());
        }
        let smaller = router.finish()?;
        fs::remove_file(&bucket.path)
            .with_context(|| format!("Failed to remove temporary file {:?}", bucket.path))?;
        Ok(smaller)
    }
}

/// Appends read sets to temporary bucket files by sort key.  Read sets are staged per bucket and
/// compressed a frame at a time by worker threads.  A single writer thread appends the finished
/// frames, opening and closing a bucket's file for each one: the number of buckets is then not
/// limited by open file handles, and no two threads can interleave writes to a file.
struct BucketRouter {
    /// Keeps the bucket files' directories alive until the worker threads are done with them.
    _temp: Arc<TempSpace>,
    partition: Partition,
    /// One per partition bucket, in sort-key order.
    buckets: Vec<Bucket>,
    /// Encoded read sets waiting to fill a frame, one buffer per bucket.
    staging: Vec<Vec<u8>>,
    frame_bytes: usize,
    /// Queue of staged frames for the workers; `None` once shut down.
    jobs: Option<SyncSender<StagedFrame>>,
    /// The compression threads, then the writer thread, so that errors from the compression
    /// threads are reported first.
    threads: Vec<JoinHandle<Result<()>>>,
}

impl BucketRouter {
    /// Creates one bucket file path per bucket of `partition` in `temp`, and starts `workers`
    /// threads that compress staged read sets at `zstd_level`, a frame of up to about
    /// `frame_bytes` at a time, and the thread that writes the frames.
    fn new(
        temp: Arc<TempSpace>,
        partition: Partition,
        frame_bytes: usize,
        zstd_level: i32,
        workers: usize,
    ) -> Self {
        let paths: Vec<PathBuf> = (0..partition.count).map(|_| temp.next_path()).collect();
        let (job_sender, job_receiver) = sync_channel::<StagedFrame>(workers);
        let (frame_sender, frame_receiver) = sync_channel::<CompressedFrame>(workers);

        let job_receiver = Arc::new(Mutex::new(job_receiver));
        let mut threads: Vec<_> = (0..workers)
            .map(|_| {
                let jobs = Arc::clone(&job_receiver);
                let frames = frame_sender.clone();
                thread::spawn(move || Self::compress_frames(&jobs, &frames, zstd_level))
            })
            .collect();
        // Only the workers hold senders now, so the writer stops once they all finish.
        drop(frame_sender);
        let writer_paths = paths.clone();
        threads.push(thread::spawn(move || Self::write_frames(&writer_paths, &frame_receiver)));

        let buckets = paths
            .into_iter()
            .map(|path| Bucket { path, raw_bytes: 0, read_sets: 0, min_key: u128::MAX, max_key: 0 })
            .collect();
        Self {
            _temp: temp,
            partition,
            buckets,
            staging: vec![Vec::new(); partition.count],
            frame_bytes,
            jobs: Some(job_sender),
            threads,
        }
    }

    /// Stages the encoded read set `read_set` for the bucket covering `key`, sending the bucket's
    /// staged read sets off to be compressed whenever they would overflow a frame.  Returns the
    /// bucket's index.
    fn route(&mut self, key: u128, read_set: &[u8]) -> Result<usize> {
        let index = self.partition.bucket_of(key);
        let bucket = &mut self.buckets[index];
        bucket.read_sets += 1;
        bucket.raw_bytes += read_set.len();
        bucket.min_key = bucket.min_key.min(key);
        bucket.max_key = bucket.max_key.max(key);

        // Flushing before a read set would overflow the frame, rather than after, keeps the
        // buffer within the capacity reserved for it instead of growing and copying it.
        let staged_len = self.staging[index].len();
        if staged_len > 0 && staged_len + read_set.len() > self.frame_bytes {
            self.flush(index)?;
        }
        let staged = &mut self.staging[index];
        if staged.capacity() == 0 {
            staged.reserve_exact(self.frame_bytes.max(read_set.len()));
        }
        staged.extend_from_slice(read_set);
        if staged.len() >= self.frame_bytes {
            self.flush(index)?;
        }
        Ok(index)
    }

    /// Writes out all staged read sets, waits for the threads to finish, and returns the
    /// non-empty buckets in sort-key order.
    fn finish(mut self) -> Result<Vec<Bucket>> {
        for index in 0..self.staging.len() {
            self.flush(index)?;
        }
        self.shutdown()?;
        Ok(std::mem::take(&mut self.buckets).into_iter().filter(|b| b.read_sets > 0).collect())
    }

    /// Sends the read sets staged for bucket `index`, if any, to be compressed and written.
    fn flush(&mut self, index: usize) -> Result<()> {
        if self.staging[index].is_empty() {
            return Ok(());
        }
        let frame = std::mem::take(&mut self.staging[index]);
        let sent = self.jobs.as_ref().is_some_and(|jobs| jobs.send((index, frame)).is_ok());
        if sent {
            Ok(())
        } else {
            // The workers only hang up after failing; report why.
            Err(self.shutdown().err().unwrap_or_else(|| anyhow!("Temporary file threads stopped.")))
        }
    }

    /// Closes the job queue and waits for the threads, returning the first error any of them
    /// hit.
    fn shutdown(&mut self) -> Result<()> {
        // Dropping the sender lets the workers drain the queue and exit.
        self.jobs = None;
        let mut result = Ok(());
        for handle in self.threads.drain(..) {
            let outcome = joined(handle.join(), "A temporary file thread");
            if result.is_ok() {
                result = outcome;
            }
        }
        result
    }

    /// Worker thread body: compresses staged frames from `jobs` and passes them to the writer.
    /// Frames carry a checksum, so corruption of the temporary files is caught when they are
    /// read back.
    fn compress_frames(
        jobs: &Mutex<Receiver<StagedFrame>>,
        frames: &SyncSender<CompressedFrame>,
        zstd_level: i32,
    ) -> Result<()> {
        let mut compressor = zstd::bulk::Compressor::new(zstd_level)?;
        compressor.set_parameter(zstd::zstd_safe::CParameter::ChecksumFlag(true))?;
        loop {
            // The lock is held while waiting for the next frame, not while compressing it.
            let job = jobs.lock().map_err(|_| anyhow!("Temporary file queue poisoned."))?.recv();
            let Ok((bucket, raw)) = job else { return Ok(()) };
            let compressed = compressor.compress(&raw)?;
            // The writer only hangs up after failing, and reports its own error.
            if frames.send((bucket, raw.len(), compressed)).is_err() {
                return Ok(());
            }
        }
    }

    /// Writer thread body: appends each compressed frame from `frames` to its bucket's file in
    /// `paths`, starting each file with [`BUCKET_MAGIC`].
    fn write_frames(paths: &[PathBuf], frames: &Receiver<CompressedFrame>) -> Result<()> {
        let mut started = vec![false; paths.len()];
        for (bucket, raw_len, compressed) in frames {
            let path = &paths[bucket];
            let mut header = Vec::with_capacity(BUCKET_MAGIC.len() + FRAME_HEADER_BYTES);
            if !started[bucket] {
                header.extend_from_slice(&BUCKET_MAGIC);
                started[bucket] = true;
            }
            header.extend_from_slice(&u32::try_from(compressed.len())?.to_le_bytes());
            header.extend_from_slice(&u32::try_from(raw_len)?.to_le_bytes());

            let mut file = OpenOptions::new()
                .create(true)
                .append(true)
                .open(path)
                .with_context(|| format!("Failed to open temporary file {path:?}"))?;
            file.write_all(&header)
                .and_then(|()| file.write_all(&compressed))
                .with_context(|| format!("Failed to write temporary file {path:?}"))?;
        }
        Ok(())
    }
}

impl Drop for BucketRouter {
    /// Stops the threads if the router is dropped without finishing, e.g. on a read error.
    fn drop(&mut self) {
        let _ = self.shutdown();
    }
}

/// A temporary file of read sets whose sort keys fall in one range.
struct Bucket {
    path: PathBuf,
    /// Total bytes of encoded read sets in the file, uncompressed.
    raw_bytes: usize,
    read_sets: u64,
    /// Smallest sort key in the file.
    min_key: u128,
    /// Largest sort key in the file.
    max_key: u128,
}

impl Bucket {
    /// Reads this bucket's read sets, hashed with `seed`, into memory unsorted, and deletes its
    /// file.  Frames are decompressed and indexed by `threads` threads as they are read, so only
    /// a few compressed frames are held at once.
    fn load(self, threads: usize, seed: u64) -> Result<SortBuffer> {
        let mut reader = FrameReader::open(&self.path)?;
        let (sender, receiver) = sync_channel::<(usize, StoredFrame)>(threads);
        // Only the workers hold the receiver, so if they all fail the sender sees a hang-up
        // rather than blocking forever.
        let receiver = Arc::new(Mutex::new(receiver));
        let (worker_results, read_result) = thread::scope(|scope| {
            let workers: Vec<_> = (0..threads)
                .map(|_| {
                    let jobs = Arc::clone(&receiver);
                    scope.spawn(move || Self::decompress_frames(&jobs, seed))
                })
                .collect();
            drop(receiver);
            let read_result = (|| -> Result<()> {
                for index in 0.. {
                    let Some(frame) = reader.next_stored()? else { break };
                    // The workers only hang up after failing, and report their own errors.
                    if sender.send((index, frame)).is_err() {
                        break;
                    }
                }
                Ok(())
            })();
            // Hanging up lets the workers finish the queued frames and exit.
            drop(sender);
            let worker_results: Vec<Result<Vec<IndexedFrame>>> = workers
                .into_iter()
                .map(|worker| joined(worker.join(), "A loader thread"))
                .collect();
            (worker_results, read_result)
        });
        let mut frames = Vec::new();
        for result in worker_results {
            frames.extend(result.with_context(|| format!("Failed to load {:?}", self.path))?);
        }
        read_result?;
        drop(reader);
        fs::remove_file(&self.path)
            .with_context(|| format!("Failed to remove temporary file {:?}", self.path))?;

        let buffer = SortBuffer::from_frames(frames);
        if buffer.len() as u64 != self.read_sets {
            bail!(
                "Temporary file {:?} held {} read sets but {} were written to it.",
                self.path,
                buffer.len(),
                self.read_sets
            );
        }
        Ok(buffer)
    }

    /// Memory needed to hold and sort this bucket's read sets.
    fn footprint(&self) -> usize {
        self.raw_bytes + self.read_sets as usize * ENTRY_BYTES
    }

    /// Loader worker thread body: decompresses and indexes numbered frames from `jobs` while
    /// they are still in cache, returning them in the order they were taken.
    fn decompress_frames(
        jobs: &Mutex<Receiver<(usize, StoredFrame)>>,
        seed: u64,
    ) -> Result<Vec<IndexedFrame>> {
        let mut decompressor = zstd::bulk::Decompressor::new()?;
        let mut frames = Vec::new();
        loop {
            // The lock is held while waiting for the next frame, not while decompressing it.
            let job = jobs.lock().map_err(|_| anyhow!("Loader queue poisoned."))?.recv();
            let Ok((index, stored)) = job else { return Ok(frames) };
            let block = stored.decompress(&mut decompressor)?;
            frames.push(IndexedFrame::new(index, block, seed)?);
        }
    }
}

/// Encoded read sets held in memory, plus an index for sorting them.  Read sets are packed into
/// blocks, and no read set spans two blocks.
struct SortBuffer {
    blocks: Vec<Vec<u8>>,
    entries: Vec<Entry>,
    /// Capacity of each new block, unless a single read set needs more.
    block_bytes: usize,
    /// Total capacity of `blocks`.
    block_capacity: usize,
}

impl SortBuffer {
    /// Creates an empty buffer that allocates blocks of `block_bytes`.
    fn new(block_bytes: usize) -> Self {
        Self { blocks: Vec::new(), entries: Vec::new(), block_bytes, block_capacity: 0 }
    }

    /// Creates a buffer from all the frames of a bucket, in any order.
    fn from_frames(mut frames: Vec<IndexedFrame>) -> Self {
        frames.sort_unstable_by_key(|frame| frame.index);
        debug_assert!(frames.iter().enumerate().all(|(i, frame)| frame.index == i));
        let mut entries = Vec::with_capacity(frames.iter().map(|f| f.entries.len()).sum());
        let mut blocks = Vec::with_capacity(frames.len());
        for frame in frames {
            entries.extend_from_slice(&frame.entries);
            blocks.push(frame.block);
        }
        let block_capacity = blocks.iter().map(Vec::capacity).sum();
        Self { blocks, entries, block_bytes: 0, block_capacity }
    }

    /// Appends the encoded read set `read_set`, whose sort key is `key`.
    fn push(&mut self, key: u128, read_set: &[u8]) -> Result<()> {
        let fits = self.blocks.last().is_some_and(|b| b.capacity() - b.len() >= read_set.len());
        if !fits {
            let block = Vec::with_capacity(self.block_bytes.max(read_set.len()));
            self.block_capacity += block.capacity();
            self.blocks.push(block);
        }
        let block_index = self.blocks.len() - 1;
        let block = &mut self.blocks[block_index];
        let offset = block.len();
        block.extend_from_slice(read_set);
        self.entries.push(Entry::new(key, block_index, offset)?);
        Ok(())
    }

    fn len(&self) -> usize {
        self.entries.len()
    }

    /// Memory used by the read sets and their index.
    fn footprint(&self) -> usize {
        self.block_capacity + self.entries.capacity() * ENTRY_BYTES
    }

    fn sort(&mut self) {
        self.entries.sort_unstable();
    }

    /// Writes the read sets, in `entries` order, to `writers`, one per input.  Returns how many
    /// read sets had the same name hash as the read set before them; `previous_hash` carries the
    /// last hash written from one buffer to the next.
    fn write_to(
        &self,
        writers: &mut [PooledWriter],
        logger: &ProgLog,
        previous_hash: &mut Option<u64>,
    ) -> Result<u64> {
        let mut shared_names = 0;
        for entry in &self.entries {
            if *previous_hash == Some(entry.hash) {
                shared_names += 1;
            }
            *previous_hash = Some(entry.hash);
            let mut reads = &self.read_set(entry)[READ_SET_HEADER_BYTES..];
            for writer in writers.iter_mut() {
                let (read, rest) = split_read(reads)?;
                writer.write_all(read)?;
                reads = rest;
            }
            logger.record();
        }
        Ok(shared_names)
    }

    /// Returns the encoded read set that `entry` points to.
    fn read_set(&self, entry: &Entry) -> &[u8] {
        let bytes = &self.blocks[entry.block as usize][entry.offset as usize..];
        let reads_len = u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]);
        &bytes[..READ_SET_HEADER_BYTES + reads_len as usize]
    }
}

/// Counts describing a completed shuffle.
#[derive(Debug, Default)]
struct ShuffleStats {
    read_sets: u64,
    /// Temporary buckets the input was split into; zero if it was shuffled in memory.
    buckets: usize,
    /// Buckets too large for memory that had to be split again.
    resplits: usize,
    /// Read sets whose name hash equals the previous read set's, almost always because they
    /// share a read name.
    shared_names: u64,
}

/// A decompressed frame of a bucket with the index entries for its read sets.
struct IndexedFrame {
    /// Position of the frame in its bucket file, and so its block index once loaded.
    index: usize,
    block: Vec<u8>,
    entries: Vec<Entry>,
}

impl IndexedFrame {
    /// Indexes every encoded read set in `block`, the frame at position `index` in its bucket
    /// file, hashing names with `seed`.
    fn new(index: usize, block: Vec<u8>, seed: u64) -> Result<Self> {
        let mut entries = Vec::new();
        let mut offset = 0;
        while offset < block.len() {
            let len = read_set_len(&block[offset..])?;
            entries.push(Entry::new(read_set_key(&block[offset..], seed)?, index, offset)?);
            offset += len;
        }
        Ok(Self { index, block, entries })
    }
}

/// Sort-index entry for one read set: its sort key and where its bytes are in a [`SortBuffer`].
/// Fields are ordered so that the derived ordering sorts by key: hash, then input index.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
struct Entry {
    hash: u64,
    input_index: u64,
    block: u32,
    offset: u32,
}

impl Entry {
    /// Creates an entry for the read set with sort key `key` at `offset` in block `block`.
    fn new(key: u128, block: usize, offset: usize) -> Result<Self> {
        Ok(Self {
            hash: (key >> 64) as u64,
            input_index: key as u64,
            block: u32::try_from(block)?,
            offset: u32::try_from(offset)?,
        })
    }

    fn key(&self) -> u128 {
        sort_key(self.hash, self.input_index)
    }
}

/// Reads the frames of a temporary bucket file back.
struct FrameReader {
    path: PathBuf,
    reader: BufReader<File>,
}

impl FrameReader {
    /// Opens the bucket file at `path` and checks that it starts with [`BUCKET_MAGIC`].
    fn open(path: &Path) -> Result<Self> {
        let file =
            File::open(path).with_context(|| format!("Failed to open temporary file {path:?}"))?;
        let mut reader = BufReader::new(file);
        let mut magic = [0u8; BUCKET_MAGIC.len()];
        reader
            .read_exact(&mut magic)
            .with_context(|| format!("Failed to read temporary file {path:?}"))?;
        if magic != BUCKET_MAGIC {
            bail!("Temporary file {path:?} is not a shuffle bucket file.");
        }
        Ok(Self { path: path.to_path_buf(), reader })
    }

    /// Returns the next frame as stored, or `None` at the end of the file.
    fn next_stored(&mut self) -> Result<Option<StoredFrame>> {
        let context = || format!("Failed to read temporary file {:?}", self.path);
        if self.reader.fill_buf().with_context(context)?.is_empty() {
            return Ok(None);
        }
        let mut header = [0u8; FRAME_HEADER_BYTES];
        self.reader.read_exact(&mut header).with_context(context)?;
        let compressed_len = u32::from_le_bytes([header[0], header[1], header[2], header[3]]);
        let raw_len = u32::from_le_bytes([header[4], header[5], header[6], header[7]]) as usize;

        let mut compressed = vec![0; compressed_len as usize];
        self.reader.read_exact(&mut compressed).with_context(context)?;
        Ok(Some(StoredFrame { compressed, raw_len }))
    }
}

/// A frame of a bucket file as stored: zstd-compressed encoded read sets.
struct StoredFrame {
    compressed: Vec<u8>,
    /// Length of the frame once decompressed.
    raw_len: usize,
}

impl StoredFrame {
    /// Decompresses the frame with `decompressor`, checking its checksum and length.
    fn decompress(&self, decompressor: &mut zstd::bulk::Decompressor) -> Result<Vec<u8>> {
        let frame = decompressor
            .decompress(&self.compressed, self.raw_len)
            .context("Corrupt frame in temporary file")?;
        if frame.len() != self.raw_len {
            bail!(
                "Corrupt frame in temporary file: {} bytes, expected {}.",
                frame.len(),
                self.raw_len
            );
        }
        Ok(frame)
    }
}

/// Divides an inclusive range of sort keys into equal-width consecutive buckets, so every key
/// in bucket `i` sorts before every key in bucket `i + 1`.
#[derive(Clone, Copy, Debug)]
struct Partition {
    low: u128,
    /// Keys per bucket.
    width: u128,
    count: usize,
}

impl Partition {
    /// Divides the keys from `low` to `high` inclusive into `count` buckets, at least two.
    fn new(low: u128, high: u128, count: usize) -> Self {
        assert!(count >= 2 && low <= high, "invalid partition: {count} buckets of {low}..={high}");
        // ceil((high - low + 1) / count), written so that the full u128 range can't overflow.
        Self { low, width: (high - low) / count as u128 + 1, count }
    }

    /// Returns the bucket holding `key`, which must lie within the partitioned range.
    fn bucket_of(&self, key: u128) -> usize {
        ((key - self.low) / self.width) as usize
    }
}

/// Temporary directories for bucket files, each deleted with its contents when dropped.  Bucket
/// files are spread across the directories in turn.
struct TempSpace {
    dirs: Vec<TempDir>,
    /// Number of the next path to hand out.
    next: AtomicUsize,
}

impl TempSpace {
    /// Creates a temporary directory inside each of `parents`, or in the system temporary
    /// directory if there are none.
    fn new(parents: &[PathBuf]) -> Result<Self> {
        let mut builder = tempfile::Builder::new();
        builder.prefix("fqtk-shuffle.");
        let dirs = if parents.is_empty() {
            vec![builder.tempdir().context("Failed to create a temporary directory")?]
        } else {
            parents
                .iter()
                .map(|parent| {
                    builder.tempdir_in(parent).with_context(|| {
                        format!("Failed to create a temporary directory in {parent:?}")
                    })
                })
                .collect::<Result<_>>()?
        };
        Ok(Self { dirs, next: AtomicUsize::new(0) })
    }

    /// Returns a new, unique path for a bucket file.
    fn next_path(&self) -> PathBuf {
        let n = self.next.fetch_add(1, Ordering::Relaxed);
        self.dirs[n % self.dirs.len()].path().join(format!("bucket-{n:06}.fqs"))
    }
}

/// Builds a progress logger for read sets, e.g. "read" or "wrote".
fn progress_logger(verb: &str) -> ProgLog {
    ProgLogBuilder::new()
        .name("fqtk")
        .noun("read sets")
        .verb(verb)
        .unit(5_000_000)
        .count_formatter(CountFormatterKind::Comma)
        .level(log::Level::Info)
        .build()
}

/// Appends a read set to `out` in the in-memory and temporary-file encoding: its input `index`,
/// the length of what follows, then each record exactly as it appeared in its input, prefixed
/// by its length.  All integers are little-endian.
fn encode_read_set(out: &mut Vec<u8>, index: u64, records: &[RefRecord]) -> Result<()> {
    let start = out.len();
    out.extend_from_slice(&index.to_le_bytes());
    out.extend_from_slice(&[0; 4]);
    for record in records {
        let len_at = out.len();
        out.extend_from_slice(&[0; READ_LEN_BYTES]);
        record.write_unchanged(&mut *out)?;
        let len =
            u32::try_from(out.len() - len_at - READ_LEN_BYTES).context("FASTQ record too long")?;
        out[len_at..len_at + READ_LEN_BYTES].copy_from_slice(&len.to_le_bytes());
    }
    let reads_len =
        u32::try_from(out.len() - start - READ_SET_HEADER_BYTES).context("Read set too long")?;
    out[start + 8..start + READ_SET_HEADER_BYTES].copy_from_slice(&reads_len.to_le_bytes());
    Ok(())
}

/// Combines a read set's hash and input index into its sort key.  The index breaks ties
/// between equal hashes, so every read set has a distinct key and the output order doesn't
/// depend on how the read sets were buffered or bucketed.
fn sort_key(hash: u64, index: u64) -> u128 {
    (u128::from(hash) << 64) | u128::from(index)
}

/// Returns the sort key of the encoded read set at the start of `read_set`, hashing its name with
/// `seed`.  `read_set` must hold at least a full header (see [`read_set_len`]).
fn read_set_key(read_set: &[u8], seed: u64) -> Result<u128> {
    let index = u64::from_le_bytes(read_set[..8].try_into().expect("an 8-byte slice"));
    Ok(sort_key(read_set_hash(read_set, seed)?, index))
}

/// Derives a seed from a read set: a hash of its first read's [`base_read_name`], then of the
/// bases and qualities of each of its reads.  Comments, `/1` and `/2` suffixes and line endings
/// don't affect it.
fn derive_seed(records: &[RefRecord]) -> u64 {
    let mut hasher = Xxh3::new();
    hasher.update(base_read_name(records[0].head()));
    for record in records {
        hasher.update(record.seq());
        hasher.update(record.qual());
    }
    hasher.digest()
}

/// Hashes the name of the encoded read set `read_set` with `seed`: the first read's header, as
/// seq_io parses it (without the `@` and any trailing `\r`), reduced by [`base_read_name`].  The
/// hash isn't stored in temporary files, which saves bytes that wouldn't compress, so both
/// passes compute it with this one function from the same bytes.  xxh3's output is the same on
/// every platform and in every version, unlike some faster hashes', so a seed always gives the
/// same order.
fn read_set_hash(read_set: &[u8], seed: u64) -> Result<u64> {
    let (first_read, _) = split_read(&read_set[READ_SET_HEADER_BYTES..])?;
    let header_end = first_read.iter().position(|&b| b == b'\n').unwrap_or(first_read.len());
    let header = first_read
        .get(1..header_end)
        .ok_or_else(|| anyhow!("A read set's first read has no header line."))?;
    let header = header.strip_suffix(b"\r").unwrap_or(header);
    Ok(xxh3_64_with_seed(base_read_name(header), seed))
}

/// Returns the length of the encoded read set at the start of `bytes`, checking it is all there.
fn read_set_len(bytes: &[u8]) -> Result<usize> {
    let truncated = || anyhow!("Truncated read set in temporary data.");
    let header = bytes.get(..READ_SET_HEADER_BYTES).ok_or_else(truncated)?;
    let reads_len = u32::from_le_bytes([header[8], header[9], header[10], header[11]]);
    let len = READ_SET_HEADER_BYTES + reads_len as usize;
    if len > bytes.len() {
        return Err(truncated());
    }
    Ok(len)
}

/// Splits the first length-prefixed FASTQ record off the encoded reads of a read set, returning
/// the record and the remaining bytes.
fn split_read(reads: &[u8]) -> Result<(&[u8], &[u8])> {
    let truncated = || anyhow!("Truncated read in temporary data.");
    let len_bytes = reads.get(..READ_LEN_BYTES).ok_or_else(truncated)?;
    let len = u32::from_le_bytes([len_bytes[0], len_bytes[1], len_bytes[2], len_bytes[3]]);
    let end = READ_LEN_BYTES + len as usize;
    let read = reads.get(READ_LEN_BYTES..end).ok_or_else(truncated)?;
    Ok((read, &reads[end..]))
}

/// Returns a joined thread's result, or an error naming the thread, `name`, if it panicked.
fn joined<T>(result: thread::Result<Result<T>>, name: &str) -> Result<T> {
    result.unwrap_or_else(|payload| {
        let cause = payload
            .downcast_ref::<&str>()
            .copied()
            .or_else(|| payload.downcast_ref::<String>().map(String::as_str))
            .unwrap_or("unknown cause");
        Err(anyhow!("{name} panicked: {cause}"))
    })
}

/// Estimates the compression ratio of the gzipped file at `path` by decompressing up to its
/// first [`RATIO_PROBE_BYTES`] bytes.
fn gzip_ratio(path: &Path) -> Result<f64> {
    let file = File::open(path).with_context(|| format!("Failed to open {path:?}"))?;
    let compressed = fs::metadata(path)?.len().min(RATIO_PROBE_BYTES);
    let mut decoder = MultiGzDecoder::new(BufReader::new(file.take(RATIO_PROBE_BYTES)));
    let mut buffer = vec![0u8; 64 * 1024];
    let mut decompressed: u64 = 0;
    // Decompression fails where the probe cuts the gzip stream off; everything decompressed
    // before that still counts.
    while let Ok(n @ 1..) = decoder.read(&mut buffer) {
        decompressed += n as u64;
    }
    Ok(if compressed == 0 { 1.0 } else { decompressed as f64 / compressed as f64 })
}

/// Parses a memory size such as `4G`, `512M`, `1.5GiB` or `1048576` into bytes.  The suffixes
/// K, M, G and T are binary units and may be followed by `B` or `iB`; case is ignored.
fn parse_memory(value: &str) -> Result<usize, String> {
    let upper = value.trim().to_ascii_uppercase();
    if upper.starts_with('-') {
        return Err(format!("Memory size {value:?} is negative."));
    }
    let unit_start = upper.find(|c: char| !(c.is_ascii_digit() || c == '.')).unwrap_or(upper.len());
    let (number, unit) = upper.split_at(unit_start);
    let multiplier: u64 = match unit {
        "" | "B" => 1,
        "K" | "KB" | "KIB" => 1 << 10,
        "M" | "MB" | "MIB" => 1 << 20,
        "G" | "GB" | "GIB" => 1 << 30,
        "T" | "TB" | "TIB" => 1 << 40,
        _ => return Err(format!("Invalid memory size {value:?}: unknown unit {unit:?}.")),
    };
    let number: f64 = number.parse().map_err(|_| format!("Invalid memory size {value:?}."))?;
    let bytes = number * multiplier as f64;
    if bytes >= usize::MAX as f64 {
        return Err(format!("Memory size {value:?} is too large."));
    }
    Ok(bytes as usize)
}

#[cfg(test)]
mod tests {
    use super::*;
    use flate2::Compression;
    use flate2::write::GzEncoder;
    use seq_io::fastq::{OwnedRecord, Reader as FastqReader};

    /// The order of `read0` .. `read9` when shuffled with seed 7.
    const GOLDEN_ORDER: [&str; 10] =
        ["read6", "read9", "read3", "read2", "read4", "read0", "read8", "read7", "read5", "read1"];

    /// Returns 40 pseudo-random bases determined by `i`, so reads are distinguishable.
    fn bases(i: usize) -> String {
        let mut state = i as u64;
        (0..40)
            .map(|_| {
                state = state
                    .wrapping_mul(6_364_136_223_846_793_005)
                    .wrapping_add(1_442_695_040_888_963_407);
                char::from(b"ACGT"[(state >> 62) as usize])
            })
            .collect()
    }

    /// Returns `count` (header, bases) reads named `read{i}` followed by `suffix`, e.g. "/1".
    fn reads(count: usize, suffix: &str) -> Vec<(String, String)> {
        (0..count).map(|i| (format!("read{i}{suffix}"), bases(i))).collect()
    }

    /// Returns (header, bases) reads as FASTQ text.
    fn fastq_text(reads: &[(String, String)]) -> String {
        reads.iter().fold(String::new(), |mut text, (head, seq)| {
            text.push_str(&format!("@{head}\n{seq}\n+\n{}\n", "I".repeat(seq.len())));
            text
        })
    }

    /// Writes (header, bases) reads as FASTQ to `path`, gzipped if it ends in `.gz`.
    fn write_fastq(path: &Path, reads: &[(String, String)]) -> PathBuf {
        let mut writer = Io::default().new_writer(path).unwrap();
        writer.write_all(fastq_text(reads).as_bytes()).unwrap();
        writer.flush().unwrap();
        path.to_path_buf()
    }

    /// Writes `count` read pairs to R1 and R2 FASTQs named `in.R1.{extension}` etc. in `dir`.
    fn write_pairs(dir: &Path, count: usize, extension: &str) -> Vec<PathBuf> {
        vec![
            write_fastq(&dir.join(format!("in.R1.{extension}")), &reads(count, "/1")),
            write_fastq(&dir.join(format!("in.R2.{extension}")), &reads(count, "/2")),
        ]
    }

    /// A shuffle of `inputs` to `{dir}/{output}.R1.fq.gz` etc. with a fixed seed and enough
    /// memory to shuffle small inputs in memory.
    fn shuffle_of(dir: &Path, output: &str, inputs: &[PathBuf]) -> Shuffle {
        Shuffle {
            inputs: inputs.to_vec(),
            outputs: (1..=inputs.len()).map(|n| dir.join(format!("{output}.R{n}.fq.gz"))).collect(),
            seed: Some(7),
            threads: 2,
            compression_level: 1,
            max_memory: 64 * 1024 * 1024,
            tmp_dirs: vec![],
            temp_compression: 1,
            disable_read_name_checking: false,
            expected_reads: None,
            buckets: None,
        }
    }

    /// A shuffle like [`shuffle_of`] but with so little memory that it splits its input into
    /// `buckets` temporary buckets.
    fn spilling_shuffle_of(
        dir: &Path,
        output: &str,
        inputs: &[PathBuf],
        buckets: usize,
    ) -> Shuffle {
        Shuffle {
            max_memory: MIN_MAX_MEMORY,
            buckets: Some(buckets),
            ..shuffle_of(dir, output, inputs)
        }
    }

    /// Reads every record from the shuffle's output for input `n` (1-based).
    fn output_records(cmd: &Shuffle, n: usize) -> Vec<OwnedRecord> {
        FastqReader::new(Io::default().new_reader(&cmd.outputs[n - 1]).unwrap())
            .into_records()
            .collect::<Result<Vec<_>, _>>()
            .unwrap()
    }

    /// Returns the decompressed text of the shuffle's output for input `n` (1-based).
    fn output_text(cmd: &Shuffle, n: usize) -> String {
        let mut text = String::new();
        Io::default().new_reader(&cmd.outputs[n - 1]).unwrap().read_to_string(&mut text).unwrap();
        text
    }

    /// Returns the records' headers as strings.
    fn heads(records: &[OwnedRecord]) -> Vec<String> {
        records.iter().map(|r| String::from_utf8(r.head.clone()).unwrap()).collect()
    }

    /// Returns a header's read name, without its comment.
    fn base_name_of(head: &str) -> String {
        head.split(' ').next().unwrap().to_owned()
    }

    /// Splits FASTQ text into its four-line records, sorted.
    fn sorted_records(text: &str) -> Vec<String> {
        let lines: Vec<&str> = text.split_inclusive('\n').collect();
        let mut records: Vec<String> = lines.chunks(4).map(<[&str]>::concat).collect();
        records.sort_unstable();
        records
    }

    /// Encodes a read set with input `index` whose reads are the given FASTQ records.
    fn encoded_read_set(index: u64, reads: &[&str]) -> Vec<u8> {
        let mut out = index.to_le_bytes().to_vec();
        let reads_len: usize = reads.iter().map(|r| 4 + r.len()).sum();
        out.extend_from_slice(&u32::try_from(reads_len).unwrap().to_le_bytes());
        for read in reads {
            out.extend_from_slice(&u32::try_from(read.len()).unwrap().to_le_bytes());
            out.extend_from_slice(read.as_bytes());
        }
        out
    }

    /// Writes `count` single-read read sets into temporary buckets with `seed` and returns the
    /// non-empty buckets, along with the temporary space holding their files.
    fn bucketed_read_sets(dir: &Path, count: u64, seed: u64) -> (Arc<TempSpace>, Vec<Bucket>) {
        let temp = Arc::new(TempSpace::new(&[dir.to_path_buf()]).unwrap());
        let mut router =
            BucketRouter::new(Arc::clone(&temp), Partition::new(0, u128::MAX, 2), 4096, 1, 1);
        for index in 0..count {
            let read_set = encoded_read_set(index, &[&format!("@read{index}\nACGT\n+\nIIII\n")]);
            router.route(read_set_key(&read_set, seed).unwrap(), &read_set).unwrap();
        }
        (temp, router.finish().unwrap())
    }

    /// Makes a named pipe at `path`.
    #[cfg(unix)]
    fn make_fifo(path: &Path) -> PathBuf {
        let status = std::process::Command::new("mkfifo").arg(path).status().unwrap();
        assert!(status.success(), "mkfifo failed for {path:?}");
        path.to_path_buf()
    }

    // ---- end-to-end behaviour ----

    #[test]
    fn shuffle_writes_every_input_record_exactly_once() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 500, "fq");
        let cmd = shuffle_of(dir.path(), "out", &inputs);
        cmd.run().unwrap();

        for (n, suffix) in [(1, "/1"), (2, "/2")] {
            let mut actual: Vec<(String, String)> = output_records(&cmd, n)
                .iter()
                .map(|r| {
                    (
                        String::from_utf8(r.head.clone()).unwrap(),
                        String::from_utf8(r.seq.clone()).unwrap(),
                    )
                })
                .collect();
            let mut expected = reads(500, suffix);
            actual.sort();
            expected.sort();
            assert_eq!(actual, expected);
        }
    }

    #[test]
    fn compression_level_0_writes_bgzf_blocks_without_compression() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 500, "fq");
        let stored = Shuffle { compression_level: 0, ..shuffle_of(dir.path(), "stored", &inputs) };
        let compressed = shuffle_of(dir.path(), "compressed", &inputs);
        stored.run().unwrap();
        compressed.run().unwrap();

        let bytes = fs::read(&stored.outputs[0]).unwrap();
        assert_eq!(bytes[..4], [0x1f, 0x8b, 8, 4], "not a gzip member with extra fields");
        assert_eq!(&bytes[12..14], b"BC", "not a BGZF block");
        assert!(bytes.len() as u64 > fs::metadata(&inputs[0]).unwrap().len(), "was compressed");
        assert_eq!(output_text(&stored, 1), output_text(&compressed, 1));
    }

    #[test]
    fn shuffle_writes_paired_reads_in_the_same_order() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 500, "fq");
        let cmd = shuffle_of(dir.path(), "out", &inputs);
        cmd.run().unwrap();

        let r1: Vec<String> =
            heads(&output_records(&cmd, 1)).iter().map(|h| h.replace("/1", "")).collect();
        let r2: Vec<String> =
            heads(&output_records(&cmd, 2)).iter().map(|h| h.replace("/2", "")).collect();
        assert_eq!(r1, r2);
    }

    #[test]
    fn shuffle_order_for_a_fixed_seed_never_changes() {
        let dir = TempDir::new().unwrap();
        let input = write_fastq(&dir.path().join("in.fq"), &reads(10, ""));
        let cmd = shuffle_of(dir.path(), "out", &[input]);
        cmd.run().unwrap();

        // Pinned literally, so that any change to the hash or to how names are read from
        // headers, which would reorder every user's output, fails here.
        let expected: Vec<String> = GOLDEN_ORDER.iter().map(|&name| name.to_owned()).collect();
        assert_eq!(heads(&output_records(&cmd, 1)), expected);
    }

    #[test]
    fn read_comments_do_not_change_the_order() {
        let dir = TempDir::new().unwrap();
        let plain = write_fastq(&dir.path().join("plain.fq"), &reads(200, ""));
        let commented = write_fastq(&dir.path().join("commented.fq"), &reads(200, " 1:N:0:ACGT"));
        let plain_cmd = shuffle_of(dir.path(), "plain", &[plain]);
        let commented_cmd = shuffle_of(dir.path(), "commented", &[commented]);
        plain_cmd.run().unwrap();
        commented_cmd.run().unwrap();

        let commented_names: Vec<String> = heads(&output_records(&commented_cmd, 1))
            .iter()
            .map(|h| h.trim_end_matches(" 1:N:0:ACGT").to_owned())
            .collect();
        assert_eq!(heads(&output_records(&plain_cmd, 1)), commented_names);
    }

    #[test]
    fn pair_suffixes_do_not_change_the_order() {
        let dir = TempDir::new().unwrap();
        let plain = write_fastq(&dir.path().join("plain.fq"), &reads(200, ""));
        let suffixed = write_fastq(&dir.path().join("suffixed.fq"), &reads(200, "/1"));
        let plain_cmd = shuffle_of(dir.path(), "plain", &[plain]);
        let suffixed_cmd = shuffle_of(dir.path(), "suffixed", &[suffixed]);
        plain_cmd.run().unwrap();
        suffixed_cmd.run().unwrap();

        let suffixed_names: Vec<String> = heads(&output_records(&suffixed_cmd, 1))
            .iter()
            .map(|h| h.trim_end_matches("/1").to_owned())
            .collect();
        assert_eq!(heads(&output_records(&plain_cmd, 1)), suffixed_names);
    }

    #[test]
    fn shuffle_with_the_same_seed_is_reproducible() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 500, "fq");
        let first = shuffle_of(dir.path(), "first", &inputs);
        let second = shuffle_of(dir.path(), "second", &inputs);
        first.run().unwrap();
        second.run().unwrap();

        assert_eq!(output_records(&first, 1), output_records(&second, 1));
        assert_eq!(output_records(&first, 2), output_records(&second, 2));
    }

    #[test]
    fn shuffle_with_different_seeds_gives_different_orders() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 500, "fq");
        let first = shuffle_of(dir.path(), "first", &inputs);
        let second = Shuffle { seed: Some(8), ..shuffle_of(dir.path(), "second", &inputs) };
        first.run().unwrap();
        second.run().unwrap();

        assert_ne!(output_records(&first, 1), output_records(&second, 1));
    }

    #[test]
    fn spilled_shuffle_sized_from_input_files_matches_in_memory_shuffle() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 8000, "fq");
        let in_memory = shuffle_of(dir.path(), "in_memory", &inputs);
        let spilled =
            Shuffle { max_memory: MIN_MAX_MEMORY, ..shuffle_of(dir.path(), "spilled", &inputs) };
        assert_eq!(in_memory.run().unwrap().buckets, 0);
        assert!(spilled.run().unwrap().buckets > 0);

        assert_eq!(output_records(&spilled, 1), output_records(&in_memory, 1));
        assert_eq!(output_records(&spilled, 2), output_records(&in_memory, 2));
    }

    #[test]
    fn spilled_shuffle_of_gzipped_inputs_matches_in_memory_shuffle() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 8000, "fq.gz");
        let in_memory = shuffle_of(dir.path(), "in_memory", &inputs);
        let spilled =
            Shuffle { max_memory: MIN_MAX_MEMORY, ..shuffle_of(dir.path(), "spilled", &inputs) };
        in_memory.run().unwrap();
        assert!(spilled.run().unwrap().buckets > 0);

        assert_eq!(output_records(&spilled, 1), output_records(&in_memory, 1));
        assert_eq!(output_records(&spilled, 2), output_records(&in_memory, 2));
    }

    #[test]
    fn spilled_shuffle_with_set_bucket_count_matches_in_memory_shuffle() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 8000, "fq");
        let in_memory = shuffle_of(dir.path(), "in_memory", &inputs);
        let spilled = spilling_shuffle_of(dir.path(), "spilled", &inputs, 8);
        in_memory.run().unwrap();
        assert!(spilled.run().unwrap().buckets > 0);

        assert_eq!(output_records(&spilled, 1), output_records(&in_memory, 1));
        assert_eq!(output_records(&spilled, 2), output_records(&in_memory, 2));
    }

    #[test]
    fn oversized_buckets_are_resplit_without_changing_the_output() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 8000, "fq");
        let in_memory = shuffle_of(dir.path(), "in_memory", &inputs);
        let resplit = spilling_shuffle_of(dir.path(), "resplit", &inputs, 2);
        in_memory.run().unwrap();
        assert!(resplit.run().unwrap().resplits > 0);

        assert_eq!(output_records(&resplit, 1), output_records(&in_memory, 1));
        assert_eq!(output_records(&resplit, 2), output_records(&in_memory, 2));
    }

    #[test]
    fn reads_sharing_a_name_stay_in_input_order_even_when_resplit() {
        let dir = TempDir::new().unwrap();
        let input_reads: Vec<(String, String)> =
            (0..8000).map(|i| ("dup".to_owned(), bases(i))).collect();
        let input = write_fastq(&dir.path().join("in.fq"), &input_reads);
        let cmd = spilling_shuffle_of(dir.path(), "out", &[input], 2);
        assert!(cmd.run().unwrap().resplits > 0);

        let seqs: Vec<String> = output_records(&cmd, 1)
            .iter()
            .map(|r| String::from_utf8(r.seq.clone()).unwrap())
            .collect();
        let expected: Vec<String> = input_reads.into_iter().map(|(_, seq)| seq).collect();
        assert_eq!(seqs, expected);
    }

    #[test]
    fn reads_sharing_a_name_are_written_together_and_counted() {
        let dir = TempDir::new().unwrap();
        let mut input_reads = reads(100, "");
        input_reads.extend((0..5).map(|i| ("read7".to_owned(), bases(1000 + i))));
        let input = write_fastq(&dir.path().join("in.fq"), &input_reads);
        let cmd = shuffle_of(dir.path(), "out", &[input]);
        assert_eq!(cmd.run().unwrap().shared_names, 5);

        let output = heads(&output_records(&cmd, 1));
        let first = output.iter().position(|h| h == "read7").unwrap();
        assert!(output[first..first + 6].iter().all(|h| h == "read7"), "{output:?}");
    }

    #[test]
    fn reads_with_unique_names_are_not_counted_as_sharing_names() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 8000, "fq");
        let cmd = spilling_shuffle_of(dir.path(), "out", &inputs, 4);
        assert_eq!(cmd.run().unwrap().shared_names, 0);
    }

    #[test]
    fn thread_count_does_not_change_the_output() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 8000, "fq");
        let two = spilling_shuffle_of(dir.path(), "two", &inputs, 2);
        // Enough threads for buckets to be loaded, and re-split, by several threads at once.
        let many = Shuffle { threads: 11, ..spilling_shuffle_of(dir.path(), "many", &inputs, 2) };
        assert!(many.load_threads() > 1);
        two.run().unwrap();
        assert!(many.run().unwrap().resplits > 0);

        assert_eq!(output_records(&two, 1), output_records(&many, 1));
        assert_eq!(output_records(&two, 2), output_records(&many, 2));
    }

    #[test]
    fn temporary_files_are_removed_after_a_shuffle() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 8000, "fq");
        let tmp_dirs = vec![dir.path().join("tmp1"), dir.path().join("tmp2")];
        for tmp_dir in &tmp_dirs {
            fs::create_dir(tmp_dir).unwrap();
        }
        let cmd = Shuffle {
            tmp_dirs: tmp_dirs.clone(),
            ..spilling_shuffle_of(dir.path(), "out", &inputs, 8)
        };
        assert!(cmd.run().unwrap().buckets > 0);

        for tmp_dir in &tmp_dirs {
            assert_eq!(fs::read_dir(tmp_dir).unwrap().count(), 0, "{tmp_dir:?} is not empty");
        }
    }

    #[test]
    fn single_input_is_shuffled() {
        let dir = TempDir::new().unwrap();
        let input = write_fastq(&dir.path().join("in.fq"), &reads(500, ""));
        let cmd = shuffle_of(dir.path(), "out", &[input]);
        cmd.run().unwrap();

        let mut output = heads(&output_records(&cmd, 1));
        let input_heads: Vec<String> = reads(500, "").into_iter().map(|(head, _)| head).collect();
        assert_ne!(output, input_heads);
        output.sort();
        let mut sorted_input = input_heads;
        sorted_input.sort();
        assert_eq!(output, sorted_input);
    }

    #[test]
    fn three_inputs_are_shuffled_in_the_same_order() {
        let dir = TempDir::new().unwrap();
        let inputs: Vec<PathBuf> = ["/1", "/2", "/3"]
            .iter()
            .enumerate()
            .map(|(n, suffix)| {
                let reads: Vec<(String, String)> =
                    reads(500, "").into_iter().map(|(h, s)| (format!("{h} {suffix}"), s)).collect();
                write_fastq(&dir.path().join(format!("in{n}.fq")), &reads)
            })
            .collect();
        let cmd = spilling_shuffle_of(dir.path(), "out", &inputs, 4);
        cmd.run().unwrap();

        let names = |n| -> Vec<String> {
            heads(&output_records(&cmd, n)).iter().map(|h| base_name_of(h)).collect()
        };
        assert_eq!(names(1), names(2));
        assert_eq!(names(1), names(3));
    }

    #[test]
    fn empty_inputs_give_empty_outputs() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 0, "fq");
        let cmd = shuffle_of(dir.path(), "out", &inputs);
        assert_eq!(cmd.run().unwrap().read_sets, 0);

        assert!(output_records(&cmd, 1).is_empty());
        assert!(output_records(&cmd, 2).is_empty());
    }

    #[test]
    fn records_are_written_exactly_as_read() {
        let dir = TempDir::new().unwrap();
        let text = "@a/1 comment one\nACGT\n+a/1 comment one\nIIII\n@b/1\tcomment\nGGCC\n+\n#III\n";
        let input = dir.path().join("in.fq");
        fs::write(&input, text).unwrap();
        let cmd = shuffle_of(dir.path(), "out", &[input]);
        cmd.run().unwrap();

        assert_eq!(sorted_records(&output_text(&cmd, 1)), sorted_records(text));
    }

    #[test]
    fn records_with_crlf_line_endings_are_written_exactly_as_read_when_spilled() {
        let dir = TempDir::new().unwrap();
        let text = (0..12_000).fold(String::new(), |mut text, i| {
            text.push_str(&format!(
                "@read{i} c{i}\r\n{}\r\n+read{i}\r\n{}\r\n",
                bases(i),
                "F".repeat(40)
            ));
            text
        });
        let input = dir.path().join("in.fq");
        fs::write(&input, &text).unwrap();
        let in_memory = shuffle_of(dir.path(), "in_memory", std::slice::from_ref(&input));
        let spilled = spilling_shuffle_of(dir.path(), "spilled", &[input], 4);
        in_memory.run().unwrap();
        spilled.run().unwrap();

        // The order is the same as in memory, so names hash the same way from both passes.
        assert_eq!(output_text(&spilled, 1), output_text(&in_memory, 1));
        assert_eq!(sorted_records(&output_text(&spilled, 1)), sorted_records(&text));
    }

    #[test]
    fn gzipped_inputs_give_the_same_output_as_plain_inputs() {
        let dir = TempDir::new().unwrap();
        let plain = shuffle_of(dir.path(), "plain", &write_pairs(dir.path(), 500, "fq"));
        let gzipped = shuffle_of(dir.path(), "gzipped", &write_pairs(dir.path(), 500, "fq.gz"));
        plain.run().unwrap();
        gzipped.run().unwrap();

        assert_eq!(output_records(&plain, 1), output_records(&gzipped, 1));
        assert_eq!(output_records(&plain, 2), output_records(&gzipped, 2));
    }

    #[cfg(unix)]
    #[test]
    fn fifo_inputs_are_shuffled_like_files() {
        let dir = TempDir::new().unwrap();
        let files = write_pairs(dir.path(), 8000, "fq");
        let fifos: Vec<PathBuf> = ["fifo.R1.fq", "fifo.R2.fq"]
            .iter()
            .map(|name| make_fifo(&dir.path().join(name)))
            .collect();
        let feeders: Vec<_> = files
            .iter()
            .zip(&fifos)
            .map(|(file, fifo)| {
                let (file, fifo) = (file.clone(), fifo.clone());
                thread::spawn(move || fs::write(&fifo, fs::read(&file).unwrap()).unwrap())
            })
            .collect();
        let from_files = Shuffle {
            seed: None,
            max_memory: MIN_MAX_MEMORY,
            ..shuffle_of(dir.path(), "files", &files)
        };
        let from_fifos = Shuffle {
            seed: None,
            max_memory: MIN_MAX_MEMORY,
            ..shuffle_of(dir.path(), "fifos", &fifos)
        };
        assert!(from_fifos.run().unwrap().buckets > 0);
        for feeder in feeders {
            feeder.join().unwrap();
        }
        from_files.run().unwrap();

        assert_eq!(output_records(&from_fifos, 1), output_records(&from_files, 1));
        assert_eq!(output_records(&from_fifos, 2), output_records(&from_files, 2));
    }

    // ---- input errors ----

    #[test]
    fn mismatched_read_names_are_an_error_naming_the_inputs() {
        let dir = TempDir::new().unwrap();
        let r1 = write_fastq(&dir.path().join("in.R1.fq"), &reads(10, "/1"));
        let mut r2_reads = reads(10, "/2");
        r2_reads[5].0 = "other/2".to_owned();
        let r2 = write_fastq(&dir.path().join("in.R2.fq"), &r2_reads);
        let cmd = shuffle_of(dir.path(), "out", &[r1, r2]);

        let error = cmd.run().unwrap_err().to_string();
        assert!(error.contains("Read name mismatch at read 6"), "{error}");
        assert!(error.contains("in.R2.fq"), "{error}");
    }

    #[test]
    fn mismatched_read_names_are_allowed_when_checking_is_disabled() {
        let dir = TempDir::new().unwrap();
        let r1 = write_fastq(&dir.path().join("in.R1.fq"), &reads(10, "/1"));
        let mut r2_reads = reads(10, "/2");
        r2_reads[5].0 = "other/2".to_owned();
        let r2 = write_fastq(&dir.path().join("in.R2.fq"), &r2_reads);
        let cmd = Shuffle {
            disable_read_name_checking: true,
            ..shuffle_of(dir.path(), "out", &[r1, r2])
        };

        assert_eq!(cmd.run().unwrap().read_sets, 10);
    }

    #[test]
    fn inputs_with_different_read_counts_are_an_error_naming_the_short_input() {
        let dir = TempDir::new().unwrap();
        let r1 = write_fastq(&dir.path().join("in.R1.fq"), &reads(10, "/1"));
        let r2 = write_fastq(&dir.path().join("in.R2.fq"), &reads(9, "/2"));
        let cmd = shuffle_of(dir.path(), "out", &[r1, r2]);

        let error = cmd.run().unwrap_err().to_string();
        assert!(error.contains("out of sync"), "{error}");
        assert!(error.contains("in.R2.fq"), "{error}");
    }

    #[test]
    fn a_malformed_record_is_an_error_naming_the_input() {
        let dir = TempDir::new().unwrap();
        let r1 = write_fastq(&dir.path().join("in.R1.fq"), &reads(10, "/1"));
        let r2 = dir.path().join("in.R2.fq");
        fs::write(&r2, "@read0/2\nACGT\n+\nII\n").unwrap();
        let cmd = shuffle_of(dir.path(), "out", &[r1, r2]);

        let error = cmd.run().unwrap_err().to_string();
        assert!(error.contains("in.R2.fq"), "{error}");
    }

    // ---- failures ----

    #[test]
    fn a_failed_shuffle_removes_its_outputs_and_temporary_files() {
        let dir = TempDir::new().unwrap();
        let tmp_dir = dir.path().join("tmp");
        fs::create_dir(&tmp_dir).unwrap();
        let r1 = write_fastq(&dir.path().join("in.R1.fq"), &reads(8000, "/1"));
        let r2 = write_fastq(&dir.path().join("in.R2.fq"), &reads(7999, "/2"));
        let cmd = Shuffle {
            tmp_dirs: vec![tmp_dir.clone()],
            ..spilling_shuffle_of(dir.path(), "out", &[r1, r2], 4)
        };

        assert!(cmd.run().is_err());
        for output in &cmd.outputs {
            assert!(!output.exists(), "{output:?} was left behind");
        }
        assert_eq!(fs::read_dir(&tmp_dir).unwrap().count(), 0);
    }

    #[cfg(unix)]
    #[test]
    fn a_failed_shuffle_leaves_a_symlinked_output_in_place() {
        let dir = TempDir::new().unwrap();
        let r1 = write_fastq(&dir.path().join("in.R1.fq"), &reads(10, "/1"));
        let r2 = write_fastq(&dir.path().join("in.R2.fq"), &reads(9, "/2"));
        let cmd = shuffle_of(dir.path(), "out", &[r1, r2]);
        let target = dir.path().join("target.fq.gz");
        fs::write(&target, b"").unwrap();
        std::os::unix::fs::symlink(&target, &cmd.outputs[0]).unwrap();

        assert!(cmd.run().is_err());
        let link = fs::symlink_metadata(&cmd.outputs[0]).expect("the symlink was removed");
        assert!(link.file_type().is_symlink());
        assert!(!cmd.outputs[1].exists(), "the regular output was left behind");
    }

    #[cfg(unix)]
    #[test]
    fn an_unwritable_tmp_dir_fails_before_creating_outputs() {
        use std::os::unix::fs::PermissionsExt;
        let dir = TempDir::new().unwrap();
        let tmp_dir = dir.path().join("tmp");
        fs::create_dir(&tmp_dir).unwrap();
        fs::set_permissions(&tmp_dir, fs::Permissions::from_mode(0o555)).unwrap();
        if fs::write(tmp_dir.join("probe"), b"").is_ok() {
            return; // Running as root, which can write anywhere.
        }
        let inputs = write_pairs(dir.path(), 10, "fq");
        let cmd =
            Shuffle { tmp_dirs: vec![tmp_dir.clone()], ..shuffle_of(dir.path(), "out", &inputs) };

        let error = cmd.run().unwrap_err();
        fs::set_permissions(&tmp_dir, fs::Permissions::from_mode(0o755)).unwrap();
        assert!(format!("{error:#}").contains("temporary directory"), "{error:#}");
        assert!(cmd.outputs.iter().all(|output| !output.exists()));
    }

    #[test]
    fn an_invalid_gzip_input_fails_before_opening_the_outputs() {
        let dir = TempDir::new().unwrap();
        let input = dir.path().join("in.R1.fq.gz");
        fs::write(&input, "@r\nACGT\n+\nIIII\n").unwrap();
        let cmd = shuffle_of(dir.path(), "out", &[input]);
        // A failed run removes the outputs it created, so only an output that was never opened
        // keeps its earlier contents.
        fs::write(&cmd.outputs[0], b"earlier contents").unwrap();

        let error = cmd.run().unwrap_err();
        assert!(format!("{error:#}").contains("in.R1.fq.gz"), "{error:#}");
        assert_eq!(fs::read(&cmd.outputs[0]).unwrap(), b"earlier contents");
    }

    #[cfg(unix)]
    #[test]
    fn an_output_write_failure_is_an_error_not_a_panic() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 8000, "fq");
        let fifo = make_fifo(&dir.path().join("out.R1.fq.gz"));
        // A reader that takes a few bytes and then goes away, so later writes fail.
        let reader = {
            let fifo = fifo.clone();
            thread::spawn(move || {
                let mut buffer = [0u8; 100];
                File::open(&fifo).unwrap().read_exact(&mut buffer).unwrap();
            })
        };
        let cmd = Shuffle {
            outputs: vec![fifo, dir.path().join("out.R2.fq.gz")],
            ..shuffle_of(dir.path(), "out", &inputs)
        };

        let error = cmd.run().unwrap_err();
        reader.join().unwrap();
        assert!(format!("{error:#}").contains("Failed to write the outputs"), "{error:#}");
    }

    #[test]
    fn loading_a_bucket_with_a_corrupted_frame_is_an_error() {
        let dir = TempDir::new().unwrap();
        let (_temp, mut buckets) = bucketed_read_sets(dir.path(), 100, 7);
        let bucket = buckets.remove(0);
        let mut bytes = fs::read(&bucket.path).unwrap();
        let middle = bytes.len() / 2;
        bytes[middle] ^= 0xFF;
        fs::write(&bucket.path, bytes).unwrap();

        assert!(bucket.load(1, 7).is_err());
    }

    #[test]
    fn loading_a_bucket_missing_read_sets_is_an_error() {
        let dir = TempDir::new().unwrap();
        let (_temp, mut buckets) = bucketed_read_sets(dir.path(), 100, 7);
        let mut bucket = buckets.remove(0);
        bucket.read_sets += 1;

        let error = bucket.load(1, 7).err().unwrap();
        assert!(error.to_string().contains("were written to it"), "{error}");
    }

    #[test]
    fn resplitting_a_bucket_missing_read_sets_is_an_error() {
        let dir = TempDir::new().unwrap();
        let (temp, mut buckets) = bucketed_read_sets(dir.path(), 100, 7);
        let mut bucket = buckets.remove(0);
        bucket.read_sets += 1;
        let loader = BucketLoader {
            temp,
            budget: MemoryBudget { total: MIN_MAX_MEMORY },
            seed: 7,
            zstd_level: 1,
            load_threads: 1,
            resplit_workers: 1,
        };

        let error = loader.resplit(bucket).err().unwrap();
        assert!(error.to_string().contains("is corrupt"), "{error}");
    }

    // ---- seed ----

    #[test]
    fn default_seed_gives_the_same_order_for_plain_and_gzipped_copies() {
        let dir = TempDir::new().unwrap();
        let plain = Shuffle {
            seed: None,
            ..shuffle_of(dir.path(), "plain", &write_pairs(dir.path(), 500, "fq"))
        };
        let gzipped = Shuffle {
            seed: None,
            ..shuffle_of(dir.path(), "gzipped", &write_pairs(dir.path(), 500, "fq.gz"))
        };
        plain.run().unwrap();
        gzipped.run().unwrap();

        assert_eq!(output_records(&plain, 1), output_records(&gzipped, 1));
    }

    #[test]
    fn default_seed_differs_between_datasets_with_the_same_read_names() {
        let dir = TempDir::new().unwrap();
        let first_reads = reads(500, "");
        let second_reads: Vec<(String, String)> = first_reads
            .iter()
            .enumerate()
            .map(|(i, (h, _))| (h.clone(), bases(i + 1000)))
            .collect();
        let first = write_fastq(&dir.path().join("first.fq"), &first_reads);
        let second = write_fastq(&dir.path().join("second.fq"), &second_reads);
        let first_cmd = Shuffle { seed: None, ..shuffle_of(dir.path(), "first", &[first]) };
        let second_cmd = Shuffle { seed: None, ..shuffle_of(dir.path(), "second", &[second]) };
        first_cmd.run().unwrap();
        second_cmd.run().unwrap();

        assert_ne!(heads(&output_records(&first_cmd, 1)), heads(&output_records(&second_cmd, 1)));
    }

    #[test]
    fn default_seed_ignores_comments_separator_lines_and_line_endings() {
        let dir = TempDir::new().unwrap();
        let plain = write_fastq(&dir.path().join("plain.fq"), &reads(500, ""));
        let decorated = dir.path().join("decorated.fq");
        let text = reads(500, "").iter().fold(String::new(), |mut text, (head, seq)| {
            let qual = "I".repeat(seq.len());
            text.push_str(&format!("@{head} 1:N:0:ACGT\r\n{seq}\r\n+{head}\r\n{qual}\r\n"));
            text
        });
        fs::write(&decorated, text).unwrap();
        let plain_cmd = Shuffle { seed: None, ..shuffle_of(dir.path(), "plain", &[plain]) };
        let decorated_cmd =
            Shuffle { seed: None, ..shuffle_of(dir.path(), "decorated", &[decorated]) };
        plain_cmd.run().unwrap();
        decorated_cmd.run().unwrap();

        let names = |cmd: &Shuffle| -> Vec<String> {
            heads(&output_records(cmd, 1)).iter().map(|head| base_name_of(head)).collect()
        };
        assert_eq!(names(&plain_cmd), names(&decorated_cmd));
    }

    // ---- bucket count ----

    #[test]
    fn bucket_count_aims_for_buckets_a_quarter_of_memory() {
        let dir = TempDir::new().unwrap();
        let input = dir.path().join("in.fq");
        File::create(&input).unwrap().set_len(64 * 1024 * 1024).unwrap();
        let cmd = shuffle_of(dir.path(), "out", &[input]);

        // 1 MiB of text took 1 MiB of memory, so all 64 MiB will take 64 MiB: four times the
        // 64 MiB budget's 16 MiB target bucket size.
        assert_eq!(cmd.bucket_count(1024 * 1024, 1000, 1024 * 1024).unwrap(), 4);
    }

    #[cfg(unix)]
    #[test]
    fn bucket_count_uses_a_fixed_default_when_an_input_is_not_a_regular_file() {
        let dir = TempDir::new().unwrap();
        let inputs = vec![PathBuf::from("/dev/null")];
        let cmd = Shuffle { max_memory: 4 << 30, ..shuffle_of(dir.path(), "out", &inputs) };

        assert_eq!(cmd.bucket_count(1024, 1, 1024).unwrap(), STREAMING_BUCKETS);
    }

    #[test]
    fn bucket_count_override_is_used_as_given() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 10, "fq");
        let cmd = Shuffle { buckets: Some(3), ..shuffle_of(dir.path(), "out", &inputs) };

        assert_eq!(cmd.bucket_count(1024, 1, 1024).unwrap(), 3);
    }

    #[cfg(unix)]
    #[test]
    fn bucket_count_extrapolates_from_expected_reads() {
        let dir = TempDir::new().unwrap();
        let inputs = vec![PathBuf::from("/dev/null")];
        let cmd =
            Shuffle { expected_reads: Some(64_000), ..shuffle_of(dir.path(), "out", &inputs) };

        // 1,000 read sets took 1 MiB, so 64,000 will take 64 MiB: four times the 64 MiB budget's
        // 16 MiB target bucket size.
        assert_eq!(cmd.bucket_count(1024 * 1024, 1000, 1024 * 1024).unwrap(), 4);
    }

    #[test]
    fn expected_reads_take_precedence_over_input_sizes() {
        let dir = TempDir::new().unwrap();
        let input = dir.path().join("in.fq");
        File::create(&input).unwrap().set_len(64 * 1024 * 1024).unwrap();
        let cmd =
            Shuffle { expected_reads: Some(128_000), ..shuffle_of(dir.path(), "out", &[input]) };

        // The file's size alone would give 4 buckets; 128,000 read sets of 1 KiB give 8.
        assert_eq!(cmd.bucket_count(1024 * 1024, 1000, 1024 * 1024).unwrap(), 8);
    }

    #[cfg(unix)]
    #[test]
    fn expected_reads_below_the_reads_already_seen_are_ignored() {
        let dir = TempDir::new().unwrap();
        let inputs = vec![PathBuf::from("/dev/null")];
        let too_few =
            Shuffle { expected_reads: Some(10), ..shuffle_of(dir.path(), "out", &inputs) };

        // 1,000 read sets have already taken 48 MiB, which needs 3 buckets of the 64 MiB
        // budget's 16 MiB target size; trusting the hint of 10 would give the minimum of 2.
        assert_eq!(too_few.bucket_count(48 * 1024 * 1024, 1000, 1).unwrap(), 3);
    }

    #[test]
    fn spilled_shuffle_with_expected_reads_matches_in_memory_shuffle() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 8000, "fq");
        let in_memory = shuffle_of(dir.path(), "in_memory", &inputs);
        let hinted = Shuffle {
            max_memory: MIN_MAX_MEMORY,
            expected_reads: Some(8000),
            ..shuffle_of(dir.path(), "hinted", &inputs)
        };
        in_memory.run().unwrap();
        assert!(hinted.run().unwrap().buckets > 0);

        assert_eq!(output_records(&hinted, 1), output_records(&in_memory, 1));
        assert_eq!(output_records(&hinted, 2), output_records(&in_memory, 2));
    }

    // ---- gzip ratio ----

    #[test]
    fn gzip_ratio_measures_compression_of_small_files() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("in.fq.gz");
        let mut encoder = GzEncoder::new(File::create(&path).unwrap(), Compression::default());
        encoder.write_all(&[b'A'; 1024 * 1024]).unwrap();
        encoder.finish().unwrap();

        let ratio = gzip_ratio(&path).unwrap();
        let expected = (1024 * 1024) as f64 / fs::metadata(&path).unwrap().len() as f64;
        assert!((ratio - expected).abs() < 1e-9, "{ratio} != {expected}");
    }

    #[test]
    fn gzip_ratio_measures_only_the_start_of_large_files() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("in.fq.gz");
        // A stored (uncompressed) member longer than the probe, then a member that compresses
        // enormously: the whole file's ratio is about 5, but the probe only sees the first.
        let mut file = File::create(&path).unwrap();
        let mut stored = GzEncoder::new(Vec::new(), Compression::none());
        stored.write_all(&vec![b'A'; 5 * 1024 * 1024]).unwrap();
        file.write_all(&stored.finish().unwrap()).unwrap();
        let mut compressed = GzEncoder::new(Vec::new(), Compression::best());
        compressed.write_all(&vec![b'A'; 20 * 1024 * 1024]).unwrap();
        file.write_all(&compressed.finish().unwrap()).unwrap();

        let ratio = gzip_ratio(&path).unwrap();
        assert!((0.98..=1.0).contains(&ratio), "{ratio}");
    }

    // ---- validation ----

    #[test]
    fn validate_rejects_a_missing_input() {
        let dir = TempDir::new().unwrap();
        let cmd = shuffle_of(dir.path(), "out", &[dir.path().join("missing.fq")]);
        let error = cmd.validate().unwrap_err();
        assert!(error.to_string().contains("does not exist"), "{error}");
    }

    #[test]
    fn validate_rejects_an_input_that_is_a_directory() {
        let dir = TempDir::new().unwrap();
        let cmd = shuffle_of(dir.path(), "out", &[dir.path().to_path_buf()]);
        let error = cmd.validate().unwrap_err();
        assert!(error.to_string().contains("is a directory"), "{error}");
    }

    #[test]
    fn validate_rejects_fewer_than_two_threads() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 1, "fq");
        let cmd = Shuffle { threads: 1, ..shuffle_of(dir.path(), "out", &inputs) };
        let error = cmd.validate().unwrap_err();
        assert!(error.to_string().contains("Threads must be at least 2"), "{error}");
    }

    #[test]
    fn validate_rejects_max_memory_below_one_mebibyte() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 1, "fq");
        let cmd = Shuffle { max_memory: 1024, ..shuffle_of(dir.path(), "out", &inputs) };
        let error = cmd.validate().unwrap_err();
        assert!(error.to_string().contains("Max memory must be at least 1M"), "{error}");
    }

    #[test]
    fn validate_rejects_a_missing_tmp_dir() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 1, "fq");
        let cmd = Shuffle {
            tmp_dirs: vec![dir.path().join("missing")],
            ..shuffle_of(dir.path(), "out", &inputs)
        };
        let error = cmd.validate().unwrap_err();
        assert!(error.to_string().contains("does not exist"), "{error}");
    }

    #[test]
    fn validate_rejects_a_tmp_dir_that_is_a_file() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 1, "fq");
        let cmd =
            Shuffle { tmp_dirs: vec![inputs[0].clone()], ..shuffle_of(dir.path(), "out", &inputs) };
        let error = cmd.validate().unwrap_err();
        assert!(error.to_string().contains("is not a directory"), "{error}");
    }

    #[test]
    fn validate_rejects_a_different_number_of_outputs_and_inputs() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 1, "fq");
        let cmd = Shuffle {
            outputs: vec![dir.path().join("out.fq.gz")],
            ..shuffle_of(dir.path(), "out", &inputs)
        };
        let error = cmd.validate().unwrap_err();
        assert!(error.to_string().contains("got 2 inputs and 1 outputs"), "{error}");
    }

    #[test]
    fn validate_rejects_an_output_that_is_also_an_input() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 1, "fq");
        // A different spelling of the first input's path.
        let same_as_r1 = dir.path().join(".").join("in.R1.fq");
        let outputs = vec![same_as_r1, dir.path().join("out.R2.fq.gz")];
        let cmd = Shuffle { outputs, ..shuffle_of(dir.path(), "out", &inputs) };
        let error = cmd.validate().unwrap_err();
        assert!(error.to_string().contains("is also an input"), "{error}");
    }

    #[test]
    fn validate_rejects_the_same_output_given_twice() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 1, "fq");
        let outputs = vec![dir.path().join("out.fq.gz"), dir.path().join("out.fq.gz")];
        let cmd = Shuffle { outputs, ..shuffle_of(dir.path(), "out", &inputs) };
        let error = cmd.validate().unwrap_err();
        assert!(error.to_string().contains("given more than once"), "{error}");
    }

    #[test]
    fn validate_rejects_an_output_in_a_missing_directory() {
        let dir = TempDir::new().unwrap();
        let inputs = write_pairs(dir.path(), 1, "fq");
        let outputs =
            vec![dir.path().join("missing/out.R1.fq.gz"), dir.path().join("out.R2.fq.gz")];
        let cmd = Shuffle { outputs, ..shuffle_of(dir.path(), "out", &inputs) };
        let error = cmd.validate().unwrap_err();
        assert!(error.to_string().contains("Output parent directory"), "{error}");
    }

    // ---- output names ----

    #[test]
    fn outputs_without_a_gzip_extension_get_a_warning() {
        let dir = TempDir::new().unwrap();
        let cmd = Shuffle {
            outputs: vec![dir.path().join("out.fq")],
            ..shuffle_of(dir.path(), "out", &[dir.path().join("in.fq")])
        };
        let warnings = cmd.output_warnings();
        assert_eq!(warnings.len(), 1);
        assert!(warnings[0].contains("will be BGZF compressed"), "{warnings:?}");
    }

    #[test]
    fn outputs_with_a_gzip_extension_get_no_warning() {
        let dir = TempDir::new().unwrap();
        let cmd = Shuffle {
            outputs: vec![dir.path().join("out.fq.gz"), dir.path().join("out.fq.bgz")],
            ..shuffle_of(dir.path(), "out", &[dir.path().join("a.fq"), dir.path().join("b.fq")])
        };
        assert!(cmd.output_warnings().is_empty());
    }

    #[test]
    fn a_dash_output_gets_a_warning_pointing_at_stdout() {
        let dir = TempDir::new().unwrap();
        let cmd = Shuffle {
            outputs: vec![PathBuf::from("-")],
            ..shuffle_of(dir.path(), "out", &[dir.path().join("in.fq")])
        };
        let warnings = cmd.output_warnings();
        assert_eq!(warnings.len(), 1);
        assert!(warnings[0].contains("/dev/stdout"), "{warnings:?}");
    }

    #[cfg(unix)]
    #[test]
    fn outputs_that_are_devices_get_no_warning() {
        let dir = TempDir::new().unwrap();
        let cmd = Shuffle {
            outputs: vec![PathBuf::from("/dev/null")],
            ..shuffle_of(dir.path(), "out", &[dir.path().join("in.fq")])
        };
        assert!(cmd.output_warnings().is_empty());
    }

    // ---- command line ----

    #[test]
    fn cli_parses_inputs_and_outputs_in_order() {
        let cmd = Shuffle::try_parse_from(["shuffle", "-i", "a", "b", "-o", "c", "d"]).unwrap();
        assert_eq!(cmd.inputs, vec![PathBuf::from("a"), PathBuf::from("b")]);
        assert_eq!(cmd.outputs, vec![PathBuf::from("c"), PathBuf::from("d")]);
    }

    #[test]
    fn cli_defaults_are_as_documented() {
        let cmd = Shuffle::try_parse_from(["shuffle", "-i", "a", "-o", "b"]).unwrap();
        assert_eq!(cmd.seed, None);
        assert_eq!(cmd.threads, 8);
        assert_eq!(cmd.compression_level, 5);
        assert_eq!(cmd.max_memory, 4 << 30);
        assert_eq!(cmd.temp_compression, 1);
        assert!(cmd.tmp_dirs.is_empty());
        assert_eq!(cmd.expected_reads, None);
    }

    #[test]
    fn cli_rejects_a_compression_level_above_12() {
        assert!(Shuffle::try_parse_from(["shuffle", "-i", "a", "-o", "b", "-c", "13"]).is_err());
    }

    #[test]
    fn cli_accepts_compression_level_0() {
        let cmd = Shuffle::try_parse_from(["shuffle", "-i", "a", "-o", "b", "-c", "0"]).unwrap();
        assert_eq!(cmd.compression_level, 0);
    }

    #[test]
    fn cli_rejects_a_temp_compression_outside_1_to_19() {
        let args = ["shuffle", "-i", "a", "-o", "b", "--temp-compression", "20"];
        assert!(Shuffle::try_parse_from(args).is_err());
    }

    #[test]
    fn cli_rejects_fewer_than_two_buckets() {
        assert!(
            Shuffle::try_parse_from(["shuffle", "-i", "a", "-o", "b", "--buckets", "1"]).is_err()
        );
    }

    // ---- partition ----

    #[test]
    fn partition_puts_the_ends_of_the_range_in_the_first_and_last_buckets() {
        let partition = Partition::new(0, u128::MAX, 7);
        assert_eq!(partition.bucket_of(0), 0);
        assert_eq!(partition.bucket_of(u128::MAX), 6);
    }

    #[test]
    fn partition_buckets_never_decrease_as_keys_increase() {
        let partition = Partition::new(1000, 1_000_000, 13);
        let buckets: Vec<usize> =
            (1000..=1_000_000).step_by(997).map(|k| partition.bucket_of(k)).collect();
        assert!(buckets.windows(2).all(|pair| pair[0] <= pair[1]));
        assert_eq!(buckets.first(), Some(&0));
        assert_eq!(partition.bucket_of(1_000_000), 12);
    }

    #[test]
    fn partition_splits_a_narrow_range_evenly() {
        let partition = Partition::new(100, 103, 2);
        let buckets: Vec<usize> = (100..=103).map(|k| partition.bucket_of(k)).collect();
        assert_eq!(buckets, vec![0, 0, 1, 1]);
    }

    // ---- temporary data ----

    #[test]
    fn temporary_files_are_spread_across_tmp_dirs_in_turn() {
        let dir = TempDir::new().unwrap();
        let parents = vec![dir.path().join("a"), dir.path().join("b")];
        for parent in &parents {
            fs::create_dir(parent).unwrap();
        }
        let temp = TempSpace::new(&parents).unwrap();

        let in_parent = |path: PathBuf, parent: &Path| path.starts_with(parent);
        assert!(in_parent(temp.next_path(), &parents[0]));
        assert!(in_parent(temp.next_path(), &parents[1]));
        assert!(in_parent(temp.next_path(), &parents[0]));
    }

    /// Returns an encoded read-set header claiming `reads_len` bytes of reads, followed by
    /// `present` bytes of them.
    fn read_set_with_reads(reads_len: u32, present: usize) -> Vec<u8> {
        let mut encoded = Vec::new();
        encoded.extend_from_slice(&2u64.to_le_bytes());
        encoded.extend_from_slice(&reads_len.to_le_bytes());
        encoded.extend_from_slice(&vec![0; present]);
        encoded
    }

    #[test]
    fn read_set_len_rejects_a_truncated_read_set() {
        assert!(read_set_len(&read_set_with_reads(10, 9)).is_err());
    }

    #[test]
    fn read_set_len_includes_the_header() {
        assert_eq!(read_set_len(&read_set_with_reads(10, 10)).unwrap(), READ_SET_HEADER_BYTES + 10);
    }

    #[test]
    fn read_set_hash_ignores_comments_pair_suffixes_and_carriage_returns() {
        let plain = encoded_read_set(0, &["@read1\nACGT\n+\nIIII\n"]);
        let decorated = encoded_read_set(0, &["@read1/1 1:N:0\r\nACGT\r\n+\r\nIIII\r\n"]);
        assert_eq!(read_set_hash(&plain, 7).unwrap(), read_set_hash(&decorated, 7).unwrap());
    }

    #[test]
    fn frame_reader_rejects_a_file_without_the_bucket_magic() {
        let dir = TempDir::new().unwrap();
        let path = dir.path().join("bucket.fqs");
        fs::write(&path, b"NOPE").unwrap();
        let error = FrameReader::open(&path).err().unwrap();
        assert!(error.to_string().contains("not a shuffle bucket file"), "{error}");
    }

    // ---- memory sizes ----

    #[test]
    fn parse_memory_reads_a_lower_case_unit() {
        assert_eq!(parse_memory("1k"), Ok(1 << 10));
    }

    #[test]
    fn parse_memory_reads_megabytes() {
        assert_eq!(parse_memory("512M"), Ok(512 << 20));
    }

    #[test]
    fn parse_memory_reads_gigabytes() {
        assert_eq!(parse_memory("4G"), Ok(4 << 30));
    }

    #[test]
    fn parse_memory_reads_terabytes() {
        assert_eq!(parse_memory("1T"), Ok(1 << 40));
    }

    #[test]
    fn parse_memory_reads_an_ib_suffix_as_binary() {
        assert_eq!(parse_memory("2GiB"), Ok(2 << 30));
    }

    #[test]
    fn parse_memory_reads_a_b_suffix_as_binary() {
        assert_eq!(parse_memory("3MB"), Ok(3 << 20));
    }

    #[test]
    fn parse_memory_reads_a_plain_byte_count() {
        assert_eq!(parse_memory("1048576"), Ok(1_048_576));
    }

    #[test]
    fn parse_memory_reads_a_byte_count_with_a_b_suffix() {
        assert_eq!(parse_memory("100B"), Ok(100));
    }

    #[test]
    fn parse_memory_reads_fractional_sizes() {
        assert_eq!(parse_memory("1.5G"), Ok(3 << 29));
    }

    #[test]
    fn parse_memory_rejects_an_unknown_unit() {
        assert!(parse_memory("4X").is_err());
    }

    #[test]
    fn parse_memory_rejects_ib_without_a_unit() {
        assert!(parse_memory("1048576ib").is_err());
    }

    #[test]
    fn parse_memory_rejects_a_missing_number() {
        assert!(parse_memory("G").is_err());
    }

    #[test]
    fn parse_memory_rejects_a_negative_size() {
        let error = parse_memory("-1G").unwrap_err();
        assert!(error.contains("negative"), "{error}");
    }

    #[test]
    fn parse_memory_rejects_a_size_too_large_for_memory() {
        assert!(parse_memory("99999999999T").is_err());
    }
}
