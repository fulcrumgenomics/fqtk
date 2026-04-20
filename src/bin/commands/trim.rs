//! The `fqtk trim` subcommand: short-read FASTQ trim + filter. Takes one or two
//! synchronized FASTQ files (single-end or paired-end), runs a fixed-order pipeline of
//! read-structure hard-trim → poly-G → adapter (PE-overlap and/or sequence-based) →
//! poly-X → quality (5' then 3' sliding window) → length / N / mean-quality / low-qual
//! fraction filters, and emits BGZF-compressed output plus optional metrics TSV and
//! fastp-shaped JSON report.
//!
//! # Threading model
//!
//! One or two background reader threads decompress each input FASTQ and feed owned
//! records through bounded channels. The main thread collects batches and sends them
//! to a pool of worker threads; each worker does the full trim+filter+serialize+BGZF-
//! compress pipeline on its batch, then hands the compressed bytes to a single writer
//! thread per output. Records are written in input order (preserved via per-batch
//! [`oneshot`] channels).
//!
//! # SIMD kernels
//!
//! The hot per-base kernels live in this module and use the [`wide`] crate for
//! portable 128/256-bit SIMD: [`observe_stats`] (Q20/Q30/N counting), [`count_mismatches_ci_bounded`]
//! (PE overlap + adapter-library comparison), [`find_polyx_tail_len`] (poly-G/X tail
//! scan), [`trim_quality_sliding_3prime`] (signed-i8 sum trick for cut-tail),
//! [`count_bases_below_q`] (low-quality fraction filter), and [`reverse_complement_acgt_into`]
//! (nibble-indexed LUT + lane reverse).
//!
//! # PE-overlap detection
//!
//! [`detect_pe_overlap`] probes overlap lengths in one of three walk modes
//! ([`OverlapWalkMode::Descending`], [`OverlapWalkMode::Ascending`], or
//! [`OverlapWalkMode::Outward`] around a known-good center) and confirms candidate
//! overlaps with a post-cut adapter-evidence check ([`post_cut_matches_library`]).
//! Each worker tracks running overlap statistics in [`OverlapStats`] and adapts its
//! walk mode to the observed insert-size distribution.
//!
//! # Output format
//!
//! The optional metrics TSV is a single-row [`TrimMetrics`] serialization; the JSON
//! report mirrors fastp's schema ([`FastpJsonReport`]) so MultiQC's existing `fastp`
//! module parses it unchanged.

use crate::commands::command::Command;
use crate::commands::utils::fmt_count;
use anyhow::{Result, anyhow};
use bgzf::{CompressionLevel, Compressor};
use clap::Parser;
use crossbeam_channel::{Receiver, Sender, bounded};
use fgoxide::io::{DelimFile, Io};
use fgoxide::iter::IntoChunkedReadAheadIterator;
use fqtk_lib::IUPAC_MASKS;
use fqtk_lib::adapter_db::{KitAdapter, expand_kit_name};
use log::info;
use read_structure::{ReadStructure, SegmentType};
use seq_io::fastq::OwnedRecord;
use seq_io::fastq::Reader as FastqReader;
use seq_io::fastq::Record;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufRead, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::thread;
use wide::{CmpLt, i8x16, u8x16, u8x32};

/// BufReader / BufWriter capacity used throughout the pipeline. 1 MiB is large enough
/// to amortize syscall cost on BGZF-compressed IO without holding an unreasonable amount
/// of per-reader memory.
const BUFFER_SIZE: usize = 1024 * 1024;

/// Emit a progress log message every N input records processed.
const LOG_EVERY: u64 = 5_000_000;

/// Separator used between multiple UMI segments concatenated into the read name.
/// Matches fgumi's `--extract-umis-from-read-names` parser (which normalizes `+` → `-`).
const UMI_JOIN: u8 = b'-';
/// Separator between the Illumina read-id's 7 canonical fields and the UMI field 8.
const UMI_ID_SEP: u8 = b':';
/// Maximum colon-separated fields allowed in the read-id before we consider the header
/// malformed. Standard Illumina ids have 7; 8 means field 8 is already a UMI and we append.
const MAX_READ_ID_FIELDS: usize = 8;

/// Pair count between walk-mode re-evaluations. Small enough to react within a few
/// seconds on the hot path, large enough to amortize the division.
const INSERT_STATS_UPDATE_INTERVAL: u64 = 1000;
/// Minimum number of detected overlaps before we trust the mean enough to switch modes.
const INSERT_STATS_MIN_DETECTIONS: u64 = 64;

/// Length of the post-cut probe used for adapter-evidence confirmation. One SIMD vector
/// on NEON/SSE2, half a register on AVX2.
const ADAPTER_EVIDENCE_PROBE_LEN: usize = 16;
/// Mismatch budget used when checking a read's post-cut bases against the adapter
/// prefix library. Generous because requiring a *both-mate* match makes false positives
/// vanishingly rare even at this rate.
const ADAPTER_EVIDENCE_MAX_MM: usize = 5;
/// Post-cut bases must be at least this long to bother with the evidence check; if the
/// proposed trim is this size or smaller, accept without checking (the risk from a tiny
/// wrong trim is also small).
const ADAPTER_EVIDENCE_MIN_TRIM: usize = 16;

/// Trim and filter short-read FASTQ files.
///
/// Accepts one or two synchronized FASTQ files (single-end or paired-end). Inputs may be
/// plain, gzip-compressed, or bgzf-compressed (auto-detected). Outputs are always
/// BGZF-compressed.
///
/// Supported operations run in a fixed order:
///
///   1. read-structure hard-trim + UMI extraction
///   2. poly-G 3' trim (on by default; `--trim-polyg 0` to disable)
///   3. adapter trim (PE overlap + user/FASTA/kit adapters)
///   4. poly-X 3' trim (opt-in)
///   5. 5'→3' then 3'→5' sliding-window quality trim (opt-in)
///   6. length filter
///   7. N-base filter (opt-in)
///   8. per-read mean-quality filter (opt-in)
///   9. low-quality-fraction filter (opt-in; runs last)
///
/// In paired-end mode, each mate is evaluated independently. If either mate fails any
/// filter, both mates are dropped.
///
/// # Example
///
/// ```bash
/// # Single-end, no trimming (pass-through)
/// fqtk trim -i r1.fq.gz -o out.fq.gz
///
/// # Paired-end with 5-base hard-trim on R1 and 8-base UMI extraction on both mates
/// fqtk trim -i r1.fq.gz r2.fq.gz -o out_r1.fq.gz out_r2.fq.gz -r 5S+T 8M+T
/// ```
#[derive(Parser, Debug)]
#[command(version)]
#[clap(verbatim_doc_comment)]
pub(crate) struct Trim {
    /// One or two input FASTQ files. A single path is single-end; two paths are paired-end.
    /// Inputs may be plain, gzip, or bgzf (auto-detected).
    #[clap(long, short = 'i', required = true, num_args = 1..=2)]
    inputs: Vec<PathBuf>,

    /// One or two output FASTQ paths. Number of outputs must equal number of inputs.
    /// Outputs are always BGZF-compressed.
    #[clap(long, short = 'o', required = true, num_args = 1..=2)]
    outputs: Vec<PathBuf>,

    /// Number of worker threads. Each worker does the full pipeline (trim + filter +
    /// serialize + BGZF compress) on a batch of records. The reader (main) and writer
    /// threads are extras not counted here; their CPU footprint is minimal. Minimum 1.
    #[clap(long, short = 't', default_value = "4")]
    threads: usize,

    /// BGZF compression level for output files (1-12).
    #[clap(long, short = 'c', default_value = "5")]
    compression_level: usize,

    /// Optional path for the trimming metrics TSV. When unset, only the stdout summary is
    /// emitted at end of run.
    #[clap(long, short = 'm')]
    metrics: Option<PathBuf>,

    /// Optional read-structure(s) for hard-trimming and UMI extraction. Either zero (leave
    /// reads untouched) or exactly one per input. Supported segment kinds:
    ///
    ///   T (template) — kept as the output sequence
    ///   M (molecular barcode) — extracted and appended to the read name as a UMI
    ///   S (skip) — discarded
    ///   B (sample barcode) — error; run `fqtk demux` first (see --discard-unsupported-segments)
    ///   C (cellular barcode) — error; no standard FASTQ convention (see --discard-unsupported-segments)
    ///
    /// When multiple M segments are present (across R1 and R2, or multiple within one read),
    /// their bases are concatenated in read order, joined with `-`, and appended to the
    /// read-id as a colon-delimited 8th field — the format fgumi's
    /// `--extract-umis-from-read-names` parser accepts.
    #[clap(long, short = 'r', num_args = 1..=2)]
    read_structures: Vec<ReadStructure>,

    /// Treat unsupported segment types (B sample barcode, C cellular barcode) in a
    /// read-structure as skip (S) segments instead of erroring out. Useful when copying a
    /// read-structure from a previous `fqtk demux` invocation that already extracted those
    /// segments.
    #[clap(long)]
    discard_unsupported_segments: bool,

    /// 3' adapter sequence(s). One value for single-end, or one or two values for paired-end
    /// (R1, R2). Adapter bases may be ACGT or IUPAC codes (e.g. N matches any read base).
    /// Sequences containing IUPAC codes bypass the SIMD ACGT fast path and fall back to a
    /// scalar matcher — plain ACGT adapters are meaningfully faster.
    /// Multiple adapters may additionally be supplied via `--adapter-fasta` and `--kit`.
    #[clap(long, short = 'a', num_args = 1..=2)]
    adapter_sequence: Vec<String>,

    /// FASTA file of adapter sequences to trim from the 3' end of reads. Each record's
    /// sequence is a candidate adapter; the best match is trimmed. Useful for kits with
    /// many index-specific adapter variants.
    #[clap(long, short = 'f')]
    adapter_fasta: Option<PathBuf>,

    /// Built-in adapter kit preset(s). Repeatable to union multiple kits. Known names:
    /// `truseq`, `nextera`, `small-rna`, `aviti`, `mgi` (alias: `dnbseq`), and `all`.
    #[clap(long, short = 'k', num_args = 1..)]
    kit: Vec<String>,

    /// Enable paired-end overlap-based adapter detection. For short inserts where R1 and
    /// R2 overlap, this is more sensitive than sequence-based matching because it does not
    /// rely on the adapter sequence being known or the adapter's early bases being
    /// error-free.
    #[clap(long)]
    detect_adapter_for_pe: bool,

    /// Minimum overlap length (in bases) required to declare R1 and R2 overlap in PE-
    /// overlap trimming. Only effective with `--detect-adapter-for-pe`.
    #[clap(long, default_value = "30")]
    overlap_min_length: usize,

    /// Maximum fraction of mismatches permitted in the PE-overlap probe window
    /// (0.0 to 1.0). E.g. 0.10 accepts up to one mismatch per 10 compared bases.
    #[clap(long, default_value = "0.10")]
    overlap_max_mismatch_rate: f64,

    /// Upper bound on the probe length used per overlap-length candidate when evaluating
    /// PE overlap. Shorter probes = less work per candidate; long enough that a probe
    /// pass reliably implies a true overlap. Default 64 aligns with common SIMD vector
    /// widths (16/32/64).
    #[clap(long, default_value = "64")]
    overlap_diagnostic_length: usize,

    /// Minimum match length (in bases) required when searching a read's 3' end for an
    /// adapter sequence (from `--adapter-sequence`, `--adapter-fasta`, or `--kit`). Short
    /// matches risk false positives; the default trades a little sensitivity for
    /// specificity.
    #[clap(long, default_value = "8")]
    adapter_min_length: usize,

    /// Maximum fraction of mismatches permitted when matching an adapter sequence against
    /// a read's 3' end (0.0 to 1.0). Default 0.125 = one mismatch per 8 bases, matching
    /// fastp's trim-by-sequence rule.
    #[clap(long, default_value = "0.125")]
    adapter_mismatch_rate: f64,

    /// 3' poly-G trim minimum run length. Always on since the SIMD kernel is essentially
    /// free; pass `--trim-polyg 0` to disable. Default 10.
    #[clap(long, default_value = "10")]
    trim_polyg: usize,

    /// Enable 3' poly-X trimming (A/C/T homopolymer tails, e.g. poly-A from RNA-seq) with
    /// the given minimum run length (default 10). Off unless the flag is provided.
    #[clap(long, num_args = 0..=1, default_missing_value = "10")]
    trim_polyx: Option<usize>,

    /// Enable sliding-window 3' quality trim as WINDOW:QUAL (e.g. `8:20` scans 8-base
    /// windows from the 3' end toward 5' and trims trailing low-quality bases until a
    /// window with mean Phred quality ≥ threshold is found). Off unless the flag is
    /// provided; passing the bare flag uses `8:20`. The window of 8 (vs the common 4)
    /// averages over more bases and is less noise-sensitive.
    #[clap(long, num_args = 0..=1, default_missing_value = "8:20")]
    quality_trim_3p: Option<QualityTrim>,

    /// Enable sliding-window 5' quality trim as WINDOW:QUAL. Mirror of
    /// `--quality-trim-3p`: scans 5'→3' from the start of the read and trims leading
    /// low-quality bases until a window with mean Phred quality ≥ threshold is reached.
    /// Off unless the flag is provided; passing the bare flag uses `8:20`. Equivalent to
    /// fastp's `--cut_front`.
    #[clap(long, num_args = 0..=1, default_missing_value = "8:20")]
    quality_trim_5p: Option<QualityTrim>,

    /// Length filter as MIN[:MAX]. Reads/pairs with post-trim length below MIN (or above
    /// MAX if given) are dropped. In paired-end mode, if either mate fails, both are
    /// dropped.
    #[clap(long = "filter-length", short = 'l', default_value = "15")]
    filter_length: LengthFilter,

    /// Drop reads/pairs whose per-mate count of ambiguous (N) bases exceeds this
    /// threshold. In paired-end mode each mate is checked independently — if either R1 or
    /// R2 exceeds the limit, the pair is dropped. Off by default. Runs after trimming.
    #[clap(long = "filter-max-ns")]
    filter_max_ns: Option<usize>,

    /// Drop reads/pairs whose post-trim mean Phred quality is below this threshold. In
    /// paired-end mode each mate is checked independently — if either fails, the pair is
    /// dropped. Off by default. Runs LAST, after every trimming stage.
    #[clap(long = "filter-mean-qual")]
    filter_mean_qual: Option<u8>,

    /// Drop reads/pairs where the fraction of bases below a Phred-quality threshold
    /// exceeds a maximum allowed fraction, specified as `Q:F` (e.g. `15:0.4` drops reads
    /// where >40% of bases are below Q15). In paired-end mode each mate is checked
    /// independently — if either fails, the pair is dropped. Off by default. Matches
    /// fastp's `--qualified_quality_phred` + `--unqualified_percent_limit`.
    #[clap(long = "filter-low-qual")]
    filter_low_qual: Option<LowQualFilter>,

    /// Optional path for a fastp-shape JSON report. The schema follows fastp's so that
    /// MultiQC's existing `fastp` module parses the output unchanged.
    #[clap(long, short = 'j')]
    json: Option<PathBuf>,

    /// Records per batch handed from the reader to each worker. Worker-side compression
    /// naturally chunks a batch's serialized output at the BGZF 64KB block boundary, so
    /// a value here of a few hundred upward keeps BGZF blocks well-filled. Exposed as a
    /// hidden tuning knob for benchmarking; rarely useful to users.
    #[clap(long, hide = true, default_value = "1024")]
    batch_size: usize,

    /// Hint at the typical insert size of the library (in bp). Seeds the PE-overlap
    /// candidate-walk order: overlaps near the expected insert are tested first, so
    /// a correct overlap is found in fewer probe iterations. Workers also update their
    /// own local estimate from observed overlaps, so the hint only needs to be
    /// roughly right — good enough to beat the default-descending start.
    #[clap(long)]
    expected_insert_size: Option<usize>,
}

impl Trim {
    /// Validates all inputs; aggregates every problem into one error.
    fn validate(&self) -> Result<()> {
        let mut errors: Vec<String> = Vec::new();

        if self.inputs.len() != self.outputs.len() {
            errors.push(format!(
                "Number of outputs ({}) must equal number of inputs ({}).",
                self.outputs.len(),
                self.inputs.len()
            ));
        }

        for path in &self.inputs {
            if !path.exists() {
                errors.push(format!("Input file {path:?} does not exist."));
            }
        }

        for path in &self.outputs {
            if let Some(parent) = path.parent()
                && !parent.as_os_str().is_empty()
                && !parent.exists()
            {
                errors.push(format!(
                    "Output parent directory {parent:?} does not exist (for {path:?})."
                ));
            }
        }

        self.check_no_output_overwrites_input(&mut errors);

        self.check_read_structures(&mut errors);

        self.check_adapter_args(&mut errors);

        self.check_filter_args(&mut errors);

        if self.threads < 1 {
            errors.push(format!("Threads must be at least 1, got {}.", self.threads));
        }

        if !(1..=12).contains(&self.compression_level) {
            errors.push(format!(
                "Compression level must be in 1..=12, got {}.",
                self.compression_level
            ));
        }

        if errors.is_empty() {
            Ok(())
        } else {
            use std::fmt::Write;
            let detail = errors.iter().fold(String::new(), |mut s, e| {
                let _ = writeln!(s, "    - {e}");
                s
            });
            Err(anyhow!("Input validation failed:\n{detail}"))
        }
    }

    /// Validates that the number of read-structures is 0 or matches the input count, and
    /// (when `--discard-unsupported-segments` is not set) that none contain B or C segments.
    fn check_read_structures(&self, errors: &mut Vec<String>) {
        if self.read_structures.is_empty() {
            return;
        }
        if self.read_structures.len() != self.inputs.len() {
            errors.push(format!(
                "Number of read-structures ({}) must be 0 or equal to number of inputs ({}).",
                self.read_structures.len(),
                self.inputs.len()
            ));
            // Skip per-segment checks when the count is wrong — the caller needs to fix the
            // count first; per-input segment errors would be misaligned with input indices.
            return;
        }
        // When the user has opted in to discard unsupported segments, skip segment-kind checks.
        if self.discard_unsupported_segments {
            return;
        }
        for (idx, rs) in self.read_structures.iter().enumerate() {
            for seg in rs.iter() {
                match seg.kind {
                    SegmentType::SampleBarcode => errors.push(format!(
                        "Read-structure for input {} ({rs}) contains a sample barcode (B) segment. \
                         Sample barcodes are not supported by `fqtk trim` — run `fqtk demux` first \
                         to assign reads to samples, or pass `--discard-unsupported-segments` to \
                         treat B as skip.",
                        idx + 1
                    )),
                    SegmentType::CellularBarcode => errors.push(format!(
                        "Read-structure for input {} ({rs}) contains a cellular barcode (C) \
                         segment. Cellular barcodes have no widely-adopted FASTQ convention; \
                         single-cell tools (CellRanger, STARsolo, alevin, kallisto-bustools) \
                         expect CB to remain in the R1 sequence. Leave C as template, or pass \
                         `--discard-unsupported-segments` to treat as skip.",
                        idx + 1
                    )),
                    SegmentType::Template | SegmentType::MolecularBarcode | SegmentType::Skip => {}
                    _ => {}
                }
            }
        }
    }

    /// Validates adapter-related arguments: mismatch rate range, min overlap non-zero, kit
    /// name validity, FASTA file existence, and explicit adapter count matching inputs.
    fn check_adapter_args(&self, errors: &mut Vec<String>) {
        if !(0.0..=1.0).contains(&self.adapter_mismatch_rate) {
            errors.push(format!(
                "--adapter-mismatch-rate must be in 0.0..=1.0, got {}.",
                self.adapter_mismatch_rate
            ));
        }
        if !(0.0..=1.0).contains(&self.overlap_max_mismatch_rate) {
            errors.push(format!(
                "--overlap-max-mismatch-rate must be in 0.0..=1.0, got {}.",
                self.overlap_max_mismatch_rate
            ));
        }
        if self.adapter_min_length == 0 {
            errors.push("--adapter-min-length must be at least 1.".to_string());
        }
        if self.overlap_min_length == 0 {
            errors.push("--overlap-min-length must be at least 1.".to_string());
        }
        if self.overlap_diagnostic_length == 0 {
            errors.push("--overlap-diagnostic-length must be at least 1.".to_string());
        }
        if self.adapter_sequence.len() > self.inputs.len() {
            errors.push(format!(
                "{} --adapter-sequence values supplied but only {} input(s); expected at most one \
                 per input.",
                self.adapter_sequence.len(),
                self.inputs.len()
            ));
        }
        for seq in &self.adapter_sequence {
            if seq.is_empty() {
                errors.push("--adapter-sequence values must not be empty.".to_string());
                continue;
            }
            if let Err(msg) = validate_adapter_bases(seq.as_bytes()) {
                errors.push(format!("--adapter-sequence {seq:?}: {msg}"));
            }
        }
        for kit in &self.kit {
            if expand_kit_name(kit).is_none() {
                errors.push(format!(
                    "--kit {kit:?} is not a recognized preset. Known: truseq, nextera, small-rna, \
                     aviti, mgi (alias dnbseq), all."
                ));
            }
        }
        if let Some(fasta) = &self.adapter_fasta
            && !fasta.exists()
        {
            errors.push(format!("--adapter-fasta path {fasta:?} does not exist."));
        }
        if self.detect_adapter_for_pe && self.inputs.len() < 2 {
            errors
                .push("--detect-adapter-for-pe requires two input files (paired-end).".to_string());
        }
    }

    /// Validates filter-related arguments: poly-X min-run > 0 (poly-G allows 0 to disable).
    fn check_filter_args(&self, errors: &mut Vec<String>) {
        if let Some(n) = self.trim_polyx
            && n == 0
        {
            errors.push("--trim-polyx min run must be >= 1.".to_string());
        }
    }

    /// Appends an error to `errors` for every output (including the metrics file) that
    /// would overwrite one of the inputs. Resolves inputs via `canonicalize` (they must
    /// already exist) and outputs via `std::path::absolute` (they typically don't).
    fn check_no_output_overwrites_input(&self, errors: &mut Vec<String>) {
        let input_abs: Vec<PathBuf> = self
            .inputs
            .iter()
            .filter(|p| p.exists())
            .filter_map(|p| std::fs::canonicalize(p).ok())
            .collect();

        let mut check = |candidate: &Path, label: &str| {
            let abs = match resolve_absolute(candidate) {
                Some(a) => a,
                None => return,
            };
            for input in &input_abs {
                if &abs == input {
                    errors.push(format!(
                        "{label} path {candidate:?} resolves to an input file ({input:?}); \
                         refusing to overwrite."
                    ));
                }
            }
        };

        for out in &self.outputs {
            check(out, "Output");
        }
        if let Some(m) = &self.metrics {
            check(m, "Metrics");
        }
    }

    /// Opens all input FASTQ readers.
    fn open_inputs(&self) -> Result<Vec<FastqReader<Box<dyn BufRead + Send>>>> {
        // fgoxide's helper pool handles gzip decompression on background threads. Size
        // to match the number of inputs (1 or 2 for SE/PE); benchmarks showed more
        // threads here don't help — decompression isn't the bottleneck.
        let fgio = Io::new(self.inputs.len().max(1) as u32, BUFFER_SIZE);
        self.inputs
            .iter()
            .map(|p| {
                fgio.new_reader(p)
                    .map(|r| FastqReader::with_capacity(r, BUFFER_SIZE))
                    .map_err(|e| anyhow!("Failed to open input {p:?}: {e}"))
            })
            .collect()
    }

    /// Emits a terse end-of-run summary via `info!`. Structured to give a quick eyeball
    /// without opening the JSON or TSV.
    fn emit_summary(&self, metrics: &TrimMetrics, elapsed: std::time::Duration) {
        let pass_pct = pct(metrics.reads_out, metrics.reads_in);
        let filt = metrics.reads_filtered_length
            + metrics.reads_filtered_n
            + metrics.reads_filtered_quality
            + metrics.reads_filtered_low_qual;
        let q20_before = metrics.q20_before_r1 + metrics.q20_before_r2;
        let q30_before = metrics.q30_before_r1 + metrics.q30_before_r2;
        let q20_after = metrics.q20_after_r1 + metrics.q20_after_r2;
        let q30_after = metrics.q30_after_r1 + metrics.q30_after_r2;
        let bases_before_qual = metrics.bases_in;
        let bases_after_qual = metrics.total_bases_after_r1 + metrics.total_bases_after_r2;
        let secs = elapsed.as_secs();
        info!("fqtk trim complete:");
        info!(
            "  input:      {} {}",
            fmt_count(metrics.reads_in),
            if self.inputs.len() == 2 { "pairs" } else { "reads" }
        );
        info!(
            "  output:     {} ({pass_pct:.2}%)  filtered: {}  (length {}, n-base {}, quality {}, low-qual {})",
            fmt_count(metrics.reads_out),
            fmt_count(filt),
            fmt_count(metrics.reads_filtered_length),
            fmt_count(metrics.reads_filtered_n),
            fmt_count(metrics.reads_filtered_quality),
            fmt_count(metrics.reads_filtered_low_qual),
        );
        info!(
            "  bases:      {} in / {} out",
            fmt_count(metrics.bases_in),
            fmt_count(metrics.bases_out)
        );
        info!(
            "  trimmed:    read-structure {}  adapter {}  poly-G {}  poly-X {}  quality {}",
            fmt_count(metrics.bases_trimmed_read_structure),
            fmt_count(metrics.bases_trimmed_adapter),
            fmt_count(metrics.bases_trimmed_polyg),
            fmt_count(metrics.bases_trimmed_polyx),
            fmt_count(metrics.bases_trimmed_quality),
        );
        info!(
            "  Q20 rate:   {:.2}% \u{2192} {:.2}%",
            pct(q20_before, bases_before_qual),
            pct(q20_after, bases_after_qual),
        );
        info!(
            "  Q30 rate:   {:.2}% \u{2192} {:.2}%",
            pct(q30_before, bases_before_qual),
            pct(q30_after, bases_after_qual),
        );
        info!("  elapsed:    {}m {}s", secs / 60, secs % 60);
    }
}

impl Command for Trim {
    /// Runs the full trim pipeline end-to-end. Validates CLI inputs, opens readers,
    /// builds the adapter set and overlap-evidence library, spawns reader threads,
    /// spawns the worker pool, spawns per-output writer threads, pumps batches through
    /// the channels in input order, writes the metrics TSV (if requested) and fastp-
    /// JSON report (if requested), and emits a terse stdout summary.
    fn execute(&self) -> Result<()> {
        let start = std::time::Instant::now();
        self.validate()?;

        info!(
            "Trimming {} input file(s) to {} output file(s)",
            self.inputs.len(),
            self.outputs.len()
        );

        let sources = self.open_inputs()?;
        let num_inputs = sources.len();
        let adapters =
            build_adapter_set(&self.adapter_sequence, &self.adapter_fasta, &self.kit, num_inputs)?;
        let batch_size = self.batch_size.max(1);

        // One background thread per input runs gzip decompression + seq_io parsing +
        // RefRecord→OwnedRecord copy, then ships chunks of owned records to the main
        // thread via a bounded channel. This gets decompression (~68% of the reader's
        // on-CPU time by profile) off the main thread while keeping each individual
        // stream's decompression serial (standard .fastq.gz can't be split).
        let read_ahead_chunk = batch_size.min(1024);
        let read_ahead_buffer = 4usize;
        let mut iters: Vec<_> = sources
            .into_iter()
            .map(|reader| {
                OwnedRecordIter { reader }.read_ahead(read_ahead_chunk, read_ahead_buffer)
            })
            .collect();

        // Grab the first batch so it can be handed to the worker pool below. The
        // reader loop picks up where this leaves off.
        let first_batch = fill_batch_from_iters(&mut iters, batch_size, num_inputs, 0)?;

        // Poly-G trimming is always on unless the user explicitly sets `--trim-polyg 0`.
        let polyg_min_run: Option<usize> =
            if self.trim_polyg == 0 { None } else { Some(self.trim_polyg) };

        let overlap_adapter_library =
            build_overlap_adapter_library(&self.adapter_sequence, &self.adapter_fasta, &adapters)?;
        let cfg = PipelineConfig {
            num_inputs,
            read_structures: self.read_structures.clone(),
            discard_unsupported_segments: self.discard_unsupported_segments,
            adapters,
            use_pe_overlap: self.detect_adapter_for_pe && num_inputs == 2,
            overlap_min_length: self.overlap_min_length,
            overlap_max_mismatch_rate: self.overlap_max_mismatch_rate,
            overlap_diagnostic_length: self.overlap_diagnostic_length,
            overlap_adapter_library,
            expected_insert_size: self.expected_insert_size,
            adapter_min_length: self.adapter_min_length,
            adapter_mismatch_rate: self.adapter_mismatch_rate,
            polyg_min_run,
            polyx_min_run: self.trim_polyx,
            quality_trim_3p: self.quality_trim_3p,
            quality_trim_5p: self.quality_trim_5p,
            filter_length: self.filter_length,
            filter_max_ns: self.filter_max_ns,
            filter_mean_qual: self.filter_mean_qual,
            filter_low_qual: self.filter_low_qual,
        };

        // `validate()` constrains compression_level to 1..=12, which fits libdeflate's
        // levels 1..=12 directly.
        let compression_level = CompressionLevel::new(
            u8::try_from(self.compression_level).expect("compression level validated in 1..=12"),
        )?;

        let n_workers = self.threads;
        let outputs = &self.outputs;

        let agg = thread::scope(|s| -> Result<WorkerAggregate> {
            // Channel topology:
            //   `batch`    — zipper (main) → worker pool; WorkPacket per batch.
            //   `order[i]` — zipper → writer[i]; oneshot receivers in submit order so each
            //                writer sees its output stream sequentially.
            // Per-input gzip decompression + seq_io parsing already runs on dedicated
            // background threads spawned by `read_ahead`, so the main thread here is a
            // light-weight zipper that pulls one OwnedRecord per source and builds batches.
            let (batch_tx, batch_rx) = bounded::<WorkPacket>(n_workers * 2);
            let (order_txs, order_rxs): (Vec<_>, Vec<_>) = (0..num_inputs)
                .map(|_| bounded::<oneshot::Receiver<Result<Vec<u8>>>>(n_workers * 4))
                .unzip();

            let mut worker_handles = Vec::with_capacity(n_workers);
            for _ in 0..n_workers {
                let rx = batch_rx.clone();
                let cfg_ref = &cfg;
                worker_handles.push(s.spawn(move || worker_loop(rx, cfg_ref, compression_level)));
            }
            drop(batch_rx); // main no longer holds a receiver

            // One writer per output file so per-file syscalls run in parallel; with PE
            // output that doubles effective write throughput.
            let mut writer_handles = Vec::with_capacity(num_inputs);
            for (idx, order_rx) in order_rxs.into_iter().enumerate() {
                let path = &outputs[idx];
                writer_handles.push(s.spawn(move || writer_loop(order_rx, path)));
            }

            let first_batch_len = first_batch.records.len() as u64;
            let first_submit = submit_batch(first_batch, &batch_tx, &order_txs);

            // Reader loop (runs on this thread). Pulls one OwnedRecord from each
            // read-ahead iterator per slot, builds batches, submits to workers.
            let mut records_read = first_batch_len;
            // Collect every error that surfaces during the run. Channel-closed errors
            // (from `submit_batch`) are often a *symptom* of a worker/writer failure
            // that closed the channel; we keep them all and pick the most-specific one
            // at the end so the user sees the root cause rather than the symptom.
            let mut errors: Vec<anyhow::Error> = Vec::new();
            match first_submit {
                Err(e) => errors.push(e),
                Ok(()) => loop {
                    let batch = match fill_batch_from_iters(
                        &mut iters,
                        batch_size,
                        num_inputs,
                        records_read,
                    ) {
                        Ok(b) => b,
                        Err(e) => {
                            errors.push(e);
                            break;
                        }
                    };
                    if batch.records.is_empty() {
                        break; // clean EOF
                    }
                    records_read += batch.records.len() as u64;
                    if let Err(e) = submit_batch(batch, &batch_tx, &order_txs) {
                        errors.push(e);
                        break;
                    }
                    if records_read.is_multiple_of(LOG_EVERY) {
                        info!(
                            "[fqtk trim] read {} {}",
                            fmt_count(records_read),
                            if num_inputs == 1 { "reads" } else { "pairs" }
                        );
                    }
                },
            }
            // Dropping the senders lets the worker pool and writer threads observe EOF
            // and exit cleanly once their queues drain.
            drop(batch_tx);
            drop(order_txs);

            // Join all spawned threads. `thread::scope` would propagate a panic at the
            // end of the scope, but we explicitly `join()` here so we can (a) collect
            // every error rather than the first panic only, and (b) convert panic
            // payloads into `anyhow::Error` so the caller sees them in context rather
            // than as a chained panic.
            let mut agg = WorkerAggregate::new(num_inputs);
            for h in worker_handles {
                match h.join() {
                    Ok(Ok(partial)) => agg.merge(partial),
                    Ok(Err(e)) => errors.push(e),
                    Err(payload) => {
                        errors.push(anyhow!("worker thread panicked: {}", panic_message(&payload)));
                    }
                }
            }
            for h in writer_handles {
                match h.join() {
                    Ok(Ok(())) => {}
                    Ok(Err(e)) => errors.push(e),
                    Err(payload) => {
                        errors.push(anyhow!("writer thread panicked: {}", panic_message(&payload)));
                    }
                }
            }
            if let Some(e) = select_most_specific_error(errors) {
                return Err(e);
            }
            Ok(agg)
        })?;

        let mut metrics = agg.metrics;
        flatten_mate_stats(&agg.mate_before, &agg.mate_after, &mut metrics);

        if let Some(path) = &self.metrics {
            DelimFile::default().write_tsv(path, std::iter::once(&metrics))?;
        }

        if let Some(path) = &self.json {
            let report = FastpJsonReport::build(self, &metrics, &agg.mate_before, &agg.mate_after);
            write_json_report(path, &report)?;
        }

        self.emit_summary(&metrics, start.elapsed());
        Ok(())
    }
}

/// Parameters for the 3' sliding-window quality trim: window size (in bases) and minimum
/// mean Phred quality within the window. Parsed from a `WINDOW:QUAL` CLI string.
#[derive(Debug, Clone, Copy)]
pub(crate) struct QualityTrim {
    pub window: usize,
    pub threshold: u8,
}

impl FromStr for QualityTrim {
    type Err = String;

    /// Parses `WINDOW:QUAL`; rejects a zero window. Returned errors are user-facing
    /// strings rendered by clap at parse time.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (w, q) = s
            .split_once(':')
            .ok_or_else(|| format!("quality trim must be WINDOW:QUAL (e.g. 4:20), got {s:?}"))?;
        let window: usize = w.parse().map_err(|e| format!("invalid window {w:?}: {e}"))?;
        let threshold: u8 = q.parse().map_err(|e| format!("invalid qual {q:?}: {e}"))?;
        if window == 0 {
            return Err("quality trim WINDOW must be >= 1".into());
        }
        Ok(QualityTrim { window, threshold })
    }
}

/// Length filter spec. Parsed from a `MIN` or `MIN:MAX` CLI string.
#[derive(Debug, Clone, Copy)]
pub(crate) struct LengthFilter {
    pub min: usize,
    pub max: Option<usize>,
}

impl FromStr for LengthFilter {
    type Err = String;

    /// Parses `MIN` or `MIN:MAX`; rejects `MAX < MIN`. Returned errors are user-facing
    /// strings rendered by clap at parse time.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if let Some((min_s, max_s)) = s.split_once(':') {
            let min: usize = min_s.parse().map_err(|e| format!("invalid min {min_s:?}: {e}"))?;
            let max: usize = max_s.parse().map_err(|e| format!("invalid max {max_s:?}: {e}"))?;
            if max < min {
                return Err(format!("length-filter max ({max}) < min ({min})"));
            }
            Ok(LengthFilter { min, max: Some(max) })
        } else {
            let min: usize = s.parse().map_err(|e| format!("invalid min {s:?}: {e}"))?;
            Ok(LengthFilter { min, max: None })
        }
    }
}

/// Low-quality fraction filter spec. Parsed from a `Q:F` CLI string — drop reads whose
/// per-mate fraction of bases below Phred quality `Q` exceeds `F` (0.0..=1.0).
#[derive(Debug, Clone, Copy)]
pub(crate) struct LowQualFilter {
    pub threshold: u8,
    pub max_fraction: f64,
}

impl FromStr for LowQualFilter {
    type Err = String;

    /// Parses `Q:F`; rejects `F` outside `0.0..=1.0`. Returned errors are user-facing
    /// strings rendered by clap at parse time.
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (q, f) = s
            .split_once(':')
            .ok_or_else(|| format!("low-qual filter must be Q:F (e.g. 15:0.4), got {s:?}"))?;
        let threshold: u8 = q.parse().map_err(|e| format!("invalid Q {q:?}: {e}"))?;
        let max_fraction: f64 = f.parse().map_err(|e| format!("invalid F {f:?}: {e}"))?;
        if !(0.0..=1.0).contains(&max_fraction) {
            return Err(format!("low-qual filter F must be in 0.0..=1.0, got {max_fraction}"));
        }
        Ok(LowQualFilter { threshold, max_fraction })
    }
}

/// Flat, serializable summary of one trim run. Written as a single-row TSV.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub(crate) struct TrimMetrics {
    /// Total input reads (single-end) or pairs (paired-end) read from input file(s).
    pub reads_in: u64,
    /// Total output reads/pairs written to output file(s).
    pub reads_out: u64,
    /// Total input bases summed across all input files.
    pub bases_in: u64,
    /// Total output bases summed across all output files.
    pub bases_out: u64,
    /// Bases removed by read-structure hard-trim (S/M/discarded-B/discarded-C segments).
    pub bases_trimmed_read_structure: u64,
    /// Reads/pairs filtered by the length filter. Currently a single bucket for both
    /// min-length and max-length failures; the fastp-JSON `too_short_reads` counter
    /// includes both, and `too_long_reads` is always 0.
    pub reads_filtered_length: u64,
    /// Reads/pairs filtered by the per-read mean-quality filter.
    pub reads_filtered_quality: u64,
    /// Reads/pairs filtered by the N-base-count filter.
    pub reads_filtered_n: u64,
    /// Reads/pairs filtered by the low-quality fraction filter (`--filter-low-qual Q:F`).
    pub reads_filtered_low_qual: u64,
    /// Bases trimmed by adapter matching.
    pub bases_trimmed_adapter: u64,
    /// Bases trimmed by sliding-window quality trim.
    pub bases_trimmed_quality: u64,
    /// Bases trimmed by poly-G 3' trim.
    pub bases_trimmed_polyg: u64,
    /// Bases trimmed by poly-X 3' trim.
    pub bases_trimmed_polyx: u64,
    /// Bases in reads/pairs that were dropped by a filter stage (length, N-base, or
    /// mean-quality). These bases were never written to output but also do not belong to
    /// any `bases_trimmed_*` counter. Included so that
    /// `bases_in == bases_out + sum(bases_trimmed_*) + bases_filtered` holds identically.
    pub bases_filtered: u64,
    /// Count of individual reads (not pairs) that had at least one base trimmed by
    /// adapter matching (PE-overlap or sequence-based). For PE runs where both mates
    /// were adapter-trimmed, this increments by 2 — matching fastp's
    /// `adapter_cutting.adapter_trimmed_reads` semantics so MultiQC's percentages line up.
    pub reads_with_adapter_trimmed: u64,
    /// R1 input Q20 base count.
    pub q20_before_r1: u64,
    /// R1 input Q30 base count.
    pub q30_before_r1: u64,
    /// R1 output Q20 base count.
    pub q20_after_r1: u64,
    /// R1 output Q30 base count.
    pub q30_after_r1: u64,
    /// R1 input total bases.
    pub total_bases_before_r1: u64,
    /// R1 output total bases.
    pub total_bases_after_r1: u64,
    /// R2 input Q20 base count (0 for single-end).
    pub q20_before_r2: u64,
    /// R2 input Q30 base count (0 for single-end).
    pub q30_before_r2: u64,
    /// R2 output Q20 base count (0 for single-end).
    pub q20_after_r2: u64,
    /// R2 output Q30 base count (0 for single-end).
    pub q30_after_r2: u64,
    /// R2 input total bases (0 for single-end).
    pub total_bases_before_r2: u64,
    /// R2 output total bases (0 for single-end).
    pub total_bases_after_r2: u64,
}

impl TrimMetrics {
    /// Adds `other`'s per-pipeline totals into `self`. The per-mate Q20/Q30/total fields
    /// are intentionally NOT merged here — they are populated from the merged per-mate
    /// `MateStats` vectors by [`flatten_mate_stats`] after all workers join.
    fn merge_totals(&mut self, other: &Self) {
        self.reads_in += other.reads_in;
        self.reads_out += other.reads_out;
        self.bases_in += other.bases_in;
        self.bases_out += other.bases_out;
        self.bases_trimmed_read_structure += other.bases_trimmed_read_structure;
        self.reads_filtered_length += other.reads_filtered_length;
        self.reads_filtered_quality += other.reads_filtered_quality;
        self.reads_filtered_n += other.reads_filtered_n;
        self.reads_filtered_low_qual += other.reads_filtered_low_qual;
        self.bases_trimmed_adapter += other.bases_trimmed_adapter;
        self.bases_trimmed_quality += other.bases_trimmed_quality;
        self.bases_trimmed_polyg += other.bases_trimmed_polyg;
        self.bases_trimmed_polyx += other.bases_trimmed_polyx;
        self.bases_filtered += other.bases_filtered;
        self.reads_with_adapter_trimmed += other.reads_with_adapter_trimmed;
    }
}

/// Configuration shared read-only across workers for the duration of one trim run.
/// Derived from the `Trim` CLI struct in `execute()` before workers are spawned.
struct PipelineConfig {
    num_inputs: usize,
    read_structures: Vec<ReadStructure>,
    discard_unsupported_segments: bool,
    adapters: AdapterSet,
    use_pe_overlap: bool,
    overlap_min_length: usize,
    overlap_max_mismatch_rate: f64,
    overlap_diagnostic_length: usize,
    /// 5' prefixes of candidate 3' adapters used by the PE-overlap evidence check.
    /// Split by mate: `r1_prefixes` are the adapters expected to appear past the 3' end
    /// of R1 (e.g. `AGATCGGAAGAGCACA` for TruSeq), `r2_prefixes` past the 3' end of R2
    /// (e.g. `AGATCGGAAGAGCGTC` for TruSeq). For Nextera and similar symmetric chemistries
    /// the two sides hold the same sequence. Assembled from `ALL_KITS`, user-supplied
    /// sequences, and FASTA-loaded adapters. Shared across workers.
    overlap_adapter_library: OverlapAdapterLibrary,
    /// Optional user-supplied insert size hint; seeds each worker's initial walk mode.
    expected_insert_size: Option<usize>,
    adapter_min_length: usize,
    adapter_mismatch_rate: f64,
    polyg_min_run: Option<usize>,
    polyx_min_run: Option<usize>,
    quality_trim_3p: Option<QualityTrim>,
    quality_trim_5p: Option<QualityTrim>,
    filter_length: LengthFilter,
    filter_max_ns: Option<usize>,
    filter_mean_qual: Option<u8>,
    filter_low_qual: Option<LowQualFilter>,
}

/// Per-worker pipeline state: borrowed read-only [`PipelineConfig`] plus every piece of
/// mutable per-worker scratch (the running [`WorkerAggregate`], the adaptive
/// [`OverlapStats`] tracker, scratch `Vec<u8>` buffers for read-structure application and
/// reverse-complement, the UMI-parts accumulator, and the per-output serialization
/// buffers). One `Pipeline` is constructed by each worker thread and driven via `run()`
/// on every record set; the owned buffers preserve capacity across records and across
/// batches.
struct Pipeline<'a> {
    cfg: &'a PipelineConfig,
    agg: WorkerAggregate,
    overlap_stats: OverlapStats,
    rs_seq_scratch: Vec<u8>,
    rs_qual_scratch: Vec<u8>,
    rc_scratch: Vec<u8>,
    umi_parts: Vec<Vec<u8>>,
    serialize_bufs: Vec<Vec<u8>>,
}

impl<'a> Pipeline<'a> {
    /// Allocates the per-worker scratch once up-front. `serialize_bufs` are sized at
    /// `2 * BGZF_BLOCK_SIZE` so the common case (one block's worth of serialized FASTQ)
    /// doesn't grow the buffer.
    fn new(cfg: &'a PipelineConfig) -> Self {
        let num_inputs = cfg.num_inputs;
        Self {
            cfg,
            agg: WorkerAggregate::new(num_inputs),
            overlap_stats: OverlapStats::new(cfg.expected_insert_size, 0),
            rs_seq_scratch: Vec::new(),
            rs_qual_scratch: Vec::new(),
            rc_scratch: Vec::new(),
            umi_parts: Vec::new(),
            serialize_bufs: (0..num_inputs)
                .map(|_| Vec::with_capacity(bgzf::BGZF_BLOCK_SIZE * 2))
                .collect(),
        }
    }

    /// Clears the per-output serialization buffers at the start of each batch so
    /// downstream compression sees only the current batch's bytes. Buffer capacity is
    /// preserved.
    fn reset_batch_bufs(&mut self) {
        for buf in &mut self.serialize_bufs {
            buf.clear();
        }
    }

    /// Runs all trim stages on one record set (one record for SE, a paired R1/R2 for
    /// PE), updates `self.agg`, and appends serialized FASTQ bytes to
    /// `self.serialize_bufs` on success. On filter drop, updates counters but emits no
    /// output bytes.
    fn run(&mut self, records: &mut [OwnedRecord]) -> Result<()> {
        let cfg = self.cfg;
        let num_inputs = cfg.num_inputs;

        // Base counts for each stage so every stage's contribution is tracked
        // independently rather than inferred algebraically.
        let mut stage_bases = sum_seq_bases(records);
        self.agg.metrics.bases_in += stage_bases;
        self.agg.metrics.reads_in += 1;

        // Pre-trim per-mate stats (single SIMD pass over seq+qual).
        for (i, rec) in records.iter().enumerate() {
            let s = observe_stats(&rec.seq, &rec.qual);
            self.agg.mate_before[i].absorb(&s);
        }

        // Stage 1: read-structure hard-trim + UMI extraction.
        if !cfg.read_structures.is_empty() {
            self.umi_parts.clear();
            for (i, rec) in records.iter_mut().enumerate() {
                apply_read_structure(
                    &cfg.read_structures[i],
                    rec,
                    cfg.discard_unsupported_segments,
                    &mut self.umi_parts,
                    &mut self.rs_seq_scratch,
                    &mut self.rs_qual_scratch,
                )?;
            }
            if !self.umi_parts.is_empty() {
                let umi_suffix = join_umi(&self.umi_parts);
                for rec in records.iter_mut() {
                    append_umi_to_head(&mut rec.head, &umi_suffix)?;
                }
            }
        }
        let after_rs = sum_seq_bases(records);
        self.agg.metrics.bases_trimmed_read_structure += stage_bases - after_rs;
        stage_bases = after_rs;

        // Stage 2: poly-G trim (3') — runs before adapter so the sharp G-transition
        // doesn't shift adapter matches past the true 3' end.
        if let Some(min_run) = cfg.polyg_min_run {
            for rec in records.iter_mut() {
                self.agg.metrics.bases_trimmed_polyg += trim_polyx_tail(rec, b'G', min_run);
            }
        }
        let _ = stage_bases;
        stage_bases = sum_seq_bases(records);

        // Stage 3: adapter trimming. First PE-overlap (if enabled), then sequence-based
        // match.
        let mut pre_adapter_lens = [0usize; 2];
        for (i, rec) in records.iter().enumerate() {
            pre_adapter_lens[i] = rec.seq.len();
        }
        let overlap_result = if cfg.use_pe_overlap {
            let result = detect_pe_overlap(
                &records[0].seq,
                &records[1].seq,
                cfg.overlap_min_length,
                cfg.overlap_max_mismatch_rate,
                cfg.overlap_diagnostic_length,
                &cfg.overlap_adapter_library,
                self.overlap_stats.mode,
                &mut self.rc_scratch,
            );
            self.overlap_stats.observe(result, records[0].seq.len());
            result
        } else {
            None
        };
        if let Some(insert_len) = overlap_result {
            if insert_len < records[0].seq.len() {
                records[0].seq.truncate(insert_len);
                records[0].qual.truncate(insert_len);
            }
            if insert_len < records[1].seq.len() {
                records[1].seq.truncate(insert_len);
                records[1].qual.truncate(insert_len);
            }
        }
        if !cfg.adapters.is_empty() {
            for (i, rec) in records.iter_mut().enumerate() {
                let mate_adapters = cfg.adapters.for_mate(i);
                if mate_adapters.is_empty() {
                    continue;
                }
                if let Some(pos) = find_best_adapter_match(
                    &rec.seq,
                    mate_adapters,
                    cfg.adapter_min_length,
                    cfg.adapter_mismatch_rate,
                ) {
                    rec.seq.truncate(pos);
                    rec.qual.truncate(pos);
                }
            }
        }
        let after_adapter = sum_seq_bases(records);
        self.agg.metrics.bases_trimmed_adapter += stage_bases - after_adapter;
        for (i, rec) in records.iter().enumerate() {
            if rec.seq.len() < pre_adapter_lens[i] {
                self.agg.metrics.reads_with_adapter_trimmed += 1;
            }
        }

        // Stage 4: poly-X trim.
        if let Some(min_run) = cfg.polyx_min_run {
            for rec in records.iter_mut() {
                let best =
                    b"ACT".iter().map(|&x| find_polyx_tail_len(&rec.seq, x)).max().unwrap_or(0);
                if best >= min_run && best > 0 {
                    let new_len = rec.seq.len() - best;
                    self.agg.metrics.bases_trimmed_polyx += best as u64;
                    rec.seq.truncate(new_len);
                    rec.qual.truncate(new_len);
                }
            }
        }

        // Stage 5: sliding-window 5' quality trim (cut_front), then 3' (cut_tail).
        // Applied in this order so the 3' scan operates on the post-5'-trim read.
        if let Some(qt) = cfg.quality_trim_5p {
            for rec in records.iter_mut() {
                self.agg.metrics.bases_trimmed_quality +=
                    trim_quality_sliding_5prime(rec, qt.window, qt.threshold);
            }
        }
        if let Some(qt) = cfg.quality_trim_3p {
            for rec in records.iter_mut() {
                self.agg.metrics.bases_trimmed_quality +=
                    trim_quality_sliding_3prime(rec, qt.window, qt.threshold);
            }
        }

        // Post-trim stats: one SIMD pass feeding both the N-filter and the mate_after
        // aggregate.
        let mut post_stats = [BaseStats::default(); 2];
        for (i, rec) in records.iter().enumerate() {
            post_stats[i] = observe_stats(&rec.seq, &rec.qual);
        }

        match evaluate_filters(
            records,
            &post_stats[..num_inputs],
            cfg.filter_length,
            cfg.filter_max_ns,
            cfg.filter_mean_qual,
            cfg.filter_low_qual,
        ) {
            None => {
                for (i, rec) in records.iter().enumerate() {
                    self.agg.metrics.bases_out += post_stats[i].total;
                    self.agg.mate_after[i].absorb(&post_stats[i]);
                    rec.write(&mut self.serialize_bufs[i])
                        .map_err(|e| anyhow!("failed to serialize record: {e}"))?;
                }
                self.agg.metrics.reads_out += 1;
            }
            Some(reason) => {
                self.agg.metrics.bases_filtered += sum_seq_bases(records);
                match reason {
                    FilterReject::Length => self.agg.metrics.reads_filtered_length += 1,
                    FilterReject::NBases => self.agg.metrics.reads_filtered_n += 1,
                    FilterReject::Quality => self.agg.metrics.reads_filtered_quality += 1,
                    FilterReject::LowQual => self.agg.metrics.reads_filtered_low_qual += 1,
                }
            }
        }
        Ok(())
    }
}

/// A chunk of synchronized records sent from the reader to a worker. `records[i]` is the
/// i-th record in the batch; each inner `Vec<OwnedRecord>` holds the paired mates for
/// that record (length 1 for SE, 2 for PE). An empty `records` vector signals EOF.
struct Batch {
    records: Vec<Vec<OwnedRecord>>,
}

/// A batch bundled with one `oneshot::Sender` per output file — workers deliver each
/// mate's compressed bytes directly to the corresponding writer thread, so per-file
/// writes proceed in parallel. One `WorkPacket` per batch.
struct WorkPacket {
    batch: Batch,
    result_txs: Vec<oneshot::Sender<Result<Vec<u8>>>>,
}

/// Per-worker running totals. Merged onto a shared `WorkerAggregate` after all workers
/// join so the main thread sees the same metrics the serial loop used to produce.
#[derive(Debug)]
struct WorkerAggregate {
    metrics: TrimMetrics,
    mate_before: Vec<MateStats>,
    mate_after: Vec<MateStats>,
}

impl WorkerAggregate {
    /// Constructs an empty aggregate with per-mate stats vectors sized for SE (1) or PE (2).
    fn new(num_inputs: usize) -> Self {
        Self {
            metrics: TrimMetrics::default(),
            mate_before: vec![MateStats::default(); num_inputs],
            mate_after: vec![MateStats::default(); num_inputs],
        }
    }

    /// Folds another aggregate (produced by a different worker thread) into this one by
    /// summing both the scalar metrics totals and the per-mate stats.
    fn merge(&mut self, other: WorkerAggregate) {
        self.metrics.merge_totals(&other.metrics);
        for (dst, src) in self.mate_before.iter_mut().zip(other.mate_before.iter()) {
            dst.merge(src);
        }
        for (dst, src) in self.mate_after.iter_mut().zip(other.mate_after.iter()) {
            dst.merge(src);
        }
    }
}

/// Wraps a `FastqReader` as an `Iterator<Item = Result<OwnedRecord>>`. Moving ownership
/// of the reader into the iterator makes it easy to hand off to a background thread
/// (e.g. via fgoxide's [`IntoChunkedReadAheadIterator::read_ahead`]), which is how trim
/// parallelizes decompression + seq_io parsing across inputs.
struct OwnedRecordIter {
    reader: FastqReader<Box<dyn BufRead + Send>>,
}

impl Iterator for OwnedRecordIter {
    type Item = Result<OwnedRecord>;

    /// Advances the underlying [`FastqReader`] and materializes each `RefRecord` as an
    /// `OwnedRecord` so it can cross thread boundaries. Parse errors are wrapped in
    /// `anyhow::Error` with context.
    fn next(&mut self) -> Option<Self::Item> {
        match self.reader.next()? {
            Ok(refrec) => Some(Ok(refrec.to_owned_record())),
            Err(e) => Some(Err(anyhow!("FASTQ read error: {e}"))),
        }
    }
}

/// Per-mate quality aggregate, internal to the execute loop. Flattened into `TrimMetrics`
/// at end-of-run for TSV output and read out of the per-mate vectors directly for the
/// fastp-shape JSON report.
#[derive(Debug, Clone, Copy, Default)]
struct MateStats {
    reads: u64,
    bases: u64,
    q20_bases: u64,
    q30_bases: u64,
}

impl MateStats {
    /// Fold a `BaseStats` produced by [`observe_stats`] into the running aggregate.
    fn absorb(&mut self, s: &BaseStats) {
        self.reads += 1;
        self.bases += s.total;
        self.q20_bases += s.q20;
        self.q30_bases += s.q30;
    }

    /// Add another `MateStats` into this one — used to merge per-worker totals after join.
    fn merge(&mut self, other: &Self) {
        self.reads += other.reads;
        self.bases += other.bases;
        self.q20_bases += other.q20_bases;
        self.q30_bases += other.q30_bases;
    }
}

/// Per-record base statistics produced by a single pass over `(seq, qual)`. Shared by the
/// per-mate aggregate ([`MateStats`]) and the per-record filter check ([`evaluate_filters`]),
/// so each record's bytes are only traversed once post-trim instead of twice (Q-count and
/// then N-count).
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
struct BaseStats {
    /// Total bases (equal to `seq.len() == qual.len()`).
    total: u64,
    /// Count of quality bytes with Phred >= 20.
    q20: u64,
    /// Count of quality bytes with Phred >= 30.
    q30: u64,
    /// Count of sequence bytes equal to `N` or `n`.
    n_bases: u64,
}

/// Root of the fastp-shape JSON report emitted when `--json` is set. Key names mirror
/// fastp's output so that MultiQC's `fastp` module parses fqtk's report unchanged.
///
/// See <https://github.com/OpenGene/fastp/blob/master/src/jsonreporter.cpp> and MultiQC's
/// fastp parser for which fields MultiQC consumes. fqtk does not attempt to populate
/// fastp-only stats we do not compute (per-position quality curves, duplication rate,
/// insert-size histograms, etc.); MultiQC simply skips missing keys.
#[derive(Debug, Serialize)]
struct FastpJsonReport<'a> {
    summary: SummarySection<'a>,
    filtering_result: FilteringResultSection,
    adapter_cutting: AdapterCuttingSection,
    #[serde(skip_serializing_if = "Option::is_none")]
    read1_before_filtering: Option<FilteringStats>,
    #[serde(skip_serializing_if = "Option::is_none")]
    read1_after_filtering: Option<FilteringStats>,
    #[serde(skip_serializing_if = "Option::is_none")]
    read2_before_filtering: Option<FilteringStats>,
    #[serde(skip_serializing_if = "Option::is_none")]
    read2_after_filtering: Option<FilteringStats>,
    command: String,
}

impl FastpJsonReport<'_> {
    /// Assembles the fastp-compatible JSON tree from the aggregated metrics and per-mate
    /// before/after stats. Handles the single-end ↔ paired-end asymmetry and the fastp
    /// convention of reporting filter counts per-mate (so PE filter counts are doubled
    /// relative to our per-pair tallies).
    fn build(cmd: &Trim, metrics: &TrimMetrics, before: &[MateStats], after: &[MateStats]) -> Self {
        let paired = cmd.inputs.len() == 2;
        let sequencing = if paired { "paired end" } else { "single end" };
        // `metrics.reads_*` tally once per pair (or once per read, for SE). fastp's JSON
        // convention is read-counts across both mates, so scale by the mate count for PE.
        let mates = if paired { 2 } else { 1 };

        // `mate_before[i].reads` is the per-mate observation count, which equals the
        // pair count for a PE run (each pair contributes one record to each mate) and the
        // read count for an SE run. fastp's convention is:
        //   * `read{1,2}_before_filtering.total_reads` = per-mate count (= pair count for PE)
        //   * `summary.before_filtering.total_reads`   = sum across mates (= 2 * pairs for PE)
        // `FilteringStats::sum` adds across mates, which reproduces fastp's convention.
        let read1_before = before.first().map(FilteringStats::from_mate);
        let read1_after = after.first().map(FilteringStats::from_mate);
        let read2_before = before.get(1).map(FilteringStats::from_mate);
        let read2_after = after.get(1).map(FilteringStats::from_mate);

        let command = std::env::args().collect::<Vec<_>>().join(" ");

        Self {
            summary: SummarySection {
                fastp_version: env!("CARGO_PKG_VERSION"),
                sequencing,
                before_filtering: FilteringStats::sum(before),
                after_filtering: FilteringStats::sum(after),
            },
            filtering_result: FilteringResultSection {
                passed_filter_reads: metrics.reads_out * mates,
                // fastp lumps mean-qual and unqualified-percent failures into one bucket.
                low_quality_reads: (metrics.reads_filtered_quality
                    + metrics.reads_filtered_low_qual)
                    * mates,
                too_many_n_reads: metrics.reads_filtered_n * mates,
                too_short_reads: metrics.reads_filtered_length * mates,
                too_long_reads: 0, // length filter uses a single "Length" reason bucket
            },
            adapter_cutting: AdapterCuttingSection {
                adapter_trimmed_reads: metrics.reads_with_adapter_trimmed,
                adapter_trimmed_bases: metrics.bases_trimmed_adapter,
            },
            read1_before_filtering: read1_before,
            read1_after_filtering: read1_after,
            read2_before_filtering: read2_before,
            read2_after_filtering: read2_after,
            command,
        }
    }
}

/// JSON schema for fastp's `summary` section. `fastp_version` is actually our own
/// crate version — MultiQC parses the value but doesn't depend on it matching fastp.
#[derive(Debug, Serialize)]
struct SummarySection<'a> {
    fastp_version: &'a str,
    sequencing: &'a str,
    before_filtering: FilteringStats,
    after_filtering: FilteringStats,
}

/// JSON schema for a per-mate or combined read-quality summary. Used both in
/// `summary.{before,after}_filtering` (combined across mates) and
/// `read{1,2}_{before,after}_filtering` (per mate).
#[derive(Debug, Serialize)]
struct FilteringStats {
    total_reads: u64,
    total_bases: u64,
    q20_bases: u64,
    q30_bases: u64,
    q20_rate: f64,
    q30_rate: f64,
}

impl FilteringStats {
    /// Builds the fastp JSON schema from a single mate's aggregate. Q20/Q30 rates are
    /// computed via [`ratio`] (safe against zero-denominator empty-input cases).
    fn from_mate(stats: &MateStats) -> Self {
        Self {
            total_reads: stats.reads,
            total_bases: stats.bases,
            q20_bases: stats.q20_bases,
            q30_bases: stats.q30_bases,
            q20_rate: ratio(stats.q20_bases, stats.bases),
            q30_rate: ratio(stats.q30_bases, stats.bases),
        }
    }

    /// Builds the combined `summary.{before,after}_filtering` entry by summing across
    /// mates (1 for SE, 2 for PE) — matches fastp's convention where the summary row is
    /// a sum and the per-mate rows are per-mate.
    fn sum(mates: &[MateStats]) -> Self {
        let mut agg = MateStats::default();
        for m in mates {
            agg.reads += m.reads;
            agg.bases += m.bases;
            agg.q20_bases += m.q20_bases;
            agg.q30_bases += m.q30_bases;
        }
        Self::from_mate(&agg)
    }
}

/// JSON schema for fastp's `filtering_result` section. Counts are *per-read* (doubled
/// for PE relative to our per-pair internal tallies) to match fastp's convention.
#[derive(Debug, Serialize)]
struct FilteringResultSection {
    passed_filter_reads: u64,
    low_quality_reads: u64,
    too_many_n_reads: u64,
    too_short_reads: u64,
    too_long_reads: u64,
}

/// JSON schema for fastp's `adapter_cutting` section. `adapter_trimmed_reads` is the
/// per-mate count — in PE, a pair where both mates had adapter trim counts as 2.
#[derive(Debug, Serialize)]
struct AdapterCuttingSection {
    adapter_trimmed_reads: u64,
    adapter_trimmed_bases: u64,
}

/// Reasons a read/pair can be dropped by the final filtering stages.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FilterReject {
    Length,
    NBases,
    Quality,
    LowQual,
}

/// Resolved set of adapter sequences for R1 and R2, compiled once up-front from the CLI
/// args (explicit `-a` sequences, FASTA records, and kit presets).
#[derive(Debug, Default, Clone)]
struct AdapterSet {
    r1: Vec<Adapter>,
    r2: Vec<Adapter>,
}

impl AdapterSet {
    /// True when neither mate has any adapters configured.
    fn is_empty(&self) -> bool {
        self.r1.is_empty() && self.r2.is_empty()
    }

    /// Returns the adapter list for mate index `i` (0 = R1, 1 = R2).
    fn for_mate(&self, i: usize) -> &[Adapter] {
        if i == 0 { &self.r1 } else { &self.r2 }
    }
}

/// A single adapter sequence along with a precomputed flag indicating whether it is
/// pure ACGT (no IUPAC ambiguity codes). Pure-ACGT adapters take the SIMD fast path
/// in [`find_adapter_3prime`]; IUPAC-bearing adapters fall back to the scalar matcher
/// that respects IUPAC semantics. Computed once at [`build_adapter_set`] time so the
/// per-read scan doesn't pay for the classification.
#[derive(Debug, Clone)]
struct Adapter {
    bytes: Vec<u8>,
    pure_acgt: bool,
}

impl Adapter {
    /// Wraps a sequence, classifying it as pure-ACGT iff every byte is `A`, `C`, `G`,
    /// or `T` (case-insensitive). `N` and every other IUPAC code count as *not* pure.
    fn new(bytes: Vec<u8>) -> Self {
        let pure_acgt = bytes
            .iter()
            .all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't'));
        Self { bytes, pure_acgt }
    }
}

/// Iteration order used by the PE-overlap search. Selected per-worker, starting from
/// `Descending` (safest default for WGS) and updated as observed overlaps accumulate.
///
/// # Context
///
/// The overlap search tries candidate overlap lengths `L ∈ [min_overlap, read_len]`.
/// Each L represents a hypothesis about the insert length (short-insert case, `I = L`,
/// with `read_len − L` bases of adapter read-through on each end). The *correct* L is
/// the one that aligns `R1[0..L]` with `revcomp(R2[0..L])` within the mismatch budget
/// AND whose post-cut bases look like known adapter sequence.
///
/// The three walk modes iterate the same set of candidates in different orders. The
/// order doesn't change correctness for high-entropy inputs (only one L passes both
/// probe and evidence checks) but it has a big impact on how many probe iterations
/// we execute before that L is found. Picking the right walk for the library's true
/// insert distribution lets us short-circuit cheaply on the common case.
///
/// # Abort logic cheat sheet
///
/// All three walks share the same inner test at each candidate L (`try_overlap`):
///
///   1. **Probe check** — compare `R1[0..probe]` vs `revcomp(R2[0..probe])` over a
///      short 5' prefix (default 64 bp). Skips most of the overlap to avoid vetoing
///      correct alignments on late-read sequencing errors.
///   2. **Evidence check** — if probe passes AND the proposed trim is ≥ 16 bp,
///      compare `R1[L..L+16]` and `R2[L..L+16]` against an adapter-prefix library.
///      Both mates must look adapter-like.
///
/// A key inference: **probe-pass at L implies L_true ≥ L** for non-pathological reads,
/// because `adapter_R1 ≠ revcomp(adapter_R2)` for every standard Illumina / Nextera /
/// MGI / AVITI library prep. A smaller true insert would force disagreement inside
/// the probe window (the post-insert bases on R1 are adapter_R1 and on R2 are
/// adapter_R2 — not reverse-complements of each other), so the probe would fail.
///
/// See each variant for how this inference interacts with walk order to enable (or not)
/// early termination.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum OverlapWalkMode {
    /// Walk candidate overlap lengths from `max_overlap` down to `min_overlap`.
    ///
    /// # When this wins
    ///
    /// Typical inserts `I ≥ R` (standard WGS, most Illumina libraries). The very first
    /// probe tested is L = R, and if the reads come from a long insert the probe fails
    /// quickly. For short-insert cases (I < R), the correct L is found after
    /// `R − L_true` probe iterations — slower than Ascending for very short inserts but
    /// still fine when L_true is in the upper half of the range.
    ///
    /// ```text
    ///                         descending walk — insert shorter than read
    ///                         ↓ tested first
    ///   L:    R=151 … 150 … L_true … 30
    ///   probe: ✗     ✗  …   ✓ ↑    (stop here, probe + evidence accepted)
    /// ```
    ///
    /// # Abort logic (the "stop on evidence fail" short-circuit)
    ///
    /// When we probe-pass at L and evidence fails, we can return `None` immediately.
    /// Reason:
    ///   * The inference `L_true ≥ L` rules out every candidate below L.
    ///   * Every candidate *above* L was already tested by this walk and rejected at
    ///     probe time (descending means we started higher and came down).
    ///   * So neither direction contains a valid L_true — we're done.
    ///
    /// This is the key speedup on tandem-repeat / low-complexity reads where many Ls
    /// probe-pass spuriously: we only pay for one such failure before giving up.
    Descending,

    /// Walk candidate overlap lengths from `min_overlap` up to `max_overlap`.
    ///
    /// # When this wins
    ///
    /// Typical inserts `I ≪ R` — short-fragment libraries (cfDNA, ChIP, FFPE, small-RNA,
    /// degraded DNA). The correct L sits near `min_overlap`, so Ascending finds it in
    /// `L_true − min_overlap` probe iterations instead of `R − L_true`.
    ///
    /// ```text
    ///                         ascending walk — short insert, deep read-through
    ///                         ↓ tested first
    ///   L:    30 … L_true … 150 = R
    ///   probe: ✗  …  ✓  ↑   (stop here, probe + evidence accepted)
    ///                 (L_true small → few iterations before the match)
    /// ```
    ///
    /// # Abort logic (more constrained than Descending)
    ///
    /// When we probe-pass at L and evidence fails, we CANNOT return `None`. Reason:
    ///   * The inference `L_true ≥ L` still holds.
    ///   * But candidates *above* L have NOT been tested yet in this walk direction.
    ///   * A larger L could be the real insert with real adapter evidence downstream.
    ///   * So we continue the scan upward; only return `None` after exhausting every
    ///     L up to `max_overlap`.
    ///
    /// Net: Ascending never short-circuits on evidence failure. It relies purely on
    /// "first probe+evidence pass wins" for early exit, which is cheap when I_true is
    /// close to `min_overlap` and more expensive otherwise.
    Ascending,

    /// Walk candidate overlap lengths outward from a center point (current estimate),
    /// alternating +k / −k at each step, clipped to `[min_overlap, max_overlap]`.
    ///
    /// # When this wins
    ///
    /// Typical inserts land somewhere in the interior of `[min_overlap, R]` — mixed
    /// libraries, short-fragment RNA, anything where `R/2 ≤ mean(I) < 0.9·R`. Outward
    /// finds the correct L in `|L_true − center|` probe iterations, which is small
    /// when the estimate is accurate.
    ///
    /// ```text
    ///                         outward walk from estimate C — insert near C
    ///                        try order: C, C+1, C−1, C+2, C−2, …
    ///                             ↓
    ///   L:    30 …   C−2  C−1   C   C+1  C+2   …  150
    ///   step:        5    3     1   2    4
    ///   probe:       …    …     ✓ ↑  (match near the estimate wins)
    /// ```
    ///
    /// # Abort logic (matches Ascending)
    ///
    /// On evidence failure at L, we CANNOT return `None`. Reason:
    ///   * The inference `L_true ≥ L` holds.
    ///   * But candidates above L on either side of the center may still be untested
    ///     (outward walk doesn't scan in monotone order).
    ///   * So we continue until all candidates in `[min_overlap, max_overlap]` are
    ///     exhausted; then return `None`.
    ///
    /// In short: both non-descending walks sacrifice the stop-on-evidence-fail
    /// short-circuit in exchange for reaching the correct L in fewer iterations when
    /// the insert is where we expect it.
    Outward(usize),
}

/// Per-worker running state used to choose a walk mode. Updated on every pair processed
/// by `detect_pe_overlap`; mode is recomputed every `INSERT_STATS_UPDATE_INTERVAL`
/// pairs from the accumulated observations.
#[derive(Debug, Clone)]
struct OverlapStats {
    /// Sum of detected overlap lengths (contributes only when we actually returned a
    /// trim). Junk reads / no-detections don't inflate this.
    sum_l: u64,
    /// Number of pairs that contributed to `sum_l`.
    count_detect: u64,
    /// Total pairs processed since the last walk-mode recomputation.
    pairs_since_update: u64,
    /// Current walk mode used for dispatch in `detect_pe_overlap`.
    mode: OverlapWalkMode,
}

impl OverlapStats {
    /// Constructs initial stats. If `hint` is supplied, it seeds the walk mode; otherwise
    /// we start in `Descending` (safest default).
    fn new(hint: Option<usize>, read_len: usize) -> Self {
        let mode = match hint {
            // Start-up: no current mode yet, so use the plain (non-hysteretic) classifier.
            Some(h) => Self::pick_mode(h, read_len, OverlapWalkMode::Descending),
            None => OverlapWalkMode::Descending,
        };
        Self { sum_l: 0, count_detect: 0, pairs_since_update: 0, mode }
    }

    /// Classifier. `mean_l` is the running mean overlap length from detected-only pairs
    /// (or the user's hint if we haven't seen enough detections yet). `current` is this
    /// worker's current mode — used to widen the retention band around each threshold
    /// so a distribution whose mean sits near a cutoff doesn't flap between modes every
    /// `INSERT_STATS_UPDATE_INTERVAL` pairs.
    ///
    /// Base thresholds: `mean_l < 50% of read_len` → Ascending; `≥ 90%` → Descending;
    /// otherwise → Outward. Hysteresis adds a 5% of read_len margin to the exit side of
    /// the current mode.
    fn pick_mode(mean_l: usize, read_len: usize, current: OverlapWalkMode) -> OverlapWalkMode {
        let margin = read_len / 20;
        let low = read_len / 2;
        let high = (read_len * 9) / 10;
        match current {
            OverlapWalkMode::Ascending => {
                // Leaving Ascending requires crossing well above the base cutoff.
                if mean_l < low + margin {
                    OverlapWalkMode::Ascending
                } else if mean_l >= high + margin {
                    OverlapWalkMode::Descending
                } else {
                    OverlapWalkMode::Outward(mean_l)
                }
            }
            OverlapWalkMode::Descending => {
                // Leaving Descending requires dropping well below the base cutoff.
                if mean_l >= high.saturating_sub(margin) {
                    OverlapWalkMode::Descending
                } else if mean_l < low.saturating_sub(margin) {
                    OverlapWalkMode::Ascending
                } else {
                    OverlapWalkMode::Outward(mean_l)
                }
            }
            OverlapWalkMode::Outward(_) => {
                // Retain Outward (but refresh the center) unless the mean has moved past
                // either boundary by the full margin.
                if mean_l < low.saturating_sub(margin) {
                    OverlapWalkMode::Ascending
                } else if mean_l >= high + margin {
                    OverlapWalkMode::Descending
                } else {
                    OverlapWalkMode::Outward(mean_l)
                }
            }
        }
    }

    /// Call on every pair. Records an observed overlap (if any), then re-selects the
    /// walk mode every `INSERT_STATS_UPDATE_INTERVAL` pairs.
    fn observe(&mut self, detected: Option<usize>, read_len: usize) {
        if let Some(l) = detected {
            self.sum_l += l as u64;
            self.count_detect += 1;
        }
        self.pairs_since_update += 1;
        if self.pairs_since_update >= INSERT_STATS_UPDATE_INTERVAL {
            self.pairs_since_update = 0;
            // Only update when we've seen enough detections since the last update;
            // otherwise keep the current mode to avoid noise-driven flapping.
            if self.count_detect >= INSERT_STATS_MIN_DETECTIONS {
                let mean_l = (self.sum_l / self.count_detect) as usize;
                self.mode = Self::pick_mode(mean_l, read_len, self.mode);
            }
        }
    }
}

/// Library of adapter 5' prefixes used by the PE-overlap evidence check, split by
/// mate. R1 post-cut bases are checked against `r1_prefixes`; R2 post-cut bases against
/// `r2_prefixes`. Splitting matters because the R1 and R2 adapters of a kit are
/// typically *different* (e.g. TruSeq: `AGATCGGAAGAGCACA` vs `AGATCGGAAGAGCGTC`) —
/// cross-side matching would both cost extra iterations and admit false positives on
/// unrelated sequence. For symmetric chemistries (Nextera) the same prefix appears in
/// both lists.
#[derive(Debug, Clone, Default)]
struct OverlapAdapterLibrary {
    r1_prefixes: Vec<Vec<u8>>,
    r2_prefixes: Vec<Vec<u8>>,
}

impl OverlapAdapterLibrary {
    fn is_empty(&self) -> bool {
        self.r1_prefixes.is_empty() && self.r2_prefixes.is_empty()
    }
}

/// Tests a single overlap length L: probe check, then evidence check if trim is large
/// enough. Returns `Some(MatchOutcome)` describing the decision; `None` when probe
/// fails.
enum ProbeOutcome {
    Accept,
    EvidenceFail,
    ProbeFail,
}

/// Pulls one record from each source iterator to assemble up to `batch_size` paired
/// record sets. Returns a `Batch` with `active` set to the fill count. An empty return
/// (empty `batch.records`) signals a clean EOF on all sources; errors out with the
/// reads-out-of-sync message if one source runs out before another at the same slot.
fn fill_batch_from_iters<I>(
    iters: &mut [I],
    batch_size: usize,
    num_inputs: usize,
    seen_before: u64,
) -> Result<Batch>
where
    I: Iterator<Item = Result<OwnedRecord>>,
{
    let mut records: Vec<Vec<OwnedRecord>> = Vec::with_capacity(batch_size);
    for slot_idx in 0..batch_size {
        let mut mates: Vec<OwnedRecord> = Vec::with_capacity(num_inputs);
        let mut eof_count = 0usize;
        for iter in iters.iter_mut() {
            match iter.next() {
                Some(Ok(rec)) => mates.push(rec),
                Some(Err(e)) => return Err(e),
                None => eof_count += 1,
            }
        }
        if eof_count == num_inputs {
            break;
        }
        anyhow::ensure!(
            mates.len() == num_inputs,
            "FASTQ files are out of sync: {}/{} files produced a record at record {}",
            mates.len(),
            num_inputs,
            seen_before + slot_idx as u64 + 1
        );
        records.push(mates);
    }
    Ok(Batch { records })
}

/// Hand a batch to the worker pool and the corresponding oneshot receiver to the writer
/// in order. Returns Err if either channel's counterpart is closed (workers or writer
/// exited early — surface the underlying error via join).
fn submit_batch(
    batch: Batch,
    batch_tx: &Sender<WorkPacket>,
    order_txs: &[Sender<oneshot::Receiver<Result<Vec<u8>>>>],
) -> Result<()> {
    let mut result_txs = Vec::with_capacity(order_txs.len());
    for order_tx in order_txs {
        let (tx, rx) = oneshot::channel::<Result<Vec<u8>>>();
        // Writers must see their receivers BEFORE we hand the batch to a worker; otherwise
        // a fast worker could deliver output before a writer knows to expect it.
        order_tx.send(rx).map_err(|_| anyhow!("writer exited before receiving batch slot"))?;
        result_txs.push(tx);
    }
    batch_tx
        .send(WorkPacket { batch, result_txs })
        .map_err(|_| anyhow!("workers exited before receiving batch"))?;
    Ok(())
}

/// Worker loop: drain `WorkPacket`s, run the per-record pipeline, compress each mate's
/// serialized output, and deliver the compressed bytes through the packet's oneshot
/// sender. The oneshot is `send`-once, so a worker that panics mid-batch naturally signals
/// the writer via an `Err` from the receiver side.
fn worker_loop(
    batch_rx: Receiver<WorkPacket>,
    cfg: &PipelineConfig,
    compression_level: CompressionLevel,
) -> Result<WorkerAggregate> {
    let num_inputs = cfg.num_inputs;

    // Each worker owns its own libdeflate-backed Compressor per output file; reused across
    // every batch to amortize the cost of creating the libdeflate context.
    let mut compressors: Vec<Compressor> =
        (0..num_inputs).map(|_| Compressor::new(compression_level)).collect();
    let mut pipeline = Pipeline::new(cfg);

    while let Ok(packet) = batch_rx.recv() {
        let WorkPacket { mut batch, result_txs } = packet;
        debug_assert_eq!(result_txs.len(), num_inputs);
        pipeline.reset_batch_bufs();

        let processed: Result<()> = (|| {
            for mates in &mut batch.records {
                pipeline.run(mates)?;
            }
            Ok(())
        })();

        // Compress each mate separately so we can dispatch per-output to the
        // corresponding writer via that mate's oneshot sender.
        let mut per_mate = match processed {
            Ok(()) => match compress_mates(&mut compressors, &pipeline.serialize_bufs) {
                Ok(v) => v.into_iter().map(Ok).collect::<Vec<_>>(),
                Err(e) => {
                    // Build an error per mate so each writer sees the failure through its
                    // own oneshot rather than a dropped sender.
                    let msg = format!("{e}");
                    (0..num_inputs).map(|_| Err(anyhow!("{msg}"))).collect()
                }
            },
            Err(e) => {
                let msg = format!("{e}");
                (0..num_inputs).map(|_| Err(anyhow!("{msg}"))).collect()
            }
        };

        for tx in result_txs.into_iter() {
            let payload = per_mate.remove(0);
            if tx.send(payload).is_err() {
                return Err(anyhow!("writer dropped before worker could deliver batch output"));
            }
        }
        // Batch is dropped here; its allocations (outer Vec + inner OwnedRecords) are
        // freed by the allocator. We no longer reuse batches across iterations — the
        // reader threads freshly allocate OwnedRecords via `to_owned_record()`.
    }
    Ok(pipeline.agg)
}

/// Compresses each mate's serialized output into BGZF blocks, returning one `Vec<u8>` per
/// output file. `Compressor::compress` emits one BGZF block per call and errors if the
/// compressed output wouldn't fit in one block (~64KB), so large batches are chunked at
/// `BGZF_BLOCK_SIZE` byte boundaries; the concatenated blocks form a valid BGZF stream.
fn compress_mates(
    compressors: &mut [Compressor],
    serialize_bufs: &[Vec<u8>],
) -> Result<Vec<Vec<u8>>> {
    let mut out: Vec<Vec<u8>> = Vec::with_capacity(serialize_bufs.len());
    let mut block_buf: Vec<u8> = Vec::with_capacity(bgzf::BGZF_BLOCK_SIZE);
    for (m, ser) in serialize_bufs.iter().enumerate() {
        let mut compressed = Vec::with_capacity(ser.len().max(1024));
        let mut offset = 0;
        while offset < ser.len() {
            let end = (offset + bgzf::BGZF_BLOCK_SIZE).min(ser.len());
            block_buf.clear();
            compressors[m]
                .compress(&ser[offset..end], &mut block_buf)
                .map_err(|e| anyhow!("BGZF compression failed: {e}"))?;
            compressed.extend_from_slice(&block_buf);
            offset = end;
        }
        out.push(compressed);
    }
    Ok(out)
}

/// Writer loop for ONE output file: pull oneshot receivers in reader-submit order,
/// block on each until its worker fills it, and write the compressed bytes. Ordering
/// is enforced structurally by the sequence in which the reader pushed receivers into
/// `order_rx`; running writers per-output lets syscalls for different files proceed
/// in parallel.
fn writer_loop(order_rx: Receiver<oneshot::Receiver<Result<Vec<u8>>>>, path: &Path) -> Result<()> {
    let file = File::create(path).map_err(|e| anyhow!("creating output {path:?}: {e}"))?;
    let mut writer = BufWriter::with_capacity(256 * 1024, file);

    while let Ok(result_rx) = order_rx.recv() {
        let bytes = result_rx
            .recv()
            .map_err(|_| anyhow!("worker dropped without delivering a batch"))??;
        writer.write_all(&bytes)?;
    }

    // BGZF spec requires an empty terminator block at EOF so readers know the stream
    // wasn't truncated.
    let mut eof = Vec::with_capacity(28);
    Compressor::append_eof(&mut eof);
    writer.write_all(&eof)?;
    writer.flush()?;
    Ok(())
}

/// Sum of `rec.seq().len()` across every record in a set, as a `u64`.
fn sum_seq_bases(records: &[OwnedRecord]) -> u64 {
    records.iter().map(|r| r.seq().len() as u64).sum()
}

/// Copies per-mate statistics into the flat `TrimMetrics` fields so the TSV row captures
/// the same numbers the JSON report exposes under `read{1,2}_{before,after}_filtering`.
///
/// Uses `.first()` / `.get(1)` rather than indexing — for single-end runs the slices
/// contain a single entry, and leaving R2 fields at their `Default::default()` zeros is
/// the intentional on-disk convention.
fn flatten_mate_stats(before: &[MateStats], after: &[MateStats], m: &mut TrimMetrics) {
    if let Some(b) = before.first() {
        m.total_bases_before_r1 = b.bases;
        m.q20_before_r1 = b.q20_bases;
        m.q30_before_r1 = b.q30_bases;
    }
    if let Some(a) = after.first() {
        m.total_bases_after_r1 = a.bases;
        m.q20_after_r1 = a.q20_bases;
        m.q30_after_r1 = a.q30_bases;
    }
    if let Some(b) = before.get(1) {
        m.total_bases_before_r2 = b.bases;
        m.q20_before_r2 = b.q20_bases;
        m.q30_before_r2 = b.q30_bases;
    }
    if let Some(a) = after.get(1) {
        m.total_bases_after_r2 = a.bases;
        m.q20_after_r2 = a.q20_bases;
        m.q30_after_r2 = a.q30_bases;
    }
}

/// Compute per-record base statistics in a single SIMD pass over `(seq, qual)`. Each loop
/// iteration processes 16 bytes of quality and 16 bytes of sequence in parallel. A scalar
/// tail handles any bytes past the last full 16-byte chunk.
///
/// Caller must guarantee `seq.len() == qual.len()` (FASTQ invariant); the debug assert
/// catches regressions in development builds.
fn observe_stats(seq: &[u8], qual: &[u8]) -> BaseStats {
    const PHRED33: u8 = 33;
    debug_assert_eq!(seq.len(), qual.len());

    let total = qual.len() as u64;
    let mut q20 = 0u64;
    let mut q30 = 0u64;
    let mut n_bases = 0u64;

    // Threshold vectors are in raw Phred+33 ASCII space (e.g. Q20 == b'5' == 53).
    let q20_thr = u8x16::splat(PHRED33 + 20);
    let q30_thr = u8x16::splat(PHRED33 + 30);
    let n_upper = u8x16::splat(b'N');
    let n_lower = u8x16::splat(b'n');

    let qual_chunks = qual.chunks_exact(16);
    let qual_tail = qual_chunks.remainder();
    let seq_chunks = seq.chunks_exact(16);
    let seq_tail = seq_chunks.remainder();

    for (qchunk, schunk) in qual_chunks.zip(seq_chunks) {
        // `try_into().unwrap()` on a slice of exactly-16 bytes is infallible and the
        // compiler elides the check; `chunks_exact` guarantees the length.
        let qv = u8x16::new(qchunk.try_into().unwrap());
        let sv = u8x16::new(schunk.try_into().unwrap());
        q20 += qv.simd_ge(q20_thr).to_bitmask().count_ones() as u64;
        q30 += qv.simd_ge(q30_thr).to_bitmask().count_ones() as u64;
        let n_mask = sv.simd_eq(n_upper).to_bitmask() | sv.simd_eq(n_lower).to_bitmask();
        n_bases += n_mask.count_ones() as u64;
    }

    for (&q, &s) in qual_tail.iter().zip(seq_tail.iter()) {
        let phred = q.saturating_sub(PHRED33);
        if phred >= 20 {
            q20 += 1;
        }
        if phred >= 30 {
            q30 += 1;
        }
        if s == b'N' || s == b'n' {
            n_bases += 1;
        }
    }

    BaseStats { total, q20, q30, n_bases }
}

/// Case-insensitive bounded mismatch counter: returns the number of positions where
/// `a[i].eq_ignore_ascii_case(&b[i])` is false, short-circuiting the moment the count
/// exceeds `limit`. The returned value is therefore `<= limit` when the slices agree
/// within budget, and `> limit` otherwise; callers use it as a pass/fail predicate.
///
/// Case folding is done via the classic `| 0x20` trick: for ASCII letters, ORing with
/// 0x20 uppercases `[A-Z]` into `[a-z]`, so two letters match iff their `| 0x20`
/// forms are equal. This is equivalent to `eq_ignore_ascii_case` for ASCII letters
/// and for bytes that are identical; byte pairs that differ only in the 0x20 bit are
/// collapsed either way, which matches the existing scalar behavior the callers rely on.
///
/// Caller must pass same-length slices (FASTQ invariant for PE overlap).
fn count_mismatches_ci_bounded(a: &[u8], b: &[u8], limit: usize) -> usize {
    debug_assert_eq!(a.len(), b.len());
    let case_mask = u8x16::splat(0x20);
    let mut count = 0usize;
    let chunks = a.len() / 16;
    for i in 0..chunks {
        let start = i * 16;
        let av = u8x16::new(a[start..start + 16].try_into().unwrap()) | case_mask;
        let bv = u8x16::new(b[start..start + 16].try_into().unwrap()) | case_mask;
        // simd_eq: 0xFF on matching lanes, 0x00 on non-matching.
        // to_bitmask: 16 bits, one per lane. count_ones gives # matches.
        let matches = av.simd_eq(bv).to_bitmask().count_ones() as usize;
        count += 16 - matches;
        if count > limit {
            return count;
        }
    }
    for j in (chunks * 16)..a.len() {
        if !a[j].eq_ignore_ascii_case(&b[j]) {
            count += 1;
            if count > limit {
                return count;
            }
        }
    }
    count
}

/// SIMD reverse-complement kernel specialized for ACGT/N input. Uses a 16-byte
/// nibble-indexed lookup table so each byte's complement is found by its low 4 bits,
/// then reverses the vector lane order via a swizzle with a descending index pattern.
///
/// **Limitation**: IUPAC ambiguity codes (R, Y, S, W, K, M, B, D, H, V) and any
/// non-ACGT/N input byte are mapped to `N` (uppercase or lowercase, following the
/// input's case bit). This is acceptable for callers that feed Illumina-style reads,
/// which contain only ACGT/N; the general-purpose [`reverse_complement_into`] keeps
/// IUPAC codes intact for other callers.
fn reverse_complement_acgt_into(seq: &[u8], out: &mut Vec<u8>) {
    // Low-nibble complement table. Only ACGT and N are mapped; all other low-nibble
    // positions return `N`.
    //   position  hex  ASCII  input base (uppercase)   complement
    //   0         0x0  NUL    (unused)                 N
    //   1         0x1  A      A                        T
    //   2         0x2         (unused; R has nib 2)    N  (IUPAC-lossy)
    //   3         0x3  C      C                        G
    //   4         0x4  T      T                        A
    //   5-6       0x5-6       (unused)                 N
    //   7         0x7  G      G                        C
    //   8-13      0x8-D       (unused)                 N
    //   14        0xE  N      N                        N
    //   15        0xF         (unused)                 N
    const COMP_LUT: [u8; 16] = [
        b'N', b'T', b'N', b'G', b'A', b'N', b'N', b'C', b'N', b'N', b'N', b'N', b'N', b'N', b'N',
        b'N',
    ];
    // Descending swizzle indices — `swizzle_relaxed(self, rhs)` uses `rhs[i]` to pick
    // lane `self[rhs[i]]`; `[15, 14, ..., 0]` reverses a 16-byte vector.
    const REV_IDX: [u8; 16] = [15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0];

    let n = seq.len();
    out.clear();
    out.resize(n, 0u8);

    let lut = u8x16::new(COMP_LUT);
    let rev_idx = u8x16::new(REV_IDX);
    let nibble_mask = u8x16::splat(0x0F);
    let case_mask = u8x16::splat(0x20);
    // Upper-case ASCII letters by clearing the case bit (bit 5). Works for letters only;
    // for non-letters the result feeds only the nibble index, so stray effects on
    // punctuation bytes are irrelevant — they all land in an `N` LUT slot anyway.
    let upper_mask = u8x16::splat(0xDF);

    let chunks = n / 16;
    for i in 0..chunks {
        let start = i * 16;
        let v = u8x16::new(seq[start..start + 16].try_into().unwrap());
        let case_bits = v & case_mask; // 0 for uppercase/non-letter, 0x20 for lowercase
        let upper = v & upper_mask;
        let nibs = upper & nibble_mask;
        // LUT lookup via swizzle: output[j] = lut[nibs[j]].
        let complemented = lut.swizzle_relaxed(nibs);
        // Re-apply the case bit: uppercase LUT output -> lowercase when input was lowercase.
        let cased = complemented | case_bits;
        // Reverse lane order within the 16-byte chunk.
        let reversed = cased.swizzle_relaxed(rev_idx);
        // Write to the mirrored position in output.
        let out_start = n - start - 16;
        out[out_start..out_start + 16].copy_from_slice(reversed.as_array());
    }

    // Scalar tail: the leftover `n % 16` bytes sit at the high end of input but write
    // to the low end of output (positions `0..n - chunks*16`).
    let tail_start = chunks * 16;
    for j in tail_start..n {
        let b = seq[j];
        let case_bit = b & 0x20;
        let upper = b & 0xDF;
        let comp = COMP_LUT[(upper & 0x0F) as usize] | case_bit;
        out[n - 1 - j] = comp;
    }
}

/// Returns the length of the 3' homopolymer run of `x` at the tail of `seq`. Strict
/// match — the run ends at the first non-`x` base. Case-insensitive.
///
/// A tolerant variant (allowing occasional sequencing errors within the run) was
/// considered but creates false-positive over-trimming on sequences like
/// `ACGTGGGGGGGG`, where the algorithm would extend the tail past the real boundary to
/// include non-G bases. Users who want stronger/weaker poly-G tolerance can dial the
/// `--trim-polyg N` min-run-length knob.
/// Returns the length of a homopolymer run of byte `x` at the 3' end of `seq`
/// (case-insensitive). Walks backward from the 3' end in u8x16 chunks: each chunk
/// is case-folded to lowercase via `| 0x20`, compared lane-wise to the target, and
/// the comparison mask's leading-ones (which correspond to the highest-indexed lanes
/// and therefore the 3'-most bytes) tell us how many bytes at the chunk's tail matched.
/// Scans stop at the first non-match; scalar fallback finishes any residual bytes
/// at the 5' end if `seq` is shorter than 16 bytes.
fn find_polyx_tail_len(seq: &[u8], x: u8) -> usize {
    // `| 0x20` maps A-Z → a-z; case-fold both the input and the target this way so the
    // SIMD compare is case-insensitive with no extra ops.
    let target = u8x16::splat(x | 0x20);
    let case = u8x16::splat(0x20);

    let mut count = 0usize;
    let mut pos = seq.len();
    while pos >= 16 {
        let chunk = u8x16::new(seq[pos - 16..pos].try_into().unwrap()) | case;
        // to_bitmask lane-to-bit mapping: lane i sets bit i. The 3'-most byte of this
        // chunk is lane 15 (bit 15); trailing matches therefore occupy the HIGH bits of
        // the 16-bit mask, so `leading_ones` on the u16 cast counts the tail.
        let bits = chunk.simd_eq(target).to_bitmask() as u16;
        let trailing = bits.leading_ones() as usize;
        count += trailing;
        if trailing < 16 {
            return count;
        }
        pos -= 16;
    }
    // Scalar tail at the 5' end. The existing loop condition (`pos` counts down to 0
    // or breaks on first non-match) reproduces the original byte-at-a-time behavior.
    for &b in seq[..pos].iter().rev() {
        if (b | 0x20) == (x | 0x20) {
            count += 1;
        } else {
            break;
        }
    }
    count
}

/// Trims a 3' poly-X tail in place if the tail length meets `min_run`. Returns the number
/// of bases removed.
fn trim_polyx_tail(rec: &mut OwnedRecord, x: u8, min_run: usize) -> u64 {
    let tail = find_polyx_tail_len(&rec.seq, x);
    if tail >= min_run && tail > 0 {
        let new_len = rec.seq.len() - tail;
        rec.seq.truncate(new_len);
        rec.qual.truncate(new_len);
        tail as u64
    } else {
        0
    }
}

/// 3' sliding-window quality trim with **cut-tail semantics**: scans from the 3' end
/// toward the 5' end with a window of `window` bases and trims only the trailing
/// bad-quality bases. When a passing window is found (mean Phred ≥ threshold), we
/// keep every base up to and including that window's 5' edge — drops the last
/// `window − 1` bases of the passing window as a conservative edge guard. This matches
/// fastp's `--cut_tail` algorithm (filter.cpp:166-194).
///
/// Returns the number of bases trimmed from the 3' end (possibly 0).
///
/// # Why cut-tail
///
/// The alternative (cut-right, our previous behavior) scans 5'→3' and trims at the
/// *first* failing window. On reads with a mid-sequence quality dip (microsatellite
/// stretches, dark-cycle artifacts in the middle), cut-right destroys perfectly good
/// 3' bases downstream of the dip. Cut-tail only trims the contiguous low-quality
/// tail, preserving any good bases before it.
///
/// # Signed-sum trick
///
/// The condition `mean(q − 33) < threshold` is equivalent to
/// `sum(q − (threshold + 33)) < 0`. Shifting each byte into signed "distance from
/// threshold" space turns the per-window decision into a sign check — no threshold
/// multiplication and no overflow concerns even at extreme window/threshold combos
/// (saturation in i8 preserves the sign of the sum).
///
/// # SIMD fast path
///
/// For `window ≤ 16` and long enough reads, process 16 candidate windows per chunk by
/// `window` overlapping `i8x16` loads at offsets `0..window-1`. Each load is
/// saturating-subtracted by `splat(threshold + 33)` and the `window` vectors are
/// saturating-added together lane-wise, yielding 16 window-sums in one vector.
/// `simd_lt(0).to_bitmask()` gives a "bad window" mask; the highest unset bit
/// identifies the 3'-most *passing* window and determines the trim position. Chunks
/// advance from the 3' end toward 5'.
///
/// # Fallback
///
/// Very short reads, unusual window sizes (> 16 or combined with thresholds pushing
/// `threshold + 33 > 127`), and the residual 5' tail of the read all go through a
/// scalar running-window scan.
fn trim_quality_sliding_3prime(rec: &mut OwnedRecord, window: usize, threshold: u8) -> u64 {
    const PHRED33: u8 = 33;
    let qual = &rec.qual;
    if qual.len() < window || window == 0 {
        return 0;
    }
    let offset_total = (threshold as u16) + (PHRED33 as u16);
    let max_s = qual.len() - window; // largest valid starting index
    let mut trim_pos = qual.len();

    // `can_simd` requires: window in 1..=16 (we do `window` u8x16 loads per chunk),
    // the signed offset fits in i8, and we have at least one chunk's worth of data
    // (`window + 16` bytes so that `qual[s_start..s_start + 16 + window - 1]` is a
    // valid range when s_start = 0).
    let can_simd = window <= 16 && offset_total <= 127 && qual.len() >= window + 16;

    // `scan_start` is the 3'-most starting position the scalar tail should handle.
    // If SIMD runs, it's the largest `s_start - 1` of an all-bad chunk; else it's
    // `max_s` (the very 3' end).
    let mut scan_start: usize = max_s;
    let mut done = false;

    if can_simd {
        let offset_vec = i8x16::splat(offset_total as i8);
        let zero_vec = i8x16::splat(0);

        let mut s_end = max_s; // inclusive 3'-most s in the current chunk
        while s_end >= 15 {
            let s_start = s_end - 15;

            // Compute 16 window-sums. For window w, lane k of `sums` is the sum of
            // (qual[s_start + k + i] - (threshold + 33)) for i in 0..w — i.e., the
            // signed-shifted window sum for start position s_start + k.
            let mut sums = i8x16::splat(0);
            for i in 0..window {
                let bytes: [u8; 16] = qual[s_start + i..s_start + i + 16].try_into().unwrap();
                let v = i8x16::new(bytemuck::cast(bytes)).saturating_sub(offset_vec);
                sums = sums.saturating_add(v);
            }

            // Lane k in `bad_bits` is 1 if window at s_start+k is below threshold.
            let bad_bits = sums.simd_lt(zero_vec).to_bitmask() as u16;

            if bad_bits == 0xFFFF {
                // Every window in this chunk is bad; record the 5'-most bad s and
                // continue 5'-ward.
                trim_pos = s_start;
                if s_start == 0 {
                    done = true;
                    break;
                }
                s_end = s_start - 1;
            } else {
                // At least one window passes. Find the *3'-most* passing lane — that's
                // the highest set bit of the inverse mask. Lane p = 15 means the
                // 3'-most window in the chunk already passes, so no bad window was
                // seen in this chunk and `trim_pos` stays at whatever earlier chunks
                // (or the initial `qual.len()`) set.
                let good_bits = !bad_bits;
                let p = 15 - good_bits.leading_zeros() as usize;
                if p < 15 {
                    // Bad windows live at lanes p+1..15 (3' of the passing one); the
                    // 5'-most of those is at s_start + p + 1 — that's the last bad s.
                    trim_pos = s_start + p + 1;
                }
                done = true;
                break;
            }
        }
        scan_start = if done { 0 } else { s_end };
    }

    if !done {
        // Scalar running-window tail: scan s from `scan_start` toward 5'. Uses a
        // sum-maintenance pattern — one add and one subtract per step — to keep the
        // scan O(n) rather than the naïve O(n · window) of the original impl.
        if scan_start + window <= qual.len() {
            let win = window as u32;
            let threshold_total = u32::from(threshold) * win;
            let mut s = scan_start;
            let mut sum: u32 =
                qual[s..s + window].iter().map(|&q| u32::from(q.saturating_sub(PHRED33))).sum();
            if sum < threshold_total {
                trim_pos = s;
                while s > 0 {
                    // Slide the window one base toward 5': add the new 5'-most byte,
                    // subtract the old 3'-most byte.
                    sum += u32::from(qual[s - 1].saturating_sub(PHRED33));
                    sum -= u32::from(qual[s - 1 + window].saturating_sub(PHRED33));
                    s -= 1;
                    if sum < threshold_total {
                        trim_pos = s;
                    } else {
                        break;
                    }
                }
            }
        }
    }

    let trimmed = (rec.seq.len() - trim_pos) as u64;
    rec.seq.truncate(trim_pos);
    rec.qual.truncate(trim_pos);
    trimmed
}

/// 5' sliding-window quality trim with **cut-front semantics**: scans from the 5' end
/// toward the 3' end with a window of `window` bases and trims leading bases up to the
/// first passing window (mean Phred ≥ threshold). Matches fastp's `--cut_front`
/// algorithm (filter.cpp:130-163). Symmetric to [`trim_quality_sliding_3prime`].
///
/// Returns the number of bases trimmed from the 5' end (possibly 0).
///
/// Uses a scalar running-window scan: O(n) add/subtract per step. A SIMD fast path is
/// possible but 5' trim is a much rarer operation than 3' trim in practice (2 of 17
/// nf-core fastp-using pipelines hardcode it), so scalar is fine for now.
fn trim_quality_sliding_5prime(rec: &mut OwnedRecord, window: usize, threshold: u8) -> u64 {
    const PHRED33: u8 = 33;
    let qual = &rec.qual;
    if qual.len() < window || window == 0 {
        return 0;
    }
    let win = window as u32;
    let threshold_total = u32::from(threshold) * win;
    let max_s = qual.len() - window;

    let mut sum: u32 = qual[..window].iter().map(|&q| u32::from(q.saturating_sub(PHRED33))).sum();
    // First passing s wins; if none, trim the entire read.
    let mut pass_s: Option<usize> = None;
    if sum >= threshold_total {
        pass_s = Some(0);
    } else {
        for s in 1..=max_s {
            // Slide 5'→3': subtract the base leaving the 5' edge, add the base entering
            // the 3' edge.
            sum -= u32::from(qual[s - 1].saturating_sub(PHRED33));
            sum += u32::from(qual[s + window - 1].saturating_sub(PHRED33));
            if sum >= threshold_total {
                pass_s = Some(s);
                break;
            }
        }
    }

    let trim_len = pass_s.unwrap_or(qual.len());
    if trim_len == 0 {
        return 0;
    }
    // `copy_within + truncate` is equivalent to `drain(..trim_len)` but skips Drain's
    // iterator machinery — the underlying memmove cost is the same but the surrounding
    // bookkeeping is nil. 5' trim is rare, but on libraries where it fires frequently
    // (smRNA-seq 5' adapter rewind, NovaSeq dark-cycle stretches) the constant factor
    // adds up.
    let new_len = rec.seq.len() - trim_len;
    rec.seq.copy_within(trim_len.., 0);
    rec.seq.truncate(new_len);
    rec.qual.copy_within(trim_len.., 0);
    rec.qual.truncate(new_len);
    trim_len as u64
}

/// SIMD count of Phred qualities strictly below `threshold` in a 33-offset FASTQ quality
/// line. Matches the scalar semantics `q.saturating_sub(33) < threshold`, i.e.
/// `q_byte < threshold + 33` (with saturating-subtract so bytes 0..32 count as Phred 0).
///
/// 32-byte chunks via `u8x32::simd_lt`; scalar tail for the last `< 32` bytes.
fn count_bases_below_q(qual: &[u8], threshold: u8) -> u64 {
    const PHRED33: u8 = 33;
    // Using `simd_lt` against a splat of `threshold + PHRED33` is equivalent to the
    // saturating-sub form: bytes < PHRED33 have saturated Phred=0, which is < any
    // threshold ≥ 1.
    let cutoff = PHRED33.saturating_add(threshold);
    let cutoff_vec = u8x32::splat(cutoff);
    let mut count = 0u64;
    let chunks = qual.chunks_exact(32);
    let tail = chunks.remainder();
    for chunk in chunks {
        let v = u8x32::new(chunk.try_into().unwrap());
        count += v.simd_lt(cutoff_vec).to_bitmask().count_ones() as u64;
    }
    for &q in tail {
        if q < cutoff {
            count += 1;
        }
    }
    count
}

/// Checks the per-read filters (length, N-base, mean quality, low-qual fraction) against
/// all mates in a PE pair. If ANY mate fails ANY filter, the reason is returned and the
/// caller should drop all mates. Filters are evaluated in order so the first failing one
/// is reported.
fn evaluate_filters(
    records: &[OwnedRecord],
    stats: &[BaseStats],
    length: LengthFilter,
    filter_max_ns: Option<usize>,
    filter_mean_qual: Option<u8>,
    filter_low_qual: Option<LowQualFilter>,
) -> Option<FilterReject> {
    debug_assert_eq!(records.len(), stats.len());
    for rec in records {
        let len = rec.seq.len();
        // Zero-length reads are always rejected as length failures, regardless of the
        // user's `--filter-length MIN` — carrying empty reads forward has no valid
        // interpretation downstream (empty mean-quality, 0-of-0 N-fraction, empty serialized
        // FASTQ records) and attaches a surprising `Quality` reject reason otherwise.
        if len == 0 || len < length.min {
            return Some(FilterReject::Length);
        }
        if let Some(max) = length.max
            && len > max
        {
            return Some(FilterReject::Length);
        }
    }
    if let Some(n_max) = filter_max_ns {
        // Reuse the N count computed during the post-trim stats pass — no second scan.
        for s in stats {
            if s.n_bases > n_max as u64 {
                return Some(FilterReject::NBases);
            }
        }
    }
    if let Some(min_q) = filter_mean_qual {
        const PHRED33: u32 = 33;
        for rec in records {
            // Empty reads are already rejected by the length check above, so qual is
            // non-empty here. Sum raw Phred+33 bytes into u32 (safe for reads up to
            // ~16 Mb: 255 × len < u32::MAX) so the hot loop is a textbook
            // widening-reduce that LLVM auto-vectorizes. Compare in u64 space against
            // `(min_q + 33) × len` — equivalent to `mean(q − 33) < min_q` for valid
            // Phred+33 input (where q ≥ 33), no saturating_sub needed.
            let raw: u32 = rec.qual.iter().map(|&q| u32::from(q)).sum();
            let threshold = (u64::from(min_q) + u64::from(PHRED33)) * rec.qual.len() as u64;
            if u64::from(raw) < threshold {
                return Some(FilterReject::Quality);
            }
        }
    }
    if let Some(lq) = filter_low_qual {
        for rec in records {
            let below = count_bases_below_q(&rec.qual, lq.threshold);
            // Integer comparison to avoid FP: below / total > max_fraction
            //   <=> below > max_fraction * total.
            if (below as f64) > lq.max_fraction * rec.qual.len() as f64 {
                return Some(FilterReject::LowQual);
            }
        }
    }
    None
}

/// Validates that adapter bases are IUPAC-compatible (ACGT + IUPAC codes including N,
/// in either case). Used to catch obvious typos up-front rather than at record-matching
/// time.
fn validate_adapter_bases(seq: &[u8]) -> Result<(), String> {
    for (i, &b) in seq.iter().enumerate() {
        if IUPAC_MASKS[b.to_ascii_uppercase() as usize] == 0 {
            return Err(format!("invalid base {:?} at position {}", b as char, i));
        }
    }
    Ok(())
}

/// Returns true if a read base and an adapter base are compatible under IUPAC semantics.
///
/// The adapter may carry IUPAC ambiguity codes (including N); the read is typically plain
/// ACGT + N. Two bases match when their IUPAC masks share at least one bit.
#[inline]
fn base_matches_iupac(read_base: u8, adapter_base: u8) -> bool {
    let r = IUPAC_MASKS[read_base.to_ascii_uppercase() as usize];
    let a = IUPAC_MASKS[adapter_base.to_ascii_uppercase() as usize];
    r != 0 && a != 0 && (r & a) != 0
}

/// Finds the best 3'-anchored match of `adapter` against the 3' end of `read`, returning
/// the trim position (index in the read where the adapter starts) or `None`. A match is
/// valid when the overlap is at least `min_length` bases and its mismatch rate is at most
/// `max_mm_rate`.
///
/// When multiple overlap lengths meet the threshold, the LONGEST overlap is returned
/// (i.e. the smallest trim position) to be maximally conservative about adapter bases.
///
/// Dispatches on `adapter.pure_acgt` to pick the per-position compare kernel:
/// * pure ACGT → [`count_mismatches_ci_bounded`] (u8x32 SIMD, bounded early-exit);
/// * any IUPAC code → the scalar IUPAC-aware counter via [`base_matches_iupac`].
fn find_adapter_3prime(
    read: &[u8],
    adapter: &Adapter,
    min_length: usize,
    max_mm_rate: f64,
) -> Option<usize> {
    if read.len() < min_length || adapter.bytes.is_empty() {
        return None;
    }
    let max_start = read.len() - min_length;
    for k in 0..=max_start {
        let alignment_len = adapter.bytes.len().min(read.len() - k);
        if alignment_len < min_length {
            continue;
        }
        let max_mm = (alignment_len as f64 * max_mm_rate).floor() as usize;
        let mismatches = if adapter.pure_acgt {
            count_mismatches_ci_bounded(
                &read[k..k + alignment_len],
                &adapter.bytes[..alignment_len],
                max_mm,
            )
        } else {
            read[k..k + alignment_len]
                .iter()
                .zip(adapter.bytes[..alignment_len].iter())
                .filter(|(r, a)| !base_matches_iupac(**r, **a))
                .count()
        };
        if mismatches <= max_mm {
            return Some(k);
        }
    }
    None
}

/// Returns the longest adapter-trim position across a set of candidate adapters, or `None`
/// if no adapter matches. "Longest" means the smallest returned index (earliest trim).
fn find_best_adapter_match(
    read: &[u8],
    adapters: &[Adapter],
    min_length: usize,
    max_mm_rate: f64,
) -> Option<usize> {
    let mut best: Option<usize> = None;
    for adapter in adapters {
        if let Some(k) = find_adapter_3prime(read, adapter, min_length, max_mm_rate)
            && best.is_none_or(|b| k < b)
        {
            best = Some(k);
        }
    }
    best
}

/// Detects paired-end read-through via overlap analysis and returns the inferred insert
/// length when both reads agree within the mismatch-rate budget.
///
/// This handles only the "short-insert / adapter read-through" scenario (insert length L
/// is less than or equal to the read length). The "long-insert / middle overlap" scenario
/// (read_length < L < 2*read_length) is a different geometry used for read merging and is
/// not implemented here — but also does not produce trimmable adapter, so it's out of
/// scope for adapter trimming.
///
/// Algorithm: compute R2_rc (reverse complement of R2). In the short-insert case, R2_rc
/// lays out as `[rc(adapter2)][insert]` where the insert is at the END (positions
/// `r2_rc_len - L .. r2_rc_len`), and R1 lays out as `[insert][adapter1]` with the insert
/// at the START. The two inserts should match, so we look for the largest L such that
/// `R1[0..L]` equals `R2_rc[r2_rc_len - L .. r2_rc_len]` within the mismatch budget.
///
/// Evidence check: when a candidate overlap implies a trim of at least
/// `ADAPTER_EVIDENCE_MIN_TRIM` bases, the post-cut bases on *both* mates are compared
/// against an adapter-prefix library. The overlap is only accepted when both mates show
/// adapter-like sequence — disambiguating tandem-repeat matches where multiple overlap
/// lengths pass the probe check.
#[allow(clippy::too_many_arguments)]
fn detect_pe_overlap(
    r1: &[u8],
    r2: &[u8],
    min_overlap: usize,
    max_mm_rate: f64,
    diagnostic_len: usize,
    adapter_library: &OverlapAdapterLibrary,
    walk_mode: OverlapWalkMode,
    rc_scratch: &mut Vec<u8>,
) -> Option<usize> {
    if r1.len() < min_overlap || r2.len() < min_overlap {
        return None;
    }
    // ACGT-specialized RC is correct for Illumina reads (the only source of R2 in
    // practice). The subsequent comparison is already case-insensitive, so case drift
    // on unknown bytes wouldn't change the match outcome anyway.
    reverse_complement_acgt_into(r2, rc_scratch);
    let max_overlap = r1.len().min(rc_scratch.len());

    match walk_mode {
        OverlapWalkMode::Descending => detect_descending(
            r1,
            r2,
            rc_scratch,
            min_overlap,
            max_overlap,
            max_mm_rate,
            diagnostic_len,
            adapter_library,
        ),
        OverlapWalkMode::Ascending => detect_ascending(
            r1,
            r2,
            rc_scratch,
            min_overlap,
            max_overlap,
            max_mm_rate,
            diagnostic_len,
            adapter_library,
        ),
        OverlapWalkMode::Outward(center) => detect_outward(
            r1,
            r2,
            rc_scratch,
            min_overlap,
            max_overlap,
            max_mm_rate,
            diagnostic_len,
            adapter_library,
            center,
        ),
    }
}

/// Tests a single candidate overlap length `L`. First compares `R1[0..probe]` against
/// `revcomp(R2)[R-L..R-L+probe]` with [`count_mismatches_ci_bounded`]; if that probe
/// passes, and the implied trim is large enough, performs the adapter-evidence check
/// via [`post_cut_matches_library`] on both mates. Returns one of three
/// [`ProbeOutcome`] variants describing the decision for the caller's walk to act on.
#[inline]
fn try_overlap(
    r1: &[u8],
    r2: &[u8],
    r2_rc: &[u8],
    overlap: usize,
    max_mm_rate: f64,
    diagnostic_len: usize,
    adapter_library: &OverlapAdapterLibrary,
) -> ProbeOutcome {
    let probe_len = overlap.min(diagnostic_len);
    let r1_probe = &r1[..probe_len];
    let r2_probe = &r2_rc[r2_rc.len() - overlap..r2_rc.len() - overlap + probe_len];
    let max_mm = (probe_len as f64 * max_mm_rate).floor() as usize;
    let mismatches = count_mismatches_ci_bounded(r1_probe, r2_probe, max_mm);
    if mismatches > max_mm {
        return ProbeOutcome::ProbeFail;
    }
    let trim = r1.len() - overlap;
    if trim < ADAPTER_EVIDENCE_MIN_TRIM || adapter_library.is_empty() {
        return ProbeOutcome::Accept;
    }
    // Check each mate's post-cut bases against only the prefixes expected on that
    // side. Cross-side checks would both cost extra iterations and weaken specificity
    // (R1 post-cut matching an R2-only adapter is not informative).
    if post_cut_matches_library(&r1[overlap..], &adapter_library.r1_prefixes)
        && post_cut_matches_library(&r2[overlap..], &adapter_library.r2_prefixes)
    {
        ProbeOutcome::Accept
    } else {
        ProbeOutcome::EvidenceFail
    }
}

/// Descending walk — `(min_overlap..=max_overlap).rev()`.
///
/// Termination rules:
/// * `ProbeOutcome::Accept` — return `Some(overlap)` immediately. Correct L found
///   (first one seen walking down, so also the largest).
/// * `ProbeOutcome::EvidenceFail` — return `None`. This is the stop-on-evidence-fail
///   short-circuit. See rationale below.
/// * `ProbeOutcome::ProbeFail` — continue to the next smaller overlap.
///
/// # Why evidence failure can terminate the search
///
/// A probe-pass at `overlap = L` means `R1[0..probe]` matched `revcomp(R2[0..probe])`
/// within the mismatch budget. For standard library preps, `adapter_R1` and
/// `revcomp(adapter_R2)` differ in their first bases (TruSeq, Nextera, MGI, AVITI all
/// verified). So if the true insert were smaller than L, the probe window would include
/// some adapter_R1 bases paired against `revcomp(adapter_R2)` bases — those would NOT
/// match, and the probe would have failed. Therefore probe-pass at L implies the true
/// insert is ≥ L.
///
/// Combine that inference with descending iteration:
/// * Every L' > L has already been tested by this walk and rejected at probe time.
/// * The inference above rules out every L' < L.
///
/// So the current L is the *only* candidate left. Evidence failure at L means one of:
/// 1. The true insert IS L but the adapter isn't in our library (custom / novel), or
/// 2. The probe match was spurious (tandem-repeat / low-complexity coincidence).
///
/// Both cases are safer to leave untrimmed than to trim incorrectly: return `None`.
#[allow(clippy::too_many_arguments)]
fn detect_descending(
    r1: &[u8],
    r2: &[u8],
    rc_scratch: &[u8],
    min_overlap: usize,
    max_overlap: usize,
    max_mm_rate: f64,
    diagnostic_len: usize,
    adapter_library: &OverlapAdapterLibrary,
) -> Option<usize> {
    for overlap in (min_overlap..=max_overlap).rev() {
        match try_overlap(r1, r2, rc_scratch, overlap, max_mm_rate, diagnostic_len, adapter_library)
        {
            ProbeOutcome::Accept => return Some(overlap),
            ProbeOutcome::EvidenceFail => return None,
            ProbeOutcome::ProbeFail => continue,
        }
    }
    None
}

/// Ascending walk — `min_overlap..=max_overlap`.
///
/// Termination rules:
/// * `ProbeOutcome::Accept` — return `Some(overlap)`. First probe+evidence pass wins.
/// * `ProbeOutcome::EvidenceFail` — continue. We can NOT short-circuit here.
/// * `ProbeOutcome::ProbeFail` — continue to the next larger overlap.
///
/// # Why evidence failure cannot terminate the search
///
/// The same inference applies (`probe-pass ⇒ L_true ≥ L`) — but in ascending order, we
/// haven't yet tested any L' > L. A larger L could well be the true insert with real
/// adapter evidence at its post-cut window; giving up after the first evidence failure
/// would miss those. So we must continue scanning upward.
///
/// The practical consequence: Ascending is fast when the true insert is small (the
/// probe+evidence match appears early). When the true insert is near or above the read
/// length, Ascending is slow — it scans the full range.
#[allow(clippy::too_many_arguments)]
fn detect_ascending(
    r1: &[u8],
    r2: &[u8],
    rc_scratch: &[u8],
    min_overlap: usize,
    max_overlap: usize,
    max_mm_rate: f64,
    diagnostic_len: usize,
    adapter_library: &OverlapAdapterLibrary,
) -> Option<usize> {
    for overlap in min_overlap..=max_overlap {
        match try_overlap(r1, r2, rc_scratch, overlap, max_mm_rate, diagnostic_len, adapter_library)
        {
            ProbeOutcome::Accept => return Some(overlap),
            ProbeOutcome::EvidenceFail | ProbeOutcome::ProbeFail => continue,
        }
    }
    None
}

/// Outward walk — iterate `[center, center+1, center-1, center+2, center-2, ...]`,
/// clipped to `[min_overlap, max_overlap]`.
///
/// Termination rules:
/// * `ProbeOutcome::Accept` — return `Some(overlap)`. First probe+evidence pass wins.
/// * `ProbeOutcome::EvidenceFail` — continue. We can NOT short-circuit here.
/// * `ProbeOutcome::ProbeFail` — continue to the next outward step.
///
/// # Why evidence failure cannot terminate the search
///
/// `probe-pass ⇒ L_true ≥ L` still holds, but outward iteration visits candidates in a
/// non-monotone order — when we fail evidence at L, there may still be untested
/// candidates both above and below L (only the ones closer to `center` have been
/// checked). Terminating would risk missing a valid L on either side.
///
/// # When this wins
///
/// When `center` is a good estimate of the typical insert, the correct L is found
/// within a few outward steps. Accuracy of the estimate dominates the iteration count.
#[allow(clippy::too_many_arguments)]
fn detect_outward(
    r1: &[u8],
    r2: &[u8],
    rc_scratch: &[u8],
    min_overlap: usize,
    max_overlap: usize,
    max_mm_rate: f64,
    diagnostic_len: usize,
    adapter_library: &OverlapAdapterLibrary,
    center: usize,
) -> Option<usize> {
    let c = center.clamp(min_overlap, max_overlap);
    let max_offset = (c - min_overlap).max(max_overlap - c);
    // Probe the center exactly once before walking outward — `c - 0 == c + 0 == c` so
    // a naïve `k in 0..=max_offset` loop would otherwise probe the center twice.
    match try_overlap(r1, r2, rc_scratch, c, max_mm_rate, diagnostic_len, adapter_library) {
        ProbeOutcome::Accept => return Some(c),
        ProbeOutcome::EvidenceFail | ProbeOutcome::ProbeFail => {}
    }
    for k in 1..=max_offset {
        for overlap in [
            c.checked_sub(k).filter(|l| *l >= min_overlap),
            (c + k <= max_overlap).then_some(c + k),
        ]
        .into_iter()
        .flatten()
        {
            match try_overlap(
                r1,
                r2,
                rc_scratch,
                overlap,
                max_mm_rate,
                diagnostic_len,
                adapter_library,
            ) {
                ProbeOutcome::Accept => return Some(overlap),
                ProbeOutcome::EvidenceFail | ProbeOutcome::ProbeFail => continue,
            }
        }
    }
    None
}

/// Returns true if `post_cut` bases, when compared against any prefix in `library`,
/// match within the mismatch budget. Compares the first `min(ADAPTER_EVIDENCE_PROBE_LEN,
/// len(post_cut), len(prefix))` bases of each pair.
fn post_cut_matches_library(post_cut: &[u8], library: &[Vec<u8>]) -> bool {
    if post_cut.len() < ADAPTER_EVIDENCE_PROBE_LEN {
        // Proportionally scale the mm budget for short probes to maintain a similar
        // false-match probability; very short probes are too noisy and fall back to
        // a permissive accept.
        if post_cut.len() < 8 {
            return true;
        }
    }
    for prefix in library {
        let n = post_cut.len().min(prefix.len()).min(ADAPTER_EVIDENCE_PROBE_LEN);
        if n < 8 {
            continue;
        }
        let budget = (n * ADAPTER_EVIDENCE_MAX_MM).div_ceil(ADAPTER_EVIDENCE_PROBE_LEN);
        let mm = count_mismatches_ci_bounded(&post_cut[..n], &prefix[..n], budget);
        if mm <= budget {
            return true;
        }
    }
    false
}

/// Compiles the effective adapter set for R1 (and R2 if paired) from the CLI arguments.
/// The caller is responsible for ensuring validation has already succeeded.
fn build_adapter_set(
    adapter_sequence: &[String],
    adapter_fasta: &Option<PathBuf>,
    kits: &[String],
    num_inputs: usize,
) -> Result<AdapterSet> {
    let mut r1: Vec<Vec<u8>> = Vec::new();
    let mut r2: Vec<Vec<u8>> = Vec::new();

    // Explicit per-mate sequences.
    if let Some(s) = adapter_sequence.first() {
        r1.push(s.as_bytes().to_ascii_uppercase());
    }
    if let Some(s) = adapter_sequence.get(1) {
        r2.push(s.as_bytes().to_ascii_uppercase());
    }

    // Kit presets contribute to both mates.
    for kit_name in kits {
        let kits = expand_kit_name(kit_name).ok_or_else(|| anyhow!("Unknown kit {kit_name:?}"))?;
        for kit in kits {
            push_kit_adapters(kit, &mut r1, &mut r2);
        }
    }

    // FASTA adapters are applied to both mates.
    if let Some(path) = adapter_fasta {
        let seqs = load_adapter_fasta(path)?;
        for seq in seqs {
            let upper = seq.to_ascii_uppercase();
            r1.push(upper.clone());
            if num_inputs >= 2 {
                r2.push(upper);
            }
        }
    }

    dedupe_sort_by_len(&mut r1);
    dedupe_sort_by_len(&mut r2);

    // Classify ACGT-vs-IUPAC once at build time so the per-read scan just reads a bool.
    Ok(AdapterSet {
        r1: r1.into_iter().map(Adapter::new).collect(),
        r2: r2.into_iter().map(Adapter::new).collect(),
    })
}

/// Loads adapter sequences from a FASTA file. Only returns the sequences themselves
/// (names are discarded). Sequences must be IUPAC-compatible.
fn load_adapter_fasta(path: &Path) -> Result<Vec<Vec<u8>>> {
    let reader = Io::new(5, BUFFER_SIZE)
        .new_reader(path)
        .map_err(|e| anyhow!("Failed to open adapter FASTA {path:?}: {e}"))?;
    let mut out: Vec<Vec<u8>> = Vec::new();
    let mut current: Vec<u8> = Vec::new();
    for line in reader.lines() {
        let line = line.map_err(|e| anyhow!("Read error in {path:?}: {e}"))?;
        let trimmed = line.trim_end_matches(&['\r', '\n'][..]);
        if trimmed.starts_with('>') {
            if !current.is_empty() {
                out.push(std::mem::take(&mut current));
            }
        } else {
            for &b in trimmed.as_bytes() {
                if !b.is_ascii_whitespace() {
                    current.push(b);
                }
            }
        }
    }
    if !current.is_empty() {
        out.push(current);
    }
    for (i, seq) in out.iter().enumerate() {
        validate_adapter_bases(seq)
            .map_err(|m| anyhow!("--adapter-fasta record {}: {m} (sequence {seq:?})", i + 1))?;
    }
    Ok(out)
}

/// Builds the adapter-prefix library used by the PE-overlap evidence check. Seeds with
/// every kit in `ALL_KITS` so common Illumina / Nextera / MGI / AVITI libraries are
/// recognized out of the box, then unions in any user-supplied adapter sequences (CLI,
/// FASTA, and kit-derived) so custom / third-party adapters also produce positive
/// evidence. Prefixes are truncated to the first 16 bp to align with
/// `ADAPTER_EVIDENCE_PROBE_LEN` and deduplicated.
fn build_overlap_adapter_library(
    adapter_sequence: &[String],
    adapter_fasta: &Option<PathBuf>,
    user_adapters: &AdapterSet,
) -> Result<OverlapAdapterLibrary> {
    // Truncate to ADAPTER_EVIDENCE_PROBE_LEN and uppercase for case-insensitive compares;
    // skip prefixes shorter than 8 bp (too short to be discriminative).
    fn push(dst: &mut Vec<Vec<u8>>, seq: &[u8]) {
        let n = seq.len().min(ADAPTER_EVIDENCE_PROBE_LEN);
        if n >= 8 {
            dst.push(seq[..n].to_ascii_uppercase());
        }
    }

    let mut r1_prefixes: Vec<Vec<u8>> = Vec::new();
    let mut r2_prefixes: Vec<Vec<u8>> = Vec::new();

    // Kit presets shipped with the crate. Each kit's R1 adapter goes in the R1 list
    // and R2 adapter (if any) in the R2 list. Kits without an R2 adapter (small-RNA)
    // contribute only to R1 — callers running PE small-RNA must pass R2 explicitly.
    for kit in fqtk_lib::adapter_db::ALL_KITS {
        push(&mut r1_prefixes, kit.seq_r1);
        if let Some(s2) = kit.seq_r2 {
            push(&mut r2_prefixes, s2);
        }
    }
    // CLI --adapter-sequence takes up to two values: first is R1, second is R2.
    if let Some(s) = adapter_sequence.first() {
        push(&mut r1_prefixes, s.as_bytes());
    }
    if let Some(s) = adapter_sequence.get(1) {
        push(&mut r2_prefixes, s.as_bytes());
    }
    // FASTA adapters don't carry mate-side annotation; apply to both sides.
    if let Some(path) = adapter_fasta {
        let seqs = load_adapter_fasta(path)?;
        for seq in &seqs {
            push(&mut r1_prefixes, seq);
            push(&mut r2_prefixes, seq);
        }
    }
    // `AdapterSet` is already split by mate.
    for ad in &user_adapters.r1 {
        push(&mut r1_prefixes, &ad.bytes);
    }
    for ad in &user_adapters.r2 {
        push(&mut r2_prefixes, &ad.bytes);
    }

    r1_prefixes.sort();
    r1_prefixes.dedup();
    r2_prefixes.sort();
    r2_prefixes.dedup();
    Ok(OverlapAdapterLibrary { r1_prefixes, r2_prefixes })
}

/// Appends a kit's R1 (and optional R2) adapter sequences to the mate-split adapter
/// accumulators used by [`build_adapter_set`]. See the note about R2 semantics inline.
fn push_kit_adapters(kit: &KitAdapter, r1: &mut Vec<Vec<u8>>, r2: &mut Vec<Vec<u8>>) {
    r1.push(kit.seq_r1.to_vec());
    // Kits that don't specify an R2 adapter (e.g. small-RNA, a single-end preset) leave
    // R2 untouched. Using the R1 adapter on R2 is wrong for small-RNA PE (which carries
    // a different adapter on R2) and it's better to require users to supply that
    // explicitly than to silently guess.
    if let Some(s2) = kit.seq_r2 {
        r2.push(s2.to_vec());
    }
}

/// In-place: deduplicate a list of adapter sequences and reorder longest-first. Used on
/// the per-mate adapter lists before they are installed on [`AdapterSet`] so scans prefer
/// the longest (most specific) match.
fn dedupe_sort_by_len(v: &mut Vec<Vec<u8>>) {
    v.sort();
    v.dedup();
    // Longer adapters first so the 3'-anchored scan prefers specific matches.
    v.sort_by_key(|x| std::cmp::Reverse(x.len()));
}

/// Applies a read-structure to one FASTQ record, replacing its `seq` and `qual` with the
/// concatenated template bases and extending `umi_parts` with any extracted M-segment
/// bases. When `discard_unsupported` is true, B and C segments are treated as Skip;
/// otherwise they have already been rejected in validation.
fn apply_read_structure(
    rs: &ReadStructure,
    rec: &mut OwnedRecord,
    discard_unsupported: bool,
    umi_parts: &mut Vec<Vec<u8>>,
    template_seq: &mut Vec<u8>,
    template_qual: &mut Vec<u8>,
) -> Result<()> {
    let min_len: usize = rs.iter().map(|s| s.length().unwrap_or(1)).sum();
    if rec.seq.len() < min_len {
        return Err(anyhow!(
            "Read {} has too few bases ({}) for read-structure {rs} (needs at least {min_len}).",
            String::from_utf8_lossy(&rec.head),
            rec.seq.len()
        ));
    }

    template_seq.clear();
    template_qual.clear();

    for seg in rs.iter() {
        let (seg_seq, seg_qual) =
            seg.extract_bases_and_quals(&rec.seq, &rec.qual).map_err(|e| {
                anyhow!(
                    "Error extracting bases/quals for segment {seg} of read-structure {rs} \
                     from read {}: {e}",
                    String::from_utf8_lossy(&rec.head)
                )
            })?;
        match seg.kind {
            SegmentType::Template => {
                template_seq.extend_from_slice(seg_seq);
                template_qual.extend_from_slice(seg_qual);
            }
            SegmentType::MolecularBarcode => umi_parts.push(seg_seq.to_vec()),
            SegmentType::Skip => {}
            SegmentType::SampleBarcode | SegmentType::CellularBarcode => {
                // `validate()` rejects these unless discard_unsupported; reaching here with
                // `discard_unsupported == false` is a logic bug.
                debug_assert!(discard_unsupported);
            }
            _ => {
                // `SegmentType` is #[non_exhaustive]; future variants default to discard so
                // trim stays safe when the read-structure crate adds new kinds. Surface a
                // warning so the user can spot this behavior change on their data.
                log::warn!(
                    "Unknown read-structure segment kind {:?} treated as skip — \
                     update fqtk to handle this kind explicitly.",
                    seg.kind
                );
            }
        }
    }

    // Swap keeps the caller's scratch Vecs alive across iterations while handing the
    // freshly-built template off into the record. After the swap, the scratch Vecs hold
    // the old seq/qual capacity, ready for reuse next call.
    std::mem::swap(&mut rec.seq, template_seq);
    std::mem::swap(&mut rec.qual, template_qual);
    Ok(())
}

/// Joins multiple M-segment bases with `-`, matching fgumi's concatenation convention.
fn join_umi(parts: &[Vec<u8>]) -> Vec<u8> {
    let total = parts.iter().map(|p| p.len()).sum::<usize>() + parts.len().saturating_sub(1);
    let mut out = Vec::with_capacity(total);
    for (i, part) in parts.iter().enumerate() {
        if i > 0 {
            out.push(UMI_JOIN);
        }
        out.extend_from_slice(part);
    }
    out
}

/// Rewrites a FASTQ head (the bytes after `@` and before the newline) so that the read-id
/// carries the given UMI as its 8th colon-delimited field.
///
/// - If the read-id has ≤ 6 colons (0–7 fields): append `:UMI` to extend it to field 8.
/// - If the read-id has exactly 7 colons (8 fields): the 8th field is presumed to already
///   hold a UMI; append `-UMI` to extend it with the new UMI. (Note: `fqtk demux` uses `+`
///   as the separator in the same situation; we use `-` for internal consistency with the
///   multi-segment M join and because fgumi normalizes `+` → `-` on read. Running
///   `fqtk demux` output through `fqtk trim` will therefore yield headers with mixed
///   `+` / `-` separators in field 8; both are parsed identically by fgumi.)
/// - If the read-id has ≥ 8 colons: return an error (malformed header).
///
/// The space-separated comment (read-num / filter-flag / control / index fields) is preserved
/// untouched.
fn append_umi_to_head(head: &mut Vec<u8>, umi: &[u8]) -> Result<()> {
    let space_idx = head.iter().position(|&b| b == b' ');
    let name_end = space_idx.unwrap_or(head.len());
    let name = &head[..name_end];
    let colons = name.iter().filter(|&&b| b == UMI_ID_SEP).count();

    if colons + 1 > MAX_READ_ID_FIELDS {
        return Err(anyhow!(
            "Cannot append UMI to read-id with more than {MAX_READ_ID_FIELDS} colon-delimited \
             fields: {}",
            String::from_utf8_lossy(name)
        ));
    }

    let joiner = if colons + 1 == MAX_READ_ID_FIELDS { UMI_JOIN } else { UMI_ID_SEP };

    // Insert `[joiner, umi...]` at `name_end` without allocating a separate copy of the
    // trailing comment. `splice` uses the ExactSizeIterator length hint to shift the
    // tail exactly once.
    head.splice(name_end..name_end, std::iter::once(joiner).chain(umi.iter().copied()));
    Ok(())
}

/// Percentage helper (0.0 when denom is 0). Built on `ratio` so the zero-denominator
/// guard lives in one place.
fn pct(num: u64, denom: u64) -> f64 {
    ratio(num, denom) * 100.0
}

/// Ratio-as-fraction (0.0 to 1.0). Returns 0.0 when the denominator is zero (avoiding
/// NaN in the JSON output for empty inputs — MultiQC prefers finite numbers).
fn ratio(num: u64, denom: u64) -> f64 {
    if denom == 0 { 0.0 } else { num as f64 / denom as f64 }
}

/// Writes the pretty-printed fastp-shape JSON report to `path`. Errors from `create`
/// or `serde_json` are wrapped with the path / context for user-facing reporting.
fn write_json_report(path: &Path, report: &FastpJsonReport) -> Result<()> {
    let file =
        File::create(path).map_err(|e| anyhow!("Failed to create JSON report at {path:?}: {e}"))?;
    let writer = BufWriter::new(file);
    serde_json::to_writer_pretty(writer, report)
        .map_err(|e| anyhow!("Failed to serialize JSON report: {e}"))?;
    Ok(())
}

/// Resolves a path to an absolute form for input-vs-output equality checks.
/// If the path exists we use `canonicalize` (resolves symlinks and `..`);
/// otherwise we fall back to `std::path::absolute` (prepends the CWD).
fn resolve_absolute(p: &Path) -> Option<PathBuf> {
    if p.exists() { std::fs::canonicalize(p).ok() } else { std::path::absolute(p).ok() }
}

/// Extracts a human-readable message from a thread-panic payload. Panics in Rust carry
/// a `Box<dyn Any + Send>` value which is most commonly `&'static str` (from
/// `panic!("literal")`) or `String` (from `panic!("{}", ...)`); we handle both and fall
/// back to a generic note for unusual payloads.
fn panic_message(payload: &Box<dyn std::any::Any + Send>) -> String {
    if let Some(s) = payload.downcast_ref::<&'static str>() {
        (*s).to_string()
    } else if let Some(s) = payload.downcast_ref::<String>() {
        s.clone()
    } else {
        "(non-string panic payload)".to_string()
    }
}

/// Picks the "best" error to surface to the user from a list collected during a run.
/// Prefers errors that are likely a *root cause* (e.g. a specific worker/writer IO
/// failure) over errors that are a *symptom* of one (e.g. "workers exited before
/// receiving batch" — those always follow an earlier failure). Returns `None` if the
/// list is empty.
fn select_most_specific_error(mut errors: Vec<anyhow::Error>) -> Option<anyhow::Error> {
    // Errors that are typical *symptoms* of a peer failure (channel-closed semantics).
    let is_symptom = |e: &anyhow::Error| {
        let m = e.to_string();
        m.contains("exited before") || m.contains("dropped")
    };
    // Prefer a non-symptom error. Otherwise return the first error seen (or None).
    if let Some(i) = errors.iter().position(|e| !is_symptom(e)) {
        Some(errors.swap_remove(i))
    } else if errors.is_empty() {
        None
    } else {
        Some(errors.swap_remove(0))
    }
}

/// Scalar reverse-complement that preserves IUPAC ambiguity codes and case. The
/// production PE-overlap hot path uses [`reverse_complement_acgt_into`] instead —
/// this general-purpose variant is kept for future callers (e.g. FASTA adapter
/// preprocessing) and as a behavioral reference in tests.
#[cfg(test)]
fn reverse_complement_into(seq: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.reserve(seq.len());
    for &b in seq.iter().rev() {
        out.push(match b {
            b'A' => b'T',
            b'T' | b'U' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            b'a' => b't',
            b't' | b'u' => b'a',
            b'c' => b'g',
            b'g' => b'c',
            b'R' => b'Y',
            b'Y' => b'R',
            b'r' => b'y',
            b'y' => b'r',
            b'S' => b'S',
            b'W' => b'W',
            b's' => b's',
            b'w' => b'w',
            b'K' => b'M',
            b'M' => b'K',
            b'k' => b'm',
            b'm' => b'k',
            b'B' => b'V',
            b'V' => b'B',
            b'b' => b'v',
            b'v' => b'b',
            b'D' => b'H',
            b'H' => b'D',
            b'd' => b'h',
            b'h' => b'd',
            b'N' | b'n' => b,
            other => other,
        });
    }
}

/// Reverse-complement a byte sequence, returning a fresh Vec. Used only in tests;
/// production callers go through [`reverse_complement_into`] with a reusable buffer.
#[cfg(test)]
fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut out = Vec::new();
    reverse_complement_into(seq, &mut out);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use seq_io::fastq::OwnedRecord;
    use tempfile::TempDir;

    /// Builds FASTQ content lines for `n` reads with names `@{prefix}_0`, `@{prefix}_1`, ...
    /// All qualities are `I` (Q40 at Phred+33).
    fn fq_lines(prefix: &str, reads: &[&str]) -> Vec<String> {
        reads
            .iter()
            .enumerate()
            .flat_map(|(i, &seq)| {
                vec![
                    format!("@{prefix}_{i}"),
                    seq.to_string(),
                    "+".to_string(),
                    "I".repeat(seq.len()),
                ]
            })
            .collect()
    }

    fn write_fastq(tmp: &TempDir, name: &str, lines: &[String]) -> PathBuf {
        let path = tmp.path().join(format!("{name}.fq"));
        Io::default().write_lines(&path, lines).unwrap();
        path
    }

    fn read_fastq(path: &Path) -> Vec<OwnedRecord> {
        let io = Io::default();
        FastqReader::new(io.new_reader(path).unwrap())
            .into_records()
            .collect::<Result<Vec<_>, seq_io::fastq::Error>>()
            .unwrap()
    }

    fn trim_cmd(inputs: Vec<PathBuf>, outputs: Vec<PathBuf>, metrics: Option<PathBuf>) -> Trim {
        Trim {
            inputs,
            outputs,
            threads: 2,
            compression_level: 1,
            metrics,
            read_structures: vec![],
            discard_unsupported_segments: false,
            adapter_sequence: vec![],
            adapter_fasta: None,
            kit: vec![],
            detect_adapter_for_pe: false,
            overlap_min_length: 5,
            overlap_max_mismatch_rate: 0.1,
            overlap_diagnostic_length: usize::MAX,
            adapter_min_length: 5,
            adapter_mismatch_rate: 0.1,
            trim_polyg: 0, // tests default poly-G off unless they opt in
            trim_polyx: None,
            quality_trim_3p: None,
            quality_trim_5p: None,
            filter_length: LengthFilter { min: 0, max: None },
            filter_max_ns: None,
            filter_mean_qual: None,
            filter_low_qual: None,
            json: None,
            batch_size: 1024,
            expected_insert_size: None,
        }
    }

    fn owned_rec(head: &str, seq: &str, qual: &str) -> OwnedRecord {
        assert_eq!(seq.len(), qual.len());
        OwnedRecord {
            head: head.as_bytes().to_vec(),
            seq: seq.as_bytes().to_vec(),
            qual: qual.as_bytes().to_vec(),
        }
    }

    fn rs(s: &str) -> ReadStructure {
        s.parse().unwrap()
    }

    /// Test wrapper that supplies ephemeral scratch Vecs; production callers pass
    /// persistent buffers that amortize across the whole run.
    fn apply_rs(
        rs_spec: &ReadStructure,
        rec: &mut OwnedRecord,
        discard_unsupported: bool,
        umi_parts: &mut Vec<Vec<u8>>,
    ) -> Result<()> {
        let mut s = Vec::new();
        let mut q = Vec::new();
        apply_read_structure(rs_spec, rec, discard_unsupported, umi_parts, &mut s, &mut q)
    }

    /// Test wrapper that supplies an ephemeral RC scratch Vec. Uses `usize::MAX` for the
    /// diagnostic-length knob and an empty adapter library so legacy tests continue to
    /// probe the full overlap without invoking the evidence check — matching the
    /// pre-diagnostic-probe semantics these tests were written for.
    fn detect_overlap(r1: &[u8], r2: &[u8], min_overlap: usize, max_mm_rate: f64) -> Option<usize> {
        let mut scratch = Vec::new();
        let empty_library = OverlapAdapterLibrary::default();
        detect_pe_overlap(
            r1,
            r2,
            min_overlap,
            max_mm_rate,
            usize::MAX,
            &empty_library,
            OverlapWalkMode::Descending,
            &mut scratch,
        )
    }

    /// Test wrapper that computes per-record stats on the fly; tests exercise filter
    /// logic, not the stats-precomputation optimization used in production.
    fn eval_filters(
        recs: &[OwnedRecord],
        length: LengthFilter,
        filter_max_ns: Option<usize>,
        filter_mean_qual: Option<u8>,
    ) -> Option<FilterReject> {
        let stats: Vec<BaseStats> = recs.iter().map(|r| observe_stats(&r.seq, &r.qual)).collect();
        evaluate_filters(recs, &stats, length, filter_max_ns, filter_mean_qual, None)
    }

    // ---- validation ----

    #[test]
    fn validation_rejects_missing_input() {
        let tmp = TempDir::new().unwrap();
        let cmd =
            trim_cmd(vec![tmp.path().join("nope.fq")], vec![tmp.path().join("out.fq.gz")], None);
        let err = cmd.validate().unwrap_err().to_string();
        assert!(err.contains("does not exist"), "{err}");
    }

    #[test]
    fn validation_rejects_output_count_mismatch() {
        let tmp = TempDir::new().unwrap();
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &["ACGT"]));
        let cmd = trim_cmd(
            vec![r1],
            vec![tmp.path().join("o1.fq.gz"), tmp.path().join("o2.fq.gz")],
            None,
        );
        let err = cmd.validate().unwrap_err().to_string();
        assert!(err.contains("must equal number of inputs"), "{err}");
    }

    #[test]
    fn validation_rejects_output_overwriting_input() {
        let tmp = TempDir::new().unwrap();
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &["ACGT"]));
        // Point output at the same path as the input
        let cmd = trim_cmd(vec![r1.clone()], vec![r1.clone()], None);
        let err = cmd.validate().unwrap_err().to_string();
        assert!(err.contains("refusing to overwrite"), "{err}");
    }

    #[test]
    fn validation_rejects_metrics_overwriting_input() {
        let tmp = TempDir::new().unwrap();
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &["ACGT"]));
        let cmd = trim_cmd(vec![r1.clone()], vec![tmp.path().join("out.fq.gz")], Some(r1));
        let err = cmd.validate().unwrap_err().to_string();
        assert!(err.contains("refusing to overwrite"), "{err}");
    }

    #[test]
    fn validation_rejects_output_parent_missing() {
        let tmp = TempDir::new().unwrap();
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &["ACGT"]));
        let cmd =
            trim_cmd(vec![r1], vec![tmp.path().join("nonexistent_dir").join("out.fq.gz")], None);
        let err = cmd.validate().unwrap_err().to_string();
        assert!(err.contains("parent directory"), "{err}");
    }

    #[test]
    fn validation_rejects_too_few_threads() {
        let tmp = TempDir::new().unwrap();
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &["ACGT"]));
        let mut cmd = trim_cmd(vec![r1], vec![tmp.path().join("out.fq.gz")], None);
        cmd.threads = 0;
        let err = cmd.validate().unwrap_err().to_string();
        assert!(err.contains("Threads must be at least 1"), "{err}");
    }

    #[test]
    fn validation_rejects_bad_compression_level() {
        let tmp = TempDir::new().unwrap();
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &["ACGT"]));
        let mut cmd = trim_cmd(vec![r1], vec![tmp.path().join("out.fq.gz")], None);
        cmd.compression_level = 0;
        let err = cmd.validate().unwrap_err().to_string();
        assert!(err.contains("Compression level"), "{err}");
    }

    #[test]
    fn validation_aggregates_multiple_errors() {
        let tmp = TempDir::new().unwrap();
        let mut cmd =
            trim_cmd(vec![tmp.path().join("missing.fq")], vec![tmp.path().join("out.fq.gz")], None);
        cmd.threads = 0;
        cmd.compression_level = 99;
        let err = cmd.validate().unwrap_err().to_string();
        assert!(err.contains("does not exist"), "{err}");
        assert!(err.contains("Threads must be at least 1"), "{err}");
        assert!(err.contains("Compression level"), "{err}");
    }

    // ---- execute: pass-through ----

    #[test]
    fn execute_single_end_passes_through() {
        let tmp = TempDir::new().unwrap();
        let reads: Vec<&str> = vec!["ACGTACGTAC", "GGGGAAAACCCC", "TTTT"];
        let r1 = write_fastq(&tmp, "r1", &fq_lines("read", &reads));
        let out = tmp.path().join("out.fq.gz");
        let metrics_path = tmp.path().join("trim-metrics.txt");

        let cmd = trim_cmd(vec![r1], vec![out.clone()], Some(metrics_path.clone()));
        cmd.execute().unwrap();

        let written = read_fastq(&out);
        assert_eq!(written.len(), reads.len());
        for (rec, &expected) in written.iter().zip(reads.iter()) {
            assert_eq!(rec.seq.as_slice(), expected.as_bytes());
        }
        // Metrics TSV exists and is non-empty
        assert!(metrics_path.exists());
    }

    #[test]
    fn execute_paired_end_passes_through() {
        let tmp = TempDir::new().unwrap();
        let r1_seqs: Vec<&str> = vec!["AAAAAAAA", "CCCCCCCC", "GGGGGGGG"];
        let r2_seqs: Vec<&str> = vec!["TTTTTTTT", "TTTTTTTT", "TTTTTTTT"];
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &r1_seqs));
        let r2 = write_fastq(&tmp, "r2", &fq_lines("r", &r2_seqs));
        let o1 = tmp.path().join("o1.fq.gz");
        let o2 = tmp.path().join("o2.fq.gz");

        let cmd = trim_cmd(vec![r1, r2], vec![o1.clone(), o2.clone()], None);
        cmd.execute().unwrap();

        let w1 = read_fastq(&o1);
        let w2 = read_fastq(&o2);
        assert_eq!(w1.len(), r1_seqs.len());
        assert_eq!(w2.len(), r2_seqs.len());
        for (rec, &seq) in w1.iter().zip(r1_seqs.iter()) {
            assert_eq!(rec.seq.as_slice(), seq.as_bytes());
        }
        for (rec, &seq) in w2.iter().zip(r2_seqs.iter()) {
            assert_eq!(rec.seq.as_slice(), seq.as_bytes());
        }
    }

    #[test]
    fn execute_empty_input_yields_empty_output() {
        let tmp = TempDir::new().unwrap();
        let r1 = write_fastq(&tmp, "r1", &[]);
        let out = tmp.path().join("out.fq.gz");
        let metrics_path = tmp.path().join("trim-metrics.txt");
        let cmd = trim_cmd(vec![r1], vec![out.clone()], Some(metrics_path.clone()));
        cmd.execute().unwrap();

        assert!(out.exists());
        let written = read_fastq(&out);
        assert_eq!(written.len(), 0);
        assert!(metrics_path.exists());
    }

    #[test]
    fn execute_errors_on_mismatched_pair_counts() {
        let tmp = TempDir::new().unwrap();
        let r1_seqs: Vec<&str> = vec!["ACGT"; 5];
        let r2_seqs: Vec<&str> = vec!["ACGT"; 3];
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &r1_seqs));
        let r2 = write_fastq(&tmp, "r2", &fq_lines("r", &r2_seqs));
        let o1 = tmp.path().join("o1.fq.gz");
        let o2 = tmp.path().join("o2.fq.gz");
        let cmd = trim_cmd(vec![r1, r2], vec![o1.clone(), o2.clone()], None);
        let err = cmd.execute().unwrap_err().to_string();
        assert!(err.contains("out of sync"), "{err}");
        // The reader batches input before workers/writers start, so a pair-mismatch
        // detected while reading aborts before any output file is created. That's the
        // desirable outcome — no partially-filled outputs to accidentally ingest
        // downstream.
        assert!(!o1.exists(), "o1 should not exist after pair-mismatch abort");
        assert!(!o2.exists(), "o2 should not exist after pair-mismatch abort");
    }

    #[test]
    fn execute_writes_metrics_counts() {
        let tmp = TempDir::new().unwrap();
        let reads: Vec<&str> = vec!["AAAA", "CCCCCC"]; // 4 + 6 = 10 bases, 2 reads
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &reads));
        let out = tmp.path().join("out.fq.gz");
        let metrics_path = tmp.path().join("metrics.txt");

        let cmd = trim_cmd(vec![r1], vec![out], Some(metrics_path.clone()));
        cmd.execute().unwrap();

        let contents = std::fs::read_to_string(&metrics_path).unwrap();
        // Header + one row
        let lines: Vec<&str> = contents.lines().collect();
        assert_eq!(lines.len(), 2);
        // Header contains our field names
        assert!(lines[0].contains("reads_in"));
        assert!(lines[0].contains("bases_in"));
        // Row has the correct counts
        let values: Vec<&str> = lines[1].split('\t').collect();
        let header: Vec<&str> = lines[0].split('\t').collect();
        let idx = |name: &str| header.iter().position(|h| *h == name).unwrap();
        assert_eq!(values[idx("reads_in")], "2");
        assert_eq!(values[idx("reads_out")], "2");
        assert_eq!(values[idx("bases_in")], "10");
        assert_eq!(values[idx("bases_out")], "10");
    }

    // ---- join_umi ----

    #[test]
    fn join_umi_single_part() {
        assert_eq!(join_umi(&[b"AACC".to_vec()]), b"AACC");
    }

    #[test]
    fn join_umi_two_parts() {
        assert_eq!(join_umi(&[b"AACC".to_vec(), b"GGTT".to_vec()]), b"AACC-GGTT");
    }

    #[test]
    fn join_umi_three_parts() {
        assert_eq!(join_umi(&[b"AA".to_vec(), b"CC".to_vec(), b"GG".to_vec()]), b"AA-CC-GG");
    }

    // ---- append_umi_to_head ----

    #[test]
    fn append_umi_to_short_head() {
        let mut head = b"readname".to_vec();
        append_umi_to_head(&mut head, b"AACCGG").unwrap();
        assert_eq!(head, b"readname:AACCGG");
    }

    #[test]
    fn append_umi_to_illumina_7_field_head() {
        let mut head = b"INSTR:123:FLOWCELL:1:1101:1000:2000".to_vec();
        append_umi_to_head(&mut head, b"AACCGG").unwrap();
        assert_eq!(head, b"INSTR:123:FLOWCELL:1:1101:1000:2000:AACCGG");
    }

    #[test]
    fn append_umi_extends_existing_field_8() {
        // 7 colons = 8 fields; field 8 is already a UMI. Append new UMI with `-`.
        let mut head = b"INSTR:123:FLOWCELL:1:1101:1000:2000:TTTT".to_vec();
        append_umi_to_head(&mut head, b"AACCGG").unwrap();
        assert_eq!(head, b"INSTR:123:FLOWCELL:1:1101:1000:2000:TTTT-AACCGG");
    }

    #[test]
    fn append_umi_preserves_comment_after_space() {
        let mut head = b"INSTR:123:FLOWCELL:1:1101:1000:2000 1:N:0:CTAG".to_vec();
        append_umi_to_head(&mut head, b"AACCGG").unwrap();
        assert_eq!(head, b"INSTR:123:FLOWCELL:1:1101:1000:2000:AACCGG 1:N:0:CTAG");
    }

    #[test]
    fn append_umi_rejects_too_many_fields() {
        let mut head = b"a:b:c:d:e:f:g:h:i".to_vec();
        let err = append_umi_to_head(&mut head, b"AACCGG").unwrap_err().to_string();
        assert!(err.contains("more than"), "{err}");
    }

    // ---- apply_read_structure ----

    #[test]
    fn apply_read_structure_template_only() {
        let mut rec = owned_rec("read1", "ACGTACGT", "IIIIIIII");
        let mut umis = vec![];
        apply_rs(&rs("+T"), &mut rec, false, &mut umis).unwrap();
        assert_eq!(rec.seq, b"ACGTACGT");
        assert_eq!(rec.qual, b"IIIIIIII");
        assert!(umis.is_empty());
    }

    #[test]
    fn apply_read_structure_hard_trim_5s_plus_t() {
        let mut rec = owned_rec("read1", "AAAAACCCC", "112233445");
        let mut umis = vec![];
        apply_rs(&rs("5S+T"), &mut rec, false, &mut umis).unwrap();
        assert_eq!(rec.seq, b"CCCC");
        assert_eq!(rec.qual, b"3445");
        assert!(umis.is_empty());
    }

    #[test]
    fn apply_read_structure_extracts_umi() {
        let mut rec = owned_rec("read1", "AACCGGTTTTAA", "111222333444");
        let mut umis = vec![];
        apply_rs(&rs("4M+T"), &mut rec, false, &mut umis).unwrap();
        assert_eq!(rec.seq, b"GGTTTTAA");
        assert_eq!(rec.qual, b"22333444");
        assert_eq!(umis, vec![b"AACC".to_vec()]);
    }

    #[test]
    fn apply_read_structure_skip_segment() {
        let mut rec = owned_rec("read1", "SSSTEMPL", "12345678");
        let mut umis = vec![];
        apply_rs(&rs("3S+T"), &mut rec, false, &mut umis).unwrap();
        assert_eq!(rec.seq, b"TEMPL");
        assert_eq!(rec.qual, b"45678");
        assert!(umis.is_empty());
    }

    #[test]
    fn apply_read_structure_too_short_errors() {
        let mut rec = owned_rec("read1", "ACG", "III");
        let mut umis = vec![];
        let err = apply_rs(&rs("5M+T"), &mut rec, false, &mut umis).unwrap_err().to_string();
        assert!(err.contains("too few bases"), "{err}");
    }

    #[test]
    fn apply_read_structure_multiple_m_segments() {
        let mut rec = owned_rec("read1", "AAABBBCCCDDD", "123456789012");
        let mut umis = vec![];
        apply_rs(&rs("3M3S3M+T"), &mut rec, false, &mut umis).unwrap();
        assert_eq!(rec.seq, b"DDD");
        assert_eq!(rec.qual, b"012");
        assert_eq!(umis, vec![b"AAA".to_vec(), b"CCC".to_vec()]);
    }

    #[test]
    fn apply_read_structure_discards_b_when_permitted() {
        let mut rec = owned_rec("read1", "BBBBTEMPL", "123456789");
        let mut umis = vec![];
        apply_rs(&rs("4B+T"), &mut rec, true, &mut umis).unwrap();
        assert_eq!(rec.seq, b"TEMPL");
        assert_eq!(rec.qual, b"56789");
    }

    // ---- validation for read-structures ----

    #[test]
    fn validation_rejects_b_segment() {
        let tmp = TempDir::new().unwrap();
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &["ACGT"]));
        let mut cmd = trim_cmd(vec![r1], vec![tmp.path().join("out.fq.gz")], None);
        cmd.read_structures = vec![rs("4B+T")];
        let err = cmd.validate().unwrap_err().to_string();
        assert!(err.contains("sample barcode (B)"), "{err}");
        assert!(err.contains("fqtk demux"), "{err}");
    }

    #[test]
    fn validation_rejects_c_segment() {
        let tmp = TempDir::new().unwrap();
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &["ACGT"]));
        let mut cmd = trim_cmd(vec![r1], vec![tmp.path().join("out.fq.gz")], None);
        cmd.read_structures = vec![rs("4C+T")];
        let err = cmd.validate().unwrap_err().to_string();
        assert!(err.contains("cellular barcode (C)"), "{err}");
    }

    #[test]
    fn validation_allows_b_with_discard_flag() {
        let tmp = TempDir::new().unwrap();
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &["ACGTACGT"]));
        let mut cmd = trim_cmd(vec![r1], vec![tmp.path().join("out.fq.gz")], None);
        cmd.read_structures = vec![rs("4B+T")];
        cmd.discard_unsupported_segments = true;
        cmd.validate().unwrap();
    }

    #[test]
    fn validation_rejects_wrong_read_structure_count() {
        let tmp = TempDir::new().unwrap();
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &["ACGT"]));
        let r2 = write_fastq(&tmp, "r2", &fq_lines("r", &["ACGT"]));
        let mut cmd = trim_cmd(
            vec![r1, r2],
            vec![tmp.path().join("o1.fq.gz"), tmp.path().join("o2.fq.gz")],
            None,
        );
        cmd.read_structures = vec![rs("+T")]; // only 1 but 2 inputs
        let err = cmd.validate().unwrap_err().to_string();
        assert!(err.contains("must be 0 or equal to number of inputs"), "{err}");
    }

    // ---- execute: read-structure end-to-end ----

    #[test]
    fn execute_se_hard_trim() {
        let tmp = TempDir::new().unwrap();
        let reads: Vec<&str> = vec!["AAAAATEMPL", "BBBBBSHORT"];
        let r1 = write_fastq(&tmp, "r1", &fq_lines("read", &reads));
        let out = tmp.path().join("out.fq.gz");
        let mut cmd = trim_cmd(vec![r1], vec![out.clone()], None);
        cmd.read_structures = vec![rs("5S+T")];
        cmd.execute().unwrap();

        let written = read_fastq(&out);
        assert_eq!(written.len(), 2);
        assert_eq!(written[0].seq.as_slice(), b"TEMPL");
        assert_eq!(written[1].seq.as_slice(), b"SHORT");
    }

    #[test]
    fn execute_se_umi_extraction_appends_to_head() {
        let tmp = TempDir::new().unwrap();
        // Write FASTQ with a standard 7-colon Illumina header
        let lines = vec![
            "@A:1:B:1:1:1:1".to_string(),
            "AACCGGTTTTAA".to_string(),
            "+".to_string(),
            "IIIIIIIIIIII".to_string(),
        ];
        let r1 = write_fastq(&tmp, "r1", &lines);
        let out = tmp.path().join("out.fq.gz");
        let mut cmd = trim_cmd(vec![r1], vec![out.clone()], None);
        cmd.read_structures = vec![rs("4M+T")];
        cmd.execute().unwrap();

        let written = read_fastq(&out);
        assert_eq!(written.len(), 1);
        assert_eq!(written[0].head.as_slice(), b"A:1:B:1:1:1:1:AACC");
        assert_eq!(written[0].seq.as_slice(), b"GGTTTTAA");
    }

    #[test]
    fn execute_pe_umi_from_both_mates_is_joined() {
        let tmp = TempDir::new().unwrap();
        let r1_lines = vec![
            "@A:1:B:1:1:1:1".to_string(),
            "AAAAGGGGGG".to_string(),
            "+".to_string(),
            "IIIIIIIIII".to_string(),
        ];
        let r2_lines = vec![
            "@A:1:B:1:1:1:1".to_string(),
            "TTTTCCCCCC".to_string(),
            "+".to_string(),
            "IIIIIIIIII".to_string(),
        ];
        let r1 = write_fastq(&tmp, "r1", &r1_lines);
        let r2 = write_fastq(&tmp, "r2", &r2_lines);
        let o1 = tmp.path().join("o1.fq.gz");
        let o2 = tmp.path().join("o2.fq.gz");

        let mut cmd = trim_cmd(vec![r1, r2], vec![o1.clone(), o2.clone()], None);
        cmd.read_structures = vec![rs("4M+T"), rs("4M+T")];
        cmd.execute().unwrap();

        let w1 = read_fastq(&o1);
        let w2 = read_fastq(&o2);
        assert_eq!(w1.len(), 1);
        assert_eq!(w2.len(), 1);
        // Both mates carry the same combined UMI `AAAA-TTTT`
        assert_eq!(w1[0].head.as_slice(), b"A:1:B:1:1:1:1:AAAA-TTTT");
        assert_eq!(w2[0].head.as_slice(), b"A:1:B:1:1:1:1:AAAA-TTTT");
        assert_eq!(w1[0].seq.as_slice(), b"GGGGGG");
        assert_eq!(w2[0].seq.as_slice(), b"CCCCCC");
    }

    #[test]
    fn apply_read_structure_rejects_read_exactly_at_fixed_length() {
        // Variable `+T` requires at least 1 base (matches demux's `min_len` convention).
        // Read of 4 bases with `4M+T` fails the length check — confirms we do not silently
        // emit a zero-length template record.
        let mut rec = owned_rec("read1", "AAAA", "IIII");
        let mut umis = vec![];
        let err = apply_rs(&rs("4M+T"), &mut rec, false, &mut umis).unwrap_err().to_string();
        assert!(err.contains("too few bases"), "{err}");
    }

    #[test]
    fn apply_read_structure_minimum_variable_length() {
        // Minimum valid case: fixed_sum + 1 bases, variable +T consumes exactly 1.
        let mut rec = owned_rec("read1", "AAAAT", "11112");
        let mut umis = vec![];
        apply_rs(&rs("4M+T"), &mut rec, false, &mut umis).unwrap();
        assert_eq!(rec.seq, b"T");
        assert_eq!(rec.qual, b"2");
        assert_eq!(umis, vec![b"AAAA".to_vec()]);
    }

    #[test]
    fn append_umi_to_head_extends_existing_field_8_idempotently() {
        // Simulates running `fqtk trim` twice: once produces `...:UMI1`, a second pass
        // appends `-UMI2` → `...:UMI1-UMI2`.
        let mut head = b"A:1:B:1:1:1:1:AAAA".to_vec();
        append_umi_to_head(&mut head, b"BBBB").unwrap();
        assert_eq!(head, b"A:1:B:1:1:1:1:AAAA-BBBB");
    }

    // ---- reverse_complement ----

    #[test]
    fn reverse_complement_basic() {
        assert_eq!(reverse_complement(b"ACGT"), b"ACGT");
        assert_eq!(reverse_complement(b"AAAA"), b"TTTT");
        assert_eq!(reverse_complement(b"ACCGTTT"), b"AAACGGT");
    }

    #[test]
    fn reverse_complement_preserves_n() {
        assert_eq!(reverse_complement(b"ACNGT"), b"ACNGT");
    }

    #[test]
    fn reverse_complement_empty() {
        assert_eq!(reverse_complement(b""), b"");
    }

    // ---- base_matches_iupac ----

    #[test]
    fn iupac_exact_match() {
        assert!(base_matches_iupac(b'A', b'A'));
        assert!(base_matches_iupac(b'C', b'C'));
        assert!(!base_matches_iupac(b'A', b'C'));
    }

    #[test]
    fn iupac_n_in_adapter_matches_any_read_base() {
        assert!(base_matches_iupac(b'A', b'N'));
        assert!(base_matches_iupac(b'C', b'N'));
        assert!(base_matches_iupac(b'G', b'N'));
        assert!(base_matches_iupac(b'T', b'N'));
    }

    #[test]
    fn iupac_ambiguity_code_matches_compatible_reads() {
        // R = A or G
        assert!(base_matches_iupac(b'A', b'R'));
        assert!(base_matches_iupac(b'G', b'R'));
        assert!(!base_matches_iupac(b'C', b'R'));
        assert!(!base_matches_iupac(b'T', b'R'));
    }

    #[test]
    fn iupac_case_insensitive() {
        assert!(base_matches_iupac(b'a', b'A'));
        assert!(base_matches_iupac(b'a', b'n'));
    }

    // ---- find_adapter_3prime ----

    // Helper: build an `Adapter` from a byte literal in tests.
    fn ad(bytes: &[u8]) -> Adapter {
        Adapter::new(bytes.to_vec())
    }

    #[test]
    fn adapter_3prime_full_adapter_at_end() {
        // read = INSERT + adapter; full adapter fits
        let read = b"TTTTTTTTTTAGATCGGAAGAG";
        let adapter = ad(b"AGATCGGAAGAG");
        assert_eq!(find_adapter_3prime(read, &adapter, 5, 0.1), Some(10));
    }

    #[test]
    fn adapter_3prime_partial_at_end() {
        let read = b"TTTTTTTTTTAGAT"; // only first 4 bp of adapter
        let adapter = ad(b"AGATCGGAAGAG");
        assert_eq!(find_adapter_3prime(read, &adapter, 4, 0.1), Some(10));
    }

    #[test]
    fn adapter_3prime_partial_below_min_overlap_no_match() {
        let read = b"TTTTTTTTTTAGA"; // only 3 bp of adapter
        let adapter = ad(b"AGATCGGAAGAG");
        assert_eq!(find_adapter_3prime(read, &adapter, 5, 0.1), None);
    }

    #[test]
    fn adapter_3prime_tolerates_one_mismatch() {
        // read has TGAT where adapter has AGAT → 1 mismatch in 4 bases = 25% rate
        let read = b"TTTTTTTTTTTGATCGGAAGAG";
        let adapter = ad(b"AGATCGGAAGAG");
        // 12-base overlap with 1 mismatch → 8.3% rate, within 10% threshold
        assert_eq!(find_adapter_3prime(read, &adapter, 5, 0.1), Some(10));
    }

    #[test]
    fn adapter_3prime_rejects_too_many_mismatches() {
        // read has CCCCCCCC at the end; adapter is AGATCGGA → all 8 bases mismatch
        let read = b"TTTTTTTTTTCCCCCCCC";
        let adapter = ad(b"AGATCGGA");
        assert_eq!(find_adapter_3prime(read, &adapter, 5, 0.1), None);
    }

    #[test]
    fn adapter_3prime_iupac_n_in_adapter() {
        // Exercises the scalar IUPAC fallback path: adapter contains `N` so `pure_acgt`
        // is false and we do NOT take the SIMD fast path.
        let read = b"TTTTTTTTTTAGATCGGAAGAG";
        let adapter = ad(b"AGANCGGAAGAG"); // N at pos 3 matches read's T
        assert!(!adapter.pure_acgt);
        assert_eq!(find_adapter_3prime(read, &adapter, 5, 0.0), Some(10));
    }

    #[test]
    fn adapter_3prime_pure_acgt_takes_simd_path() {
        // Sanity: a plain-ACGT adapter classifies as pure_acgt; test indirectly by
        // observing correctness on a match that the SIMD kernel must handle.
        let adapter = ad(b"AGATCGGAAGAG");
        assert!(adapter.pure_acgt);
    }

    #[test]
    fn adapter_3prime_empty_adapter_is_none() {
        let adapter = ad(b"");
        assert_eq!(find_adapter_3prime(b"ACGT", &adapter, 5, 0.1), None);
    }

    #[test]
    fn adapter_3prime_read_shorter_than_min_overlap() {
        let adapter = ad(b"AGATCGG");
        assert_eq!(find_adapter_3prime(b"ACGT", &adapter, 5, 0.1), None);
    }

    #[test]
    fn adapter_3prime_adapter_longer_than_read() {
        // Read has 5 bases matching the adapter's first 5 bases (partial adapter covers entire read).
        let read = b"AGATC";
        let adapter = ad(b"AGATCGGAAGAG");
        assert_eq!(find_adapter_3prime(read, &adapter, 5, 0.1), Some(0));
    }

    #[test]
    fn adapter_3prime_read_entirely_adapter() {
        let read = b"AGATCGGAAGAG";
        let adapter = ad(b"AGATCGGAAGAG");
        assert_eq!(find_adapter_3prime(read, &adapter, 5, 0.0), Some(0));
    }

    // ---- find_best_adapter_match ----

    #[test]
    fn best_adapter_match_picks_earliest_trim() {
        let read = b"INSERTAGAT"; // adapter at position 6 ("AGAT")
        let adapters = vec![ad(b"CCCC"), ad(b"AGAT"), ad(b"AGATCGG")];
        // AGAT at pos 6 matches (4bp overlap); AGATCGG at pos 6 matches first 4 bp (4bp overlap).
        // Both yield trim position 6. No other adapter matches.
        assert_eq!(find_best_adapter_match(read, &adapters, 4, 0.1), Some(6));
    }

    #[test]
    fn best_adapter_match_none_match() {
        let adapters = vec![ad(b"TTTT")];
        assert_eq!(find_best_adapter_match(b"AAAACCCC", &adapters, 4, 0.1), None);
    }

    // ---- detect_pe_overlap ----

    #[test]
    fn pe_overlap_short_insert() {
        // Insert = TTTT (4bp), R1 adapter = "AA", R2 adapter = "CC"
        // R1 = TTTTAA (6bp), R2 = AAAACC (6bp; reverse of insert + R2 adapter)
        let r1 = b"TTTTAA";
        let r2 = b"AAAACC";
        // R2_rc = reverse_complement("AAAACC") = "GGTTTT"
        // R1[0..4] = "TTTT", R2_rc[2..6] = "TTTT" → match!
        assert_eq!(detect_overlap(r1, r2, 4, 0.0), Some(4));
    }

    #[test]
    fn pe_overlap_no_overlap_when_insert_longer_than_read() {
        // Reads don't overlap at all in the sense of adapter read-through — but the
        // algorithm will still try. Use reads with incompatible sequences to confirm
        // no match is found.
        let r1 = b"AAAAAAAA";
        let r2 = b"CCCCCCCC";
        assert_eq!(detect_overlap(r1, r2, 4, 0.0), None);
    }

    #[test]
    fn pe_overlap_tolerates_mismatches() {
        // Insert = AAAAT (5), R1 = AAAATCC (7), R2 = ATTTTGG (7)
        // R2_rc = CCAAAAT, R1[0..5] = "AAAAT", R2_rc[2..7] = "AAAAT" → exact
        // Now inject a mismatch:
        let r1 = b"AAAATCCC";
        let r2 = b"AGTTTTGG"; // R2_rc = CCAAAACT; compare to R1[0..5]=AAAAT vs R2_rc[3..8]=AAACT
        // 5-bp overlap with 1 mismatch (T vs C at position 3) → 20% rate
        assert_eq!(detect_overlap(r1, r2, 5, 0.25), Some(5));
        assert_eq!(detect_overlap(r1, r2, 5, 0.1), None);
    }

    #[test]
    fn pe_overlap_respects_min_overlap() {
        let r1 = b"TTTTAA";
        let r2 = b"AAAACC";
        // insert = 4bp; min_overlap 5 rejects it
        assert_eq!(detect_overlap(r1, r2, 5, 0.0), None);
    }

    // ---- AdapterSet build ----

    #[test]
    fn build_adapter_set_from_kit() {
        let set = build_adapter_set(&[], &None, &["truseq".to_string()], 2).unwrap();
        assert!(!set.r1.is_empty());
        assert!(!set.r2.is_empty());
        assert!(set.r1.iter().any(|a| a.bytes == b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"));
        assert!(set.r2.iter().any(|a| a.bytes == b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"));
    }

    #[test]
    fn build_adapter_set_explicit_sequences_lowercase_normalized() {
        let set =
            build_adapter_set(&["agat".to_string(), "cgat".to_string()], &None, &[], 2).unwrap();
        let r1_bytes: Vec<Vec<u8>> = set.r1.iter().map(|a| a.bytes.clone()).collect();
        let r2_bytes: Vec<Vec<u8>> = set.r2.iter().map(|a| a.bytes.clone()).collect();
        assert_eq!(r1_bytes, vec![b"AGAT".to_vec()]);
        assert_eq!(r2_bytes, vec![b"CGAT".to_vec()]);
    }

    #[test]
    fn build_adapter_set_dedupes() {
        // truseq + explicit same seq → one entry each
        let set = build_adapter_set(
            &["AGATCGGAAGAGCACACGTCTGAACTCCAGTCA".to_string()],
            &None,
            &["truseq".to_string()],
            1,
        )
        .unwrap();
        let count =
            set.r1.iter().filter(|a| a.bytes == b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA").count();
        assert_eq!(count, 1);
    }

    #[test]
    fn build_adapter_set_classifies_pure_acgt() {
        // Plain ACGT adapters should be flagged as pure_acgt = true (SIMD path-eligible).
        let set = build_adapter_set(&["AGATCGG".to_string()], &None, &[], 1).unwrap();
        assert!(set.r1.iter().all(|a| a.pure_acgt));
        // IUPAC adapters should be flagged as pure_acgt = false (scalar fallback).
        let set = build_adapter_set(&["AGANCGG".to_string()], &None, &[], 1).unwrap();
        assert!(!set.r1[0].pure_acgt);
    }

    #[test]
    fn build_adapter_set_is_empty_without_configuration() {
        let set = build_adapter_set(&[], &None, &[], 2).unwrap();
        assert!(set.is_empty());
    }

    #[test]
    fn build_adapter_set_small_rna_leaves_r2_empty() {
        // small-rna kit has no R2 adapter; leaking R1 into R2 would mis-trim PE small-RNA
        // libraries that use a distinct R2 adapter.
        let set = build_adapter_set(&[], &None, &["small-rna".to_string()], 2).unwrap();
        assert_eq!(set.r1.len(), 1);
        assert!(set.r2.is_empty());
    }

    #[test]
    fn build_adapter_set_all_unions_every_kit() {
        let set = build_adapter_set(&[], &None, &["all".to_string()], 2).unwrap();
        // TruSeq + Nextera + small-rna + AVITI + MGI gives 5 R1 entries (all distinct).
        assert_eq!(set.r1.len(), 5);
        // R2: TruSeq + Nextera + AVITI + MGI (small-rna contributes nothing to R2).
        assert_eq!(set.r2.len(), 4);
    }

    #[test]
    fn reverse_complement_iupac_codes() {
        // R(A/G) ↔ Y(C/T); S/W are self-complementary; M(A/C) ↔ K(G/T)
        assert_eq!(reverse_complement(b"RMKY"), b"RMKY");
        assert_eq!(reverse_complement(b"ARGC"), b"GCYT");
    }

    // ---- validation ----

    #[test]
    fn validation_rejects_bad_adapter_mismatch_rate() {
        let tmp = TempDir::new().unwrap();
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &["ACGT"]));
        let mut cmd = trim_cmd(vec![r1], vec![tmp.path().join("o.fq.gz")], None);
        cmd.adapter_mismatch_rate = 1.5;
        let err = cmd.validate().unwrap_err().to_string();
        assert!(err.contains("adapter-mismatch-rate"), "{err}");
    }

    #[test]
    fn validation_rejects_unknown_kit() {
        let tmp = TempDir::new().unwrap();
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &["ACGT"]));
        let mut cmd = trim_cmd(vec![r1], vec![tmp.path().join("o.fq.gz")], None);
        cmd.kit = vec!["unicorn".to_string()];
        let err = cmd.validate().unwrap_err().to_string();
        assert!(err.contains("not a recognized preset"), "{err}");
    }

    #[test]
    fn validation_rejects_invalid_adapter_base() {
        let tmp = TempDir::new().unwrap();
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &["ACGT"]));
        let mut cmd = trim_cmd(vec![r1], vec![tmp.path().join("o.fq.gz")], None);
        cmd.adapter_sequence = vec!["ACGZ".to_string()];
        let err = cmd.validate().unwrap_err().to_string();
        assert!(err.contains("invalid base"), "{err}");
    }

    #[test]
    fn validation_rejects_detect_pe_on_single_end() {
        let tmp = TempDir::new().unwrap();
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &["ACGT"]));
        let mut cmd = trim_cmd(vec![r1], vec![tmp.path().join("o.fq.gz")], None);
        cmd.detect_adapter_for_pe = true;
        let err = cmd.validate().unwrap_err().to_string();
        assert!(err.contains("detect-adapter-for-pe"), "{err}");
    }

    // ---- execute: adapter trimming ----

    #[test]
    fn execute_se_trims_explicit_adapter() {
        let tmp = TempDir::new().unwrap();
        let reads: Vec<&str> = vec!["AAAAAAAAAAAGATCGGAAGAG", "CCCCCCCCCCAGATCGGAAGAG"];
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &reads));
        let out = tmp.path().join("out.fq.gz");
        let mut cmd = trim_cmd(vec![r1], vec![out.clone()], None);
        cmd.adapter_sequence = vec!["AGATCGGAAGAG".to_string()];
        cmd.execute().unwrap();

        let written = read_fastq(&out);
        assert_eq!(written.len(), 2);
        assert_eq!(written[0].seq.as_slice(), b"AAAAAAAAAA");
        assert_eq!(written[1].seq.as_slice(), b"CCCCCCCCCC");
    }

    #[test]
    fn execute_se_trims_truseq_kit() {
        let tmp = TempDir::new().unwrap();
        let truseq = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
        let insert = b"TTTTTTTTTT";
        let mut read = Vec::from(&insert[..]);
        read.extend_from_slice(truseq);
        let reads_str = String::from_utf8(read).unwrap();
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &[reads_str.as_str()]));
        let out = tmp.path().join("out.fq.gz");
        let mut cmd = trim_cmd(vec![r1], vec![out.clone()], None);
        cmd.kit = vec!["truseq".to_string()];
        cmd.execute().unwrap();

        let written = read_fastq(&out);
        assert_eq!(written[0].seq.as_slice(), insert);
    }

    #[test]
    fn execute_pe_overlap_trims_both_mates() {
        let tmp = TempDir::new().unwrap();
        // Short insert of 8bp; 4bp adapters on each mate. Both reads are 12bp.
        // Insert = ACGTACGT, R1 adapter = NNNN, R2 adapter = NNNN (use arbitrary letters)
        // R1 = ACGTACGT + AAAA (12bp); R2 = reverse_complement(ACGTACGT) + CCCC = ACGTACGT + CCCC (12bp)
        // Wait: reverse_complement(ACGTACGT) = ACGTACGT (palindrome)
        // R2_rc = reverse_complement(ACGTACGTCCCC) = GGGGACGTACGT
        // R1[0..8] = "ACGTACGT", R2_rc[4..12] = "ACGTACGT" → match!
        let r1 = write_fastq(&tmp, "r1", &fq_lines("p", &["ACGTACGTAAAA"]));
        let r2 = write_fastq(&tmp, "r2", &fq_lines("p", &["ACGTACGTCCCC"]));
        let o1 = tmp.path().join("o1.fq.gz");
        let o2 = tmp.path().join("o2.fq.gz");
        let mut cmd = trim_cmd(vec![r1, r2], vec![o1.clone(), o2.clone()], None);
        cmd.detect_adapter_for_pe = true;
        cmd.overlap_min_length = 5;
        cmd.execute().unwrap();

        let w1 = read_fastq(&o1);
        let w2 = read_fastq(&o2);
        assert_eq!(w1[0].seq.as_slice(), b"ACGTACGT");
        assert_eq!(w2[0].seq.as_slice(), b"ACGTACGT");
    }

    #[test]
    fn execute_pe_overlap_evidence_check_fires_for_custom_adapter() {
        // Exercises the post-cut adapter-evidence confirmation path (which was not
        // covered by the plain PE-overlap test). Supplies a unique R1 / R2 adapter via
        // `--adapter-sequence` so both mates' post-cut bases match a library entry and
        // the Descending walk's evidence check passes.
        //
        // Geometry: insert = 40bp `ACGTACGTACGT...` (length 12, repeated). Adapter
        // read-through = 20 bp beyond the insert on each mate. Reads = 40bp total.
        //   R1 = insert[..20] + adapter_r1[..20]
        //   R2 = revcomp(insert[20..]) + adapter_r2[..20]
        let tmp = TempDir::new().unwrap();
        let insert_r1_half = b"AAAACCCCGGGGTTTTAAAA"; // 20bp
        let insert_r2_half = b"TTTTAAAACCCCGGGGTTTT"; // 20bp, revcomp of insert[20..40] is this padded
        let adapter_r1 = b"CAGCAGATCTCGGTGG"; // 16bp unique, distinct from kit adapters
        let adapter_r2 = b"GGAGATCAGCAGTCGC"; // 16bp unique, distinct from adapter_r1
        let r1_seq: Vec<u8> = [&insert_r1_half[..], &adapter_r1[..], &[b'N'; 4][..]].concat();
        let r2_seq: Vec<u8> = [&insert_r2_half[..], &adapter_r2[..], &[b'N'; 4][..]].concat();
        let r1_path =
            write_fastq(&tmp, "r1", &fq_lines("p", &[std::str::from_utf8(&r1_seq).unwrap()]));
        let r2_path =
            write_fastq(&tmp, "r2", &fq_lines("p", &[std::str::from_utf8(&r2_seq).unwrap()]));
        let o1 = tmp.path().join("o1.fq.gz");
        let o2 = tmp.path().join("o2.fq.gz");
        let mut cmd = trim_cmd(vec![r1_path, r2_path], vec![o1.clone(), o2.clone()], None);
        cmd.detect_adapter_for_pe = true;
        cmd.overlap_min_length = 10;
        cmd.adapter_sequence = vec![
            std::str::from_utf8(adapter_r1).unwrap().to_string(),
            std::str::from_utf8(adapter_r2).unwrap().to_string(),
        ];
        cmd.execute().unwrap();

        // This test does NOT assert exact trim positions — insert content is chosen so
        // overlap is at least detected; the key behavior under test is that the
        // evidence check does not reject a legitimate custom-adapter overlap.
        let w1 = read_fastq(&o1);
        let w2 = read_fastq(&o2);
        assert_eq!(w1.len(), 1);
        assert_eq!(w2.len(), 1);
        // Post-trim length must be strictly less than the 40 input bases (overlap + evidence fired).
        assert!(w1[0].seq.len() < 40, "R1 should be trimmed, got {}", w1[0].seq.len());
        assert!(w2[0].seq.len() < 40, "R2 should be trimmed, got {}", w2[0].seq.len());
    }

    #[test]
    fn execute_adapter_trim_tracks_bases_metric() {
        let tmp = TempDir::new().unwrap();
        // 10bp insert + 12bp adapter = 22bp read; after trim, 10bp remains, 12bp trimmed.
        let reads: Vec<&str> = vec!["AAAAAAAAAAAGATCGGAAGAG"];
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &reads));
        let out = tmp.path().join("out.fq.gz");
        let metrics_path = tmp.path().join("m.txt");
        let mut cmd = trim_cmd(vec![r1], vec![out], Some(metrics_path.clone()));
        cmd.adapter_sequence = vec!["AGATCGGAAGAG".to_string()];
        cmd.execute().unwrap();

        let contents = std::fs::read_to_string(&metrics_path).unwrap();
        let lines: Vec<&str> = contents.lines().collect();
        let header: Vec<&str> = lines[0].split('\t').collect();
        let values: Vec<&str> = lines[1].split('\t').collect();
        let idx = |name: &str| header.iter().position(|h| *h == name).unwrap();
        assert_eq!(values[idx("bases_in")], "22");
        assert_eq!(values[idx("bases_out")], "10");
        assert_eq!(values[idx("bases_trimmed_adapter")], "12");
    }

    #[test]
    fn execute_fasta_adapter_trimming() {
        let tmp = TempDir::new().unwrap();
        let fasta_path = tmp.path().join("adapters.fa");
        let fasta_lines = vec![
            ">adapter1".to_string(),
            "AGATCGGAAGAG".to_string(),
            ">adapter2".to_string(),
            "TGTCTCTTATAC".to_string(),
        ];
        Io::default().write_lines(&fasta_path, &fasta_lines).unwrap();

        // insert + adapter for each read; chosen so the adapter's first base does not
        // match the last base of the insert (avoids the well-known "adapter eats insert"
        // aliasing effect).
        let reads: Vec<&str> = vec!["AAAAAAAAAAAGATCGGAAGAG", "CCCCCCCCCCTGTCTCTTATAC"];
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &reads));
        let out = tmp.path().join("out.fq.gz");
        let mut cmd = trim_cmd(vec![r1], vec![out.clone()], None);
        cmd.adapter_fasta = Some(fasta_path);
        cmd.execute().unwrap();

        let written = read_fastq(&out);
        assert_eq!(written[0].seq.as_slice(), b"AAAAAAAAAA");
        assert_eq!(written[1].seq.as_slice(), b"CCCCCCCCCC");
    }

    #[test]
    fn execute_hard_trim_tracks_bases() {
        let tmp = TempDir::new().unwrap();
        let reads: Vec<&str> = vec!["AAAAATEMPL"]; // 10 bases in, 5 out after 5S+T
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &reads));
        let out = tmp.path().join("out.fq.gz");
        let metrics_path = tmp.path().join("m.txt");
        let mut cmd = trim_cmd(vec![r1], vec![out], Some(metrics_path.clone()));
        cmd.read_structures = vec![rs("5S+T")];
        cmd.execute().unwrap();

        let contents = std::fs::read_to_string(&metrics_path).unwrap();
        let lines: Vec<&str> = contents.lines().collect();
        let header: Vec<&str> = lines[0].split('\t').collect();
        let values: Vec<&str> = lines[1].split('\t').collect();
        let idx = |name: &str| header.iter().position(|h| *h == name).unwrap();
        assert_eq!(values[idx("bases_in")], "10");
        assert_eq!(values[idx("bases_out")], "5");
        assert_eq!(values[idx("bases_trimmed_read_structure")], "5");
    }

    // ---- QualityTrim / LengthFilter parsing ----

    #[test]
    fn quality_trim_parse_ok() {
        let qt: QualityTrim = "4:20".parse().unwrap();
        assert_eq!(qt.window, 4);
        assert_eq!(qt.threshold, 20);
    }

    #[test]
    fn quality_trim_parse_rejects_missing_colon() {
        let err: String = "4".parse::<QualityTrim>().unwrap_err();
        assert!(err.contains("WINDOW:QUAL"), "{err}");
    }

    #[test]
    fn quality_trim_parse_rejects_zero_window() {
        assert!("0:20".parse::<QualityTrim>().is_err());
    }

    #[test]
    fn length_filter_parse_min_only() {
        let lf: LengthFilter = "15".parse().unwrap();
        assert_eq!(lf.min, 15);
        assert!(lf.max.is_none());
    }

    #[test]
    fn length_filter_parse_min_max() {
        let lf: LengthFilter = "15:150".parse().unwrap();
        assert_eq!(lf.min, 15);
        assert_eq!(lf.max, Some(150));
    }

    #[test]
    fn length_filter_parse_rejects_max_less_than_min() {
        let err: String = "50:10".parse::<LengthFilter>().unwrap_err();
        assert!(err.contains("max"), "{err}");
    }

    // ---- find_polyx_tail_len ----

    #[test]
    fn polyx_tail_all_g() {
        assert_eq!(find_polyx_tail_len(b"ACGTGGGGGGGG", b'G'), 8);
    }

    #[test]
    fn polyx_tail_no_g() {
        assert_eq!(find_polyx_tail_len(b"ACGTACGTACGT", b'G'), 0);
    }

    #[test]
    fn polyx_tail_strict_stops_at_first_non_g() {
        // Tail is 7 G's; the T at position -7 halts the run even though more G's appear
        // further 5' in the read.
        assert_eq!(find_polyx_tail_len(b"ACGTGGGGTGGGGGGG", b'G'), 7);
    }

    #[test]
    fn polyx_tail_single_trailing_base() {
        assert_eq!(find_polyx_tail_len(b"AAAAAAG", b'G'), 1);
    }

    #[test]
    fn polyx_tail_empty() {
        assert_eq!(find_polyx_tail_len(b"", b'G'), 0);
    }

    // ---- trim_polyx_tail ----

    #[test]
    fn trim_polyg_tail_applied() {
        let mut rec = owned_rec("r", "ACGTGGGGGGGGGGGGGGGGG", "IIIIIIIIIIIIIIIIIIIII");
        let trimmed = trim_polyx_tail(&mut rec, b'G', 10);
        assert_eq!(trimmed, 17);
        assert_eq!(rec.seq, b"ACGT");
    }

    #[test]
    fn trim_polyg_tail_below_min_run_untouched() {
        let mut rec = owned_rec("r", "ACGTGGGG", "IIIIIIII");
        let trimmed = trim_polyx_tail(&mut rec, b'G', 10);
        assert_eq!(trimmed, 0);
        assert_eq!(rec.seq, b"ACGTGGGG");
    }

    // ---- trim_quality_sliding_3prime ----

    #[test]
    fn quality_sliding_trims_tail() {
        // First 10 bases Q40 ("I"), last 10 Q5 ("&"). Window size 4 at threshold 20.
        // The first "bad" window starts at position 9 ([Q40, Q5, Q5, Q5] mean = 13.75 < 20),
        // so trim occurs at position 9 — consistent with fastp's cut_right / Trimmomatic's
        // SLIDINGWINDOW (they cut at the window's start, not past its end).
        let mut rec = owned_rec("r", "AAAAAAAAAACCCCCCCCCC", "IIIIIIIIII&&&&&&&&&&");
        let trimmed = trim_quality_sliding_3prime(&mut rec, 4, 20);
        assert_eq!(trimmed, 11);
        assert_eq!(rec.seq, b"AAAAAAAAA");
    }

    #[test]
    fn quality_sliding_no_trim_for_high_quality_read() {
        let mut rec = owned_rec("r", "AAAAAAAA", "IIIIIIII");
        let trimmed = trim_quality_sliding_3prime(&mut rec, 4, 20);
        assert_eq!(trimmed, 0);
    }

    #[test]
    fn quality_sliding_shorter_than_window_is_noop() {
        let mut rec = owned_rec("r", "AAA", "III");
        let trimmed = trim_quality_sliding_3prime(&mut rec, 4, 20);
        assert_eq!(trimmed, 0);
    }

    #[test]
    fn quality_sliding_window_equals_read_length() {
        let mut rec = owned_rec("r", "AAAA", "IIII");
        let trimmed = trim_quality_sliding_3prime(&mut rec, 4, 20);
        assert_eq!(trimmed, 0);
    }

    #[test]
    fn quality_sliding_preserves_good_bases_past_mid_dip() {
        // cut_tail semantics: a short low-quality dip IN THE MIDDLE of an otherwise
        // good read must NOT cause the read to be trimmed — only trailing bad bases
        // get removed. (Under the old cut_right semantics we used to have, this
        // read would have been truncated at the dip.)
        //
        // 30 bases total: Q40 (I) × 10, Q5 (&) × 4 dip, Q40 × 16 tail.
        let mut rec =
            owned_rec("r", "AAAAAAAAAACCCCAAAAAAAAAAAAAAAA", "IIIIIIIIII&&&&IIIIIIIIIIIIIIII");
        let trimmed = trim_quality_sliding_3prime(&mut rec, 4, 20);
        assert_eq!(trimmed, 0, "mid-dip must not trim trailing good bases");
        assert_eq!(rec.seq.len(), 30);
    }

    #[test]
    fn quality_sliding_trims_only_trailing_bad_with_mixed_interior() {
        // Good start, bad dip, good middle, bad tail. Cut-tail scans from the 3'
        // end, hits the bad tail first, walks back through it, stops at the good
        // middle — preserving the start and middle and dropping only the bad tail.
        //
        // Layout (20 bases): Q40 × 6, Q5 × 3, Q40 × 7, Q5 × 4.
        //   qual: "IIIIII&&&IIIIIII&&&&"
        //                      ^ bad tail starts at position 16
        //
        // Scanning 3'→5' with window=4 threshold=20:
        //   s=16 [16..20) = Q5*4,              sum=152. FAIL.
        //   s=15 [15..19) = Q40, Q5*3,         sum=187. FAIL.
        //   s=14 [14..18) = Q40*2, Q5*2,       sum=222. PASSES (>= 212).
        //
        // Last bad s = 15. trim_pos = 15. We keep 15 bases and drop the bad tail.
        // The interior Q5 dip at positions 6-8 is NEVER scanned (we stop at the
        // first good window encountered from the 3' end), so it doesn't influence
        // the trim — exactly the cut_tail property we want.
        let mut rec = owned_rec("r", "AAAAAACCCAAAAAAACCCC", "IIIIII&&&IIIIIII&&&&");
        let trimmed = trim_quality_sliding_3prime(&mut rec, 4, 20);
        assert_eq!(trimmed, 5);
        assert_eq!(rec.seq.len(), 15);
    }

    #[test]
    fn quality_sliding_window_8_trims_tail() {
        // New default window is 8. Q40 × 40, Q5 × 10 tail. Threshold 20.
        // Scanning 3'→5' with window=8: the tail windows all fail until we reach a
        // window spanning enough Q40 bases. Specifically, window [32..40) is Q40*8
        // mean 40, passing. So trim is somewhere between s=32 and s=40.
        //   s=42 window [42..50) all Q5 fails.
        //   ...
        //   s=33 [33..41) Q40*7,Q5 sum=7*73+38=549; 549 >= 8*53=424 YES, passes.
        // Last bad s = 34 (window [34..42) has 6 Q40 + 2 Q5; sum=6*73+2*38=514 >= 424, actually passes).
        // Let me compute: at s=34, [34..42) = Q40 at 34..39 and Q5 at 40,41.
        //    Q40*6 + Q5*2 = 6*73+2*38 = 438+76 = 514. 514 >= 424? Yes. PASSES.
        // At s=35 [35..43) = Q40 at 35..39, Q5 at 40..42: 5*73+3*38=365+114=479. >= 424. PASSES.
        // s=36 [36..44): 4*73+4*38=292+152=444. PASSES.
        // s=37 [37..45): 3*73+5*38=219+190=409. FAILS.
        // So last bad s = 37. trim_pos = 37. Keep 37 bases.
        let seq = "A".repeat(50);
        let qual = format!("{}{}", "I".repeat(40), "&".repeat(10));
        let mut rec = owned_rec("r", &seq, &qual);
        let trimmed = trim_quality_sliding_3prime(&mut rec, 8, 20);
        assert_eq!(trimmed, 13);
        assert_eq!(rec.seq.len(), 37);
    }

    #[test]
    fn quality_sliding_all_bad_drops_everything() {
        // Whole read under threshold: trim_pos should end up at 0 (empty output).
        let mut rec = owned_rec("r", "AAAAAAAAAAAAAAAAAAAA", "&&&&&&&&&&&&&&&&&&&&");
        let trimmed = trim_quality_sliding_3prime(&mut rec, 4, 20);
        assert_eq!(trimmed, 20);
        assert_eq!(rec.seq.len(), 0);
    }

    #[test]
    fn quality_sliding_scalar_matches_simd_on_random_data() {
        // SIMD fast path requires >= window + 16 bytes. Generate reads long enough to
        // cross that threshold and compare against a reference scalar implementation
        // that matches cut-tail semantics explicitly. Runs across a window sweep so
        // a bug specific to a single window size (e.g. the CLI default of 8) doesn't
        // slip past. Uses xorshift to avoid bringing in an rng dep.
        fn reference(qual: &[u8], window: usize, threshold: u8) -> usize {
            let mut trim_pos = qual.len();
            if qual.len() < window {
                return trim_pos;
            }
            let max_s = qual.len() - window;
            let threshold_total = (threshold as u32) * (window as u32);
            let mut s = max_s as i64;
            while s >= 0 {
                let start = s as usize;
                let sum: u32 = qual[start..start + window]
                    .iter()
                    .map(|&q| u32::from(q.saturating_sub(33)))
                    .sum();
                if sum < threshold_total {
                    trim_pos = start;
                    s -= 1;
                } else {
                    break;
                }
            }
            trim_pos
        }

        let mut rng: u64 = 0xC0FFEE1234CAFE99;
        // Exercise window=1 (edge), 4 (fastp default), 8 (our default), 16 (SIMD cap).
        for &window in &[1usize, 4, 8, 16] {
            for _ in 0..100 {
                // xorshift
                rng ^= rng << 13;
                rng ^= rng >> 7;
                rng ^= rng << 17;
                // Length well above window + 16 so the SIMD path exercises at least one
                // chunk plus a scalar-tail trip.
                let len = window + 20 + (rng as usize % 25);
                let mut qual: Vec<u8> = Vec::with_capacity(len);
                for _ in 0..len {
                    rng ^= rng << 13;
                    rng ^= rng >> 7;
                    rng ^= rng << 17;
                    // Phred+33 range: 33..=73 (Q0..Q40). Bias toward Q40 slightly.
                    qual.push(33 + (rng as u8 % 41));
                }
                let seq = vec![b'A'; len];
                let mut rec = owned_rec(
                    "r",
                    std::str::from_utf8(&seq).unwrap(),
                    std::str::from_utf8(&qual).unwrap(),
                );
                let expected_trim_pos = reference(&qual, window, 20);
                let _trimmed = trim_quality_sliding_3prime(&mut rec, window, 20);
                assert_eq!(
                    rec.seq.len(),
                    expected_trim_pos,
                    "trim_pos mismatch on window={window} len={len} qual={qual:?}"
                );
            }
        }
    }

    // ---- observe_stats (SIMD stats kernel) ----

    /// Scalar reference implementation matching the semantics of `observe_stats`, used
    /// to cross-check the SIMD kernel across boundaries and random inputs.
    fn observe_stats_scalar(seq: &[u8], qual: &[u8]) -> BaseStats {
        const PHRED33: u8 = 33;
        assert_eq!(seq.len(), qual.len());
        let mut s = BaseStats { total: qual.len() as u64, ..Default::default() };
        for (&q, &b) in qual.iter().zip(seq.iter()) {
            let phred = q.saturating_sub(PHRED33);
            if phred >= 20 {
                s.q20 += 1;
            }
            if phred >= 30 {
                s.q30 += 1;
            }
            if b == b'N' || b == b'n' {
                s.n_bases += 1;
            }
        }
        s
    }

    #[test]
    fn observe_stats_empty_is_default() {
        assert_eq!(observe_stats(b"", b""), BaseStats::default());
    }

    #[test]
    fn observe_stats_single_byte() {
        // b'5' = 53 = Q20 exactly; b'?' = 63 = Q30 exactly.
        assert_eq!(observe_stats(b"A", b"5"), BaseStats { total: 1, q20: 1, q30: 0, n_bases: 0 });
        assert_eq!(observe_stats(b"N", b"?"), BaseStats { total: 1, q20: 1, q30: 1, n_bases: 1 });
        assert_eq!(observe_stats(b"n", b"!"), BaseStats { total: 1, q20: 0, q30: 0, n_bases: 1 });
    }

    #[test]
    fn observe_stats_exact_chunk_boundary() {
        // Exactly 16 bytes: no scalar tail. All Q40 (I = 73 = Q40).
        let seq = b"ACGTACGTNNNNACGT";
        let qual = b"IIIIIIIIIIIIIIII";
        let s = observe_stats(seq, qual);
        assert_eq!(s.total, 16);
        assert_eq!(s.q20, 16);
        assert_eq!(s.q30, 16);
        assert_eq!(s.n_bases, 4);
    }

    #[test]
    fn observe_stats_one_past_chunk() {
        // 17 bytes: 16 chunk + 1 scalar tail.
        let seq = b"ACGTACGTNNNNACGTN";
        let qual = b"IIIIIIIIIIIIIIII!"; // last byte = Q0
        let s = observe_stats(seq, qual);
        assert_eq!(s.total, 17);
        assert_eq!(s.q20, 16);
        assert_eq!(s.q30, 16);
        assert_eq!(s.n_bases, 5);
    }

    #[test]
    fn observe_stats_sub_chunk_length() {
        // 5 bytes: all in scalar tail.
        let seq = b"NNACG";
        let qual = b"!!III";
        let s = observe_stats(seq, qual);
        assert_eq!(s, BaseStats { total: 5, q20: 3, q30: 3, n_bases: 2 });
    }

    #[test]
    fn observe_stats_pe150_length_matches_scalar() {
        // PE150-ish synthetic read with a mix of bases and quality values spanning Q0–Q40.
        let seq: Vec<u8> =
            (0..150).map(|i| if i % 13 == 0 { b'N' } else { b"ACGT"[i % 4] }).collect();
        let qual: Vec<u8> = (0..150).map(|i| 33 + ((i * 7 + 3) % 41) as u8).collect();
        assert_eq!(observe_stats(&seq, &qual), observe_stats_scalar(&seq, &qual));
    }

    #[test]
    fn observe_stats_random_lengths_match_scalar() {
        // Smoke test across lengths that straddle the chunk boundary in all alignments.
        // Uses a simple LCG so the sequence is reproducible without bringing in `rand`.
        let mut state: u64 = 0x1234_5678_9abc_def0;
        let mut rng = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            state
        };
        for len in 0..=200 {
            let seq: Vec<u8> = (0..len).map(|_| b"ACGTNn"[(rng() % 6) as usize]).collect();
            let qual: Vec<u8> = (0..len).map(|_| 33 + (rng() % 45) as u8).collect();
            assert_eq!(
                observe_stats(&seq, &qual),
                observe_stats_scalar(&seq, &qual),
                "mismatch at len={len}"
            );
        }
    }

    // ---- count_mismatches_ci_bounded ----

    /// Scalar reference for cross-checking the SIMD kernel.
    fn ci_mismatches_scalar(a: &[u8], b: &[u8], limit: usize) -> usize {
        assert_eq!(a.len(), b.len());
        let mut count = 0;
        for (&x, &y) in a.iter().zip(b.iter()) {
            if !x.eq_ignore_ascii_case(&y) {
                count += 1;
                if count > limit {
                    return count;
                }
            }
        }
        count
    }

    #[test]
    fn ci_bounded_identical_returns_zero() {
        let a = b"ACGTACGTACGTACGT";
        assert_eq!(count_mismatches_ci_bounded(a, a, 0), 0);
        assert_eq!(count_mismatches_ci_bounded(a, a, 10), 0);
    }

    #[test]
    fn ci_bounded_case_insensitive() {
        let a = b"ACGTACGTACGTACGT";
        let b = b"acgtacgtacgtacgt";
        assert_eq!(count_mismatches_ci_bounded(a, b, 0), 0);
    }

    #[test]
    fn ci_bounded_all_mismatch_stops_early() {
        // 16 all-mismatches, limit 3 — kernel should return immediately with >3.
        let a = b"AAAAAAAAAAAAAAAA";
        let b = b"TTTTTTTTTTTTTTTT";
        let mm = count_mismatches_ci_bounded(a, b, 3);
        assert!(mm > 3, "got {mm}");
    }

    #[test]
    fn ci_bounded_empty_is_zero() {
        assert_eq!(count_mismatches_ci_bounded(b"", b"", 0), 0);
    }

    #[test]
    fn ci_bounded_sub_chunk_tail() {
        // Length 5 — only scalar tail runs.
        let a = b"AAGCT";
        let b = b"ATGCT"; // one mismatch at position 1
        assert_eq!(count_mismatches_ci_bounded(a, b, 5), 1);
        let mm = count_mismatches_ci_bounded(a, b, 0);
        assert!(mm > 0);
    }

    #[test]
    fn ci_bounded_boundary_lengths_match_scalar() {
        // Cover 0..=33 so we hit pre-chunk, full-chunk, and chunk+tail configurations.
        let base_a = b"ACGTACGTACGTACGTNNNNAAAAACCCCTGG";
        let base_b = b"ACGTACGTacgTACGTnnnnAAAATCCCCTGG"; // differs at pos 8 (T vs T) hmm let me think
        // Actually just use two similar-but-differing strings:
        for len in 0..=base_a.len() {
            let a = &base_a[..len];
            let b = &base_b[..len];
            for limit in [0, 1, 3, 5, 16, 64] {
                assert_eq!(
                    count_mismatches_ci_bounded(a, b, limit),
                    ci_mismatches_scalar(a, b, limit),
                    "len={len} limit={limit}"
                );
            }
        }
    }

    // ---- reverse_complement_acgt_into ----

    fn rc_acgt(seq: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        reverse_complement_acgt_into(seq, &mut out);
        out
    }

    #[test]
    fn rc_acgt_empty() {
        assert_eq!(rc_acgt(b""), b"");
    }

    #[test]
    fn rc_acgt_palindrome() {
        assert_eq!(rc_acgt(b"ACGT"), b"ACGT");
    }

    #[test]
    fn rc_acgt_simple_uppercase() {
        assert_eq!(rc_acgt(b"AAAA"), b"TTTT");
        assert_eq!(rc_acgt(b"CCCC"), b"GGGG");
        assert_eq!(rc_acgt(b"ACCGTTT"), b"AAACGGT");
    }

    #[test]
    fn rc_acgt_preserves_n() {
        assert_eq!(rc_acgt(b"ACNGT"), b"ACNGT");
    }

    #[test]
    fn rc_acgt_preserves_lowercase() {
        assert_eq!(rc_acgt(b"aaaa"), b"tttt");
        // aCgT is NOT a palindrome under RC with case preserved:
        //   aCgT reversed = TgCa; complementing each base gives AcGt.
        assert_eq!(rc_acgt(b"aCgT"), b"AcGt");
        assert_eq!(rc_acgt(b"nACG"), b"CGTn");
    }

    #[test]
    fn rc_acgt_exact_chunk_boundary() {
        // 16 bytes, full chunk, no tail.
        let seq = b"AAAACCCCGGGGTTTT";
        assert_eq!(rc_acgt(seq), b"AAAACCCCGGGGTTTT");
    }

    #[test]
    fn rc_acgt_one_past_chunk() {
        // 17 bytes: chunk + 1-byte tail.
        let seq = b"AAAACCCCGGGGTTTTN";
        assert_eq!(rc_acgt(seq), b"NAAAACCCCGGGGTTTT");
    }

    #[test]
    fn rc_acgt_thirty_three_bytes() {
        // Two full chunks + 1-byte tail, ensures chunk ordering in output is correct.
        let seq = b"AAAACCCCGGGGTTTTACGTACGTACGTACGTN";
        // Expected: reverse-complement of the whole string.
        // rc(AAAACCCCGGGGTTTTACGTACGTACGTACGTN)
        // = rc(N) + rc(ACGTACGTACGTACGT) + rc(AAAACCCCGGGGTTTT)
        // = N + ACGTACGTACGTACGT + AAAACCCCGGGGTTTT
        assert_eq!(rc_acgt(seq), b"NACGTACGTACGTACGTAAAACCCCGGGGTTTT");
    }

    #[test]
    fn rc_acgt_iupac_lossy_as_documented() {
        // IUPAC ambiguity codes are intentionally mapped to N in this specialized
        // kernel. Callers that need IUPAC fidelity use reverse_complement_into.
        assert_eq!(rc_acgt(b"R"), b"N");
    }

    #[test]
    fn rc_acgt_matches_scalar_on_acgt_n_inputs() {
        // On the supported alphabet, the SIMD kernel and the full scalar implementation
        // produce identical output (the scalar also keeps ACGT/N stable under RC).
        let inputs: &[&[u8]] = &[
            b"",
            b"A",
            b"AC",
            b"ACG",
            b"ACGT",
            b"ACGTN",
            b"ACGTACGTAC",
            b"ACGTACGTACGTACGT",
            b"ACGTACGTACGTACGTA",
            b"NNNNAAAACCCCGGGGTTTTacgt",
        ];
        for input in inputs {
            assert_eq!(rc_acgt(input), reverse_complement(input), "input = {:?}", input);
        }
    }

    // ---- evaluate_filters ----

    #[test]
    fn filters_length_reject_min() {
        let recs = vec![owned_rec("r", "AAAA", "IIII")];
        let res = eval_filters(&recs, LengthFilter { min: 10, max: None }, None, None);
        assert_eq!(res, Some(FilterReject::Length));
    }

    #[test]
    fn filters_length_reject_max() {
        let recs = vec![owned_rec("r", "AAAAAAAAAA", "IIIIIIIIII")];
        let res = eval_filters(&recs, LengthFilter { min: 0, max: Some(5) }, None, None);
        assert_eq!(res, Some(FilterReject::Length));
    }

    #[test]
    fn filters_n_base_rejects_over_limit() {
        let recs = vec![owned_rec("r", "ANNNNN", "IIIIII")];
        let res = eval_filters(&recs, LengthFilter { min: 0, max: None }, Some(2), None);
        assert_eq!(res, Some(FilterReject::NBases));
    }

    #[test]
    fn filters_mean_quality_rejects_low() {
        // "!" = 33 = Q0; mean is 0 which is below any threshold >= 1
        let recs = vec![owned_rec("r", "AAAA", "!!!!")];
        let res = eval_filters(&recs, LengthFilter { min: 0, max: None }, None, Some(20));
        assert_eq!(res, Some(FilterReject::Quality));
    }

    #[test]
    fn filters_none_passes() {
        let recs = vec![owned_rec("r", "AAAAAAAAAA", "IIIIIIIIII")];
        let res = eval_filters(&recs, LengthFilter { min: 5, max: Some(20) }, Some(5), Some(20));
        assert!(res.is_none());
    }

    #[test]
    fn filters_pe_any_mate_fails_rejects() {
        let recs = vec![owned_rec("r", "AAAAAAAAAA", "IIIIIIIIII"), owned_rec("r", "AAA", "III")];
        let res = eval_filters(&recs, LengthFilter { min: 5, max: None }, None, None);
        assert_eq!(res, Some(FilterReject::Length));
    }

    // ---- execute end-to-end for filter stages ----

    #[test]
    fn execute_polyg_trims_tail_when_enabled() {
        let tmp = TempDir::new().unwrap();
        let lines = vec![
            "@NB501234:1:FLOWCELL:1:1101:1000:2000".to_string(),
            "ACGTGGGGGGGGGGGGGGGGG".to_string(), // 4 real bases + 17 G's
            "+".to_string(),
            "I".repeat(21),
        ];
        let r1 = write_fastq(&tmp, "r1", &lines);
        let out = tmp.path().join("out.fq.gz");
        let cmd = trim_cmd(vec![r1], vec![out.clone()], None);
        // trim_cmd defaults trim_polyg = 0 (off); enable with the standard min-run of 10.
        let mut cmd = cmd;
        cmd.trim_polyg = 10;
        cmd.execute().unwrap();

        let written = read_fastq(&out);
        assert_eq!(written[0].seq.as_slice(), b"ACGT");
    }

    #[test]
    fn execute_length_filter_drops_short_reads() {
        let tmp = TempDir::new().unwrap();
        let reads: Vec<&str> = vec!["ACGTACGTACGT", "AAA", "CCCCCCCCCC"];
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &reads));
        let out = tmp.path().join("out.fq.gz");
        let metrics = tmp.path().join("m.txt");
        let mut cmd = trim_cmd(vec![r1], vec![out.clone()], Some(metrics.clone()));
        cmd.filter_length = LengthFilter { min: 5, max: None };
        cmd.execute().unwrap();

        let written = read_fastq(&out);
        assert_eq!(written.len(), 2);
        assert_eq!(written[0].seq.as_slice(), b"ACGTACGTACGT");
        assert_eq!(written[1].seq.as_slice(), b"CCCCCCCCCC");

        let contents = std::fs::read_to_string(&metrics).unwrap();
        let lines: Vec<&str> = contents.lines().collect();
        let header: Vec<&str> = lines[0].split('\t').collect();
        let values: Vec<&str> = lines[1].split('\t').collect();
        let idx = |name: &str| header.iter().position(|h| *h == name).unwrap();
        assert_eq!(values[idx("reads_in")], "3");
        assert_eq!(values[idx("reads_out")], "2");
        assert_eq!(values[idx("reads_filtered_length")], "1");
    }

    #[test]
    fn execute_pe_one_mate_short_drops_both() {
        let tmp = TempDir::new().unwrap();
        // Two pairs: first pair R2 is short; second pair both OK.
        let r1_lines = vec![
            "@pair1".to_string(),
            "AAAAAAAAAA".to_string(),
            "+".to_string(),
            "IIIIIIIIII".to_string(),
            "@pair2".to_string(),
            "CCCCCCCCCC".to_string(),
            "+".to_string(),
            "IIIIIIIIII".to_string(),
        ];
        let r2_lines = vec![
            "@pair1".to_string(),
            "TTT".to_string(),
            "+".to_string(),
            "III".to_string(),
            "@pair2".to_string(),
            "GGGGGGGGGG".to_string(),
            "+".to_string(),
            "IIIIIIIIII".to_string(),
        ];
        let r1 = write_fastq(&tmp, "r1", &r1_lines);
        let r2 = write_fastq(&tmp, "r2", &r2_lines);
        let o1 = tmp.path().join("o1.fq.gz");
        let o2 = tmp.path().join("o2.fq.gz");
        let mut cmd = trim_cmd(vec![r1, r2], vec![o1.clone(), o2.clone()], None);
        cmd.filter_length = LengthFilter { min: 5, max: None };
        cmd.execute().unwrap();

        // pair1 is dropped entirely; only pair2 survives.
        let w1 = read_fastq(&o1);
        let w2 = read_fastq(&o2);
        assert_eq!(w1.len(), 1);
        assert_eq!(w2.len(), 1);
        assert_eq!(w1[0].seq.as_slice(), b"CCCCCCCCCC");
        assert_eq!(w2[0].seq.as_slice(), b"GGGGGGGGGG");
    }

    #[test]
    fn execute_quality_trim_sliding_window() {
        let tmp = TempDir::new().unwrap();
        let lines = vec![
            "@read1".to_string(),
            "AAAAAAAAAACCCCCCCCCC".to_string(),
            "+".to_string(),
            "IIIIIIIIII&&&&&&&&&&".to_string(), // first 10 Q40, last 10 Q5
        ];
        let r1 = write_fastq(&tmp, "r1", &lines);
        let out = tmp.path().join("out.fq.gz");
        let mut cmd = trim_cmd(vec![r1], vec![out.clone()], None);
        cmd.quality_trim_3p = Some(QualityTrim { window: 4, threshold: 20 });
        cmd.execute().unwrap();

        let written = read_fastq(&out);
        // Sliding window cuts at the first bad window's start (position 9).
        assert_eq!(written[0].seq.as_slice(), b"AAAAAAAAA");
    }

    #[test]
    fn execute_mean_qual_filter_drops_low_quality_reads() {
        let tmp = TempDir::new().unwrap();
        let lines = vec![
            "@good".to_string(),
            "AAAAAAAAAA".to_string(),
            "+".to_string(),
            "IIIIIIIIII".to_string(),
            "@bad".to_string(),
            "CCCCCCCCCC".to_string(),
            "+".to_string(),
            "!!!!!!!!!!".to_string(), // Q0
        ];
        let r1 = write_fastq(&tmp, "r1", &lines);
        let out = tmp.path().join("out.fq.gz");
        let mut cmd = trim_cmd(vec![r1], vec![out.clone()], None);
        cmd.filter_mean_qual = Some(20);
        cmd.execute().unwrap();

        let written = read_fastq(&out);
        assert_eq!(written.len(), 1);
        assert_eq!(written[0].seq.as_slice(), b"AAAAAAAAAA");
    }

    #[test]
    fn execute_base_accounting_identity_holds() {
        let tmp = TempDir::new().unwrap();
        let lines = vec![
            "@A:1:B:1:1:1:1".to_string(),
            "AAAAAAAAAAACGTACGTAGATCGGAAGAG".to_string(),
            "+".to_string(),
            "I".repeat(30),
            "@A:1:B:1:1:1:2".to_string(),
            "CCCCCCCC".to_string(),
            "+".to_string(),
            "IIIIIIII".to_string(),
        ];
        let r1 = write_fastq(&tmp, "r1", &lines);
        let out = tmp.path().join("out.fq.gz");
        let metrics = tmp.path().join("m.txt");
        let mut cmd = trim_cmd(vec![r1], vec![out.clone()], Some(metrics.clone()));
        cmd.read_structures = vec![rs("5S+T")];
        cmd.adapter_sequence = vec!["AGATCGGAAGAG".to_string()];
        cmd.filter_length = LengthFilter { min: 5, max: None };
        cmd.execute().unwrap();

        let contents = std::fs::read_to_string(&metrics).unwrap();
        let lines: Vec<&str> = contents.lines().collect();
        let header: Vec<&str> = lines[0].split('\t').collect();
        let values: Vec<&str> = lines[1].split('\t').collect();
        let get = |name: &str| -> u64 {
            let idx = header.iter().position(|h| *h == name).unwrap();
            values[idx].parse().unwrap()
        };
        let bases_in = get("bases_in");
        let bases_out = get("bases_out");
        let sum_trimmed = get("bases_trimmed_read_structure")
            + get("bases_trimmed_adapter")
            + get("bases_trimmed_polyg")
            + get("bases_trimmed_polyx")
            + get("bases_trimmed_quality");
        let filtered = get("bases_filtered");
        assert_eq!(bases_in, bases_out + sum_trimmed + filtered);
    }

    #[test]
    fn execute_polyx_picks_longest_single_base_tail() {
        // `NNNNAAATTT`: pure-A tail has 0 trailing A's; pure-T tail has 3 trailing T's.
        // Single-pass trimming removes only the T's and avoids cascading.
        let tmp = TempDir::new().unwrap();
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &["NNNNAAATTT"]));
        let out = tmp.path().join("out.fq.gz");
        let mut cmd = trim_cmd(vec![r1], vec![out.clone()], None);
        cmd.trim_polyx = Some(3);
        cmd.execute().unwrap();

        let written = read_fastq(&out);
        assert_eq!(written[0].seq.as_slice(), b"NNNNAAA");
    }

    #[test]
    fn execute_n_base_filter() {
        let tmp = TempDir::new().unwrap();
        let reads: Vec<&str> = vec!["ACGTACGT", "NNNNNNNN"];
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &reads));
        let out = tmp.path().join("out.fq.gz");
        let mut cmd = trim_cmd(vec![r1], vec![out.clone()], None);
        cmd.filter_max_ns = Some(2);
        cmd.execute().unwrap();

        let written = read_fastq(&out);
        assert_eq!(written.len(), 1);
        assert_eq!(written[0].seq.as_slice(), b"ACGTACGT");
    }

    // ---- JSON report ----

    #[test]
    fn execute_writes_json_report_with_fastp_schema_keys() {
        let tmp = TempDir::new().unwrap();
        let lines = vec![
            "@NB501234:1:FC:1:1:1:1".to_string(),
            "AAAAAAAAAA".to_string(),
            "+".to_string(),
            "IIIIIIIIII".to_string(),
            "@NB501234:1:FC:1:1:1:2".to_string(),
            "CCCCCCCCCC".to_string(),
            "+".to_string(),
            "IIIIIIIIII".to_string(),
        ];
        let r1 = write_fastq(&tmp, "r1", &lines);
        let out = tmp.path().join("out.fq.gz");
        let json_path = tmp.path().join("trim.json");
        let mut cmd = trim_cmd(vec![r1], vec![out], None);
        cmd.json = Some(json_path.clone());
        cmd.execute().unwrap();

        let json: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&json_path).unwrap()).unwrap();

        // Top-level keys MultiQC's fastp module looks for
        assert!(json.get("summary").is_some());
        assert!(json.get("filtering_result").is_some());
        assert!(json.get("adapter_cutting").is_some());
        assert!(json.get("read1_before_filtering").is_some());
        assert!(json.get("read1_after_filtering").is_some());
        // Single-end: R2 keys should be absent (serialized with skip_serializing_if).
        assert!(json.get("read2_before_filtering").is_none());
        assert!(json.get("read2_after_filtering").is_none());

        let summary = &json["summary"];
        assert_eq!(summary["sequencing"], "single end");
        assert_eq!(summary["before_filtering"]["total_reads"], 2);
        assert_eq!(summary["before_filtering"]["total_bases"], 20);
        // All Q40 inputs → q20 and q30 rates are 100%.
        assert_eq!(summary["before_filtering"]["q20_rate"], 1.0);
        assert_eq!(summary["before_filtering"]["q30_rate"], 1.0);
    }

    #[test]
    fn execute_pe_json_includes_read2() {
        let tmp = TempDir::new().unwrap();
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &["AAAAAAAA"]));
        let r2 = write_fastq(&tmp, "r2", &fq_lines("r", &["TTTTTTTT"]));
        let o1 = tmp.path().join("o1.fq.gz");
        let o2 = tmp.path().join("o2.fq.gz");
        let json_path = tmp.path().join("trim.json");
        let mut cmd = trim_cmd(vec![r1, r2], vec![o1, o2], None);
        cmd.json = Some(json_path.clone());
        cmd.execute().unwrap();

        let json: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&json_path).unwrap()).unwrap();
        assert_eq!(json["summary"]["sequencing"], "paired end");
        assert_eq!(json["read1_before_filtering"]["total_bases"], 8);
        assert_eq!(json["read2_before_filtering"]["total_bases"], 8);
        // Summary totals should sum across both mates.
        assert_eq!(json["summary"]["before_filtering"]["total_bases"], 16);
    }

    #[test]
    fn execute_zero_reads_produces_finite_json() {
        // Empty FASTQ: the JSON must still be valid (no NaN, no division by zero) so
        // MultiQC can parse it.
        let tmp = TempDir::new().unwrap();
        let r1 = write_fastq(&tmp, "r1", &[]);
        let out = tmp.path().join("out.fq.gz");
        let json_path = tmp.path().join("trim.json");
        let mut cmd = trim_cmd(vec![r1], vec![out], None);
        cmd.json = Some(json_path.clone());
        cmd.execute().unwrap();

        let raw = std::fs::read_to_string(&json_path).unwrap();
        // Trivial NaN sanity: JSON does not represent NaN so a properly finite output
        // parses cleanly as `serde_json::Value`.
        assert!(!raw.contains("NaN"));
        let json: serde_json::Value = serde_json::from_str(&raw).unwrap();
        assert_eq!(json["summary"]["before_filtering"]["total_reads"], 0);
        assert_eq!(json["summary"]["before_filtering"]["q20_rate"], 0.0);
    }

    #[test]
    fn execute_adapter_trimming_populates_adapter_cutting() {
        let tmp = TempDir::new().unwrap();
        let reads: Vec<&str> = vec!["AAAAAAAAAAAGATCGGAAGAG", "CCCCCCCCCC"];
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &reads));
        let out = tmp.path().join("out.fq.gz");
        let json_path = tmp.path().join("trim.json");
        let mut cmd = trim_cmd(vec![r1], vec![out], None);
        cmd.adapter_sequence = vec!["AGATCGGAAGAG".to_string()];
        cmd.json = Some(json_path.clone());
        cmd.execute().unwrap();

        let json: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&json_path).unwrap()).unwrap();
        // Only the first read has adapter trimmed (12 bases).
        assert_eq!(json["adapter_cutting"]["adapter_trimmed_reads"], 1);
        assert_eq!(json["adapter_cutting"]["adapter_trimmed_bases"], 12);
    }

    #[test]
    fn execute_json_rates_are_fractions_not_percentages() {
        // Quality string "I?" = Q40, Q30 → both pass Q20, only one passes Q30.
        // Actually "?" is ASCII 63, Phred = 30 → both pass Q30.
        // Use "I" (Q40), "5" (Q20 exactly ASCII 53 Phred 20 → Q30 fails), for an easier mix.
        // "I" = 73 → Q40, "5" = 53 → Q20, "!" = 33 → Q0. With sequence "AAAA" and quals
        // "I5I5", we get two Q40 bases and two Q20 bases → all 4 Q20, 2 Q30 → q30_rate=0.5.
        let tmp = TempDir::new().unwrap();
        let lines =
            vec!["@read".to_string(), "AAAA".to_string(), "+".to_string(), "I5I5".to_string()];
        let r1 = write_fastq(&tmp, "r1", &lines);
        let out = tmp.path().join("out.fq.gz");
        let json_path = tmp.path().join("trim.json");
        let mut cmd = trim_cmd(vec![r1], vec![out], None);
        cmd.json = Some(json_path.clone());
        cmd.execute().unwrap();

        let json: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&json_path).unwrap()).unwrap();
        let q20 = json["summary"]["before_filtering"]["q20_rate"].as_f64().unwrap();
        let q30 = json["summary"]["before_filtering"]["q30_rate"].as_f64().unwrap();
        assert!((0.0..=1.0).contains(&q20), "q20_rate {q20} out of range");
        assert!((0.0..=1.0).contains(&q30), "q30_rate {q30} out of range");
        assert!((q20 - 1.0).abs() < 1e-9, "q20_rate should be 1.0, got {q20}");
        assert!((q30 - 0.5).abs() < 1e-9, "q30_rate should be 0.5, got {q30}");
    }

    #[test]
    fn execute_json_has_fastp_version_and_valid_json() {
        let tmp = TempDir::new().unwrap();
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &["ACGT"]));
        let out = tmp.path().join("out.fq.gz");
        let json_path = tmp.path().join("trim.json");
        let mut cmd = trim_cmd(vec![r1], vec![out], None);
        cmd.json = Some(json_path.clone());
        cmd.execute().unwrap();

        let raw = std::fs::read_to_string(&json_path).unwrap();
        assert!(!raw.contains("NaN"));
        assert!(!raw.contains("Infinity"));
        let json: serde_json::Value = serde_json::from_str(&raw).unwrap();
        let version = json["summary"]["fastp_version"].as_str().unwrap();
        assert!(!version.is_empty(), "fastp_version should be populated");
    }

    #[test]
    fn execute_pe_adapter_trimmed_counts_per_mate() {
        // Both mates have a trimmable adapter. The counter must increment by 2, matching
        // fastp's convention of counting individual reads.
        let tmp = TempDir::new().unwrap();
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &["AAAAAAAAAAAGATCGGAAGAG"]));
        let r2 = write_fastq(&tmp, "r2", &fq_lines("r", &["CCCCCCCCCCAGATCGGAAGAG"]));
        let o1 = tmp.path().join("o1.fq.gz");
        let o2 = tmp.path().join("o2.fq.gz");
        let json_path = tmp.path().join("trim.json");
        let mut cmd = trim_cmd(vec![r1, r2], vec![o1, o2], None);
        cmd.adapter_sequence = vec!["AGATCGGAAGAG".to_string(), "AGATCGGAAGAG".to_string()];
        cmd.json = Some(json_path.clone());
        cmd.execute().unwrap();

        let json: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&json_path).unwrap()).unwrap();
        assert_eq!(json["adapter_cutting"]["adapter_trimmed_reads"], 2);
    }

    #[test]
    fn execute_filter_populates_filtering_result() {
        let tmp = TempDir::new().unwrap();
        let reads: Vec<&str> = vec!["ACGTACGTACGT", "AAA", "CCCCCCCCCC"];
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &reads));
        let out = tmp.path().join("out.fq.gz");
        let json_path = tmp.path().join("trim.json");
        let mut cmd = trim_cmd(vec![r1], vec![out], None);
        cmd.filter_length = LengthFilter { min: 5, max: None };
        cmd.json = Some(json_path.clone());
        cmd.execute().unwrap();

        let json: serde_json::Value =
            serde_json::from_str(&std::fs::read_to_string(&json_path).unwrap()).unwrap();
        assert_eq!(json["filtering_result"]["passed_filter_reads"], 2);
        assert_eq!(json["filtering_result"]["too_short_reads"], 1);
        assert_eq!(json["filtering_result"]["low_quality_reads"], 0);
    }

    // ---- flatten_mate_stats ----

    #[test]
    fn flatten_mate_stats_paired() {
        let before = vec![
            MateStats { reads: 10, bases: 100, q20_bases: 90, q30_bases: 80 },
            MateStats { reads: 10, bases: 100, q20_bases: 85, q30_bases: 70 },
        ];
        let after = vec![
            MateStats { reads: 9, bases: 90, q20_bases: 88, q30_bases: 80 },
            MateStats { reads: 9, bases: 90, q20_bases: 82, q30_bases: 70 },
        ];
        let mut m = TrimMetrics::default();
        flatten_mate_stats(&before, &after, &mut m);
        assert_eq!(m.q20_before_r1, 90);
        assert_eq!(m.q30_before_r1, 80);
        assert_eq!(m.total_bases_before_r1, 100);
        assert_eq!(m.q20_after_r1, 88);
        assert_eq!(m.q30_after_r1, 80);
        assert_eq!(m.q20_before_r2, 85);
        assert_eq!(m.q30_before_r2, 70);
        assert_eq!(m.q20_after_r2, 82);
        assert_eq!(m.total_bases_after_r2, 90);
    }

    #[test]
    fn flatten_mate_stats_single_end_leaves_r2_zero() {
        let before = vec![MateStats { reads: 5, bases: 50, q20_bases: 50, q30_bases: 45 }];
        let after = vec![MateStats { reads: 4, bases: 40, q20_bases: 40, q30_bases: 35 }];
        let mut m = TrimMetrics::default();
        flatten_mate_stats(&before, &after, &mut m);
        assert_eq!(m.q20_before_r1, 50);
        assert_eq!(m.q20_before_r2, 0);
        assert_eq!(m.total_bases_before_r2, 0);
    }

    // ---- ratio / pct helpers ----

    #[test]
    fn ratio_zero_denom_returns_zero() {
        assert_eq!(ratio(5, 0), 0.0);
    }

    #[test]
    fn pct_zero_denom_returns_zero() {
        assert_eq!(pct(5, 0), 0.0);
    }

    // ---- LowQualFilter parser ----

    #[test]
    fn low_qual_filter_parses_q_f() {
        let f: LowQualFilter = "15:0.4".parse().unwrap();
        assert_eq!(f.threshold, 15);
        assert!((f.max_fraction - 0.4).abs() < 1e-9);
    }

    #[test]
    fn low_qual_filter_rejects_fraction_out_of_range() {
        assert!("15:1.5".parse::<LowQualFilter>().unwrap_err().contains("0.0..=1.0"));
        assert!("15:-0.1".parse::<LowQualFilter>().unwrap_err().contains("0.0..=1.0"));
    }

    #[test]
    fn low_qual_filter_rejects_missing_colon() {
        assert!("15".parse::<LowQualFilter>().unwrap_err().contains("Q:F"));
    }

    // ---- count_bases_below_q ----

    #[test]
    fn count_bases_below_q_simple() {
        // '!' = Phred 0, '5' = Phred 20, 'I' = Phred 40.
        assert_eq!(count_bases_below_q(b"!!!!", 1), 4);
        assert_eq!(count_bases_below_q(b"IIII", 40), 0);
        assert_eq!(count_bases_below_q(b"IIII", 41), 4);
        assert_eq!(count_bases_below_q(b"5555", 20), 0);
        assert_eq!(count_bases_below_q(b"5555", 21), 4);
    }

    #[test]
    fn count_bases_below_q_mixed_crosses_simd_boundary() {
        // 20 low + 20 high = 40 bytes total (crosses the 32-byte SIMD chunk boundary).
        let mut qual = vec![b'!'; 20];
        qual.extend(std::iter::repeat_n(b'I', 20));
        assert_eq!(count_bases_below_q(&qual, 10), 20);
    }

    // ---- 5' quality trim ----

    #[test]
    fn quality_trim_5p_trims_leading_bad_bases() {
        let mut rec =
            OwnedRecord { head: vec![], seq: b"ACGTACGTAC".to_vec(), qual: b"!!!IIIIIII".to_vec() };
        // Window 4, threshold 20: window sum >= 80 passes.
        //   s=0: Q0+Q0+Q0+Q40 = 40, fail
        //   s=1: Q0+Q0+Q40+Q40 = 80, pass
        // cut_front keeps bases from s=1 onward → trims 1 base, 9 remain.
        let trimmed = trim_quality_sliding_5prime(&mut rec, 4, 20);
        assert_eq!(trimmed, 1);
        assert_eq!(rec.seq.as_slice(), b"CGTACGTAC");
        assert_eq!(rec.qual.as_slice(), b"!!IIIIIII");
    }

    #[test]
    fn quality_trim_5p_empties_all_bad_read() {
        let mut rec = OwnedRecord { head: vec![], seq: b"ACGT".to_vec(), qual: b"!!!!".to_vec() };
        let trimmed = trim_quality_sliding_5prime(&mut rec, 2, 20);
        assert_eq!(trimmed, 4);
        assert!(rec.seq.is_empty());
        assert!(rec.qual.is_empty());
    }

    #[test]
    fn quality_trim_5p_noop_when_all_good() {
        let mut rec =
            OwnedRecord { head: vec![], seq: b"ACGTACGT".to_vec(), qual: b"IIIIIIII".to_vec() };
        let trimmed = trim_quality_sliding_5prime(&mut rec, 4, 20);
        assert_eq!(trimmed, 0);
        assert_eq!(rec.seq.as_slice(), b"ACGTACGT");
    }

    // ---- --filter-low-qual end-to-end ----

    #[test]
    fn execute_filter_low_qual_drops_reads_with_too_many_lowq_bases() {
        let tmp = TempDir::new().unwrap();
        let lines = vec![
            "@good".to_string(),
            "AAAAAAAAAA".to_string(),
            "+".to_string(),
            "IIIIIIIIII".to_string(), // 0% below Q15
            "@bad".to_string(),
            "CCCCCCCCCC".to_string(),
            "+".to_string(),
            "!!!!!!!!!!".to_string(), // 100% below Q15
        ];
        let r1 = write_fastq(&tmp, "r1", &lines);
        let out = tmp.path().join("o.fq.gz");
        let mut cmd = trim_cmd(vec![r1], vec![out.clone()], None);
        cmd.filter_low_qual = Some(LowQualFilter { threshold: 15, max_fraction: 0.4 });
        cmd.execute().unwrap();
        let written = read_fastq(&out);
        assert_eq!(written.len(), 1);
        assert_eq!(written[0].seq.as_slice(), b"AAAAAAAAAA");
    }

    #[test]
    fn polyg_zero_disables_trim() {
        // `--trim-polyg 0` is the user-facing way to disable poly-G. The read has 12 G's
        // on the tail that would be removed if trimming were on; they should survive.
        let tmp = TempDir::new().unwrap();
        let lines =
            vec!["@r".to_string(), "ACGTGGGGGGGGGGGG".to_string(), "+".to_string(), "I".repeat(16)];
        let r1 = write_fastq(&tmp, "r1", &lines);
        let out = tmp.path().join("o.fq.gz");
        let mut cmd = trim_cmd(vec![r1], vec![out.clone()], None);
        cmd.trim_polyg = 0;
        cmd.execute().unwrap();
        let written = read_fastq(&out);
        assert_eq!(written[0].seq.as_slice(), b"ACGTGGGGGGGGGGGG");
    }

    // ---- panic_message ----

    #[test]
    fn panic_message_downcasts_static_str() {
        let payload = std::panic::catch_unwind(|| panic!("static str panic")).unwrap_err();
        assert_eq!(panic_message(&payload), "static str panic");
    }

    #[test]
    fn panic_message_downcasts_string() {
        // `format!(...)` in panic!() produces a `String` payload rather than `&'static str`.
        let payload = std::panic::catch_unwind(|| panic!("string panic: {}", 42)).unwrap_err();
        assert_eq!(panic_message(&payload), "string panic: 42");
    }

    #[test]
    fn panic_message_falls_back_on_unknown_payload() {
        // Custom non-string panic payload.
        #[derive(Debug)]
        struct Weird;
        let payload = std::panic::catch_unwind(|| std::panic::panic_any(Weird)).unwrap_err();
        assert_eq!(panic_message(&payload), "(non-string panic payload)");
    }

    // ---- select_most_specific_error ----

    #[test]
    fn select_most_specific_error_empty_returns_none() {
        assert!(select_most_specific_error(vec![]).is_none());
    }

    #[test]
    fn select_most_specific_error_prefers_non_symptom() {
        // First error is a channel-closed *symptom*; second is the real cause.
        let errors = vec![
            anyhow!("writer exited before receiving batch slot"),
            anyhow!("disk full: ENOSPC"),
        ];
        let e = select_most_specific_error(errors).unwrap();
        assert!(e.to_string().contains("disk full"));
    }

    #[test]
    fn select_most_specific_error_falls_back_to_first_when_all_symptoms() {
        // If every error is a symptom, we return the first-seen one rather than nothing.
        let errors = vec![
            anyhow!("writer exited before receiving batch slot"),
            anyhow!("worker dropped before handoff"),
        ];
        let e = select_most_specific_error(errors).unwrap();
        assert!(e.to_string().contains("writer exited"));
    }

    // ---- OverlapStats::pick_mode hysteresis ----

    #[test]
    fn overlap_stats_pick_mode_does_not_flap_at_thresholds() {
        // Read length 100 → base thresholds: Ascending for mean<50, Descending for
        // mean>=90, Outward otherwise. Margin = 100/20 = 5.
        // Feeding means that straddle the Ascending/Outward boundary (48 vs 52) should
        // NOT flip out of Ascending because 52 < 50+5.
        let read_len = 100;
        let mut mode = OverlapWalkMode::Ascending;
        for &mean in &[48usize, 52, 49, 53, 47, 54] {
            mode = OverlapStats::pick_mode(mean, read_len, mode);
            assert!(
                matches!(mode, OverlapWalkMode::Ascending),
                "flapped out of Ascending at mean={mean}"
            );
        }
        // A mean clearly past the retention threshold should move us into Outward.
        mode = OverlapStats::pick_mode(60, read_len, mode);
        assert!(matches!(mode, OverlapWalkMode::Outward(60)));

        // Symmetric: starting in Descending, straddling the Descending/Outward boundary
        // (88 vs 92) should NOT flip out because 88 >= 90-5.
        let mut mode = OverlapWalkMode::Descending;
        for &mean in &[88usize, 92, 86, 93, 85, 94] {
            mode = OverlapStats::pick_mode(mean, read_len, mode);
            assert!(
                matches!(mode, OverlapWalkMode::Descending),
                "flapped out of Descending at mean={mean}"
            );
        }
        // A mean clearly past retention moves us into Outward.
        mode = OverlapStats::pick_mode(80, read_len, mode);
        assert!(matches!(mode, OverlapWalkMode::Outward(80)));
    }

    // ---- asymmetric EOF ----

    #[test]
    fn execute_errors_when_r1_and_r2_have_different_record_counts() {
        let tmp = TempDir::new().unwrap();
        let r1 = write_fastq(&tmp, "r1", &fq_lines("p", &["ACGT", "ACGT", "ACGT"]));
        let r2 = write_fastq(&tmp, "r2", &fq_lines("p", &["ACGT", "ACGT"]));
        let o1 = tmp.path().join("o1.fq.gz");
        let o2 = tmp.path().join("o2.fq.gz");
        let cmd = trim_cmd(vec![r1, r2], vec![o1, o2], None);
        let err = cmd.execute().unwrap_err().to_string();
        assert!(err.contains("out of sync"), "expected out-of-sync error, got: {err}");
    }

    // ---- filter drops every read ----

    #[test]
    fn execute_filter_drops_all_reads_accounting_holds() {
        // Every read is Q0 so the mean-quality filter drops them all. bases_in must
        // equal bases_out (0) + sum(trimmed=0) + bases_filtered; reads_out must be 0
        // and reads_filtered_quality must equal reads_in.
        let tmp = TempDir::new().unwrap();
        let lines = vec![
            "@a".to_string(),
            "AAAA".to_string(),
            "+".to_string(),
            "!!!!".to_string(),
            "@b".to_string(),
            "CCCC".to_string(),
            "+".to_string(),
            "!!!!".to_string(),
        ];
        let r1 = write_fastq(&tmp, "r1", &lines);
        let out = tmp.path().join("o.fq.gz");
        let metrics = tmp.path().join("m.txt");
        let mut cmd = trim_cmd(vec![r1], vec![out], Some(metrics.clone()));
        cmd.filter_mean_qual = Some(20);
        cmd.execute().unwrap();

        let contents = std::fs::read_to_string(&metrics).unwrap();
        let lines: Vec<&str> = contents.lines().collect();
        let header: Vec<&str> = lines[0].split('\t').collect();
        let values: Vec<&str> = lines[1].split('\t').collect();
        let get = |name: &str| -> u64 {
            let idx = header.iter().position(|h| *h == name).unwrap();
            values[idx].parse().unwrap()
        };

        assert_eq!(get("reads_in"), 2);
        assert_eq!(get("reads_out"), 0);
        assert_eq!(get("reads_filtered_quality"), 2);
        assert_eq!(get("bases_out"), 0);
        // Accounting: bases_in = bases_out + sum(trimmed_*) + bases_filtered.
        let sum_trimmed = get("bases_trimmed_read_structure")
            + get("bases_trimmed_adapter")
            + get("bases_trimmed_polyg")
            + get("bases_trimmed_polyx")
            + get("bases_trimmed_quality");
        assert_eq!(get("bases_in"), get("bases_out") + sum_trimmed + get("bases_filtered"));
    }

    // ---- IUPAC adapter end-to-end ----

    #[test]
    fn execute_iupac_adapter_trims_via_scalar_fallback() {
        // Adapter contains an `N` so the SIMD ACGT fast-path is bypassed and the scalar
        // IUPAC matcher is exercised. Read = 10bp insert + 12bp adapter-with-N.
        // Insert = `AAAAAAAAAA` (10bp). Adapter pattern = `AGATCGGNAGAG` (12bp, N at pos 7).
        // The read contains `AGATCGGAAGAG` at its 3' end (real base `A` at the N slot);
        // `N` in the adapter matches any read base, so the scalar scanner should match
        // starting at position 10 and trim 12 bp.
        let tmp = TempDir::new().unwrap();
        let reads: Vec<&str> = vec!["AAAAAAAAAAAGATCGGAAGAG"];
        let r1 = write_fastq(&tmp, "r1", &fq_lines("r", &reads));
        let out = tmp.path().join("out.fq.gz");
        let mut cmd = trim_cmd(vec![r1], vec![out.clone()], None);
        cmd.adapter_sequence = vec!["AGATCGGNAGAG".to_string()];
        cmd.execute().unwrap();

        let written = read_fastq(&out);
        assert_eq!(written.len(), 1);
        assert_eq!(written[0].seq.as_slice(), b"AAAAAAAAAA");
    }
}
