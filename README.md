# fqtk

<p align="center">
  <a href="https://github.com/fulcrumgenomics/fqtk/actions?query=workflow%3ACheck"><img src="https://github.com/fulcrumgenomics/fqtk/actions/workflows/build_and_test.yml/badge.svg" alt="Build Status"></a>
  <img src="https://img.shields.io/crates/l/fqtk.svg" alt="license">
  <a href="https://crates.io/crates/fqtk"><img src="https://img.shields.io/crates/v/fqtk.svg?colorB=319e8c" alt="Version info"></a>
  <a href="http://bioconda.github.io/recipes/fqtk/README.html"><img src="https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat" alt="Install with bioconda"></a>
  <a href="https://doi.org/10.5281/zenodo.13345413"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.13345413.svg" alt="DOI"></a>
  <br>
</p>

A toolkit for working with FASTQ files, written in Rust.

<p>
<a href="https://fulcrumgenomics.com">
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/fulcrumgenomics/fqtk/main/.github/logos/fulcrumgenomics-dark.svg">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/fulcrumgenomics/fqtk/main/.github/logos/fulcrumgenomics-light.svg">
  <img alt="Fulcrum Genomics" src="https://raw.githubusercontent.com/fulcrumgenomics/fqtk/main/.github/logos/fulcrumgenomics-light.svg" height="100">
</picture>
</a>
</p>

[Visit us at Fulcrum Genomics](https://www.fulcrumgenomics.com) to learn more about how we can power your Bioinformatics with fqtk and beyond.

<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img alt="Email Fulcrum Genomics" src="https://img.shields.io/badge/Email_us-%2338b44a.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
<a href="https://www.fulcrumgenomics.com"><img alt="Visit Fulcrum Genomics" src="https://img.shields.io/badge/Visit_Us-%2326a8e0.svg?&style=for-the-badge&logo=wordpress&logoColor=white"/></a>

`fqtk` provides several tools for working with FASTQ files:

- **`demux`** — demultiplex one or more FASTQ files into per-sample FASTQs using sample barcodes at fixed positions within the reads.
- **`shard`** — split one or more matched FASTQs (e.g. R1/R2) into N shards, assigning reads round-robin so each input read ends up in exactly one output FASTQ.
- **`subsample`** — randomly subsample reads from one or more synchronized FASTQs.

All tools are highly efficient and multi-threaded for high performance.

## `fqtk demux`

`fqtk demux` demultiplexes one or more FASTQ files (e.g. a set of R1, R2 and I1 FASTQ files) with any number of sample barcodes at fixed locations within the reads.

Usage for `fqtk demux` follows:

<!-- start usage:demux -->
````console

Performs sample demultiplexing on FASTQs.

The sample barcode for each sample in the metadata TSV will be compared against the sample
barcode bases extracted from the FASTQs, to assign each read to a sample.  Reads that do not
match any sample within the given error tolerance will be placed in the ``unmatched_prefix``
file.

FASTQs and associated read structures for each sub-read should be given:

- a single fragment read (with inline index) should have one FASTQ and one read structure
- paired end reads should have two FASTQs and two read structures
- a dual-index sample with paired end reads should have four FASTQs and four read structures
  given: two for the two index reads, and two for the template reads.

If multiple FASTQs are present for each sub-read, then the FASTQs for each sub-read should be
concatenated together prior to running this tool
(e.g. `zcat s_R1_L001.fq.gz s_R1_L002.fq.gz | bgzip -c > s_R1.fq.gz`).

(Read structures)[<https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures>] are made up of
`<number><operator>` pairs much like the `CIGAR` string in BAM files.
Five kinds of operators are recognized:

1. `T` identifies a template read
2. `B` identifies a sample barcode read
3. `M` identifies a unique molecular index read
4. `C` identifies a unique cellular barcode read
5. `S` identifies a set of bases that should be skipped or ignored

The last `<number><operator>` pair may be specified using a `+` sign instead of number to
denote "all remaining bases". This is useful if, e.g., fastqs have been trimmed and contain
reads of varying length. Both reads must have template bases.

Metadata about the samples should be given as a headered metadata TSV file with at least the
following two columns present:

1. `sample_id` - the id of the sample or library.
2. `barcode` - the expected barcode sequence associated with the `sample_id`.

For reads containing multiple barcodes (such as dual-indexed reads), all barcodes should be
concatenated together in the order they are read and stored in the `barcode` field.

IUPAC bases are supported in the (expected) `barcode` column.  An observed IUPAC base must be
at least as specific as the corresponding base in the expected sample barcode.  E.g. If the
observed base is an N, it will only match expected sample barcrods with an N.  And if the
observed base is an R, it will match R, V, D, and N, since the latter IUPAC codes allow both
A and G (R/V/D/N are a superset of the bases compare to R).

The read structures will be used to extract the observed sample barcode, template bases,
molecular identifiers, and cellular barcodes from each read.  The observed sample barcode will
be matched to the sample barcodes extracted from the bases in the sample metadata and associated
read structures.

An observed barcode matches an expected barcode if all the following are true:
1. The number of mismatches (edits/substitutions) is less than or equal to the maximum
   mismatches (see `--max-mismatches`).
2. The difference between number of mismatches in the best and second best barcodes is greater
   than or equal to the minimum mismatch delta (`--min-mismatch-delta`).

The expected barcode sequence may contains Ns, which are not counted as mismatches regardless
of the observed base (e.g. the expected barcode `AAN` will have zero mismatches relative to
both the observed barcodes `AAA` and `AAN`).

## Per-sample read structures

In addition to the global `--read-structures`, the metadata TSV may include optional columns
`read_structure_1`, `read_structure_2`, ..., `read_structure_<N>` (where `N` is the number of
input FASTQs).  When any cell is non-empty for a sample, the per-sample read structures
replace the global `--read-structures` for that sample — both for matching and for output
extraction.

Per-sample structures support per-cell fall-back to the global `--read-structures`:

- A blank `read_structure_<i>` cell uses `--read-structures[i-1]` for that sample's input *i*.
- A row whose `read_structure_<n>` cells are all blank uses `--read-structures` entirely for
  that sample (equivalent to omitting the columns for that sample only).
- The unmatched pseudo-sample always uses the global `--read-structures`.

Constraints:

1. The number of `read_structure_<n>` columns must equal `--read-structures.len()`.
2. The concatenated `B`-segment length must equal the `barcode` column length for every
   sample (computed from each sample's *effective* read structures, i.e. with fall-backs
   applied).

Different samples may have different per-input `(T, B, M, C)` segment counts (and hence
produce different sets of output files).  This supports protocols with sample-dependent
read structures (e.g. CODEC, where each sample may include a stagger spacer of varying
length so that the position of the constant ligation base shifts per sample).

During matching, each sample's expected pattern is constructed from its effective read
structure by filling `B`-segment positions from the `barcode` column and treating
`M`/`S`/`C` segment positions as `N` wildcards.  All patterns are padded with trailing `N`s
to a per-input matching window equal to the longest pre-template prefix across samples.

Example metadata TSV (CODEC stagger):

```text
sample_id  barcode         read_structure_1   read_structure_2
S1         GATTACAGATTACA  3M7B1S+T           3M7B1S+T
S2         TTTTTTTTTTTTTT  3M1S7B1S+T         3M1S7B1S+T
```

## Outputs

All outputs are generated in the provided `--output` directory.  For each sample plus the
unmatched reads, FASTQ files are written for each read segment (specified in the read
structures) of one of the types supplied to `--output-types`.  FASTQ files have names
of the format:

```bash
{sample_id}.{segment_type}{read_num}.fq.gz
```

where `segment_type` is one of `R`, `I`, `U`, and `C` (for template, sample barcode/index,
molecular barcode/UMI, and cellular barcode reads, respectively) and `read_num` is a number
starting at 1 for each segment type.

In addition a `demux-metrics.txt` file is written that is a tab-delimited file with counts
of how many reads were assigned to each sample and derived metrics.

## Example Command Line

As an example, if the sequencing run was 2x100bp (paired end) with two 8bp index reads both
reading a sample barcode, as well as an in-line 8bp sample barcode in read one, the command
line would be:

```bash
fqtk demux \
    --inputs r1.fq.gz i1.fq.gz i2.fq.gz r2.fq.gz \
    --read-structures 8B92T 8B 8B 100T \
    --sample-metadata metadata.tsv \
    --output output_folder
```

Usage: fqtk demux [OPTIONS] --inputs <INPUTS>... --read-structures <READ_STRUCTURES>... --sample-metadata <SAMPLE_METADATA> --output <OUTPUT>

Options:
  -i, --inputs <INPUTS>...
          One or more input FASTQ files each corresponding to a sequencing read (e.g. R1, I1)

  -r, --read-structures <READ_STRUCTURES>...
          The read structures, one per input FASTQ in the same order.

          Per-sample read structures (see the `read_structure_<n>` metadata columns) take precedence for each matched sample, and a blank cell falls back to the corresponding `--read-structures` entry.  The unmatched pseudo-sample always uses `--read-structures` for its output extraction.  The number of `read_structure_<n>` columns must equal `--read-structures.len()`.

  -b, --output-types <OUTPUT_TYPES>...
          The read structure types to write to their own files (Must be one of T, B, M, or C for template reads, sample barcode reads, molecular barcode reads, or cellular barcode reads).

          Multiple output types may be specified as a space-delimited list.

          [default: T]

  -s, --sample-metadata <SAMPLE_METADATA>
          A file containing the metadata about the samples

  -o, --output <OUTPUT>
          The output directory into which to write per-sample FASTQs

  -u, --unmatched-prefix <UNMATCHED_PREFIX>
          Output prefix for FASTQ file(s) for reads that cannot be matched to a sample

          [default: unmatched]

      --max-mismatches <MAX_MISMATCHES>
          Maximum mismatches for a barcode to be considered a match

          [default: 1]

  -d, --min-mismatch-delta <MIN_MISMATCH_DELTA>
          Minimum difference between number of mismatches in the best and second best barcodes for a barcode to be considered a match

          [default: 2]

  -t, --threads <THREADS>
          The number of threads to use. Cannot be less than 5

          [default: 8]

  -c, --compression-level <COMPRESSION_LEVEL>
          The level of compression to use to compress outputs

          [default: 5]

  -S, --skip-reasons <SKIP_REASONS>
          Skip demultiplexing reads for any of the following reasons, otherwise panic.

          1. `too-few-bases`: there are too few bases or qualities to extract given the read structures.  For example, if a read is 8bp long but the read structure is `10B`, or if a read is empty and the read structure is `+T`.

      --template-types <TEMPLATE_TYPES>...
          The read structure types to include in the template FASTQ output files.

          By default, only template (T) segments are included. To include additional segment types (e.g. to preserve UMIs in the output read bases), specify them here. For example, `--template-types M T` will concatenate the molecular barcode and template segments.

          To output the full original reads (all segments), specify all segment types present in your read structure (e.g. `--template-types B M T`).

          Segments are only merged *within the same physical read*: a non-`T` segment is folded into the template bases only when it is co-located with a `T` in the same read structure (e.g. `8M84T`). A segment on a separate read (e.g. a UMI on its own index read) is never merged into a template on another read; route it via `--output-types` instead, or leave it in the read header. When a UMI (M) is included here it is written into the template bases and is therefore omitted from the read header (it is not written in both places).

          Note: If `--template-types` includes any non-`T` type, `T` must be included in `--output-types`; each requested non-`T` type must be co-located with a `T` in every read structure where it appears; and each read structure must contain at most one `T` segment.

          [default: T]

      --header-format <HEADER_FORMAT>
          The FASTQ header "shape" for each output record: whether the comment is rewritten (`illumina`), kept verbatim (`unmodified`), or dropped entirely (`name-only`).  See each value's own description below for exact per-format behavior and examples.

          With `unmodified` and `name-only`, sample barcode bases (and, unless `--umi-in-name true` is given, UMI bases) are retained only if routed to their own FASTQs via `--output-types` or folded into the template bases via `--template-types`; otherwise they are not present in any output.  A warning is emitted when this would silently discard bases.

          Possible values:
          - illumina:   Casava ≥1.8 / bcl-convert-style header: the read-number field in the comment is rewritten to match the output file and the sample barcode(s) are appended to the comment.  UMI(s) are folded into the read name by default (override with `--umi-in-name`).  For example `@inst:1:FC:1:1101:5:7 1:N:0:0` becomes `@inst:1:FC:1:1101:5:7:AACCGGTT 1:N:0:ACGTACGT`
          - unmodified: Emit each read's original header verbatim: the read-number field is not rewritten and no sample barcode is appended.  UMI(s) are NOT folded into the name by default (override with `--umi-in-name true`).  Nothing is parsed unless a UMI is folded in, so with the default `--umi-in-name` this (and `name-only`, below) is one of the formats that passes through headers which do not follow Illumina conventions unchanged (more than eight `:`-delimited name fields, or a comment that is not four `:`-delimited fields), such as those produced by MGI, Element, and ONT instruments
          - name-only:  Emit only the read name: any pre-existing comment is always dropped, regardless of `--umi-in-name`.  UMI(s) are NOT folded into the name by default (override with `--umi-in-name true`), in which case the same non-Illumina-header tolerance described for `unmodified` applies here too (the name is only parsed when a UMI is actually folded into it).  Intended for pipelines (e.g. `bwa mem -C`) where a Casava-style comment would break downstream parsing.  With no comment and (by default) no UMI in the name, R1/R2 headers become byte-identical (no read number anywhere) — the same tradeoff `samtools fastq -n` makes for file-paired workflows

          [default: illumina]

      --umi-in-name <UMI_IN_NAME>
          Whether to fold UMI(s) into the read name.

          If not given, this is decided by `--header-format`: `true` for `illumina`, `false` for `unmodified` and `name-only`.  An explicit value here overrides the format's default -- e.g. `--header-format unmodified --umi-in-name true` keeps the original comment but still appends the UMI(s) to the name.  Has no effect when the UMI is instead folded into the template bases via `--template-types M ...`, since there is then nothing left to fold into the name.

          [possible values: true, false]

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
````
<!-- end usage:demux -->

## `fqtk shard`

`fqtk shard` splits one or more matched FASTQs (e.g. R1 and R2) into N output shards, assigning reads to shards on a round-robin basis so that each input read ends up in exactly one output FASTQ.  This is useful for splitting large FASTQs into evenly sized pieces for parallel downstream processing.

Usage for `fqtk shard` follows:

<!-- start usage:shard -->
````console

Shards a set of FASTQs into N output shards.

Shards a set of matched FASTQs (e.g. R1 and R2) into one or more set of FASTQs where each
input read ends up in exactly one output FASTQ. Reads are assigned to shards on a round-robin
basis, so e.g. if using `--shards 10` the first read in the input files will end up in the
first shard, the second read in the second shard ... and the tenth read in the tenth shard.

Each shard will contain one output FASTQ file per input FASTQ files.  Output files are named
as follows:

```
{output_prefix}.{shard_prefix}{shard_num}.{read_number_prefix}{read_num}.fq.gz
```

where `shard_num` is n for the nth shard (starting at 1), `read_num` corresponds to the nth
file in the `inputs` list (starting at 1), and all other values in `{}` are named command
line parameters.  The `output_prefix` may contain an absolute path, or a relative path, with
relative paths interpreted relative to the working directory where the command is run.

Inputs may be uncompressed, gzipped, or block-gzipped.  Output files are _always_ block gzipped.

Usage: fqtk shard [OPTIONS] --inputs <INPUTS>... --output-prefix <OUTPUT_PREFIX> --shards <SHARDS>

Options:
  -i, --inputs <INPUTS>...
          One or more input FASTQ files each corresponding to a sequencing read (e.g. R1, R2)

  -o, --output-prefix <OUTPUT_PREFIX>
          Output prefix for sharded FASTQ file(s)

  -S, --shard-prefix <SHARD_PREFIX>
          Prefix to place before the shard number in the generated output file names

          [default: s]

  -R, --read-number-prefix <READ_NUMBER_PREFIX>
          Prefix to place before the read number in the generated output file names

          [default: r]

  -s, --shards <SHARDS>
          Number of shards to generate

  -t, --threads <THREADS>
          The number of threads to use for compressing output files.  Minimum 2

          [default: 8]

  -c, --compression-level <COMPRESSION_LEVEL>
          The level of compression to use to compress outputs.  Defaults to 1 because sharded FASTQs are typically short-lived intermediates, where write throughput matters more than squeezing out the last few percent of file size

          [default: 1]

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
````
<!-- end usage:shard -->

## `fqtk subsample`

`fqtk subsample` reads one or more synchronized FASTQs (e.g. R1 and R2) and writes a randomly chosen subset of the reads, keeping or discarding each read across all files together so paired reads stay in sync.

Usage for `fqtk subsample` follows:

<!-- start usage:subsample -->
````console

Subsamples reads from one or more synchronized FASTQ files.

Reads one or more FASTQ files (e.g. paired-end R1 and R2) and writes a
random subset of reads to output files. All input files must contain the
same number of reads in the same order; each read is either kept or
discarded across all files simultaneously.

Output files are named `{output}.R1.fq.gz`, `{output}.R2.fq.gz`, etc.
and are always BGZF compressed.

Each read is independently retained with probability equal to `--fraction`,
giving an approximate subsample without needing to know the total read count
upfront. When no explicit `--seed` is provided, a deterministic seed is
derived from all input parameters, so identical inputs and parameters always
produce identical output.

# Example

```bash
fqtk subsample \
    --input r1.fq.gz r2.fq.gz \
    --output subsampled \
    --fraction 0.1
```

Usage: fqtk subsample [OPTIONS] --inputs <INPUTS>... --output <OUTPUT> --fraction <FRACTION>

Options:
  -i, --inputs <INPUTS>...
          One or more input FASTQ files (may be gzipped). All files must have the same number of reads in the same order

  -o, --output <OUTPUT>
          Output path prefix. Files will be named {output}.R1.fq.gz, etc

  -f, --fraction <FRACTION>
          Fraction of reads to retain, in the range [0.0, 1.0]

  -t, --threads <THREADS>
          Number of threads for compression. Minimum 2

          [default: 8]

  -c, --compression-level <COMPRESSION_LEVEL>
          BGZF compression level for output files

          [default: 5]

  -s, --seed <SEED>
          Explicit RNG seed for reproducibility. When omitted, a deterministic seed is derived from all other parameters

      --disable-read-name-checking
          Disable checking that read names are in sync across input files

  -h, --help
          Print help (see a summary with '-h')

  -V, --version
          Print version
````
<!-- end usage:subsample -->

## Installing

### Installing with `conda`
To install with conda you must first [install conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html#installation).
Then, in your command line (and with the environment you wish to install fqtk into active) run:

```console
conda install -c bioconda fqtk
```

### Installing with `cargo`
To install with cargo you must first [install rust](https://doc.rust-lang.org/cargo/getting-started/installation.html).
Which (On Mac OS and Linux) can be done with the command:

```console
curl https://sh.rustup.rs -sSf | sh
```

Then, to install `fqtk` run:

```console
cargo install fqtk
```

### Building From Source

First, clone the git repo:

```console
git clone https://github.com/fulcrumgenomics/fqtk.git
```

Secondly, if you do not already have rust development tools installed, install via [rustup](https://rustup.rs/):

```console
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

Then build the toolkit in release mode:

```console
cd fqtk
cargo build --release
./target/release/fqtk --help
```

## Developing

fqtk is developed in Rust and follows the conventions of using `rustfmt` and `clippy` to ensure both code quality and standardized formatting.
When working on fqtk, before pushing any commits, please first run `./ci/check.sh` and resolve any issues that are reported.
Note that `./ci/check.sh` only checks formatting; to auto-fix formatting issues, run `cargo fmt --all`.

## Releasing a New Version

Releases are automated with [`release-plz`][release-plz-link], driven by the
[`.github/workflows/release.yml`](.github/workflows/release.yml) workflow that runs on every push to `main`.

### Semantic Versioning

This tool follows [Semantic Versioning](https://semver.org/):

* MAJOR version when you make incompatible API changes,
* MINOR version when you add functionality in a backwards-compatible manner, and
* PATCH version when you make backwards-compatible bug fixes.

`release-plz` derives the next version automatically from the [Conventional Commits][conventional-commits-link] merged since the last release.

### Release Flow

1. As commits land on `main`, `release-plz` opens (or updates) a **release PR** labelled `release`. This PR bumps the version in `Cargo.toml`, refreshes `Cargo.lock`, and updates `CHANGELOG.md` (rendered via [`cliff.toml`](cliff.toml)).
2. Review the release PR. Adjust the version bump if needed by amending commit messages or the changelog, then merge it.
3. Merging the release PR re-triggers the workflow. The `publish` job detects the new version, runs `cargo publish` to push the crate to [crates.io](https://crates.io/crates/fqtk), and creates the matching `vX.Y.Z` git tag and GitHub release.

The publish job is idempotent: if the version in `Cargo.toml` has already been published to crates.io, it skips publishing, and tag/release creation is skipped if they already exist.

### Pre-requisites (maintainers)

The workflow requires two repository secrets:

* `RELEASE_PLZ_TOKEN` — a token (PAT or GitHub App token) with `contents: write` and `pull-requests: write` so `release-plz` can open release PRs and push tags.
* `CARGO_REGISTRY_TOKEN` — a [crates.io API token](https://crates.io/settings/tokens) used by `cargo publish`.

[release-plz-link]:         https://release-plz.dev
[conventional-commits-link]: https://www.conventionalcommits.org
