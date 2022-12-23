# fqtk

<p align="center">
  <a href="https://github.com/fulcrumgenomics/fqtk/actions?query=workflow%3ACheck"><img src="https://github.com/fulcrumgenomics/fqtk/actions/workflows/build_and_test.yml/badge.svg" alt="Build Status"></a>
</p>

A toolkit for working with FASTQ files, written in Rust.

## Installing

Until the first public release, the only way to install the tool is by checking out the git repository and building with cargo.

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