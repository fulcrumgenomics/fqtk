use crate::commands::command::Command;
use anyhow::Result;
use clap::Parser;
use ext_sort::{ExternalSorter, ExternalSorterBuilder, LimitedBufferBuilder};
use seq_io::fastq::Reader as FastqReader;
use seq_io::fastq::{OwnedRecord, Record};
use serde::{Deserialize, Serialize};
use std::cmp::Ordering;
use std::fs::File;
use std::io::BufWriter;
use std::path::{Path, PathBuf};

#[derive(Deserialize, Serialize)]
struct SortableFastqRecord(OwnedRecord);

impl PartialEq for SortableFastqRecord {
    fn eq(&self, other: &Self) -> bool {
        self.0.id().unwrap() == other.0.id().unwrap()
    }
}

impl Eq for SortableFastqRecord {}

impl PartialOrd for SortableFastqRecord {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(&other))
    }
}

impl Ord for SortableFastqRecord {
    fn cmp(&self, other: &Self) -> Ordering {
        self.0.id().unwrap().cmp(&other.0.id().unwrap())
    }
}

////////////////////////////////////////////////////////////////////////////////
// SortFastq (main class) and it's impls
////////////////////////////////////////////////////////////////////////////////

/// Sorts the input FASTQ file by lexicographic ordering of its read names.
///
/// ## Example Command Line
///
/// ```
/// fqtk sort-fastq \
///     --input test.fq  \
///     --output test.sorted.fq \
///     --max-records 1000000
/// ```
///
#[derive(Parser, Debug)]
#[command(version)]
pub(crate) struct SortFastq {
    /// Input fastq file
    #[clap(long, short = 'i', required = true)]
    input: PathBuf,

    /// Output fastq path
    #[clap(long, short = 'o', required = true)]
    output: PathBuf,

    /// Maximum number of records to be kept in buffer
    #[clap(long, short = 'm', default_value = "1000000")]
    max_records: usize,
}

impl Command for SortFastq {
    /// Executes the sort_fastq command
    fn execute(&self) -> Result<()> {
        let mut fq_reader = FastqReader::from_path(&self.input)?;
        let mut fq_writer = BufWriter::new(File::create(&self.output)?);

        let sorter: ExternalSorter<
            SortableFastqRecord,
            seq_io::fastq::Error,
            LimitedBufferBuilder,
        > = ExternalSorterBuilder::new()
            .with_tmp_dir(Path::new("./"))
            .with_buffer(LimitedBufferBuilder::new(self.max_records, true))
            .build()
            .unwrap();

        let sorted = sorter
            .sort(
                fq_reader.records().map(|record| record.map(|record| SortableFastqRecord(record))),
            )
            .unwrap();

        for item in sorted.map(Result::unwrap) {
            let _ =
                seq_io::fastq::write_to(&mut fq_writer, &item.0.head, &item.0.seq, &item.0.qual);
        }

        Ok(())
    }
}
