// MIT License
//
// Copyright (c) 2025 Fulcrum Genomics LLC
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

use crate::commands::command::Command;
use anyhow::{bail, ensure, Result};
use bstr::{BString, ByteSlice};
use clap::Parser;
use fgoxide::io::Io;
use fqtk_lib::read_set::FastqSegment;
use fqtk_lib::read_set::ReadSet;
use fqtk_lib::read_set::ReadSetIterator;
use noodles::sam;
use noodles::sam::alignment::io::Write;
use noodles::sam::alignment::record::data::field::Tag;
use noodles::sam::alignment::record::data::field::Value;
use noodles::sam::alignment::record::Flags;
use noodles::sam::alignment::record_buf::QualityScores;
use noodles::sam::alignment::RecordBuf;
use noodles::sam::header::record::value::map::read_group::{self, tag as rg_tag};
use noodles::sam::header::record::value::map::ReadGroup;
use noodles::sam::header::record::value::Map;
use noodles::sam::header::Header;
use read_structure::{ReadStructure, SegmentType};
use seq_io::fastq::Reader as FastqReader;
use seq_io::fastq::Record;
use std::fs::File;
use std::io::{BufRead, BufWriter};
use std::path::PathBuf;
use std::str::FromStr;

const BUFFER_SIZE: usize = 1024 * 1024;
const QUALITY_DETECTION_SAMPLE_SIZE: usize = 400;

/// Quality encoding type
#[derive(Debug, Clone, Copy)]
enum QualityEncoding {
    Standard, // Phred+33 (Sanger)
}

impl QualityEncoding {
    /// Convert quality scores to standard numeric format (Phred+33)
    fn to_standard_numeric(self, quals: &[u8]) -> Vec<u8> {
        match self {
            QualityEncoding::Standard => quals.iter().map(|&q| q.saturating_sub(33)).collect(),
        }
    }

    /// Detect encoding from sample records
    fn detect(records: &[Vec<u8>]) -> Result<Self> {
        // Simple detection: if all quality scores are >= 33, assume Sanger/Phred+33
        for qual in records {
            for &q in qual {
                if q < 33 {
                    bail!("Quality scores below 33 detected - unsupported encoding");
                }
            }
        }
        Ok(QualityEncoding::Standard)
    }
}

/// Generates an unmapped BAM (or SAM or CRAM) file from FASTQ files.
///
/// Takes in one or more FASTQ files (optionally gzipped), each representing a different
/// sequencing read (e.g. R1, R2, I1 or I2) and can use a set of read structures to
/// allocate bases in those reads to template reads, sample indices, unique molecular
/// indices, cell barcodes, or to designate bases to be skipped over.
#[derive(Parser, Debug)]
#[command(version)]
#[clap(verbatim_doc_comment)]
pub(crate) struct Extract {
    /// Input FASTQ files corresponding to each sequencing read (e.g. R1, I1, etc.)
    #[clap(long, short = 'i', required = true, num_args = 1..)]
    inputs: Vec<PathBuf>,

    /// Output SAM or BAM file to be written
    #[clap(long, short = 'o', required = true)]
    output: PathBuf,

    /// Read structures, one for each of the FASTQs (optional if 1-2 template-only FASTQs)
    #[clap(long, short = 'r', num_args = 0..)]
    read_structures: Vec<ReadStructure>,

    /// Tag in which to store molecular barcodes/UMIs
    #[clap(long, short = 'u', default_value = "RX")]
    umi_tag: String,

    /// Tag in which to store molecular barcode/UMI qualities
    #[clap(long, short = 'q')]
    umi_qual_tag: Option<String>,

    /// Tag in which to store the cellular barcodes
    #[clap(long, short = 'c', default_value = "CB")]
    cell_tag: String,

    /// Tag in which to store the cellular barcode qualities
    #[clap(long, short = 'C')]
    cell_qual_tag: Option<String>,

    /// Store the sample barcode qualities in the QT Tag
    #[clap(long, short = 'Q')]
    store_sample_barcode_qualities: bool,

    /// Extract UMI(s) from read names and prepend to UMIs from reads
    #[clap(long, short = 'n')]
    #[allow(clippy::struct_field_names)]
    extract_umis_from_read_names: bool,

    /// Read group ID to use in the file header
    #[clap(long, default_value = "A")]
    read_group_id: String,

    /// The name of the sequenced sample
    #[clap(long, required = true)]
    sample: String,

    /// The name/ID of the sequenced library
    #[clap(long, required = true)]
    library: String,

    /// Library or Sample barcode sequence
    #[clap(long, short = 'b')]
    barcode: Option<String>,

    /// Sequencing Platform
    #[clap(long, default_value = "illumina")]
    platform: String,

    /// Platform unit (e.g. '<flowcell-barcode>.<lane>.<sample-barcode>')
    #[clap(long)]
    platform_unit: Option<String>,

    /// Platform model to insert into the group header (ex. miseq, hiseq2500, hiseqX)
    #[clap(long)]
    platform_model: Option<String>,

    /// The sequencing center from which the data originated
    #[clap(long)]
    sequencing_center: Option<String>,

    /// Predicted median insert size, to insert into the read group header
    #[clap(long)]
    predicted_insert_size: Option<u32>,

    /// Description of the read group
    #[clap(long)]
    description: Option<String>,

    /// Comment(s) to include in the output file's header
    #[clap(long, num_args = 0..)]
    comment: Vec<String>,

    /// Date the run was produced, to insert into the read group header
    #[clap(long)]
    run_date: Option<String>,
}

impl Extract {
    /// Get actual read structures (default to +T if none provided for 1-2 FASTQs)
    fn get_read_structures(&self) -> Result<Vec<ReadStructure>> {
        if self.read_structures.is_empty() && (1..=2).contains(&self.inputs.len()) {
            Ok(vec![ReadStructure::from_str("+T")?; self.inputs.len()])
        } else {
            Ok(self.read_structures.clone())
        }
    }

    /// Validate inputs
    fn validate(&self) -> Result<()> {
        let read_structures = self.get_read_structures()?;

        ensure!(
            self.inputs.len() == read_structures.len(),
            "input and read-structure must be supplied the same number of times."
        );

        let template_count: usize = read_structures
            .iter()
            .map(|rs| rs.segments_by_type(SegmentType::Template).count())
            .sum();

        ensure!(
            (1..=2).contains(&template_count),
            "read structures must contain 1-2 template reads total."
        );

        ensure!(
            !self.extract_umis_from_read_names || self.umi_qual_tag.is_none(),
            "Cannot extract UMI qualities when also extracting UMI from read names."
        );

        Ok(())
    }

    /// Create SAM header
    fn create_header(&self) -> Result<Header> {
        let mut header = Header::builder();

        // Add comments
        for comment in &self.comment {
            header = header.add_comment(comment.clone());
        }

        // Create read group

        let mut rg = Map::<ReadGroup>::builder();

        rg = rg.insert(rg_tag::SAMPLE, self.sample.clone());
        rg = rg.insert(rg_tag::LIBRARY, self.library.clone());

        if let Some(ref bc) = self.barcode {
            rg = rg.insert(rg_tag::BARCODE, bc.clone());
        }

        rg = rg.insert(rg_tag::PLATFORM, self.platform.clone());

        if let Some(ref pu) = self.platform_unit {
            rg = rg.insert(rg_tag::PLATFORM_UNIT, pu.clone());
        }

        if let Some(ref pm) = self.platform_model {
            rg = rg.insert(rg_tag::PLATFORM_MODEL, pm.clone());
        }

        if let Some(ref sc) = self.sequencing_center {
            rg = rg.insert(rg_tag::SEQUENCING_CENTER, sc.clone());
        }

        if let Some(isize) = self.predicted_insert_size {
            rg = rg.insert(rg_tag::PREDICTED_MEDIAN_INSERT_SIZE, isize.to_string());
        }

        if let Some(ref desc) = self.description {
            rg = rg.insert(rg_tag::DESCRIPTION, desc.clone());
        }

        if let Some(ref date) = self.run_date {
            // Note: noodles expects ISO 8601 format
            rg = rg.insert(read_group::tag::PRODUCED_AT, date.clone());
        }

        header = header.add_read_group(self.read_group_id.clone(), rg.build()?);

        Ok(header.build())
    }

    fn add_tag(
        data: &mut noodles::sam::alignment::record_buf::Data,
        tag: &[u8],
        value: &BString,
    ) -> Result<(), Box<dyn std::error::Error>> {
        // Create the tag
        let bytes = tag.as_bytes();
        let tag = Tag::from([bytes[0], bytes[1]]);
        // Create the value
        let bstr_value = value.as_bstr();
        // Insert it
        data.insert(tag, Value::String(bstr_value).try_into()?);
        Ok(())
    }

    /// Extracts read name, optionally extracting UMI from the 8th colon-delimited field
    fn extract_read_name_and_umi(
        header: &Vec<u8>,
        extract_umis: bool,
    ) -> (Vec<u8>, Option<Vec<u8>>) {
        // Remove @ prefix if present
        let name_bytes = if header.starts_with(b"@") { &header[1..] } else { header.as_slice() };

        // Split on space to get just the name part
        let name_part = name_bytes.find_byte(b' ').map_or(name_bytes, |pos| &name_bytes[..pos]);

        if !extract_umis {
            return (name_part.to_vec(), None);
        }

        // Count colons to find 8th field (UMI)
        let parts: Vec<&[u8]> = name_part.split(|&b| b == b':').collect();

        if parts.len() >= 8 {
            let umi = parts[7].to_vec();
            if !umi.is_empty() {
                return (name_part.to_vec(), Some(umi));
            }
        }

        (name_part.to_vec(), None)
    }

    /// Create SAM records from a read set
    #[allow(clippy::too_many_lines)]
    fn make_sam_records(
        &self,
        read_set: &ReadSet,
        encoding: QualityEncoding,
    ) -> Result<Vec<RecordBuf>> {
        let templates: Vec<&FastqSegment> = read_set.template_segments().collect();

        ensure!(!templates.is_empty(), "No template segments found");

        // Extract various barcode types as BString to ensure proper lifetimes
        let cell_barcode_bs: BString = read_set
            .cell_barcode_segments()
            .map(|s| String::from_utf8_lossy(&s.seq))
            .collect::<Vec<_>>()
            .join("-")
            .into();

        let cell_quals_bs: BString = read_set
            .cell_barcode_segments()
            .map(|s| String::from_utf8_lossy(&s.quals))
            .collect::<Vec<_>>()
            .join(" ")
            .into();

        let sample_barcode_bs: BString = read_set
            .sample_barcode_segments()
            .map(|s| String::from_utf8_lossy(&s.seq))
            .collect::<Vec<_>>()
            .join("-")
            .into();

        let sample_quals_bs: BString = read_set
            .sample_barcode_segments()
            .map(|s| String::from_utf8_lossy(&s.quals))
            .collect::<Vec<_>>()
            .join(" ")
            .into();

        let umi_bs: BString = read_set
            .molecular_barcode_segments()
            .map(|s| String::from_utf8_lossy(&s.seq))
            .collect::<Vec<_>>()
            .join("-")
            .into();

        let umi_qual_bs: BString = read_set
            .molecular_barcode_segments()
            .map(|s| String::from_utf8_lossy(&s.quals))
            .collect::<Vec<_>>()
            .join(" ")
            .into();

        // Extract UMI from read name if requested
        let (read_name, umi_from_name) =
            Self::extract_read_name_and_umi(&read_set.header, self.extract_umis_from_read_names);

        // Prepare final UMI as BString
        let final_umi_bs: BString = match (umi_bs.is_empty(), &umi_from_name) {
            (true, Some(from_name)) => String::from_utf8_lossy(from_name).to_string().into(),
            (true, None) => BString::from(""),
            (false, Some(from_name)) => {
                format!("{}-{}", String::from_utf8_lossy(from_name), umi_bs.as_bstr()).into()
            }
            (false, None) => umi_bs.clone(),
        };

        // Read group ID as BString
        let rgid_bs: BString = self.read_group_id.as_str().into();

        let mut records = Vec::new();

        for (index, template) in templates.iter().enumerate() {
            // Create record with proper API
            let mut rec = RecordBuf::default();

            // Set read name
            *rec.name_mut() = Some(read_name.clone().into());

            // Set bases and qualities - if empty, substitute with single N @ Q2
            if template.seq.is_empty() {
                *rec.sequence_mut() = "N".as_bytes().to_vec().into();
                *rec.quality_scores_mut() = QualityScores::from(vec![2u8]);
            } else {
                *rec.sequence_mut() = template.seq.clone().into();
                let numeric_quals = encoding.to_standard_numeric(&template.quals);
                *rec.quality_scores_mut() = QualityScores::from(numeric_quals);
            }

            // Set flags for unmapped reads
            let mut flags = Flags::UNMAPPED;
            if templates.len() == 2 {
                flags |= Flags::SEGMENTED | Flags::MATE_UNMAPPED;
                if index == 0 {
                    flags |= Flags::FIRST_SEGMENT;
                } else {
                    flags |= Flags::LAST_SEGMENT;
                }
            }
            *rec.flags_mut() = flags;

            // Add tags using the mutable data accessor
            let data = rec.data_mut();

            // Read group
            data.insert(Tag::READ_GROUP, Value::String(rgid_bs.as_bstr()).try_into()?);

            // Cell barcode
            if !cell_barcode_bs.is_empty() {
                Self::add_tag(data, self.cell_tag.as_bytes(), &cell_barcode_bs.clone()).unwrap();
            }

            if !cell_quals_bs.is_empty() {
                if let Some(ref qual_tag) = self.cell_qual_tag {
                    Self::add_tag(data, qual_tag.as_bytes(), &cell_quals_bs.clone()).unwrap();
                }
            }

            // Sample barcode
            if !sample_barcode_bs.is_empty() {
                data.insert(
                    Tag::SAMPLE_BARCODE_SEQUENCE,
                    Value::String(sample_barcode_bs.as_bstr()).try_into()?,
                );
            }

            if self.store_sample_barcode_qualities && !sample_quals_bs.is_empty() {
                data.insert(
                    Tag::SAMPLE_BARCODE_QUALITY_SCORES,
                    Value::String(sample_quals_bs.as_bstr()).try_into()?,
                );
            }

            // UMI
            if !final_umi_bs.is_empty() {
                Self::add_tag(data, self.umi_tag.as_bytes(), &final_umi_bs.clone()).unwrap();

                // Only add UMI qualities if not extracted from read names
                if umi_from_name.is_none() && !umi_qual_bs.is_empty() {
                    if let Some(ref qual_tag) = self.umi_qual_tag {
                        Self::add_tag(data, qual_tag.as_bytes(), &umi_qual_bs.clone()).unwrap();
                    }
                }
            }

            records.push(rec);
        }

        Ok(records)
    }
}

impl Command for Extract {
    fn execute(&self) -> Result<()> {
        // Validate inputs
        self.validate()?;

        let read_structures = self.get_read_structures()?;

        // Open input FASTQs
        let fgio = Io::new(5, BUFFER_SIZE);
        let fq_readers: Vec<Box<dyn BufRead + Send>> =
            self.inputs.iter().map(|p| fgio.new_reader(p)).collect::<Result<Vec<_>, _>>()?;

        let fq_sources: Vec<FastqReader<Box<dyn BufRead + Send>>> =
            fq_readers.into_iter().map(|fq| FastqReader::with_capacity(fq, BUFFER_SIZE)).collect();

        // Detect quality encoding from first 400 records
        let mut sample_quals = Vec::new();
        let temp_reader =
            FastqReader::with_capacity(fgio.new_reader(&self.inputs[0])?, BUFFER_SIZE);
        for (i, result) in temp_reader.into_records().enumerate() {
            if i >= QUALITY_DETECTION_SAMPLE_SIZE {
                break;
            }
            if let Ok(rec) = result {
                sample_quals.push(rec.qual().to_vec());
            }
        }

        let encoding = QualityEncoding::detect(&sample_quals)?;

        // Create iterators
        let mut fq_iterators: Vec<ReadSetIterator> = fq_sources
            .into_iter()
            .zip(read_structures.iter())
            .map(|(source, rs)| ReadSetIterator::new(rs.clone(), source, Vec::new()))
            .collect();

        // Create header and output writer
        let header = self.create_header()?;
        let output_file = File::create(&self.output)?;
        let mut writer = sam::io::Writer::new(BufWriter::new(output_file));
        writer.write_header(&header)?;

        // Process records
        loop {
            let mut next_read_sets = Vec::with_capacity(fq_iterators.len());
            for iter in &mut fq_iterators {
                if let Some(rec) = iter.next() {
                    next_read_sets.push(rec);
                } else {
                    break;
                }
            }

            if next_read_sets.is_empty() {
                break;
            }

            ensure!(next_read_sets.len() == fq_iterators.len(), "FASTQ sources out of sync");

            let read_set = ReadSet::combine_readsets(next_read_sets);
            let sam_records = self.make_sam_records(&read_set, encoding)?;

            for record in sam_records {
                writer.write_alignment_record(&header, &record)?;
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::TempDir;

    fn create_fastq(dir: &TempDir, name: &str, records: &[(&str, &str, &str)]) -> PathBuf {
        let path = dir.path().join(name);
        let mut file = File::create(&path).unwrap();
        for (name, seq, qual) in records {
            writeln!(file, "@{name}").unwrap();
            writeln!(file, "{seq}").unwrap();
            writeln!(file, "+").unwrap();
            writeln!(file, "{qual}").unwrap();
        }
        path
    }

    #[test]
    fn test_single_fastq_no_read_structure() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(
            &tmp,
            "r1.fq",
            &[("q1", "AAAAAAAAAA", "=========="), ("q2", "CCCCCCCCCC", "##########")],
        );
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1],
            output,
            read_structures: vec![],
            umi_tag: "RX".to_string(),
            umi_qual_tag: None,
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            read_group_id: "A".to_string(),
            sample: "foo".to_string(),
            library: "bar".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
        };

        extract.execute().unwrap();
    }

    #[test]
    fn test_paired_end_with_inline_umi() {
        let tmp = TempDir::new().unwrap();
        let r1 = create_fastq(&tmp, "r1.fq", &[("q1", "ACGTAAAAAA", "IIII======")]);
        let r2 = create_fastq(&tmp, "r2.fq", &[("q1", "CCCCCCCCCC", "##########")]);
        let output = tmp.path().join("output.bam");

        let extract = Extract {
            inputs: vec![r1, r2],
            output,
            read_structures: vec![
                ReadStructure::from_str("4M+T").unwrap(),
                ReadStructure::from_str("+T").unwrap(),
            ],
            umi_tag: "RX".to_string(),
            umi_qual_tag: Some("QX".to_string()),
            cell_tag: "CB".to_string(),
            cell_qual_tag: None,
            store_sample_barcode_qualities: false,
            extract_umis_from_read_names: false,
            read_group_id: "A".to_string(),
            sample: "s".to_string(),
            library: "l".to_string(),
            barcode: None,
            platform: "illumina".to_string(),
            platform_unit: None,
            platform_model: None,
            sequencing_center: None,
            predicted_insert_size: None,
            description: None,
            comment: vec![],
            run_date: None,
        };

        extract.execute().unwrap();
    }
}
