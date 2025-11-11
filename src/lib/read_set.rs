use anyhow::{anyhow, ensure, Result};
use bstr::ByteSlice;
use read_structure::ReadStructure;
use read_structure::SegmentType;
use seq_io::fastq::Reader as FastqReader;
use seq_io::fastq::Record;
use std::fmt::Display;
use std::io::{BufRead, Write};
use std::iter::Filter;
use std::slice::Iter;
use std::str::FromStr;

/// Type alias for segment type iter functions, which iterate over the segments of a ``ReadSet``
/// filtering for a specific type.
type SegmentIter<'a> = Filter<Iter<'a, FastqSegment>, fn(&&FastqSegment) -> bool>;

/// The bases and qualities associated with a segment of a FASTQ record.
#[derive(PartialEq, Eq, Debug, Clone)]
pub struct FastqSegment {
    /// bases of the FASTQ subsection
    pub seq: Vec<u8>,
    /// qualities of the FASTQ subsection
    pub quals: Vec<u8>,
    /// the type of segment being stored
    pub segment_type: SegmentType,
}

////////////////////////////////////////////////////////////////////////////////
// ReadSet and it's impls
////////////////////////////////////////////////////////////////////////////////

#[derive(Eq, Hash, PartialEq, Debug, Clone, Copy)]
pub enum SkipReason {
    /// The read had too few bases for the segment.
    TooFewBases,
}

impl Display for SkipReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SkipReason::TooFewBases => write!(f, "Too few bases"),
        }
    }
}

impl FromStr for SkipReason {
    type Err = anyhow::Error;
    fn from_str(string: &str) -> Result<Self, Self::Err> {
        match string {
            "too few bases" | "too-few-bases" | "toofewbases" => Ok(SkipReason::TooFewBases),
            _ => Err(anyhow!("Invalid skip reason: {string}")),
        }
    }
}

/// One unit of FASTQ records separated into their component read segments.
#[derive(PartialEq, Debug, Clone)]
pub struct ReadSet {
    /// Header of the FASTQ record
    pub header: Vec<u8>,
    /// Segments of reads
    pub segments: Vec<FastqSegment>,
    /// The reason the read should be skipped for demultiplexing, or None if it should be
    /// demultiplexed.
    pub skip_reason: Option<SkipReason>,
}

impl ReadSet {
    const PREFIX: u8 = b'@';
    const SPACE: u8 = b' ';
    const COLON: u8 = b':';
    const PLUS: u8 = b'+';

    /// Produces an iterator over references to the template segments stored in this ``ReadSet``.
    pub fn template_segments(&self) -> SegmentIter {
        self.segments.iter().filter(|s| s.segment_type == SegmentType::Template)
    }

    /// Produces an iterator over references to the sample barcode segments stored in this
    /// ``ReadSet``.
    pub fn sample_barcode_segments(&self) -> SegmentIter {
        self.segments.iter().filter(|s| s.segment_type == SegmentType::SampleBarcode)
    }

    /// Produces an iterator over references to the molecular barcode segments stored in this
    /// ``ReadSet``.
    pub fn molecular_barcode_segments(&self) -> SegmentIter {
        self.segments.iter().filter(|s| s.segment_type == SegmentType::MolecularBarcode)
    }

    /// Produces an iterator over references to the cell barcode segments stored in this
    /// ``ReadSet``.
    pub fn cell_barcode_segments(&self) -> impl Iterator<Item = &FastqSegment> {
        self.segments.iter().filter(|s| s.segment_type == SegmentType::CellularBarcode)
    }

    /// Generates the sample barcode sequence for this read set and returns it as a Vec of bytes.
    #[must_use]
    pub fn sample_barcode_sequence(&self) -> Vec<u8> {
        self.sample_barcode_segments().flat_map(|s| &s.seq).copied().collect()
    }

    /// Combines ``ReadSet`` structs together into a single ``ReadSet``
    ///
    /// # Panics
    /// If the total # of segments is zero.
    #[must_use]
    pub fn combine_readsets(readsets: Vec<Self>) -> Self {
        let total_segments: usize = readsets.iter().map(|s| s.segments.len()).sum();
        assert!(total_segments > 0, "Cannot call combine readsets on an empty vec!");

        let mut readset_iter = readsets.into_iter();
        let mut first = readset_iter.next().expect("Cannot call combine readsets on an empty vec!");
        first.segments.reserve_exact(total_segments - first.segments.len());

        for next_readset in readset_iter {
            first.segments.extend(next_readset.segments);
        }

        first
    }

    /// Writes the FASTQ header to the given writer.  Substitutes in the given read number into
    /// the comment section.  Also adds in the UMI(s) from the UMI segments and the sample barcodes
    /// into the appropriate places in the header.
    ///
    /// UMI and sample barcode segments and (separately) concatenated with `+`s in between. If
    /// there is an existing UMI or sample barcode in the header, the segments are appended after
    /// adding a `+`, else the segments are inserted.
    ///
    /// Supports headers that are just the `name` segment, or `name comment`.  The name must have
    /// at most 8 colon-separated parts.  If there are seven or fewer parts in the name, the
    /// UMI is appended as the last part.  If there are eight, the eighth is assumed to be an
    /// existing UMI, and any UMI is appended.
    ///
    /// If the comment is present is must have exactly four colon-separated parts.
    ///
    /// Format of the header is:
    ///   @name comment
    /// Where
    ///   name = @<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI>
    ///   comment = <read>:<is filtered>:<control number>:<index>
    ///
    /// # Errors
    /// - a read name has more than 8 segments
    pub fn write_header<W: Write>(&self, writer: &mut W, read_num: usize) -> Result<()> {
        Self::write_header_internal(
            writer,
            read_num,
            self.header.as_slice(),
            self.sample_barcode_segments(),
            self.molecular_barcode_segments(),
        )
    }

    fn write_header_internal<W: Write>(
        writer: &mut W,
        read_num: usize,
        header: &[u8],
        sample_barcode_segments: SegmentIter,
        mut molecular_barcode_segments: SegmentIter,
    ) -> Result<()> {
        // Extract the name and optionally the comment
        let (name, comment) = match header.find_byte(Self::SPACE) {
            Some(x) => (&header[0..x], Some(&header[x + 1..])),
            None => (header, None),
        };

        writer.write_all(&[Self::PREFIX])?;

        // Handle the 'name' component of the header.  If we don't have any UMI segments
        // we can emit the name part as is.  Otherwise we need to append the UMIs to the name.
        if let Some(first_seg) = molecular_barcode_segments.next() {
            let sep_count = bytecount::count(name, Self::COLON);
            ensure!(
                sep_count <= 7,
                "Can't handle read name with more than 8 segments: {}",
                String::from_utf8(header.to_vec())?
            );

            writer.write_all(name)?;
            if sep_count == 7 {
                // UMI already present, append to it with a UMI separator first
                writer.write_all(&[Self::PLUS])?;
            } else {
                // UMI not present yet, insert a minor separator
                writer.write_all(&[Self::COLON])?;
            }

            writer.write_all(first_seg.seq.as_slice())?;
            // Append all the UMI segments with pluses in between.
            for seg in molecular_barcode_segments {
                writer.write_all(&[Self::PLUS])?;
                writer.write_all(seg.seq.as_slice())?;
            }
        } else {
            writer.write_all(name)?;
        }

        writer.write_all(&[Self::SPACE])?;

        // Then the 'comment' section
        match comment {
            None => {
                // If no pre-existing comment, assume the read is a passing filter, non-control
                // read and generate a comment for it (sample barcode is added below).
                write!(writer, "{read_num}:N:0:")?;
            }
            Some(chars) => {
                // Else check it's a 4-part name... fix the read number at the front and
                // check to see if there's a real sample barcode on the back
                let sep_count = bytecount::count(chars, Self::COLON);
                if sep_count < 3 {
                    writer.write_all(chars)?;
                    if *chars.last().unwrap() != Self::COLON {
                        writer.write_all(&[Self::COLON])?;
                    }
                } else {
                    ensure!(
                        sep_count == 3,
                        "Comment in did not have 4 segments: {}",
                        String::from_utf8(header.to_vec())?
                    );
                    let first_colon_idx = chars.iter().position(|ch| *ch == Self::COLON).unwrap();

                    // Illumina, in the unmatched FASTQs, can place a "0" in the index position, sigh
                    let remainder = if chars.last().unwrap().is_ascii_digit() {
                        &chars[first_colon_idx + 1..chars.len() - 1]
                    } else {
                        &chars[first_colon_idx + 1..chars.len()]
                    };

                    write!(writer, "{read_num}:")?;
                    writer.write_all(remainder)?;

                    if *remainder.last().unwrap() != Self::COLON {
                        writer.write_all(&[Self::PLUS])?;
                    }
                }
            }
        }

        // Append all the sample barcode segments to the new comment
        for (idx, seg) in sample_barcode_segments.enumerate() {
            if idx > 0 {
                writer.write_all(&[Self::PLUS])?;
            }
            writer.write_all(seg.seq.as_slice())?;
        }

        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////
// ReadSetIterator and it's impls
////////////////////////////////////////////////////////////////////////////////

/// A struct for iterating over the records in multiple FASTQ files simultaneously, destructuring
/// them according to the provided read structures, yielding ``ReadSet`` objects on each iteration.
pub struct ReadSetIterator {
    /// Read structure object describing the structure of the reads in this file.
    read_structure: ReadStructure,
    /// The file containing FASTQ records.
    source: FastqReader<Box<dyn BufRead + Send>>,
    /// Valid reasons for skipping reads, otherwise panic!
    skip_reasons: Vec<SkipReason>,
}

impl Iterator for ReadSetIterator {
    type Item = ReadSet;

    fn next(&mut self) -> Option<Self::Item> {
        if let Some(rec) = self.source.next() {
            let mut segments = Vec::with_capacity(self.read_structure.number_of_segments());
            let next_fq_rec = rec.expect("Unexpected error parsing FASTQs");
            let read_name = next_fq_rec.head();
            let bases = next_fq_rec.seq();
            let quals = next_fq_rec.qual();

            // Check that we have enough bases for this segment. For variable length segments,
            // we must have at least one base.
            let min_len = self.read_structure.iter().map(|s| s.length().unwrap_or(1)).sum();
            if bases.len() < min_len {
                if self.skip_reasons.contains(&SkipReason::TooFewBases) {
                    return Some(ReadSet {
                        header: read_name.to_vec(),
                        segments: vec![],
                        skip_reason: Some(SkipReason::TooFewBases),
                    });
                }
                panic!(
                    "Read {} had too few bases to demux {} vs. {} needed in read structure {}.",
                    String::from_utf8(read_name.to_vec()).unwrap(),
                    bases.len(),
                    min_len,
                    self.read_structure
                );
            }

            for (read_segment_index, read_segment) in self.read_structure.iter().enumerate() {
                // Extract the bases and qualities for this segment
                let (seq, quals) =
                    read_segment.extract_bases_and_quals(bases, quals).unwrap_or_else(|e| {
                        panic!(
                            "Error extracting bases (len: {}) or quals (len: {}) for the {}th read segment ({}) in read structure ({}) from FASTQ record with name {}; {}",
                            bases.len(),
                            quals.len(),
                            read_segment_index,
                            read_segment,
                            self.read_structure,
                            String::from_utf8(read_name.to_vec()).unwrap(),
                            e
                        )
                    });
                // Add a new FastqSegment
                segments.push(FastqSegment {
                    seq: seq.to_vec(),
                    quals: quals.to_vec(),
                    segment_type: read_segment.kind,
                });
            }
            Some(ReadSet { header: read_name.to_vec(), segments, skip_reason: None })
        } else {
            None
        }
    }
}

impl ReadSetIterator {
    /// Instantiates a new iterator over the read sets for a set of FASTQs with defined read
    /// structures
    #[must_use]
    pub fn new(
        read_structure: ReadStructure,
        source: FastqReader<Box<dyn BufRead + Send>>,
        skip_reasons: Vec<SkipReason>,
    ) -> Self {
        Self { read_structure, source, skip_reasons }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn read_set(segments: Vec<FastqSegment>) -> ReadSet {
        ReadSet { header: "NOT_IMPORTANT".as_bytes().to_owned(), segments, skip_reason: None }
    }

    ////////////////////////////////////////////////////////////////////////////
    // Tests for ReadSet::write_header_internal
    ////////////////////////////////////////////////////////////////////////////

    fn seg(bases: &[u8], segment_type: SegmentType) -> FastqSegment {
        let quals = vec![b'#'; bases.len()];
        FastqSegment { seq: bases.to_vec(), quals, segment_type }
    }

    #[test]
    fn test_write_header_standard_no_umi() {
        let mut out = Vec::new();
        let header = b"inst:123:ABCDE:1:204:1022:2108 1:N:0:0";
        let barcode_segs =
            [seg(b"ACGT", SegmentType::SampleBarcode), seg(b"GGTT", SegmentType::SampleBarcode)];
        let umi_segs = [];
        let expected = "@inst:123:ABCDE:1:204:1022:2108 1:N:0:ACGT+GGTT".to_string();
        ReadSet::write_header_internal(
            &mut out,
            1,
            header,
            barcode_segs.iter().filter(|_| true),
            umi_segs.iter().filter(|_| true),
        )
        .unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    #[test]
    fn test_write_header_standard_with_umi() {
        let mut out = Vec::new();
        let header = b"inst:123:ABCDE:1:204:1022:2108 1:Y:0:0";
        let barcode_segs =
            [seg(b"ACGT", SegmentType::SampleBarcode), seg(b"GGTT", SegmentType::SampleBarcode)];
        let umi_segs = [seg(b"AACCGGTT", SegmentType::MolecularBarcode)];
        let expected = "@inst:123:ABCDE:1:204:1022:2108:AACCGGTT 2:Y:0:ACGT+GGTT".to_string();
        ReadSet::write_header_internal(
            &mut out,
            2,
            header,
            barcode_segs.iter().filter(|_| true),
            umi_segs.iter().filter(|_| true),
        )
        .unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    #[test]
    fn test_write_header_append_barcode_and_umi() {
        let mut out = Vec::new();
        let header = b"inst:123:ABCDE:1:204:1022:2108:AAAA 1:Y:0:TTTT";
        let barcode_segs =
            [seg(b"ACGT", SegmentType::SampleBarcode), seg(b"GGTT", SegmentType::SampleBarcode)];
        let umi_segs = [seg(b"AACCGGTT", SegmentType::MolecularBarcode)];
        let expected =
            "@inst:123:ABCDE:1:204:1022:2108:AAAA+AACCGGTT 2:Y:0:TTTT+ACGT+GGTT".to_string();
        ReadSet::write_header_internal(
            &mut out,
            2,
            header,
            barcode_segs.iter().filter(|_| true),
            umi_segs.iter().filter(|_| true),
        )
        .unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    #[test]
    fn test_write_header_short_name_no_comment() {
        let mut out = Vec::new();
        let header = b"q1";
        let barcode_segs =
            [seg(b"ACGT", SegmentType::SampleBarcode), seg(b"GGTT", SegmentType::SampleBarcode)];
        let umi_segs = [seg(b"AACCGGTT", SegmentType::MolecularBarcode)];
        let expected = "@q1:AACCGGTT 1:N:0:ACGT+GGTT".to_string();
        ReadSet::write_header_internal(
            &mut out,
            1,
            header,
            barcode_segs.iter().filter(|_| true),
            umi_segs.iter().filter(|_| true),
        )
        .unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    #[test]
    #[should_panic(expected = "8 segments")]
    fn test_write_header_name_too_many_parts() {
        let mut out = Vec::new();
        let header = b"q1:1:2:3:4:5:6:7:8:9:10";
        let barcode_segs =
            [seg(b"ACGT", SegmentType::SampleBarcode), seg(b"GGTT", SegmentType::SampleBarcode)];
        let umi_segs = [seg(b"AACCGGTT", SegmentType::MolecularBarcode)];
        ReadSet::write_header_internal(
            &mut out,
            1,
            header,
            barcode_segs.iter().filter(|_| true),
            umi_segs.iter().filter(|_| true),
        )
        .unwrap();
    }

    #[test]
    fn test_write_header_comment_too_few_parts() {
        let mut out = Vec::new();
        let header = b"q1 0:0";
        let barcode_segs =
            [seg(b"ACGT", SegmentType::SampleBarcode), seg(b"GGTT", SegmentType::SampleBarcode)];
        let umi_segs = [seg(b"AACCGGTT", SegmentType::MolecularBarcode)];
        let expected = "@q1:AACCGGTT 0:0:ACGT+GGTT".to_string();
        ReadSet::write_header_internal(
            &mut out,
            1,
            header,
            barcode_segs.iter().filter(|_| true),
            umi_segs.iter().filter(|_| true),
        )
        .unwrap();
        assert_eq!(String::from_utf8(out).unwrap(), expected);
    }

    // ############################################################################################
    // Test other ``ReadSet`` functions
    // ############################################################################################

    #[test]
    fn test_sample_barcode_sequence() {
        let segments = vec![
            seg("AGCT".as_bytes(), SegmentType::Template),
            seg("GATA".as_bytes(), SegmentType::SampleBarcode),
            seg("CAC".as_bytes(), SegmentType::SampleBarcode),
            seg("GACCCC".as_bytes(), SegmentType::MolecularBarcode),
        ];

        let read_set = read_set(segments);

        assert_eq!(read_set.sample_barcode_sequence(), "GATACAC".as_bytes().to_owned());
    }
    #[test]
    fn test_template_segments() {
        let segments = vec![
            seg("AGCT".as_bytes(), SegmentType::SampleBarcode),
            seg("GATA".as_bytes(), SegmentType::Template),
            seg("CAC".as_bytes(), SegmentType::Template),
            seg("GACCCC".as_bytes(), SegmentType::MolecularBarcode),
        ];
        let expected = vec![
            seg("GATA".as_bytes(), SegmentType::Template),
            seg("CAC".as_bytes(), SegmentType::Template),
        ];

        let read_set = read_set(segments);

        assert_eq!(expected, read_set.template_segments().cloned().collect::<Vec<_>>());
    }
    #[test]
    fn test_sample_barcode_segments() {
        let segments = vec![
            seg("AGCT".as_bytes(), SegmentType::Template),
            seg("GATA".as_bytes(), SegmentType::SampleBarcode),
            seg("CAC".as_bytes(), SegmentType::SampleBarcode),
            seg("GACCCC".as_bytes(), SegmentType::MolecularBarcode),
        ];
        let expected = vec![
            seg("GATA".as_bytes(), SegmentType::SampleBarcode),
            seg("CAC".as_bytes(), SegmentType::SampleBarcode),
        ];

        let read_set = read_set(segments);

        assert_eq!(expected, read_set.sample_barcode_segments().cloned().collect::<Vec<_>>());
    }
    #[test]
    fn test_molecular_barcode_segments() {
        let segments = vec![
            seg("AGCT".as_bytes(), SegmentType::Template),
            seg("GATA".as_bytes(), SegmentType::MolecularBarcode),
            seg("CAC".as_bytes(), SegmentType::MolecularBarcode),
            seg("GACCCC".as_bytes(), SegmentType::SampleBarcode),
        ];
        let expected = vec![
            seg("GATA".as_bytes(), SegmentType::MolecularBarcode),
            seg("CAC".as_bytes(), SegmentType::MolecularBarcode),
        ];

        let read_set = read_set(segments);

        assert_eq!(expected, read_set.molecular_barcode_segments().cloned().collect::<Vec<_>>());
    }

    #[test]
    fn test_combine_readsets() {
        let segments1 = vec![
            seg("A".as_bytes(), SegmentType::Template),
            seg("G".as_bytes(), SegmentType::Template),
            seg("C".as_bytes(), SegmentType::MolecularBarcode),
            seg("T".as_bytes(), SegmentType::SampleBarcode),
        ];
        let read_set1 = read_set(segments1.clone());
        let segments2 = vec![
            seg("AA".as_bytes(), SegmentType::Template),
            seg("AG".as_bytes(), SegmentType::Template),
            seg("AC".as_bytes(), SegmentType::MolecularBarcode),
            seg("AT".as_bytes(), SegmentType::SampleBarcode),
        ];
        let read_set2 = read_set(segments2.clone());
        let segments3 = vec![
            seg("AAA".as_bytes(), SegmentType::Template),
            seg("AAG".as_bytes(), SegmentType::Template),
            seg("AAC".as_bytes(), SegmentType::MolecularBarcode),
            seg("AAT".as_bytes(), SegmentType::SampleBarcode),
        ];
        let read_set3 = read_set(segments3.clone());

        let mut expected_segments = Vec::new();
        expected_segments.extend(segments1);
        expected_segments.extend(segments2);
        expected_segments.extend(segments3);
        let expected = read_set(expected_segments);

        let result = ReadSet::combine_readsets(vec![read_set1, read_set2, read_set3]);

        assert_eq!(result, expected);
    }

    #[test]
    #[should_panic(expected = "Cannot call combine readsets on an empty vec!")]
    fn test_combine_readsets_fails_on_empty_vector() {
        let _result = ReadSet::combine_readsets(Vec::new());
    }
}
