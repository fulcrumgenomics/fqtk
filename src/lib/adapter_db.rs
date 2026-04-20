//! Built-in 3' adapter sequence presets for common short-read sequencing kits.
//!
//! Sources (consulted April 2026):
//! - Illumina TruSeq — adapter sequence document (rev H); widely used across Illumina WGS,
//!   WES, ChIP-seq, RNA-seq.
//! - Illumina Nextera / Tn5 — same sequence on both mates (ligated via transposition).
//! - Illumina small-RNA / TruSeq Small-RNA — 3' adapter only; single-end by convention.
//! - Element Biosciences AVITI Elevate — Library Sequences Technical Note LT-00014 (Mar 2024).
//! - MGI / DNBSEQ / BGI — documented in MGI's oligo primer docs and cross-referenced in
//!   `OpenGene/fastp` issue #259 and `linsalrob/mgi-adapters`.
//!
//! Element Freestyle and Ultima Illumina-compatible libraries typically use the `truseq`
//! preset; Ultima's native adapters and Singular G4's S1/S2 adapters are not publicly
//! disclosed and are therefore not included.

/// A named adapter preset. When `seq_r2` is `None` the preset covers R1-only libraries
/// (historically single-end such as small-RNA-seq) or the R1 adapter is also used on R2.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KitAdapter {
    pub name: &'static str,
    pub seq_r1: &'static [u8],
    pub seq_r2: Option<&'static [u8]>,
}

/// Illumina TruSeq (standard and TruSeq LT/HT): used across most Illumina WGS/WES/RNA-seq.
pub const TRUSEQ: KitAdapter = KitAdapter {
    name: "truseq",
    seq_r1: b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    seq_r2: Some(b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"),
};

/// Illumina Nextera / Tn5: used for Nextera DNA, ATAC-seq, Nextera XT.
/// Same adapter ligated on both mates by the transposition reaction.
pub const NEXTERA: KitAdapter = KitAdapter {
    name: "nextera",
    seq_r1: b"CTGTCTCTTATACACATCT",
    seq_r2: Some(b"CTGTCTCTTATACACATCT"),
};

/// Illumina TruSeq Small-RNA / miRNA — 3' adapter. Single-end by convention;
/// no canonical R2 adapter.
pub const SMALL_RNA: KitAdapter =
    KitAdapter { name: "small-rna", seq_r1: b"TGGAATTCTCGGGTGCCAAGG", seq_r2: None };

/// Element Biosciences AVITI "Elevate" (native) adapters.
/// AVITI Adept / Cloudbreak Freestyle libraries use standard Illumina (TruSeq) adapters
/// and should therefore use the `truseq` preset instead of this one.
pub const AVITI: KitAdapter = KitAdapter {
    name: "aviti",
    seq_r1: b"ATGTCGGAAGGTGTGCAGGCTACCGCTTGTCAACT",
    seq_r2: Some(b"ATGTCGGAAGGTGTCTGGTGAGCCAATC"),
};

/// MGI / DNBSEQ / BGI native adapters. Note that many MGI library preps sold in the
/// US/EU actually use Illumina-style adapters — `truseq` may be the correct preset
/// for such libraries despite the instrument.
pub const MGI: KitAdapter = KitAdapter {
    name: "mgi",
    seq_r1: b"AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA",
    seq_r2: Some(b"AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG"),
};

/// Every published kit preset, in stable order.
pub const ALL_KITS: &[KitAdapter] = &[TRUSEQ, NEXTERA, SMALL_RNA, AVITI, MGI];

/// Returns the named kit or `None` if `name` is not recognized. The special alias
/// `"all"` is not a single kit and returns `None` here; callers wanting to union every
/// kit should use [`expand_kit_name`] instead.
pub fn kit_by_name(name: &str) -> Option<&'static KitAdapter> {
    ALL_KITS.iter().find(|k| k.name.eq_ignore_ascii_case(name))
}

/// Expands a single user-supplied kit spelling into zero or more `KitAdapter` slices.
/// Recognizes kit names and the special `"all"` alias.
pub fn expand_kit_name(name: &str) -> Option<&'static [KitAdapter]> {
    if name.eq_ignore_ascii_case("all") {
        return Some(ALL_KITS);
    }
    if name.eq_ignore_ascii_case("dnbseq") {
        return Some(std::slice::from_ref(&MGI));
    }
    kit_by_name(name).map(std::slice::from_ref)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kit_by_name_known() {
        assert_eq!(kit_by_name("truseq"), Some(&TRUSEQ));
        assert_eq!(kit_by_name("TRUSEQ"), Some(&TRUSEQ));
        assert_eq!(kit_by_name("nextera"), Some(&NEXTERA));
        assert_eq!(kit_by_name("small-rna"), Some(&SMALL_RNA));
        assert_eq!(kit_by_name("aviti"), Some(&AVITI));
        assert_eq!(kit_by_name("mgi"), Some(&MGI));
    }

    #[test]
    fn kit_by_name_unknown() {
        assert_eq!(kit_by_name("notakit"), None);
        // "all" is a collective alias; not a single kit.
        assert_eq!(kit_by_name("all"), None);
    }

    #[test]
    fn expand_kit_name_all_unions_every_kit() {
        let expanded = expand_kit_name("all").unwrap();
        assert_eq!(expanded.len(), ALL_KITS.len());
    }

    #[test]
    fn expand_kit_name_dnbseq_alias() {
        let expanded = expand_kit_name("dnbseq").unwrap();
        assert_eq!(expanded, std::slice::from_ref(&MGI));
    }

    #[test]
    fn expand_kit_name_unknown_returns_none() {
        assert_eq!(expand_kit_name("mystery"), None);
    }

    #[test]
    fn adapters_are_uppercase_acgt_only() {
        for kit in ALL_KITS {
            for seq in [Some(kit.seq_r1), kit.seq_r2].into_iter().flatten() {
                for &b in seq {
                    assert!(
                        matches!(b, b'A' | b'C' | b'G' | b'T'),
                        "kit {} has non-ACGT byte {:?} in adapter",
                        kit.name,
                        b as char,
                    );
                }
            }
        }
    }
}
