pub mod barcode_matching;
pub mod samples;

use bio_seq::prelude::*;
use lazy_static::lazy_static;

pub const DNA_BASES: [u8; 5] = *b"ACGTN";
pub const IUPAC_BASES: [u8; 15] = *b"ACGTMRWSYKVHDBN";

lazy_static! {
    pub static ref BASE_A: usize = 1;
    pub static ref BASE_C: usize = 2;
    pub static ref BASE_G: usize = 4;
    pub static ref BASE_T: usize = 8;
    pub static ref BASE_N: usize = 15;
    pub static ref DNA_MASKS: [u8; 256] = {
        let mut masks = [0; 256];
        let (a, c, g, t) = (1, 2, 4, 8);
        masks['A' as usize] = a;
        masks['C' as usize] = c;
        masks['G' as usize] = g;
        masks['T' as usize] = t;
        masks['U' as usize] = t;
        masks['N' as usize] = a | c | g | t;
        masks
    };
    pub static ref IUPAC_MASKS: [u8; 256] = {
        let mut masks = [0; 256];
        let (a, c, g, t) = (1, 2, 4, 8);
        masks['A' as usize] = a;
        masks['C' as usize] = c;
        masks['G' as usize] = g;
        masks['T' as usize] = t;
        masks['U' as usize] = t;
        masks['M' as usize] = a | c;
        masks['R' as usize] = a | g;
        masks['W' as usize] = a | t;
        masks['S' as usize] = c | g;
        masks['Y' as usize] = c | t;
        masks['K' as usize] = g | t;
        masks['V' as usize] = a | c | g;
        masks['H' as usize] = a | c | t;
        masks['D' as usize] = a | g | t;
        masks['B' as usize] = c | g | t;
        masks['N' as usize] = a | c | g | t;
        masks
    };
}

#[must_use]
fn encode_observed(bases: &[u8]) -> Seq<Iupac> {
    bases
        .iter()
        .map(|b| {
            if byte_is_nocall(*b) || !is_valid_base(*b) {
                b'-' // so that Ns in the observed read match nothing, even other Ns
            } else if *b == b'U' {
                b'T'
            } else {
                *b
            }
        })
        .collect::<Vec<_>>()
        .try_into()
        .unwrap()
}

#[must_use]
fn encode_expected(bases: &[u8]) -> Seq<Iupac> {
    bases
        .iter()
        .map(|b| {
            if byte_is_nocall(*b) || !is_valid_iupac(*b) {
                b'N'
            } else if *b == b'U' {
                b'T'
            } else {
                *b
            }
        })
        .collect::<Vec<_>>()
        .try_into()
        .unwrap()
}

/// Checks whether a given u8 byte is a "No-call"-ed base, signified by the bytes 'N', 'n' and '.'
fn byte_is_nocall(byte: u8) -> bool {
    byte == b'N' || byte == b'n' || byte == b'.'
}

/// Checks whether a provided byte is an IUPAC or nocall.
fn is_valid_base(byte: u8) -> bool {
    DNA_MASKS[byte as usize] != 0 || byte_is_nocall(byte)
}

/// Checks whether a provided byte is an IUPAC or nocall.
fn is_valid_iupac(byte: u8) -> bool {
    IUPAC_MASKS[byte as usize] != 0 || byte_is_nocall(byte)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ############################################################################################
    // Test byte_is_no_call
    // ############################################################################################
    #[test]
    fn test_byte_is_no_call() {
        assert!(byte_is_nocall(b'N'));
        assert!(byte_is_nocall(b'n'));
        assert!(byte_is_nocall(b'.'));
        assert!(!byte_is_nocall(b'A'));
        assert!(!byte_is_nocall(b'C'));
        assert!(!byte_is_nocall(b'G'));
        assert!(!byte_is_nocall(b'T'));
        assert!(!byte_is_nocall(b'a'));
        assert!(!byte_is_nocall(b'c'));
        assert!(!byte_is_nocall(b'g'));
        assert!(!byte_is_nocall(b't'));
    }

    #[test]
    fn test_is_valid_iupac() {
        assert!(is_valid_iupac(b'N'));
        assert!(is_valid_iupac(b'n'));
        assert!(is_valid_iupac(b'.'));
        assert!(!is_valid_iupac(b'a'));
        assert!(!is_valid_iupac(b'c'));
        assert!(!is_valid_iupac(b'g'));
        assert!(!is_valid_iupac(b't'));
        for base in IUPAC_BASES {
            assert!(is_valid_iupac(base));
        }
    }

    #[test]
    fn test_encode_observed() {
        assert_eq!(encode_observed(b"A").get(0).unwrap(), Iupac::A);
        assert_eq!(encode_observed(b"C").get(0).unwrap(), Iupac::C);
        assert_eq!(encode_observed(b"G").get(0).unwrap(), Iupac::G);
        assert_eq!(encode_observed(b"T").get(0).unwrap(), Iupac::T);
        assert_eq!(encode_observed(b"U").get(0).unwrap(), Iupac::T);
        assert_eq!(encode_observed(b"M").get(0).unwrap(), Iupac::X);
        assert_eq!(encode_observed(b"R").get(0).unwrap(), Iupac::X);
        assert_eq!(encode_observed(b"W").get(0).unwrap(), Iupac::X);
        assert_eq!(encode_observed(b"S").get(0).unwrap(), Iupac::X);
        assert_eq!(encode_observed(b"Y").get(0).unwrap(), Iupac::X);
        assert_eq!(encode_observed(b"K").get(0).unwrap(), Iupac::X);
        assert_eq!(encode_observed(b"V").get(0).unwrap(), Iupac::X);
        assert_eq!(encode_observed(b"H").get(0).unwrap(), Iupac::X);
        assert_eq!(encode_observed(b"D").get(0).unwrap(), Iupac::X);
        assert_eq!(encode_observed(b"B").get(0).unwrap(), Iupac::X);
        assert_eq!(encode_observed(b"N").get(0).unwrap(), Iupac::X);
        assert_eq!(encode_observed(b".").get(0).unwrap(), Iupac::X);
        assert_eq!(encode_observed(b"h").get(0).unwrap(), Iupac::X);
    }

    #[test]
    fn test_encode_expected() {
        assert_eq!(encode_expected(b"A").get(0).unwrap(), Iupac::A);
        assert_eq!(encode_expected(b"C").get(0).unwrap(), Iupac::C);
        assert_eq!(encode_expected(b"G").get(0).unwrap(), Iupac::G);
        assert_eq!(encode_expected(b"T").get(0).unwrap(), Iupac::T);
        assert_eq!(encode_expected(b"U").get(0).unwrap(), Iupac::T);
        assert_eq!(encode_expected(b"M").get(0).unwrap(), Iupac::M);
        assert_eq!(encode_expected(b"R").get(0).unwrap(), Iupac::R);
        assert_eq!(encode_expected(b"W").get(0).unwrap(), Iupac::W);
        assert_eq!(encode_expected(b"S").get(0).unwrap(), Iupac::S);
        assert_eq!(encode_expected(b"Y").get(0).unwrap(), Iupac::Y);
        assert_eq!(encode_expected(b"K").get(0).unwrap(), Iupac::K);
        assert_eq!(encode_expected(b"V").get(0).unwrap(), Iupac::V);
        assert_eq!(encode_expected(b"H").get(0).unwrap(), Iupac::H);
        assert_eq!(encode_expected(b"D").get(0).unwrap(), Iupac::D);
        assert_eq!(encode_expected(b"B").get(0).unwrap(), Iupac::B);
        assert_eq!(encode_expected(b"N").get(0).unwrap(), Iupac::N);
        assert_eq!(encode_expected(b".").get(0).unwrap(), Iupac::N);
        assert_eq!(encode_expected(b"h").get(0).unwrap(), Iupac::N);
    }
}
