pub mod barcode_matching;
pub mod bitenc;
pub mod samples;

use crate::bitenc::BitEnc;
use std::sync::LazyLock;

pub const DNA_BASES: [u8; 5] = *b"ACGTN";
pub const IUPAC_BASES: [u8; 15] = *b"ACGTMRWSYKVHDBN";

pub static BASE_A: LazyLock<usize> = LazyLock::new(|| 1);
pub static BASE_C: LazyLock<usize> = LazyLock::new(|| 2);
pub static BASE_G: LazyLock<usize> = LazyLock::new(|| 4);
pub static BASE_T: LazyLock<usize> = LazyLock::new(|| 8);
pub static BASE_N: LazyLock<usize> = LazyLock::new(|| 15);
pub static DNA_MASKS: LazyLock<[u8; 256]> = LazyLock::new(|| {
    let mut masks = [0; 256];
    let (a, c, g, t) = (1, 2, 4, 8);
    masks['A' as usize] = a;
    masks['C' as usize] = c;
    masks['G' as usize] = g;
    masks['T' as usize] = t;
    masks['U' as usize] = t;
    masks['N' as usize] = a | c | g | t;
    masks
});
pub static IUPAC_MASKS: LazyLock<[u8; 256]> = LazyLock::new(|| {
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
});

#[must_use]
pub fn encode(bases: &[u8]) -> BitEnc {
    let mut vec = BitEnc::with_capacity(4, bases.len());
    for base in bases {
        let bit: u8 = if byte_is_nocall(*base) {
            IUPAC_MASKS[b'N' as usize]
        } else {
            let value = base.to_ascii_uppercase() as usize;
            if value < 256 { IUPAC_MASKS[value] } else { 0 }
        };
        vec.push(bit);
    }
    vec
}

/// Decodes a DNA/IUPAC encoded squence.
///
/// # Panics
/// when an invalid encoding is provided.
#[must_use]
pub fn decode(bases: &BitEnc) -> String {
    let mut result = String::new();
    for bit in bases.iter() {
        let mut found = false;
        for base in &IUPAC_BASES {
            if IUPAC_MASKS[*base as usize] == bit {
                result.push(*base as char);
                found = true;
                break;
            }
        }
        assert!(found, "Invalid bit mask for base: {bit}");
    }
    result
}

/// Checks whether a given u8 byte is a "No-call"-ed base, signified by the bytes 'N', 'n' and '.'
fn byte_is_nocall(byte: u8) -> bool {
    byte == b'N' || byte == b'n' || byte == b'.'
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
    fn test_encode_dna_bases() {
        for base in DNA_BASES {
            let actual: u8 = encode(&[base]).get(0).unwrap();
            assert_eq!(actual, IUPAC_MASKS[base as usize]);
        }
    }

    #[test]
    fn test_decode_dna_bases() {
        assert_eq!(DNA_BASES, decode(&encode(&DNA_BASES)).as_bytes());
    }

    #[test]
    fn test_encode_iupac_bases() {
        for base in IUPAC_BASES {
            let actual: u8 = encode(&[base]).get(0).unwrap();
            assert_eq!(actual, IUPAC_MASKS[base as usize]);
        }
    }

    #[test]
    fn test_decode_iupac_bases() {
        assert_eq!(IUPAC_BASES, decode(&encode(&IUPAC_BASES)).as_bytes());
    }

    #[test]
    fn test_encode_no_calls() {
        for base in [b'N', b'n', b'.'] {
            let actual: u8 = encode(&[base]).get(0).unwrap();
            assert_eq!(actual as usize, 15);
        }
    }

    #[test]
    fn test_decode_no_calls() {
        let bases = [b'N', b'n', b'.'];
        assert_eq!(decode(&encode(&bases)), "NNN".to_string());
    }
}
