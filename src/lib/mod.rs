pub mod barcode_matching;
pub mod base_quality;
pub mod fastq_stats;
pub mod moving_average;
pub mod pair_overlap;
pub mod samples;

/// Checks whether a given u8 byte is a "No-call"-ed base, signified by the bytes 'N', 'n' and '.'
fn byte_is_nocall(byte: u8) -> bool {
    byte == b'N' || byte == b'n' || byte == b'.'
}

/// Checks whether a provided byte is an A, G, C, or T.
fn is_valid_base(byte: u8) -> bool {
    byte == b'A' || byte == b'C' || byte == b'G' || byte == b'T' || byte_is_nocall(byte)
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

    // ############################################################################################
    // Test is_valid_base
    // ############################################################################################
    #[test]
    fn test_is_valid_base() {
        assert!(is_valid_base(b'N'));
        assert!(is_valid_base(b'n'));
        assert!(is_valid_base(b'.'));
        assert!(!is_valid_base(b'a'));
        assert!(!is_valid_base(b'c'));
        assert!(!is_valid_base(b'g'));
        assert!(!is_valid_base(b't'));
        assert!(is_valid_base(b'A'));
        assert!(is_valid_base(b'C'));
        assert!(is_valid_base(b'G'));
        assert!(is_valid_base(b'T'));
    }
}
