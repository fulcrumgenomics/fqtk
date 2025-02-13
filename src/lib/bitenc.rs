// Copyright 2014-2016 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.
// Adapted by Nils Homer (2025)

//! A fixed-width bit encoding implementation. This allows to store a sequence of values over
//! a reduced alphabet by packing them bit-encoded into a sequence of bytes.
//!
//! Similar behaviour can be achieved using a `PackedVec` from the [packedvec](https://docs.rs/packedvec) crate.
//!
//! # Example
//!
//! ```
//! use fqtk_lib::bitenc::BitEnc;
//! let mut bitenc = BitEnc::new(2);
//! bitenc.push(0);
//! bitenc.push(2);
//! bitenc.push(1);
//! let values: Vec<u8> = bitenc.iter().collect();
//! assert_eq!(values, [0, 2, 1]);
//! ```
//!

use serde::Deserialize;
use serde::Serialize;

/// A sequence of bitencoded values.
///
/// Space complexity: O(⌈(n * width) / k⌉) * 32 bit, where n is the length of the input
/// sequence and `k = 32 - (32 % width)`  is the number of bits in each
/// 32-bit block that can be used to store values.
/// For values that are not a divider of 32, some bits will remain unused.
/// For example for `width = 7` only `4 * 7 = 28` bits are used.
/// Five 7-bit values are stored in 2 blocks.
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct BitEnc {
    storage: Vec<u32>,
    width: usize,
    mask: u32,
    len: usize,
    usable_bits_per_block: usize,
}

/// Create a mask with `width` 1-bits.
fn mask(width: usize) -> u32 {
    (1 << width) - 1
}

impl BitEnc {
    /// Create a new instance with a given encoding width (e.g. width=2 for using two bits per value).
    /// Supports widths up to 8 bits per character, i.e. `1 <= width <= 8`.
    ///
    /// Complexity: O(1)
    ///
    /// # Example
    ///
    /// ```
    /// use fqtk_lib::bitenc::BitEnc;
    /// let bitenc = BitEnc::new(3);
    /// ```
    pub fn new(width: usize) -> Self {
        assert!(width <= 8, "Only encoding widths up to 8 supported");
        BitEnc {
            storage: Vec::new(),
            width,
            mask: mask(width),
            len: 0,
            usable_bits_per_block: 32 - 32 % width,
        }
    }

    /// Create a new instance with a given capacity and encoding width
    /// (e.g. width=2 for using two bits per value).
    /// Supports widths up to 8 bits per character, i.e. `1 <= width <= 8`.
    ///
    /// Complexity: O(1)
    ///
    /// # Example
    ///
    /// ```
    /// use fqtk_lib::bitenc::BitEnc;
    ///
    /// let bitenc = BitEnc::with_capacity(3, 42);
    /// ```
    pub fn with_capacity(width: usize, n: usize) -> Self {
        assert!(width <= 8, "Only encoding widths up to 8 supported");
        BitEnc {
            storage: Vec::with_capacity(n * width / 32),
            width,
            mask: mask(width),
            len: 0,
            usable_bits_per_block: 32 - 32 % width,
        }
    }

    /// Append a character to the current bit-encoding.
    ///
    /// Complexity: O(1)
    ///
    /// # Example
    ///
    /// ```
    /// use fqtk_lib::bitenc::BitEnc;
    ///
    /// let mut bitenc = BitEnc::new(4);
    /// bitenc.push(0b0000);
    /// bitenc.push(0b1000);
    /// bitenc.push(0b1010);
    /// // The three characters added above are encoded into one u32 entry.
    /// let values: Vec<u8> = bitenc.iter().collect();
    /// assert_eq!(values, [0b0000, 0b1000, 0b1010]);
    /// ```
    pub fn push(&mut self, value: u8) {
        let (block, bit) = self.addr(self.len);
        if bit == 0 {
            self.storage.push(0);
        }
        self.set_by_addr(block, bit, value);
        self.len += 1;
    }

    /// Append the given `value` to the encoding `n` times.
    ///
    /// The added values comprise 0 to 1 blocks that need to be filled up
    /// from previous steps, 0 to m blocks that are
    /// completely filled with the value and 0 to 1 blocks
    /// that are only partially filled.
    ///
    /// Complexity: O(n)
    ///
    /// # Example
    ///
    /// ```
    /// use fqtk_lib::bitenc::BitEnc;
    ///
    /// let mut bitenc = BitEnc::new(8);
    /// // Width: 8 → 4 values per block
    /// // | __ __ __ __ | Denotes one block with 4 empty slots
    ///
    /// bitenc.push_values(5, 0b101010);
    /// // This adds one full and one partial block.
    /// // | 42 42 42 42 | __ __ __ 42 |
    ///
    /// let values: Vec<u8> = bitenc.iter().collect();
    /// assert_eq!(values, [42, 42, 42, 42, 42]);
    ///
    /// bitenc.push_values(1, 23);
    /// // This only fills up an existing block;
    /// // | 42 42 42 42 | __ __ 23 42 |
    ///
    /// let values: Vec<u8> = bitenc.iter().collect();
    /// assert_eq!(values, [42, 42, 42, 42, 42, 23]);
    ///
    /// bitenc.push_values(6, 17);
    /// // Fills up the current block, adds a whole new one but does not create a partial block.
    /// // | 42 42 42 42 | 17 17 23 42 | 17 17 17 17 |
    ///
    /// let values: Vec<u8> = bitenc.iter().collect();
    /// assert_eq!(values, [42, 42, 42, 42, 42, 23, 17, 17, 17, 17, 17, 17]);
    /// ```
    pub fn push_values(&mut self, mut n: usize, value: u8) {
        // Fill up the previous block.
        // Example: After adding 3 values with a width
        // of 8, 8 out out 32 bits in the first block
        // can still be used.
        {
            // Check if there are remaining free slots
            // in the current block, i.e. if the bit offset
            // within the block is non-zero
            let (block, bit) = self.addr(self.len);
            if bit > 0 {
                // insert as many values as required to fill
                // up the previous block and decrease the number
                // left to insert accordingly
                // The take(n) assures this iterator stops before
                // an overflow of n can occur, if less symbols are
                // added than open slots in the block.
                for bit in (bit..32).step_by(self.width).take(n) {
                    self.set_by_addr(block, bit, value);
                    n -= 1;
                    self.len += 1;
                }
            }
        }

        // If symbols remain to be inserted after
        // filling up the current block divide them into
        // full and partial blocks to add them.
        if n > 0 {
            // Create a full value block containing
            // as many copies of value as possible.
            let mut value_block = 0;
            {
                let mut v = u32::from(value);
                for _ in 0..32 / self.width {
                    value_block |= v;
                    v <<= self.width;
                }
            }

            // Append as many full value blocks as needed
            // to reach the last block
            let i = self.len + n;
            let (block, bit) = self.addr(i);
            self.storage.resize(block, value_block);

            if bit > 0 {
                // add the remaining values to a final block
                // let shifted_block = value_block >> (32 - bit);
                let shifted_block = value_block >> (self.usable_bits_per_block - bit);
                self.storage.push(shifted_block);
            }

            self.len = i;
        }
    }

    /// Replace the current value as position `i` with the given value.
    ///
    /// Complexity: O(1)
    ///
    /// ```
    /// use fqtk_lib::bitenc::BitEnc;
    ///
    /// let mut bitenc = BitEnc::new(4);
    /// bitenc.push_values(4, 0b1111);
    /// bitenc.set(2, 0b0000);
    ///
    /// let values: Vec<u8> = bitenc.iter().collect();
    /// assert_eq!(values, [0b1111, 0b1111, 0b0000, 0b1111]);
    /// ```
    pub fn set(&mut self, i: usize, value: u8) {
        let (block, bit) = self.addr(i);
        self.set_by_addr(block, bit, value);
    }

    /// Get the value at position `i`.
    ///
    /// Complexity: O(1)
    ///
    /// ```
    /// use fqtk_lib::bitenc::BitEnc;
    ///
    /// let mut bitenc = BitEnc::new(4);
    /// for value in 1..=4 {
    ///     bitenc.push(value);
    /// }
    ///
    /// let values: Vec<u8> = bitenc.iter().collect();
    /// assert_eq!(values, [0b0001, 0b0010, 0b0011, 0b0100]);
    /// ```
    pub fn get(&self, i: usize) -> Option<u8> {
        if i >= self.len {
            None
        } else {
            let (block, bit) = self.addr(i);
            Some(self.get_by_addr(block, bit))
        }
    }

    /// Iterate over stored values (values will be unpacked into bytes).
    ///
    /// Complexity: O(n), where n is the number of encoded values
    ///
    /// # Example
    ///
    /// ```
    /// use fqtk_lib::bitenc::BitEnc;
    ///
    /// // Fill bitenc with 1, 2, 3, and 4.
    /// let mut bitenc = BitEnc::new(4);
    /// for value in 1..=4 {
    ///     bitenc.push(value);
    /// }
    ///
    /// // Collect iterator for comparison
    /// let values: Vec<u8> = bitenc.iter().collect();
    /// assert_eq!(values, [0b0001, 0b0010, 0b0011, 0b0100]);
    /// ```
    pub fn iter(&self) -> BitEncIter<'_> {
        BitEncIter { bitenc: self, i: 0 }
    }

    /// Clear the sequence.
    ///
    /// Complexity: O(1)
    ///
    /// # Example
    ///
    /// ```
    /// use fqtk_lib::bitenc::BitEnc;
    ///
    /// let mut bitenc = BitEnc::new(2);
    /// bitenc.push(2);
    /// assert_eq!(bitenc.len(), 1);
    /// bitenc.clear();
    /// assert_eq!(bitenc.len(), 0);
    /// ```
    pub fn clear(&mut self) {
        self.storage.clear();
        self.len = 0;
    }

    /// Get the value stored in the given `block` at `bit`.
    fn get_by_addr(&self, block: usize, bit: usize) -> u8 {
        ((self.storage[block] >> bit) & self.mask) as u8
    }

    /// Replace the value in the given `block` at `bit` with the given `value`.
    fn set_by_addr(&mut self, block: usize, bit: usize, value: u8) {
        let mask = self.mask << bit;
        self.storage[block] |= mask;
        self.storage[block] ^= mask;
        self.storage[block] |= (u32::from(value) & self.mask) << bit;
    }

    /// Get the block and start bit for the `i`th encoded value.
    fn addr(&self, i: usize) -> (usize, usize) {
        let k = i * self.width;
        (k / self.usable_bits_per_block, k % self.usable_bits_per_block)
    }

    /// Get the number of symbols encoded.
    ///
    /// Complexity: O(1)
    ///
    /// # Example
    ///
    /// ```
    /// use fqtk_lib::bitenc::BitEnc;
    ///
    /// let mut bitenc = BitEnc::new(8);
    /// bitenc.push(2);
    /// assert_eq!(bitenc.len(), 1);
    /// bitenc.push(2);
    /// bitenc.push(2);
    /// bitenc.push(2);
    /// assert_eq!(bitenc.len(), 4);
    /// // Add another 2 to create a second block
    /// bitenc.push(2);
    /// assert_eq!(bitenc.len(), 5);
    /// ```
    #[deprecated(
        since = "0.33.0",
        note = "Please use the more specific `nr_blocks` and `nr_symbols` functions instead."
    )]
    pub fn len(&self) -> usize {
        self.len
    }

    /// Get the number of blocks used by the encoding.
    ///
    /// Complexity: O(1)
    ///
    /// # Example
    ///
    /// ```
    /// use fqtk_lib::bitenc::BitEnc;
    ///
    /// let mut bitenc = BitEnc::new(8);
    /// bitenc.push(2);
    /// assert_eq!(bitenc.nr_blocks(), 1);
    /// // Add enough 2s to completely fill the first block
    /// bitenc.push(2);
    /// bitenc.push(2);
    /// bitenc.push(2);
    /// assert_eq!(bitenc.nr_blocks(), 1);
    /// // Add another 2 to create a second block
    /// bitenc.push(2);
    /// assert_eq!(bitenc.nr_blocks(), 2);
    /// ```
    pub fn nr_blocks(&self) -> usize {
        self.storage.len()
    }

    /// Get the number of symbols encoded.
    ///
    /// Complexity: O(1)
    ///
    /// # Example
    ///
    /// ```
    /// use fqtk_lib::bitenc::BitEnc;
    ///
    /// let mut bitenc = BitEnc::new(8);
    /// bitenc.push(2);
    /// assert_eq!(bitenc.nr_symbols(), 1);
    /// bitenc.push(2);
    /// bitenc.push(2);
    /// bitenc.push(2);
    /// assert_eq!(bitenc.nr_symbols(), 4);
    /// bitenc.push(2);
    /// assert_eq!(bitenc.nr_symbols(), 5);
    /// ```
    pub fn nr_symbols(&self) -> usize {
        self.len
    }

    /// Is the encoded sequence empty?
    ///
    /// Complexity: O(1)
    ///
    /// # Example
    ///
    /// ```
    /// use fqtk_lib::bitenc::BitEnc;
    ///
    /// let mut bitenc = BitEnc::new(2);
    /// assert!(bitenc.is_empty());
    /// bitenc.push(2);
    /// assert!(!bitenc.is_empty());
    /// bitenc.clear();
    /// assert!(bitenc.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// Calculate the Hamming distance between this and another bitencoded sequence.
    ///
    /// Note: this allows IUPAC fuzzy matching with IUPAC bases in this sequence.
    /// An IUPAC base in the other sequence matches if it is at least as specific as the
    /// corresponding (IUPAC) base in this sequence. E.g. If the other sequence is an
    /// N, it will not match anything but an N, and if the other base is an R, it
    /// will match R, V, D, and N, since the latter IUPAC codes allow both A and G.
    pub fn hamming(&self, other: &BitEnc, max_mismatches: u32) -> u32 {
        assert!(self.len == other.len, "Both bitenc sequences must have the same length");
        assert!(self.width == other.width, "Both bitenc sequences must have the same width");
        let mut count: u32 = 0;
        let values_per_block = self.usable_bits_per_block / self.width;
        for block_index in 0..self.nr_blocks() {
            let intersection = self.storage[block_index] & other.storage[block_index];
            if intersection != self.storage[block_index] {
                let mut shift_i = 0;
                for _ in 0..values_per_block {
                    let intersection_sub = (intersection >> shift_i) & self.mask;
                    let self_sub = (self.storage[block_index] >> shift_i) & self.mask;
                    if intersection_sub != self_sub {
                        count += 1;
                    }
                    shift_i += self.width;
                }
                if count >= max_mismatches {
                    return max_mismatches;
                }
            }
        }
        count
    }
}

/// Iterator over values of a bitencoded sequence (values will be unpacked into bytes).
/// Used to implement the `iter` method of `BitEnc`.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize)]
pub struct BitEncIter<'a> {
    bitenc: &'a BitEnc,
    i: usize,
}

impl Iterator for BitEncIter<'_> {
    type Item = u8;

    fn next(&mut self) -> Option<u8> {
        let value = self.bitenc.get(self.i);
        self.i += 1;
        value
    }
}

#[cfg(test)]
mod tests {
    use super::BitEnc;

    #[test]
    fn test_bitenc() {
        let mut bitenc = BitEnc::new(2);
        bitenc.push(0);
        bitenc.push(2);
        bitenc.push(1);
        let mut values: Vec<u8> = bitenc.iter().collect();
        assert_eq!(values, [0, 2, 1]);
        bitenc.set(1, 3);
        values = bitenc.iter().collect();
        assert_eq!(values, [0, 3, 1]);
    }

    #[test]
    fn test_push_values() {
        let mut bitenc = BitEnc::new(2);
        bitenc.push_values(32, 0);
        assert_eq!(bitenc.storage, [0, 0]);
    }

    #[test]
    fn test_push_values_edge_cases() {
        //! This is a slight variation of the methods doc test,
        //! which also creates a new partial block in the last step
        //! and an additonal full block of 17s.
        //! As this bloats the result vectors, this was not included
        //! in the doc test.
        //! Additionally, this test uses a width of 7 bits, which leave
        //! some bits unused.

        let mut bitenc = BitEnc::new(7);
        // Width: 7 → 4 values per block
        // | __ __ __ __ | Denotes one block with 4 empty slots and 4 leftover bits

        bitenc.push_values(5, 0b101010);
        // This adds one full and one partial block.
        // | 42 42 42 42 | __ __ __ 42 |

        let values: Vec<u8> = bitenc.iter().collect();
        assert_eq!(values, [42, 42, 42, 42, 42]);
        assert_eq!(bitenc.nr_blocks(), 2);
        assert_eq!(bitenc.nr_symbols(), 5);

        bitenc.push_values(1, 23);
        // This only fills up an existing block;
        // | 42 42 42 42 | __ __ 23 42 |
        let values: Vec<u8> = bitenc.iter().collect();
        assert_eq!(values, [42, 42, 42, 42, 42, 23]);
        assert_eq!(bitenc.nr_blocks(), 2);
        assert_eq!(bitenc.nr_symbols(), 6);

        bitenc.push_values(12, 17);
        // Fills up the current block, adds a whole new one AND create a partial block.
        // | 42 42 42 42 | 17 17 23 42 | 17 17 17 17 | 17 17 17 17 | __ __ 17 17 |

        let values: Vec<u8> = bitenc.iter().collect();
        assert_eq!(values, [
            42, 42, 42, 42, 42, 23, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17
        ]);
        assert_eq!(bitenc.nr_blocks(), 5);
        assert_eq!(bitenc.nr_symbols(), 18);
    }

    #[test]
    fn test_issue29() {
        for w in 2..9 {
            let mut vec = BitEnc::with_capacity(w, 1000);
            for _ in 0..1000 {
                vec.push(1);
            }
        }
    }

    #[test]
    fn test_hamming() {
        let mut left = BitEnc::new(4);
        left.push_values(10, 0b1010); // Y (C or T)
        let mut right = BitEnc::new(4);
        right.push_values(10, 0b0010); // C
        let mut mismatches = BitEnc::new(4);
        mismatches.push_values(2, 0b1010); // Y
        mismatches.push_values(1, 0b0100); // G
        mismatches.push_values(3, 0b1010); // Y
        mismatches.push_values(2, 0b0001); // A
        mismatches.push_values(2, 0b1010); // Y

        assert_eq!(left.hamming(&left, 100), 0);
        assert_eq!(right.hamming(&right, 100), 0);
        assert_eq!(right.hamming(&left, 100), 0); // since left is less restrictive than right
        assert_eq!(left.hamming(&right, 100), 10); // since right is more restrictive than left
        assert_eq!(left.hamming(&mismatches, 100), 3);
        assert_eq!(mismatches.hamming(&left, 100), 3);
        assert_eq!(mismatches.hamming(&left, 1), 1);
        assert_eq!(mismatches.hamming(&left, 2), 2);
        assert_eq!(mismatches.hamming(&left, 3), 3);
    }
}
