/// 2-bit encoding/decoding of DNA bases.
///
/// Encoding: A=0b00, C=0b01, G=0b10, T=0b11
/// This matches jellyfish's encoding for bit-exact compatibility.

use super::Kmer;

/// Encode a single base character to its 2-bit representation.
#[inline]
pub fn encode_base(base: u8) -> Option<u8> {
    match base {
        b'A' | b'a' => Some(0b00),
        b'C' | b'c' => Some(0b01),
        b'G' | b'g' => Some(0b10),
        b'T' | b't' => Some(0b11),
        _ => None,
    }
}

/// Decode a 2-bit value back to a base character.
#[inline]
pub fn decode_base(bits: u8) -> u8 {
    match bits & 0b11 {
        0b00 => b'A',
        0b01 => b'C',
        0b10 => b'G',
        0b11 => b'T',
        _ => unreachable!(),
    }
}

/// Encode a DNA string into a packed Kmer.
pub fn encode(s: &str) -> Option<Kmer> {
    let bytes = s.as_bytes();
    let k = bytes.len();
    if k == 0 || k > 32 {
        return None;
    }

    let mut data: u64 = 0;
    for &b in bytes {
        let bits = encode_base(b)?;
        data = (data << 2) | bits as u64;
    }

    Some(Kmer::new(data, k as u8))
}

/// Decode a 2-bit packed value into a DNA string of length k.
pub fn decode(data: u64, k: u8) -> String {
    let mut s = Vec::with_capacity(k as usize);
    for i in (0..k).rev() {
        let bits = ((data >> (2 * i)) & 0b11) as u8;
        s.push(decode_base(bits));
    }
    String::from_utf8(s).expect("valid ASCII bases")
}
