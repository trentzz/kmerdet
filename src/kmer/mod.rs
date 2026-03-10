pub mod canonical;
pub mod encoding;

/// A k-mer stored as a 2-bit packed u64 (supports up to k=32).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct Kmer {
    /// 2-bit encoded bases, right-aligned.
    pub data: u64,
    /// Length of the k-mer.
    pub k: u8,
}

impl Kmer {
    /// Create a new Kmer from a 2-bit packed value and length.
    pub fn new(data: u64, k: u8) -> Self {
        debug_assert!(k <= 32);
        let mask = if k == 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
        Self { data: data & mask, k }
    }

    /// Encode a DNA string (ACGT) into a Kmer.
    pub fn from_str(s: &str) -> Option<Self> {
        encoding::encode(s)
    }

    /// Decode this Kmer back to a DNA string.
    pub fn to_string(&self) -> String {
        encoding::decode(self.data, self.k)
    }

    /// Return the canonical form (lexicographically smaller of forward and reverse complement).
    pub fn canonical(&self) -> Self {
        canonical::canonical(*self)
    }

    /// Return the reverse complement.
    pub fn reverse_complement(&self) -> Self {
        canonical::reverse_complement(*self)
    }

    /// Extend this k-mer to the right with one base, keeping length k.
    pub fn extend_right(&self, base: u8) -> Self {
        let shifted = (self.data << 2) | (base as u64 & 0b11);
        Self::new(shifted, self.k)
    }

    /// Extend this k-mer to the left with one base, keeping length k.
    pub fn extend_left(&self, base: u8) -> Self {
        let shifted = self.data >> 2;
        let placed = shifted | ((base as u64 & 0b11) << (2 * (self.k - 1)));
        Self::new(placed, self.k)
    }
}

impl std::fmt::Display for Kmer {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.to_string())
    }
}
