/// Canonical k-mer form: min(kmer, reverse_complement(kmer)).
///
/// Jellyfish stores canonical k-mers to count both strands together.

use super::Kmer;

/// Compute the reverse complement of a 2-bit packed k-mer.
pub fn reverse_complement(kmer: Kmer) -> Kmer {
    let mut rc: u64 = 0;
    let mut data = kmer.data;
    for _ in 0..kmer.k {
        let base = data & 0b11;
        // Complement: A(00)↔T(11), C(01)↔G(10) — just flip both bits
        let comp = base ^ 0b11;
        rc = (rc << 2) | comp;
        data >>= 2;
    }
    Kmer::new(rc, kmer.k)
}

/// Return the canonical form (lexicographically smaller of fwd and revcomp).
pub fn canonical(kmer: Kmer) -> Kmer {
    let rc = reverse_complement(kmer);
    if kmer.data <= rc.data {
        kmer
    } else {
        rc
    }
}
