/// Reference and alternative path representations.
///
/// A path is a sequence of k-mers that forms a contiguous DNA sequence
/// through the k-mer graph.

/// A path through the k-mer graph, represented as an ordered sequence of k-mers.
#[derive(Debug, Clone)]
pub struct KmerPath {
    /// Ordered k-mer strings forming this path.
    pub kmers: Vec<String>,
    /// Whether this is the reference path.
    pub is_reference: bool,
}

impl KmerPath {
    /// Reconstruct the full DNA sequence from the ordered k-mers.
    ///
    /// Each successive k-mer overlaps by (k-1) bases, so we take the first
    /// k-mer in full and append the last base of each subsequent k-mer.
    pub fn to_sequence(&self) -> String {
        if self.kmers.is_empty() {
            return String::new();
        }

        let mut seq = self.kmers[0].clone();
        for kmer in &self.kmers[1..] {
            if let Some(last_base) = kmer.bytes().last() {
                seq.push(last_base as char);
            }
        }
        seq
    }

    /// Number of k-mers in this path.
    pub fn len(&self) -> usize {
        self.kmers.len()
    }

    /// Whether this path is empty.
    pub fn is_empty(&self) -> bool {
        self.kmers.is_empty()
    }
}
