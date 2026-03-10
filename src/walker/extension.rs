/// Child k-mer extension with count thresholding.
///
/// Given a k-mer, try all 4 possible single-base extensions and keep
/// those whose count passes both the absolute (n_cutoff) and ratio (cutoff)
/// thresholds — matching km's `Jellyfish.get_child` logic.

use crate::jellyfish::KmerDatabase;

/// A candidate child k-mer with its count.
#[derive(Debug, Clone)]
pub struct ChildKmer {
    pub sequence: String,
    pub count: u64,
}

const BASES: [u8; 4] = [b'A', b'C', b'G', b'T'];

/// Extend a k-mer forward (to the right) by one base.
///
/// Returns all children whose count meets the threshold criteria:
/// - count >= n_cutoff (absolute minimum)
/// - count >= sum_of_siblings * ratio (relative minimum)
pub fn extend_forward(
    db: &dyn KmerDatabase,
    kmer: &str,
    n_cutoff: u32,
    ratio: f64,
) -> Vec<ChildKmer> {
    let k = kmer.len();
    let suffix = &kmer[1..]; // drop first base, keep last k-1

    let mut candidates = Vec::with_capacity(4);
    let mut total_count: u64 = 0;

    for &base in &BASES {
        let mut child = String::with_capacity(k);
        child.push_str(suffix);
        child.push(base as char);

        let count = db.query(&child);
        total_count += count;
        candidates.push(ChildKmer {
            sequence: child,
            count,
        });
    }

    let ratio_threshold = (total_count as f64 * ratio) as u64;
    let threshold = std::cmp::max(n_cutoff as u64, ratio_threshold);

    candidates
        .into_iter()
        .filter(|c| c.count >= threshold)
        .collect()
}

/// Extend a k-mer backward (to the left) by one base.
pub fn extend_backward(
    db: &dyn KmerDatabase,
    kmer: &str,
    n_cutoff: u32,
    ratio: f64,
) -> Vec<ChildKmer> {
    let k = kmer.len();
    let prefix = &kmer[..k - 1]; // drop last base, keep first k-1

    let mut candidates = Vec::with_capacity(4);
    let mut total_count: u64 = 0;

    for &base in &BASES {
        let mut child = String::with_capacity(k);
        child.push(base as char);
        child.push_str(prefix);

        let count = db.query(&child);
        total_count += count;
        candidates.push(ChildKmer {
            sequence: child,
            count,
        });
    }

    let ratio_threshold = (total_count as f64 * ratio) as u64;
    let threshold = std::cmp::max(n_cutoff as u64, ratio_threshold);

    candidates
        .into_iter()
        .filter(|c| c.count >= threshold)
        .collect()
}
