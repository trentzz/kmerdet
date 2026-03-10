use std::collections::HashMap;

use kmerdet::jellyfish::KmerDatabase;

/// Mock k-mer database for testing without a real jellyfish file.
pub struct MockDb {
    counts: HashMap<String, u64>,
    k: u8,
}

impl MockDb {
    pub fn new(k: u8) -> Self {
        Self {
            counts: HashMap::new(),
            k,
        }
    }

    /// Set the count for a specific k-mer.
    pub fn set(&mut self, kmer: &str, count: u64) {
        self.counts.insert(kmer.to_uppercase(), count);
    }

    /// Set counts for all k-mers in a sequence (sliding window of size k).
    pub fn set_sequence(&mut self, sequence: &str, count: u64) {
        let seq = sequence.to_uppercase();
        let k = self.k as usize;
        for window in seq.as_bytes().windows(k) {
            let kmer = std::str::from_utf8(window).unwrap().to_string();
            self.counts.insert(kmer, count);
        }
    }
}

impl KmerDatabase for MockDb {
    fn query(&self, kmer: &str) -> u64 {
        *self.counts.get(&kmer.to_uppercase()).unwrap_or(&0)
    }

    fn kmer_length(&self) -> u8 {
        self.k
    }
}
