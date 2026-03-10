pub mod ffi;
pub mod header;

/// Trait abstracting k-mer database lookups.
///
/// Enables both the real jellyfish FFI backend and a mock for testing.
pub trait KmerDatabase: Send + Sync {
    /// Query the count of a k-mer (as a DNA string).
    fn query(&self, kmer: &str) -> u64;

    /// Get the k-mer length used by this database.
    fn kmer_length(&self) -> u8;

    /// Query the count of a k-mer and its canonical form.
    /// Default implementation canonicalizes then queries.
    fn query_canonical(&self, kmer: &str) -> u64 {
        self.query(kmer)
    }
}
