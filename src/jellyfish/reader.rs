//! Pure Rust jellyfish database reader using the `jellyfish-reader` crate.
//!
//! Replaces the C++ FFI wrapper with memory-mapped random-access lookups.

use std::path::Path;

use anyhow::{Context, Result};
use jellyfish_reader::QueryMerFile;

use super::KmerDatabase;

/// A jellyfish database backed by the pure Rust `jellyfish-reader` crate.
///
/// Uses memory-mapped I/O for efficient random-access k-mer lookups.
pub struct JellyfishReader {
    db: QueryMerFile,
}

impl JellyfishReader {
    /// Open a jellyfish `.jf` binary database file.
    pub fn open(path: &Path) -> Result<Self> {
        let db = QueryMerFile::open(path)
            .with_context(|| format!("opening jellyfish database: {}", path.display()))?;
        Ok(Self { db })
    }
}

impl KmerDatabase for JellyfishReader {
    fn query(&self, kmer: &str) -> u64 {
        self.db.query(kmer).unwrap_or(None).unwrap_or(0)
    }

    fn kmer_length(&self) -> u8 {
        self.db.k() as u8
    }
}
