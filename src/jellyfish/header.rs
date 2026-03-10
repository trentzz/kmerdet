/// Pure Rust parser for the jellyfish .jf file header.
///
/// The header is a JSON object at the start of the file, terminated by a null byte.
/// It contains metadata including k-mer length, counter length, hash size, etc.

use std::path::Path;

use anyhow::{Context, Result};
use serde::Deserialize;

#[derive(Debug, Deserialize)]
pub struct JfHeader {
    #[serde(rename = "cmdline")]
    pub command_line: Option<Vec<String>>,
    #[serde(rename = "canonical")]
    pub canonical: Option<bool>,
    #[serde(rename = "key_len")]
    pub key_len: Option<u64>,
    #[serde(rename = "val_len")]
    pub val_len: Option<u64>,
    #[serde(rename = "counter_len")]
    pub counter_len: Option<u64>,
    #[serde(rename = "size")]
    pub size: Option<u64>,
}

impl JfHeader {
    /// Read and parse the JSON header from a .jf file.
    pub fn from_file(path: &Path) -> Result<Self> {
        let data = std::fs::read(path)
            .with_context(|| format!("reading jellyfish file: {}", path.display()))?;

        // Find the null terminator that ends the JSON header
        let null_pos = data
            .iter()
            .position(|&b| b == 0)
            .context("no null terminator found in .jf header")?;

        let json_str = std::str::from_utf8(&data[..null_pos])
            .context("invalid UTF-8 in .jf header")?;

        let header: JfHeader =
            serde_json::from_str(json_str).context("parsing .jf JSON header")?;

        Ok(header)
    }

    /// Derive the k-mer length from key_len (which is in bits: k * 2).
    pub fn kmer_length(&self) -> Option<u8> {
        self.key_len.map(|bits| (bits / 2) as u8)
    }
}
