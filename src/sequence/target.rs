/// Target FASTA loading and reference sequence decomposition.
///
/// A target is a short FASTA sequence (~100-300 bp) around a region of interest.
/// It gets decomposed into overlapping k-mers for the walking algorithm.

use std::path::Path;

use anyhow::{Context, Result};

/// A loaded target sequence with its k-mer decomposition.
#[derive(Debug, Clone)]
pub struct Target {
    /// Name from the FASTA header (e.g., "NPM1_4ins_exons_10-11utr").
    pub name: String,
    /// Full DNA sequence.
    pub sequence: String,
    /// Source file path.
    pub source: std::path::PathBuf,
}

/// Reference sequence decomposed into k-mers.
#[derive(Debug)]
pub struct RefSeq {
    /// The target this was derived from.
    pub target: Target,
    /// Ordered list of reference k-mer strings.
    pub kmers: Vec<String>,
}

impl RefSeq {
    /// Decompose a target sequence into overlapping k-mers of length k.
    pub fn from_target(target: Target, k: u8) -> Result<Self> {
        let k = k as usize;
        anyhow::ensure!(
            target.sequence.len() >= k,
            "target '{}' is shorter ({} bp) than k-mer length ({})",
            target.name,
            target.sequence.len(),
            k,
        );

        let kmers: Vec<String> = target
            .sequence
            .as_bytes()
            .windows(k)
            .map(|w| std::str::from_utf8(w).unwrap().to_uppercase())
            .collect();

        Ok(Self { target, kmers })
    }
}

/// Load all target FASTA files from a path (file or directory, recursive).
pub fn load_targets(path: &Path) -> Result<Vec<Target>> {
    let mut targets = Vec::new();

    if path.is_dir() {
        for entry in walkdir(path)? {
            if is_fasta(&entry) {
                targets.extend(load_fasta(&entry)?);
            }
        }
    } else {
        targets.extend(
            load_fasta(path).with_context(|| format!("loading target: {}", path.display()))?,
        );
    }

    Ok(targets)
}

fn is_fasta(path: &Path) -> bool {
    matches!(
        path.extension().and_then(|e| e.to_str()),
        Some("fa" | "fasta" | "fna")
    )
}

fn walkdir(dir: &Path) -> Result<Vec<std::path::PathBuf>> {
    let mut files = Vec::new();
    for entry in std::fs::read_dir(dir)
        .with_context(|| format!("reading directory: {}", dir.display()))?
    {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir() {
            files.extend(walkdir(&path)?);
        } else {
            files.push(path);
        }
    }
    files.sort();
    Ok(files)
}

fn load_fasta(path: &Path) -> Result<Vec<Target>> {
    use needletail::parse_fastx_file;

    let mut reader =
        parse_fastx_file(path).with_context(|| format!("parsing FASTA: {}", path.display()))?;

    let mut targets = Vec::new();
    while let Some(record) = reader.next() {
        let record = record.with_context(|| format!("reading record from {}", path.display()))?;
        let name = String::from_utf8_lossy(record.id()).to_string();
        let sequence = String::from_utf8(record.seq().to_vec())
            .with_context(|| format!("invalid UTF-8 in sequence '{}'", name))?
            .to_uppercase();

        targets.push(Target {
            name,
            sequence,
            source: path.to_owned(),
        });
    }

    Ok(targets)
}
