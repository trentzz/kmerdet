pub mod classifier;
pub mod cluster;
pub mod normalize;
pub mod quantifier;

/// The type of variant detected.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub enum VariantType {
    /// No change from reference.
    Reference,
    /// Single or multi-nucleotide substitution.
    Substitution,
    /// Insertion of bases.
    Insertion,
    /// Deletion of bases.
    Deletion,
    /// Internal tandem duplication.
    Itd,
    /// Complex insertion+deletion.
    Indel,
}

impl std::fmt::Display for VariantType {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Reference => write!(f, "Reference"),
            Self::Substitution => write!(f, "Substitution"),
            Self::Insertion => write!(f, "Insertion"),
            Self::Deletion => write!(f, "Deletion"),
            Self::Itd => write!(f, "ITD"),
            Self::Indel => write!(f, "Indel"),
        }
    }
}

impl std::str::FromStr for VariantType {
    type Err = anyhow::Error;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "reference" => Ok(Self::Reference),
            "substitution" => Ok(Self::Substitution),
            "insertion" => Ok(Self::Insertion),
            "deletion" => Ok(Self::Deletion),
            "itd" => Ok(Self::Itd),
            "indel" => Ok(Self::Indel),
            _ => anyhow::bail!("unknown variant type: '{}'", s),
        }
    }
}

/// A complete variant call with all associated data.
#[derive(Debug, Clone, serde::Serialize)]
pub struct VariantCall {
    /// Sample/database identifier.
    pub sample: String,
    /// Target name from FASTA header.
    pub target: String,
    /// Variant type.
    pub variant_type: VariantType,
    /// Human-readable variant description (e.g., "41:acgt/ACGTACGT:45").
    pub variant_name: String,
    /// Relative variant allele frequency (0.0–1.0).
    pub rvaf: f64,
    /// Total k-mer expression.
    pub expression: f64,
    /// Minimum k-mer count along the path.
    pub min_coverage: u64,
    /// Count of the first k-mer in the path.
    pub start_kmer_count: u64,
    /// Reference path DNA sequence.
    pub ref_sequence: String,
    /// Alternative path DNA sequence.
    pub alt_sequence: String,
    /// Status info (e.g., "vs_ref").
    pub info: String,
    /// Chromosome (extracted from target name if available).
    pub chrom: Option<String>,
    /// Genomic position.
    pub pos: Option<u64>,
    /// Reference allele.
    pub ref_allele: Option<String>,
    /// Alternative allele.
    pub alt_allele: Option<String>,
}
