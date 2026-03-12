//! Panel-of-Normals (PoN) for filtering systematic false positives.
//!
//! A Panel-of-Normals records variants found in normal (non-tumor) samples.
//! Variants appearing in many normals are likely systematic errors (sequencing
//! artifacts, reference errors, CHIP variants) rather than true somatic
//! mutations. The PoN acts as a blacklist filter.
//!
//! ## Workflow
//!
//! 1. **Build**: Run `kmerdet pon build --normals n1.tsv n2.tsv ... --output pon.json`
//!    to construct a PoN from detection results on normal samples.
//! 2. **Filter**: Run `kmerdet pon filter --input results.tsv --pon pon.json`
//!    to remove variants present in the PoN above a frequency threshold.

use std::collections::HashMap;
use std::path::Path;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

use crate::variant::normalize::{normalize_indel, NormalizedVariant};
use crate::variant::{VariantCall, VariantType};

/// Entry in the panel-of-normals.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PonEntry {
    /// Number of normal samples where this variant was observed.
    pub sample_count: u32,
    /// Total number of normal samples in the panel.
    pub total_samples: u32,
    /// Frequency = sample_count / total_samples.
    pub frequency: f64,
    /// Maximum VAF observed across normal samples.
    pub max_vaf: f64,
    /// Mean VAF across normal samples where detected.
    pub mean_vaf: f64,
}

/// Key for PoN lookup: normalized variant identity.
#[derive(Debug, Clone, Hash, PartialEq, Eq, Serialize, Deserialize)]
pub struct PonKey {
    pub target: String,
    pub chrom: String,
    pub pos: u64,
    pub ref_allele: String,
    pub alt_allele: String,
}

/// Intermediate accumulator used during PoN building.
#[derive(Debug)]
struct PonAccumulator {
    /// Set of sample names that have seen this variant.
    samples: std::collections::HashSet<String>,
    /// VAFs observed across samples.
    vafs: Vec<f64>,
}

/// A key-value pair for JSON serialization of the PoN entries.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct PonRecord {
    key: PonKey,
    entry: PonEntry,
}

/// Serialization wrapper: convert HashMap to Vec for JSON compatibility.
/// JSON requires string keys; PonKey is a struct, so we serialize as an array
/// of {key, entry} records.
mod entries_serde {
    use super::*;

    pub fn serialize<S>(
        map: &HashMap<PonKey, PonEntry>,
        serializer: S,
    ) -> std::result::Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        let records: Vec<PonRecord> = map
            .iter()
            .map(|(k, v)| PonRecord {
                key: k.clone(),
                entry: v.clone(),
            })
            .collect();
        records.serialize(serializer)
    }

    pub fn deserialize<'de, D>(
        deserializer: D,
    ) -> std::result::Result<HashMap<PonKey, PonEntry>, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        let records: Vec<PonRecord> = Vec::deserialize(deserializer)?;
        Ok(records.into_iter().map(|r| (r.key, r.entry)).collect())
    }
}

/// The panel-of-normals database.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PanelOfNormals {
    /// Map from normalized variant key to PoN entry statistics.
    #[serde(with = "entries_serde")]
    pub entries: HashMap<PonKey, PonEntry>,
    /// Total number of normal samples used to build this panel.
    pub total_samples: u32,
    /// ISO 8601 date string when this PoN was created.
    pub creation_date: String,
}

impl PanelOfNormals {
    /// Build a PoN from multiple detection result files.
    ///
    /// Each file is a TSV produced by `kmerdet detect`. Non-Reference calls
    /// with rVAF >= `min_vaf` are normalized and accumulated. The frequency
    /// of each variant across all normal samples is computed.
    pub fn build(result_files: &[&Path], min_vaf: f64) -> Result<Self> {
        let mut accumulators: HashMap<PonKey, PonAccumulator> = HashMap::new();
        let mut all_samples: std::collections::HashSet<String> = std::collections::HashSet::new();

        for path in result_files {
            let calls = crate::output::parse_detection_tsv(path)
                .with_context(|| format!("reading normal sample results: {}", path.display()))?;

            for call in &calls {
                // Skip Reference calls -- they are not variants.
                if call.variant_type == VariantType::Reference {
                    continue;
                }

                // Skip calls below the minimum VAF threshold.
                if call.rvaf < min_vaf {
                    continue;
                }

                let key = make_pon_key(call);
                all_samples.insert(call.sample.clone());

                let acc = accumulators.entry(key).or_insert_with(|| PonAccumulator {
                    samples: std::collections::HashSet::new(),
                    vafs: Vec::new(),
                });

                acc.samples.insert(call.sample.clone());
                acc.vafs.push(call.rvaf);
            }
        }

        let total_samples = all_samples.len().max(result_files.len()) as u32;

        let entries = accumulators
            .into_iter()
            .map(|(key, acc)| {
                let sample_count = acc.samples.len() as u32;
                let frequency = if total_samples > 0 {
                    sample_count as f64 / total_samples as f64
                } else {
                    0.0
                };
                let max_vaf = acc
                    .vafs
                    .iter()
                    .copied()
                    .fold(0.0_f64, f64::max);
                let mean_vaf = if acc.vafs.is_empty() {
                    0.0
                } else {
                    acc.vafs.iter().sum::<f64>() / acc.vafs.len() as f64
                };

                let entry = PonEntry {
                    sample_count,
                    total_samples,
                    frequency,
                    max_vaf,
                    mean_vaf,
                };
                (key, entry)
            })
            .collect();

        // Use a simple date string without requiring chrono.
        let creation_date = current_date_string();

        Ok(Self {
            entries,
            total_samples,
            creation_date,
        })
    }

    /// Save PoN to a JSON file.
    pub fn save(&self, path: &Path) -> Result<()> {
        let json = serde_json::to_string_pretty(self)
            .context("serializing PoN to JSON")?;
        std::fs::write(path, json)
            .with_context(|| format!("writing PoN file: {}", path.display()))?;
        Ok(())
    }

    /// Load PoN from a JSON file.
    pub fn load(path: &Path) -> Result<Self> {
        let json = std::fs::read_to_string(path)
            .with_context(|| format!("reading PoN file: {}", path.display()))?;
        let pon: Self = serde_json::from_str(&json)
            .with_context(|| format!("parsing PoN JSON: {}", path.display()))?;
        Ok(pon)
    }

    /// Check if a variant is in the PoN above a frequency threshold.
    ///
    /// Returns `Some(&PonEntry)` if the variant is found with frequency >=
    /// `frequency_threshold`, or `None` if the variant is absent or below
    /// the threshold.
    pub fn check(&self, call: &VariantCall, frequency_threshold: f64) -> Option<&PonEntry> {
        let key = make_pon_key(call);
        self.entries.get(&key).filter(|entry| entry.frequency >= frequency_threshold)
    }

    /// Filter a list of variant calls against the PoN.
    ///
    /// Returns `(passed, filtered)` where `passed` contains calls that are
    /// NOT in the PoN (or below the frequency threshold) and `filtered`
    /// contains calls that ARE in the PoN above the threshold.
    ///
    /// Reference calls always pass through (they are not variants).
    pub fn filter_calls(
        &self,
        calls: Vec<VariantCall>,
        frequency_threshold: f64,
    ) -> (Vec<VariantCall>, Vec<VariantCall>) {
        let mut passed = Vec::new();
        let mut filtered = Vec::new();

        for call in calls {
            if call.variant_type == VariantType::Reference {
                passed.push(call);
                continue;
            }

            if self.check(&call, frequency_threshold).is_some() {
                filtered.push(call);
            } else {
                passed.push(call);
            }
        }

        (passed, filtered)
    }

    /// Return the number of entries in the panel.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Return whether the panel is empty.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }
}

/// Build a `PonKey` from a `VariantCall` by normalizing the variant.
fn make_pon_key(call: &VariantCall) -> PonKey {
    let chrom = call.chrom.as_deref().unwrap_or("");
    let ref_allele = call.ref_allele.as_deref().unwrap_or("");
    let alt_allele = call.alt_allele.as_deref().unwrap_or("");
    let pos = call.pos.unwrap_or(0);

    // Use reference context for normalization when available.
    let ref_context = if call.ref_sequence.is_empty() {
        None
    } else {
        Some(call.ref_sequence.as_str())
    };

    let normalized: NormalizedVariant =
        normalize_indel(chrom, pos, ref_allele, alt_allele, ref_context);

    PonKey {
        target: call.target.clone(),
        chrom: normalized.chrom,
        pos: normalized.pos,
        ref_allele: normalized.ref_allele,
        alt_allele: normalized.alt_allele,
    }
}

/// Get the current date as an ISO 8601 string (YYYY-MM-DD).
///
/// Uses `std::time::SystemTime` to avoid requiring the `chrono` crate.
fn current_date_string() -> String {
    // Compute days since the Unix epoch.
    let secs = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap_or_default()
        .as_secs();

    let days = secs / 86400;

    // Civil calendar conversion from days since epoch.
    // Algorithm: https://howardhinnant.github.io/date_algorithms.html
    let z = days as i64 + 719468;
    let era = if z >= 0 { z } else { z - 146096 } / 146097;
    let doe = (z - era * 146097) as u64;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let d = doy - (153 * mp + 2) / 5 + 1;
    let m = if mp < 10 { mp + 3 } else { mp - 9 };
    let y = if m <= 2 { y + 1 } else { y };

    format!("{:04}-{:02}-{:02}", y, m, d)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Helper to build a VariantCall for testing.
    fn make_call(
        sample: &str,
        target: &str,
        variant_type: VariantType,
        rvaf: f64,
        chrom: &str,
        pos: u64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> VariantCall {
        VariantCall {
            sample: sample.to_string(),
            target: target.to_string(),
            variant_type,
            variant_name: format!("{}:{}/{}:{}", pos, ref_allele, alt_allele, pos),
            rvaf,
            expression: 10.0,
            min_coverage: 20,
            start_kmer_count: 50,
            ref_sequence: String::new(),
            alt_sequence: String::new(),
            info: "vs_ref".to_string(),
            chrom: if chrom.is_empty() {
                None
            } else {
                Some(chrom.to_string())
            },
            pos: if pos == 0 && chrom.is_empty() {
                None
            } else {
                Some(pos)
            },
            ref_allele: if ref_allele.is_empty() {
                None
            } else {
                Some(ref_allele.to_string())
            },
            alt_allele: if alt_allele.is_empty() {
                None
            } else {
                Some(alt_allele.to_string())
            },
            pvalue: None,
            qual: None,
            ci_lower: None,
            ci_upper: None,
        }
    }

    // ── PonEntry creation and frequency ─────────────────────────────────

    #[test]
    fn test_pon_entry_frequency_computation() {
        let entry = PonEntry {
            sample_count: 3,
            total_samples: 10,
            frequency: 3.0 / 10.0,
            max_vaf: 0.05,
            mean_vaf: 0.03,
        };
        assert!((entry.frequency - 0.3).abs() < 1e-10);
        assert_eq!(entry.sample_count, 3);
        assert_eq!(entry.total_samples, 10);
    }

    #[test]
    fn test_pon_entry_max_and_mean_vaf() {
        let entry = PonEntry {
            sample_count: 2,
            total_samples: 5,
            frequency: 0.4,
            max_vaf: 0.08,
            mean_vaf: 0.05,
        };
        assert!((entry.max_vaf - 0.08).abs() < 1e-10);
        assert!((entry.mean_vaf - 0.05).abs() < 1e-10);
    }

    // ── PonKey normalization ────────────────────────────────────────────

    #[test]
    fn test_make_pon_key_snv() {
        let call = make_call("s1", "target1", VariantType::Substitution, 0.05, "chr1", 100, "A", "T");
        let key = make_pon_key(&call);
        assert_eq!(key.target, "target1");
        assert_eq!(key.chrom, "chr1");
        assert_eq!(key.pos, 100);
        assert_eq!(key.ref_allele, "A");
        assert_eq!(key.alt_allele, "T");
    }

    #[test]
    fn test_make_pon_key_with_ref_context_normalization() {
        // An INDEL that normalizes via trimming. Position must be within the
        // ref_context so left_align can access context bytes.
        let mut call = make_call(
            "s1",
            "target1",
            VariantType::Deletion,
            0.05,
            "chr1",
            4,
            "AGT",
            "AT",
        );
        // Provide a ref context so normalization can left-align.
        // ref_context: "AAAAGTTTT" (9 bases), pos=4 points to 'G' in "GTTTT".
        call.ref_sequence = "AAAAGTTTT".to_string();
        let key = make_pon_key(&call);
        // After suffix trimming: AGT/AT -> AG/A at pos 4.
        assert_eq!(key.ref_allele, "AG");
        assert_eq!(key.alt_allele, "A");
    }

    // ── PoN check ───────────────────────────────────────────────────────

    #[test]
    fn test_check_finds_known_variant() {
        let key = PonKey {
            target: "target1".to_string(),
            chrom: "chr1".to_string(),
            pos: 100,
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
        };
        let entry = PonEntry {
            sample_count: 5,
            total_samples: 10,
            frequency: 0.5,
            max_vaf: 0.05,
            mean_vaf: 0.03,
        };
        let mut entries = HashMap::new();
        entries.insert(key, entry);

        let pon = PanelOfNormals {
            entries,
            total_samples: 10,
            creation_date: "2026-01-01".to_string(),
        };

        let call = make_call("tumor", "target1", VariantType::Substitution, 0.05, "chr1", 100, "A", "T");
        let result = pon.check(&call, 0.05);
        assert!(result.is_some());
        assert_eq!(result.unwrap().sample_count, 5);
    }

    #[test]
    fn test_check_returns_none_for_novel_variant() {
        let key = PonKey {
            target: "target1".to_string(),
            chrom: "chr1".to_string(),
            pos: 100,
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
        };
        let entry = PonEntry {
            sample_count: 5,
            total_samples: 10,
            frequency: 0.5,
            max_vaf: 0.05,
            mean_vaf: 0.03,
        };
        let mut entries = HashMap::new();
        entries.insert(key, entry);

        let pon = PanelOfNormals {
            entries,
            total_samples: 10,
            creation_date: "2026-01-01".to_string(),
        };

        // Different variant (G>C at pos 200).
        let call = make_call("tumor", "target1", VariantType::Substitution, 0.05, "chr1", 200, "G", "C");
        assert!(pon.check(&call, 0.05).is_none());
    }

    #[test]
    fn test_check_returns_none_when_below_threshold() {
        let key = PonKey {
            target: "target1".to_string(),
            chrom: "chr1".to_string(),
            pos: 100,
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
        };
        let entry = PonEntry {
            sample_count: 1,
            total_samples: 100,
            frequency: 0.01,
            max_vaf: 0.02,
            mean_vaf: 0.02,
        };
        let mut entries = HashMap::new();
        entries.insert(key, entry);

        let pon = PanelOfNormals {
            entries,
            total_samples: 100,
            creation_date: "2026-01-01".to_string(),
        };

        // The variant is in PoN but at frequency 0.01, below threshold 0.05.
        let call = make_call("tumor", "target1", VariantType::Substitution, 0.05, "chr1", 100, "A", "T");
        assert!(pon.check(&call, 0.05).is_none());
    }

    // ── filter_calls ────────────────────────────────────────────────────

    #[test]
    fn test_filter_calls_separates_passed_and_filtered() {
        let key = PonKey {
            target: "target1".to_string(),
            chrom: "chr1".to_string(),
            pos: 100,
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
        };
        let entry = PonEntry {
            sample_count: 5,
            total_samples: 10,
            frequency: 0.5,
            max_vaf: 0.05,
            mean_vaf: 0.03,
        };
        let mut entries = HashMap::new();
        entries.insert(key, entry);

        let pon = PanelOfNormals {
            entries,
            total_samples: 10,
            creation_date: "2026-01-01".to_string(),
        };

        let calls = vec![
            // This one is in the PoN.
            make_call("tumor", "target1", VariantType::Substitution, 0.05, "chr1", 100, "A", "T"),
            // This one is novel.
            make_call("tumor", "target1", VariantType::Substitution, 0.03, "chr2", 200, "G", "C"),
            // Reference calls pass through.
            make_call("tumor", "target1", VariantType::Reference, 0.0, "chr1", 100, "", ""),
        ];

        let (passed, filtered) = pon.filter_calls(calls, 0.05);
        assert_eq!(passed.len(), 2); // novel + reference
        assert_eq!(filtered.len(), 1); // PoN match
        assert_eq!(filtered[0].chrom, Some("chr1".to_string()));
        assert_eq!(filtered[0].pos, Some(100));
    }

    #[test]
    fn test_filter_calls_empty_input() {
        let pon = PanelOfNormals {
            entries: HashMap::new(),
            total_samples: 0,
            creation_date: "2026-01-01".to_string(),
        };
        let (passed, filtered) = pon.filter_calls(vec![], 0.05);
        assert!(passed.is_empty());
        assert!(filtered.is_empty());
    }

    #[test]
    fn test_filter_calls_empty_pon_passes_everything() {
        let pon = PanelOfNormals {
            entries: HashMap::new(),
            total_samples: 0,
            creation_date: "2026-01-01".to_string(),
        };
        let calls = vec![
            make_call("tumor", "t1", VariantType::Substitution, 0.05, "chr1", 100, "A", "T"),
            make_call("tumor", "t1", VariantType::Insertion, 0.10, "chr2", 200, "G", "GA"),
        ];
        let (passed, filtered) = pon.filter_calls(calls, 0.05);
        assert_eq!(passed.len(), 2);
        assert!(filtered.is_empty());
    }

    // ── PoN save/load roundtrip ─────────────────────────────────────────

    #[test]
    fn test_pon_save_load_roundtrip() {
        let key = PonKey {
            target: "target1".to_string(),
            chrom: "chr1".to_string(),
            pos: 100,
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
        };
        let entry = PonEntry {
            sample_count: 5,
            total_samples: 10,
            frequency: 0.5,
            max_vaf: 0.05,
            mean_vaf: 0.03,
        };
        let mut entries = HashMap::new();
        entries.insert(key.clone(), entry);

        let pon = PanelOfNormals {
            entries,
            total_samples: 10,
            creation_date: "2026-03-13".to_string(),
        };

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test_pon.json");

        pon.save(&path).unwrap();

        let loaded = PanelOfNormals::load(&path).unwrap();
        assert_eq!(loaded.total_samples, 10);
        assert_eq!(loaded.creation_date, "2026-03-13");
        assert_eq!(loaded.entries.len(), 1);

        let loaded_entry = loaded.entries.get(&key).unwrap();
        assert_eq!(loaded_entry.sample_count, 5);
        assert!((loaded_entry.frequency - 0.5).abs() < 1e-10);
        assert!((loaded_entry.max_vaf - 0.05).abs() < 1e-10);
        assert!((loaded_entry.mean_vaf - 0.03).abs() < 1e-10);
    }

    // ── PoN build from TSV files ────────────────────────────────────────

    /// Helper: write a detection TSV file for building PoN tests.
    fn write_detection_tsv(records: &[&str]) -> tempfile::NamedTempFile {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            f,
            "sample\ttarget\ttype\tvariant_name\trVAF\texpression\t\
             min_coverage\tstart_kmer_count\tref_sequence\talt_sequence\t\
             info\tchrom\tpos\tref_allele\talt_allele\tpvalue\tqual\tci_lower\tci_upper"
        )
        .unwrap();
        for record in records {
            writeln!(f, "{}", record).unwrap();
        }
        f.flush().unwrap();
        f
    }

    #[test]
    fn test_build_from_multiple_files() {
        // Normal sample 1: has SNV at chr1:100 A>T
        let f1 = write_detection_tsv(&[
            "normal1\ttarget1\tSubstitution\t100:A/T:100\t0.03\t10.0\t20\t50\t\t\tvs_ref\tchr1\t100\tA\tT\t\t\t\t",
            "normal1\ttarget1\tReference\tref\t0.0\t5.0\t10\t30\t\t\tvs_ref\t\t\t\t\t\t\t\t",
        ]);

        // Normal sample 2: has same SNV at chr1:100 A>T
        let f2 = write_detection_tsv(&[
            "normal2\ttarget1\tSubstitution\t100:A/T:100\t0.05\t15.0\t25\t60\t\t\tvs_ref\tchr1\t100\tA\tT\t\t\t\t",
        ]);

        // Normal sample 3: has different variant at chr2:200 G>C
        let f3 = write_detection_tsv(&[
            "normal3\ttarget2\tSubstitution\t200:G/C:200\t0.02\t8.0\t15\t40\t\t\tvs_ref\tchr2\t200\tG\tC\t\t\t\t",
        ]);

        let files: Vec<&Path> = vec![f1.path(), f2.path(), f3.path()];
        let pon = PanelOfNormals::build(&files, 0.0).unwrap();

        assert_eq!(pon.total_samples, 3);
        assert_eq!(pon.entries.len(), 2);

        // chr1:100 A>T should be in 2/3 normals.
        let key1 = PonKey {
            target: "target1".to_string(),
            chrom: "chr1".to_string(),
            pos: 100,
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
        };
        let entry1 = pon.entries.get(&key1).unwrap();
        assert_eq!(entry1.sample_count, 2);
        assert!((entry1.frequency - 2.0 / 3.0).abs() < 1e-10);
        assert!((entry1.max_vaf - 0.05).abs() < 1e-10);
        assert!((entry1.mean_vaf - 0.04).abs() < 1e-10);

        // chr2:200 G>C should be in 1/3 normals.
        let key2 = PonKey {
            target: "target2".to_string(),
            chrom: "chr2".to_string(),
            pos: 200,
            ref_allele: "G".to_string(),
            alt_allele: "C".to_string(),
        };
        let entry2 = pon.entries.get(&key2).unwrap();
        assert_eq!(entry2.sample_count, 1);
        assert!((entry2.frequency - 1.0 / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_build_skips_reference_calls() {
        let f1 = write_detection_tsv(&[
            "normal1\ttarget1\tReference\tref\t0.0\t5.0\t10\t30\t\t\tvs_ref\t\t\t\t\t\t\t\t",
        ]);

        let files: Vec<&Path> = vec![f1.path()];
        let pon = PanelOfNormals::build(&files, 0.0).unwrap();

        assert!(pon.entries.is_empty());
    }

    #[test]
    fn test_build_respects_min_vaf() {
        let f1 = write_detection_tsv(&[
            "normal1\ttarget1\tSubstitution\t100:A/T:100\t0.001\t10.0\t20\t50\t\t\tvs_ref\tchr1\t100\tA\tT\t\t\t\t",
            "normal1\ttarget1\tSubstitution\t200:G/C:200\t0.05\t15.0\t25\t60\t\t\tvs_ref\tchr2\t200\tG\tC\t\t\t\t",
        ]);

        let files: Vec<&Path> = vec![f1.path()];
        // min_vaf = 0.01 should exclude the first variant (VAF=0.001).
        let pon = PanelOfNormals::build(&files, 0.01).unwrap();

        assert_eq!(pon.entries.len(), 1);
        let key = PonKey {
            target: "target1".to_string(),
            chrom: "chr2".to_string(),
            pos: 200,
            ref_allele: "G".to_string(),
            alt_allele: "C".to_string(),
        };
        assert!(pon.entries.contains_key(&key));
    }

    #[test]
    fn test_build_same_variant_same_sample_counts_once() {
        // Same sample reports same variant twice (e.g., from overlapping targets).
        let f1 = write_detection_tsv(&[
            "normal1\ttarget1\tSubstitution\t100:A/T:100\t0.03\t10.0\t20\t50\t\t\tvs_ref\tchr1\t100\tA\tT\t\t\t\t",
            "normal1\ttarget1\tSubstitution\t100:A/T:100\t0.05\t15.0\t25\t60\t\t\tvs_ref\tchr1\t100\tA\tT\t\t\t\t",
        ]);

        let files: Vec<&Path> = vec![f1.path()];
        let pon = PanelOfNormals::build(&files, 0.0).unwrap();

        let key = PonKey {
            target: "target1".to_string(),
            chrom: "chr1".to_string(),
            pos: 100,
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
        };
        let entry = pon.entries.get(&key).unwrap();
        // Only one unique sample.
        assert_eq!(entry.sample_count, 1);
        // But both VAFs are recorded for statistics.
        assert!((entry.max_vaf - 0.05).abs() < 1e-10);
        assert!((entry.mean_vaf - 0.04).abs() < 1e-10);
    }

    // ── INDEL normalization in PoN ──────────────────────────────────────

    #[test]
    fn test_indel_normalization_matches_equivalent_variants() {
        // Two calls representing the same INDEL at different positions.
        // Without ref_context, normalization falls back to suffix/prefix trimming.
        // Call 1: chr1:100 AGT/AT (deletion of G, suffix T trimmed -> AG/A at 100)
        let call1 = make_call(
            "tumor",
            "target1",
            VariantType::Deletion,
            0.05,
            "chr1",
            100,
            "AGT",
            "AT",
        );
        // Call 2: chr1:100 AG/A (already normalized)
        let call2 = make_call(
            "normal",
            "target1",
            VariantType::Deletion,
            0.03,
            "chr1",
            100,
            "AG",
            "A",
        );

        let key1 = make_pon_key(&call1);
        let key2 = make_pon_key(&call2);
        assert_eq!(key1, key2);
    }

    // ── Misc ────────────────────────────────────────────────────────────

    #[test]
    fn test_len_and_is_empty() {
        let empty_pon = PanelOfNormals {
            entries: HashMap::new(),
            total_samples: 0,
            creation_date: "2026-01-01".to_string(),
        };
        assert!(empty_pon.is_empty());
        assert_eq!(empty_pon.len(), 0);

        let mut entries = HashMap::new();
        entries.insert(
            PonKey {
                target: "t".to_string(),
                chrom: "chr1".to_string(),
                pos: 1,
                ref_allele: "A".to_string(),
                alt_allele: "T".to_string(),
            },
            PonEntry {
                sample_count: 1,
                total_samples: 1,
                frequency: 1.0,
                max_vaf: 0.05,
                mean_vaf: 0.05,
            },
        );
        let pon = PanelOfNormals {
            entries,
            total_samples: 1,
            creation_date: "2026-01-01".to_string(),
        };
        assert!(!pon.is_empty());
        assert_eq!(pon.len(), 1);
    }

    #[test]
    fn test_current_date_string_format() {
        let date = current_date_string();
        // Should be YYYY-MM-DD format.
        assert_eq!(date.len(), 10);
        assert_eq!(date.as_bytes()[4], b'-');
        assert_eq!(date.as_bytes()[7], b'-');
        // Year should be >= 2024.
        let year: i32 = date[0..4].parse().unwrap();
        assert!(year >= 2024);
    }
}
