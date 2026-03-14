use std::path::PathBuf;

use anyhow::Result;
use serde::Serialize;

/// Verbosity level for diagnostic reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, clap::ValueEnum)]
pub enum ReportLevel {
    /// One file per run: run_summary.json + detection_summary.tsv
    Summary,
    /// Per-target decision summary with result.json and decision.txt
    Standard,
    /// Adds walking and graph diagnostics (walk_log, pruning_log, paths, etc.)
    Detailed,
    /// Adds k-mer-level database queries (ref/alt/missing k-mers, coverage profile)
    Full,
}

impl Default for ReportLevel {
    fn default() -> Self {
        Self::Summary
    }
}

/// Manages writing diagnostic reports for a detection run.
#[derive(Debug)]
pub struct ReportWriter {
    /// Root directory for all report output.
    pub root_dir: PathBuf,
    /// Verbosity level controlling what gets written.
    pub level: ReportLevel,
}

impl ReportWriter {
    /// Create a new ReportWriter, setting up the directory structure.
    pub fn new(root_dir: PathBuf, level: ReportLevel) -> Result<Self> {
        std::fs::create_dir_all(&root_dir)?;
        Ok(Self { root_dir, level })
    }

    /// Get or create the directory for a specific target.
    pub fn target_dir(&self, target_name: &str) -> Result<PathBuf> {
        let safe_name: String = target_name
            .chars()
            .map(|c| {
                if c.is_alphanumeric() || c == '-' || c == '_' {
                    c
                } else {
                    '_'
                }
            })
            .collect();
        let dir = self.root_dir.join("targets").join(&safe_name);
        std::fs::create_dir_all(&dir)?;
        Ok(dir)
    }

    /// Write the run summary (always written at any level).
    pub fn write_run_summary(&self, summary: &RunSummary) -> Result<()> {
        let path = self.root_dir.join("run_summary.json");
        let json = serde_json::to_string_pretty(summary)?;
        std::fs::write(path, json)?;
        Ok(())
    }

    /// Write the detection summary TSV (one line per target).
    pub fn write_detection_summary(&self, entries: &[DetectionSummaryEntry]) -> Result<()> {
        let path = self.root_dir.join("detection_summary.tsv");
        let mut wtr = std::fs::File::create(path)?;
        use std::io::Write;
        writeln!(
            wtr,
            "target\tresult_type\trvaf\tcoverage\tn_paths\tn_variants"
        )?;
        for e in entries {
            writeln!(
                wtr,
                "{}\t{}\t{:.6}\t{}\t{}\t{}",
                e.target, e.result_type, e.rvaf, e.coverage, e.n_paths, e.n_variants
            )?;
        }
        Ok(())
    }

    /// Write per-target result and decision (standard+ levels).
    pub fn write_target_result(&self, target_name: &str, result: &TargetResult) -> Result<()> {
        if self.level < ReportLevel::Standard {
            return Ok(());
        }
        let dir = self.target_dir(target_name)?;

        // result.json
        let json = serde_json::to_string_pretty(result)?;
        std::fs::write(dir.join("result.json"), json)?;

        // decision.txt
        std::fs::write(dir.join("decision.txt"), &result.decision_text)?;

        Ok(())
    }

    /// Write walking diagnostics (detailed+ levels).
    pub fn write_walking_report(&self, target_name: &str, report: &WalkingReport) -> Result<()> {
        if self.level < ReportLevel::Detailed {
            return Ok(());
        }
        let dir = self.target_dir(target_name)?;
        let walk_dir = dir.join("walking");
        std::fs::create_dir_all(&walk_dir)?;

        // walk_log.tsv
        {
            let mut wtr = std::fs::File::create(walk_dir.join("walk_log.tsv"))?;
            use std::io::Write;
            writeln!(wtr, "kmer\tdirection\tcount\tthreshold\taccepted")?;
            for entry in &report.extensions {
                writeln!(
                    wtr,
                    "{}\t{}\t{}\t{}\t{}",
                    entry.kmer, entry.direction, entry.count, entry.threshold, entry.accepted
                )?;
            }
        }

        // branches_explored.tsv
        {
            let mut wtr = std::fs::File::create(walk_dir.join("branches_explored.tsv"))?;
            use std::io::Write;
            writeln!(wtr, "position\tref_count\talt_counts\tchosen")?;
            for b in &report.branches {
                let alt_str = b
                    .alt_counts
                    .iter()
                    .map(|c| c.to_string())
                    .collect::<Vec<_>>()
                    .join(",");
                writeln!(
                    wtr,
                    "{}\t{}\t{}\t{}",
                    b.position, b.ref_count, alt_str, b.chosen
                )?;
            }
        }

        // walk_stats.json
        let stats_json = serde_json::to_string_pretty(&report.stats)?;
        std::fs::write(walk_dir.join("walk_stats.json"), stats_json)?;

        Ok(())
    }

    /// Write graph diagnostics (detailed+ levels).
    pub fn write_graph_report(&self, target_name: &str, report: &GraphReport) -> Result<()> {
        if self.level < ReportLevel::Detailed {
            return Ok(());
        }
        let dir = self.target_dir(target_name)?;
        let graph_dir = dir.join("graph");
        std::fs::create_dir_all(&graph_dir)?;

        // graph_stats.json
        let stats_json = serde_json::to_string_pretty(&report.stats)?;
        std::fs::write(graph_dir.join("graph_stats.json"), stats_json)?;

        // paths.tsv
        {
            let mut wtr = std::fs::File::create(graph_dir.join("paths.tsv"))?;
            use std::io::Write;
            writeln!(wtr, "path_id\tsequence\tweight\tis_ref")?;
            for p in &report.paths {
                writeln!(
                    wtr,
                    "{}\t{}\t{:.4}\t{}",
                    p.path_id, p.sequence, p.weight, p.is_ref
                )?;
            }
        }

        // pruning_log.tsv
        {
            let mut wtr = std::fs::File::create(graph_dir.join("pruning_log.tsv"))?;
            use std::io::Write;
            writeln!(wtr, "node\treason")?;
            for entry in &report.pruning_log {
                writeln!(wtr, "{}\t{}", entry.node, entry.reason)?;
            }
        }

        Ok(())
    }

    /// Write quantification diagnostics (detailed+ levels).
    pub fn write_quantification_report(
        &self,
        target_name: &str,
        report: &QuantificationReport,
    ) -> Result<()> {
        if self.level < ReportLevel::Detailed {
            return Ok(());
        }
        let dir = self.target_dir(target_name)?;
        let quant_dir = dir.join("quantification");
        std::fs::create_dir_all(&quant_dir)?;

        // nnls_matrix.tsv
        {
            let mut wtr = std::fs::File::create(quant_dir.join("nnls_matrix.tsv"))?;
            use std::io::Write;
            // Write column headers (path indices)
            let header: Vec<String> = (0..report.n_paths).map(|i| format!("path_{}", i)).collect();
            writeln!(wtr, "kmer\t{}", header.join("\t"))?;
            for row in &report.matrix_rows {
                let vals: Vec<String> = row.values.iter().map(|v| format!("{:.4}", v)).collect();
                writeln!(wtr, "{}\t{}", row.kmer, vals.join("\t"))?;
            }
        }

        // coefficients.json
        let coeff_json = serde_json::to_string_pretty(&report.coefficients)?;
        std::fs::write(quant_dir.join("coefficients.json"), coeff_json)?;

        // rvaf_calculation.json
        let rvaf_json = serde_json::to_string_pretty(&report.rvaf_details)?;
        std::fs::write(quant_dir.join("rvaf_calculation.json"), rvaf_json)?;

        Ok(())
    }

    /// Write k-mer diagnostics (full level only).
    pub fn write_kmer_report(&self, target_name: &str, report: &KmerReport) -> Result<()> {
        if self.level < ReportLevel::Full {
            return Ok(());
        }
        let dir = self.target_dir(target_name)?;
        let kmer_dir = dir.join("kmers");
        std::fs::create_dir_all(&kmer_dir)?;

        // ref_kmers.tsv
        {
            let mut wtr = std::fs::File::create(kmer_dir.join("ref_kmers.tsv"))?;
            use std::io::Write;
            writeln!(wtr, "kmer\tcanonical\tdb_count\texpected_count")?;
            for k in &report.ref_kmers {
                writeln!(
                    wtr,
                    "{}\t{}\t{}\t{}",
                    k.kmer, k.canonical, k.db_count, k.expected_count
                )?;
            }
        }

        // alt_kmers.tsv
        {
            let mut wtr = std::fs::File::create(kmer_dir.join("alt_kmers.tsv"))?;
            use std::io::Write;
            writeln!(wtr, "kmer\tdb_count\tpath_id")?;
            for k in &report.alt_kmers {
                writeln!(wtr, "{}\t{}\t{}", k.kmer, k.db_count, k.path_id)?;
            }
        }

        // missing_kmers.tsv
        {
            let mut wtr = std::fs::File::create(kmer_dir.join("missing_kmers.tsv"))?;
            use std::io::Write;
            writeln!(wtr, "kmer\texpected_in\treason")?;
            for k in &report.missing_kmers {
                writeln!(wtr, "{}\t{}\t{}", k.kmer, k.expected_in, k.reason)?;
            }
        }

        // coverage_profile.tsv
        {
            let mut wtr = std::fs::File::create(kmer_dir.join("coverage_profile.tsv"))?;
            use std::io::Write;
            writeln!(wtr, "position\tcoverage")?;
            for (pos, cov) in report.coverage_profile.iter().enumerate() {
                writeln!(wtr, "{}\t{}", pos, cov)?;
            }
        }

        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Serializable data structs
// ---------------------------------------------------------------------------

#[derive(Debug, Serialize)]
pub struct RunSummary {
    pub targets_processed: usize,
    pub variants_found: usize,
    pub targets_with_variants: usize,
    pub targets_reference_only: usize,
    pub targets_failed: usize,
    pub total_time_ms: u64,
    pub parameters: RunParameters,
}

#[derive(Debug, Serialize)]
pub struct RunParameters {
    pub kmer_length: u8,
    pub count: u32,
    pub ratio: f64,
    pub adaptive: bool,
    pub bidirectional: bool,
    pub prune_enabled: bool,
}

#[derive(Debug, Serialize)]
pub struct DetectionSummaryEntry {
    pub target: String,
    pub result_type: String,
    pub rvaf: f64,
    pub coverage: u64,
    pub n_paths: usize,
    pub n_variants: usize,
}

#[derive(Debug, Serialize)]
pub struct TargetResult {
    pub target: String,
    pub result_type: String,
    pub variants: Vec<TargetVariant>,
    pub n_ref_kmers: usize,
    pub n_walked_kmers: usize,
    pub n_graph_nodes: usize,
    pub n_paths: usize,
    pub decision_text: String,
}

#[derive(Debug, Serialize)]
pub struct TargetVariant {
    pub variant_type: String,
    pub variant_name: String,
    pub rvaf: f64,
    pub min_coverage: u64,
    pub ref_allele: String,
    pub alt_allele: String,
}

// Walking report types

#[derive(Debug, Serialize)]
pub struct WalkingReport {
    pub extensions: Vec<WalkExtension>,
    pub branches: Vec<WalkBranch>,
    pub stats: WalkStats,
}

#[derive(Debug, Serialize)]
pub struct WalkExtension {
    pub kmer: String,
    pub direction: String,
    pub count: u64,
    pub threshold: u64,
    pub accepted: bool,
}

#[derive(Debug, Serialize)]
pub struct WalkBranch {
    pub position: usize,
    pub ref_count: u64,
    pub alt_counts: Vec<u64>,
    pub chosen: String,
}

#[derive(Debug, Serialize)]
pub struct WalkStats {
    pub total_extensions: usize,
    pub accepted: usize,
    pub rejected: usize,
    pub max_depth: usize,
    pub dead_ends: usize,
    pub branch_points: usize,
}

// Graph report types

#[derive(Debug, Serialize)]
pub struct GraphReport {
    pub stats: GraphStats,
    pub paths: Vec<PathEntry>,
    pub pruning_log: Vec<PruningEntry>,
}

#[derive(Debug, Serialize)]
pub struct GraphStats {
    pub n_nodes: usize,
    pub n_edges: usize,
    pub n_ref_nodes: usize,
    pub n_alt_nodes: usize,
    pub n_pruned: usize,
}

#[derive(Debug, Serialize)]
pub struct PathEntry {
    pub path_id: usize,
    pub sequence: String,
    pub weight: f64,
    pub is_ref: bool,
}

#[derive(Debug, Serialize)]
pub struct PruningEntry {
    pub node: String,
    pub reason: String,
}

// Quantification report types

#[derive(Debug, Serialize)]
pub struct QuantificationReport {
    pub n_paths: usize,
    pub matrix_rows: Vec<MatrixRow>,
    pub coefficients: Vec<f64>,
    pub rvaf_details: RvafDetails,
}

#[derive(Debug, Serialize)]
pub struct MatrixRow {
    pub kmer: String,
    pub values: Vec<f64>,
}

#[derive(Debug, Serialize)]
pub struct RvafDetails {
    pub coefficients: Vec<f64>,
    pub sum_coefficients: f64,
    pub rvafs: Vec<f64>,
    pub min_coverages: Vec<u64>,
}

// K-mer report types

#[derive(Debug, Serialize)]
pub struct KmerReport {
    pub ref_kmers: Vec<RefKmerEntry>,
    pub alt_kmers: Vec<AltKmerEntry>,
    pub missing_kmers: Vec<MissingKmerEntry>,
    pub coverage_profile: Vec<u64>,
}

#[derive(Debug, Serialize)]
pub struct RefKmerEntry {
    pub kmer: String,
    pub canonical: String,
    pub db_count: u64,
    pub expected_count: u64,
}

#[derive(Debug, Serialize)]
pub struct AltKmerEntry {
    pub kmer: String,
    pub db_count: u64,
    pub path_id: usize,
}

#[derive(Debug, Serialize)]
pub struct MissingKmerEntry {
    pub kmer: String,
    pub expected_in: String,
    pub reason: String,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_report_level_ordering() {
        assert!(ReportLevel::Summary < ReportLevel::Standard);
        assert!(ReportLevel::Standard < ReportLevel::Detailed);
        assert!(ReportLevel::Detailed < ReportLevel::Full);
    }

    #[test]
    fn test_report_level_default() {
        assert_eq!(ReportLevel::default(), ReportLevel::Summary);
    }

    #[test]
    fn test_report_writer_creates_root_dir() {
        let tmp = tempfile::tempdir().unwrap();
        let root = tmp.path().join("reports");
        let writer = ReportWriter::new(root.clone(), ReportLevel::Summary).unwrap();
        assert!(writer.root_dir.exists());
    }

    #[test]
    fn test_target_dir_sanitizes_name() {
        let tmp = tempfile::tempdir().unwrap();
        let writer = ReportWriter::new(tmp.path().to_path_buf(), ReportLevel::Standard).unwrap();
        let dir = writer.target_dir("chr1:12345:A/T").unwrap();
        // Colons and slash should be replaced with underscores.
        let dir_name = dir.file_name().unwrap().to_str().unwrap();
        assert_eq!(dir_name, "chr1_12345_A_T");
        assert!(dir.exists());
    }

    #[test]
    fn test_write_run_summary() {
        let tmp = tempfile::tempdir().unwrap();
        let writer = ReportWriter::new(tmp.path().to_path_buf(), ReportLevel::Summary).unwrap();
        let summary = RunSummary {
            targets_processed: 100,
            variants_found: 5,
            targets_with_variants: 3,
            targets_reference_only: 95,
            targets_failed: 2,
            total_time_ms: 12345,
            parameters: RunParameters {
                kmer_length: 31,
                count: 2,
                ratio: 0.05,
                adaptive: false,
                bidirectional: false,
                prune_enabled: true,
            },
        };
        writer.write_run_summary(&summary).unwrap();

        let path = tmp.path().join("run_summary.json");
        assert!(path.exists());
        let content = std::fs::read_to_string(path).unwrap();
        let parsed: serde_json::Value = serde_json::from_str(&content).unwrap();
        assert_eq!(parsed["targets_processed"], 100);
        assert_eq!(parsed["variants_found"], 5);
        assert_eq!(parsed["parameters"]["kmer_length"], 31);
    }

    #[test]
    fn test_write_detection_summary() {
        let tmp = tempfile::tempdir().unwrap();
        let writer = ReportWriter::new(tmp.path().to_path_buf(), ReportLevel::Summary).unwrap();
        let entries = vec![
            DetectionSummaryEntry {
                target: "BRAF_V600E".to_string(),
                result_type: "Substitution".to_string(),
                rvaf: 0.032,
                coverage: 500,
                n_paths: 2,
                n_variants: 1,
            },
            DetectionSummaryEntry {
                target: "KRAS_G12D".to_string(),
                result_type: "Reference".to_string(),
                rvaf: 1.0,
                coverage: 800,
                n_paths: 1,
                n_variants: 0,
            },
        ];
        writer.write_detection_summary(&entries).unwrap();

        let path = tmp.path().join("detection_summary.tsv");
        assert!(path.exists());
        let content = std::fs::read_to_string(path).unwrap();
        assert!(content.starts_with("target\tresult_type\trvaf\tcoverage\tn_paths\tn_variants\n"));
        assert!(content.contains("BRAF_V600E\tSubstitution\t0.032000\t500\t2\t1\n"));
        assert!(content.contains("KRAS_G12D\tReference\t1.000000\t800\t1\t0\n"));
    }

    #[test]
    fn test_write_target_result_skipped_at_summary_level() {
        let tmp = tempfile::tempdir().unwrap();
        let writer = ReportWriter::new(tmp.path().to_path_buf(), ReportLevel::Summary).unwrap();
        let result = TargetResult {
            target: "test".to_string(),
            result_type: "Reference".to_string(),
            variants: vec![],
            n_ref_kmers: 10,
            n_walked_kmers: 15,
            n_graph_nodes: 20,
            n_paths: 1,
            decision_text: "Reference only".to_string(),
        };
        writer.write_target_result("test", &result).unwrap();
        // At summary level, no target directory should be created.
        assert!(!tmp.path().join("targets").exists());
    }

    #[test]
    fn test_write_target_result_at_standard_level() {
        let tmp = tempfile::tempdir().unwrap();
        let writer = ReportWriter::new(tmp.path().to_path_buf(), ReportLevel::Standard).unwrap();
        let result = TargetResult {
            target: "BRAF_V600E".to_string(),
            result_type: "Substitution".to_string(),
            variants: vec![TargetVariant {
                variant_type: "Substitution".to_string(),
                variant_name: "600:T/A".to_string(),
                rvaf: 0.032,
                min_coverage: 50,
                ref_allele: "T".to_string(),
                alt_allele: "A".to_string(),
            }],
            n_ref_kmers: 10,
            n_walked_kmers: 15,
            n_graph_nodes: 20,
            n_paths: 2,
            decision_text: "Detected Substitution at pos 600, rVAF=0.032".to_string(),
        };
        writer.write_target_result("BRAF_V600E", &result).unwrap();

        let target_dir = tmp.path().join("targets").join("BRAF_V600E");
        assert!(target_dir.join("result.json").exists());
        assert!(target_dir.join("decision.txt").exists());

        let decision = std::fs::read_to_string(target_dir.join("decision.txt")).unwrap();
        assert_eq!(decision, "Detected Substitution at pos 600, rVAF=0.032");
    }

    #[test]
    fn test_walking_report_skipped_below_detailed() {
        let tmp = tempfile::tempdir().unwrap();
        let writer = ReportWriter::new(tmp.path().to_path_buf(), ReportLevel::Standard).unwrap();
        let report = WalkingReport {
            extensions: vec![],
            branches: vec![],
            stats: WalkStats {
                total_extensions: 0,
                accepted: 0,
                rejected: 0,
                max_depth: 0,
                dead_ends: 0,
                branch_points: 0,
            },
        };
        writer.write_walking_report("test", &report).unwrap();
        // No walking directory should exist at standard level.
        assert!(!tmp.path().join("targets").exists());
    }

    #[test]
    fn test_kmer_report_skipped_below_full() {
        let tmp = tempfile::tempdir().unwrap();
        let writer = ReportWriter::new(tmp.path().to_path_buf(), ReportLevel::Detailed).unwrap();
        let report = KmerReport {
            ref_kmers: vec![],
            alt_kmers: vec![],
            missing_kmers: vec![],
            coverage_profile: vec![],
        };
        writer.write_kmer_report("test", &report).unwrap();
        // No kmer directory should exist at detailed level.
        assert!(!tmp.path().join("targets").exists());
    }
}
