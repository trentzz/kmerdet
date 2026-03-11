use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

use crate::variant::{VariantCall, VariantType};

#[derive(clap::Args, Debug)]
pub struct StatsArgs {
    /// Input result file(s) to summarize
    #[arg(short, long, num_args = 1..)]
    pub input: Vec<PathBuf>,

    /// Group statistics by this column
    #[arg(long)]
    pub group_by: Option<String>,
}

/// Summary statistics computed from detection results.
#[derive(Debug, serde::Serialize)]
struct Stats {
    total_targets: usize,
    detected_variants: usize,
    detection_rate: f64,
    vaf_mean: f64,
    vaf_median: f64,
    vaf_min: f64,
    vaf_max: f64,
    vaf_p5: f64,
    vaf_p25: f64,
    vaf_p75: f64,
    vaf_p95: f64,
    variant_types: HashMap<String, usize>,
    expression_mean: f64,
    expression_median: f64,
    expression_min: f64,
    expression_max: f64,
}

pub fn run(args: StatsArgs, global: &super::GlobalOptions) -> Result<()> {
    anyhow::ensure!(!args.input.is_empty(), "at least one input file required");

    // 1. Read all input files, parse into VariantCalls
    let mut all_calls = Vec::new();
    for path in &args.input {
        let calls = parse_detection_results(path)?;
        all_calls.extend(calls);
    }

    // 2. Prepare output
    let mut output: Box<dyn Write> = match &global.output {
        Some(path) => Box::new(
            File::create(path)
                .with_context(|| format!("creating output file: {}", path.display()))?,
        ),
        None => Box::new(std::io::stdout().lock()),
    };

    // 3. Group-by support or overall stats
    match &args.group_by {
        Some(column) => {
            let groups = group_calls(&all_calls, column)?;
            let mut group_names: Vec<&String> = groups.keys().collect();
            group_names.sort();

            match global.format {
                super::OutputFormat::Json => {
                    let mut map: HashMap<&String, Stats> = HashMap::new();
                    for name in &group_names {
                        map.insert(name, compute_stats(&groups[*name]));
                    }
                    writeln!(output, "{}", serde_json::to_string_pretty(&map)?)?;
                }
                super::OutputFormat::Tsv => {
                    write_stats_tsv_grouped(&groups, &group_names, column, &mut output)?;
                }
                _ => {
                    for name in &group_names {
                        writeln!(output, "\n=== {} = {} ===", column, name)?;
                        let stats = compute_stats(&groups[*name]);
                        write_stats_text(&stats, &mut output)?;
                    }
                }
            }
        }
        None => {
            let stats = compute_stats(&all_calls);
            match global.format {
                super::OutputFormat::Json => {
                    writeln!(output, "{}", serde_json::to_string_pretty(&stats)?)?;
                }
                super::OutputFormat::Tsv => {
                    write_stats_tsv(&stats, &mut output)?;
                }
                _ => {
                    writeln!(output, "=== kmerdet Detection Statistics ===")?;
                    write_stats_text(&stats, &mut output)?;
                }
            }
        }
    }

    Ok(())
}

/// Parse detection results from a TSV file into VariantCall structs.
fn parse_detection_results(path: &Path) -> Result<Vec<VariantCall>> {
    crate::output::parse_detection_tsv(path)
}

/// Group calls by a specified column name.
fn group_calls(
    calls: &[VariantCall],
    column: &str,
) -> Result<HashMap<String, Vec<VariantCall>>> {
    let mut groups: HashMap<String, Vec<VariantCall>> = HashMap::new();

    for call in calls {
        let key = match column {
            "sample" => call.sample.clone(),
            "target" => call.target.clone(),
            "type" => call.variant_type.to_string(),
            other => anyhow::bail!(
                "unsupported group-by column '{}'; supported: sample, target, type",
                other
            ),
        };
        groups.entry(key).or_default().push(call.clone());
    }

    Ok(groups)
}

/// Compute summary statistics from a set of variant calls.
fn compute_stats(calls: &[VariantCall]) -> Stats {
    let total = calls.len();
    let non_ref: Vec<&VariantCall> = calls
        .iter()
        .filter(|c| c.variant_type != VariantType::Reference)
        .collect();
    let detected = non_ref.len();
    let detection_rate = if total > 0 {
        detected as f64 / total as f64
    } else {
        0.0
    };

    // VAF statistics (from non-reference calls only)
    let mut vafs: Vec<f64> = non_ref.iter().map(|c| c.rvaf).collect();
    vafs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let vaf_mean = if vafs.is_empty() {
        0.0
    } else {
        vafs.iter().sum::<f64>() / vafs.len() as f64
    };
    let vaf_median = percentile(&vafs, 50.0);
    let vaf_min = vafs.first().copied().unwrap_or(0.0);
    let vaf_max = vafs.last().copied().unwrap_or(0.0);
    let vaf_p5 = percentile(&vafs, 5.0);
    let vaf_p25 = percentile(&vafs, 25.0);
    let vaf_p75 = percentile(&vafs, 75.0);
    let vaf_p95 = percentile(&vafs, 95.0);

    // Per-type counts
    let mut variant_types: HashMap<String, usize> = HashMap::new();
    for call in &non_ref {
        *variant_types
            .entry(call.variant_type.to_string())
            .or_insert(0) += 1;
    }

    // Expression statistics (from non-reference calls only)
    let mut expressions: Vec<f64> = non_ref.iter().map(|c| c.expression).collect();
    expressions.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));

    let expression_mean = if expressions.is_empty() {
        0.0
    } else {
        expressions.iter().sum::<f64>() / expressions.len() as f64
    };
    let expression_median = percentile(&expressions, 50.0);
    let expression_min = expressions.first().copied().unwrap_or(0.0);
    let expression_max = expressions.last().copied().unwrap_or(0.0);

    Stats {
        total_targets: total,
        detected_variants: detected,
        detection_rate,
        vaf_mean,
        vaf_median,
        vaf_min,
        vaf_max,
        vaf_p5,
        vaf_p25,
        vaf_p75,
        vaf_p95,
        variant_types,
        expression_mean,
        expression_median,
        expression_min,
        expression_max,
    }
}

/// Compute the p-th percentile of a sorted slice using linear interpolation.
fn percentile(sorted: &[f64], p: f64) -> f64 {
    if sorted.is_empty() {
        return 0.0;
    }
    if sorted.len() == 1 {
        return sorted[0];
    }
    let idx = (p / 100.0) * (sorted.len() - 1) as f64;
    let lower = idx.floor() as usize;
    let upper = idx.ceil() as usize;
    if lower == upper {
        return sorted[lower];
    }
    let frac = idx - lower as f64;
    sorted[lower] * (1.0 - frac) + sorted[upper] * frac
}

/// Write statistics in human-readable text format.
fn write_stats_text(stats: &Stats, writer: &mut dyn Write) -> Result<()> {
    writeln!(writer, "Total targets:      {}", stats.total_targets)?;
    writeln!(writer, "Detected variants:  {}", stats.detected_variants)?;
    writeln!(
        writer,
        "Detection rate:     {:.1}%",
        stats.detection_rate * 100.0
    )?;
    writeln!(writer)?;
    writeln!(writer, "VAF Summary:")?;
    writeln!(writer, "  Mean:    {:.4}", stats.vaf_mean)?;
    writeln!(writer, "  Median:  {:.4}", stats.vaf_median)?;
    writeln!(writer, "  Min:     {:.4}", stats.vaf_min)?;
    writeln!(writer, "  Max:     {:.4}", stats.vaf_max)?;
    writeln!(writer, "  P5:      {:.4}", stats.vaf_p5)?;
    writeln!(writer, "  P25:     {:.4}", stats.vaf_p25)?;
    writeln!(writer, "  P75:     {:.4}", stats.vaf_p75)?;
    writeln!(writer, "  P95:     {:.4}", stats.vaf_p95)?;
    writeln!(writer)?;

    if !stats.variant_types.is_empty() {
        writeln!(writer, "Variant Types:")?;
        let mut types: Vec<(&String, &usize)> = stats.variant_types.iter().collect();
        types.sort_by(|a, b| b.1.cmp(a.1)); // Sort by count descending
        for (vtype, count) in types {
            writeln!(writer, "  {:15}{}", format!("{}:", vtype), count)?;
        }
        writeln!(writer)?;
    }

    writeln!(writer, "Expression Summary:")?;
    writeln!(writer, "  Mean:    {:.2}", stats.expression_mean)?;
    writeln!(writer, "  Median:  {:.2}", stats.expression_median)?;
    writeln!(writer, "  Min:     {:.2}", stats.expression_min)?;
    writeln!(writer, "  Max:     {:.2}", stats.expression_max)?;

    Ok(())
}

/// Write statistics as TSV.
fn write_stats_tsv(stats: &Stats, writer: &mut dyn Write) -> Result<()> {
    writeln!(writer, "metric\tvalue")?;
    writeln!(writer, "total_targets\t{}", stats.total_targets)?;
    writeln!(writer, "detected_variants\t{}", stats.detected_variants)?;
    writeln!(writer, "detection_rate\t{:.4}", stats.detection_rate)?;
    writeln!(writer, "vaf_mean\t{:.6}", stats.vaf_mean)?;
    writeln!(writer, "vaf_median\t{:.6}", stats.vaf_median)?;
    writeln!(writer, "vaf_min\t{:.6}", stats.vaf_min)?;
    writeln!(writer, "vaf_max\t{:.6}", stats.vaf_max)?;
    writeln!(writer, "vaf_p5\t{:.6}", stats.vaf_p5)?;
    writeln!(writer, "vaf_p25\t{:.6}", stats.vaf_p25)?;
    writeln!(writer, "vaf_p75\t{:.6}", stats.vaf_p75)?;
    writeln!(writer, "vaf_p95\t{:.6}", stats.vaf_p95)?;

    let mut types: Vec<(&String, &usize)> = stats.variant_types.iter().collect();
    types.sort_by_key(|&(name, _)| name.clone());
    for (vtype, count) in types {
        writeln!(writer, "type_{}\t{}", vtype.to_lowercase(), count)?;
    }

    writeln!(writer, "expression_mean\t{:.2}", stats.expression_mean)?;
    writeln!(writer, "expression_median\t{:.2}", stats.expression_median)?;
    writeln!(writer, "expression_min\t{:.2}", stats.expression_min)?;
    writeln!(writer, "expression_max\t{:.2}", stats.expression_max)?;

    Ok(())
}

/// Write grouped statistics as TSV.
fn write_stats_tsv_grouped(
    groups: &HashMap<String, Vec<VariantCall>>,
    group_names: &[&String],
    column: &str,
    writer: &mut dyn Write,
) -> Result<()> {
    writeln!(
        writer,
        "{}\ttotal\tdetected\tdetection_rate\tvaf_mean\tvaf_median",
        column
    )?;
    for name in group_names {
        let stats = compute_stats(&groups[*name]);
        writeln!(
            writer,
            "{}\t{}\t{}\t{:.4}\t{:.6}\t{:.6}",
            name,
            stats.total_targets,
            stats.detected_variants,
            stats.detection_rate,
            stats.vaf_mean,
            stats.vaf_median
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_call(
        sample: &str,
        target: &str,
        vtype: VariantType,
        rvaf: f64,
        expression: f64,
    ) -> VariantCall {
        VariantCall {
            sample: sample.to_string(),
            target: target.to_string(),
            variant_type: vtype,
            variant_name: "test".to_string(),
            rvaf,
            expression,
            min_coverage: 10,
            start_kmer_count: 20,
            ref_sequence: "ACGT".to_string(),
            alt_sequence: "TCGT".to_string(),
            info: "vs_ref".to_string(),
            chrom: Some("chr1".to_string()),
            pos: Some(100),
            ref_allele: Some("A".to_string()),
            alt_allele: Some("T".to_string()),
        }
    }

    #[test]
    fn test_percentile_empty() {
        assert_eq!(percentile(&[], 50.0), 0.0);
    }

    #[test]
    fn test_percentile_single() {
        assert_eq!(percentile(&[5.0], 50.0), 5.0);
    }

    #[test]
    fn test_percentile_even() {
        let data = vec![1.0, 2.0, 3.0, 4.0];
        assert!((percentile(&data, 50.0) - 2.5).abs() < 1e-10);
        assert_eq!(percentile(&data, 0.0), 1.0);
        assert_eq!(percentile(&data, 100.0), 4.0);
    }

    #[test]
    fn test_compute_stats_basic() {
        let calls = vec![
            make_call("s1", "TP53", VariantType::Substitution, 0.05, 1.0),
            make_call("s1", "BRCA1", VariantType::Reference, 0.0, 0.0),
            make_call("s1", "EGFR", VariantType::Insertion, 0.10, 2.0),
            make_call("s1", "FLT3", VariantType::Itd, 0.20, 3.0),
        ];

        let stats = compute_stats(&calls);
        assert_eq!(stats.total_targets, 4);
        assert_eq!(stats.detected_variants, 3);
        assert!((stats.detection_rate - 0.75).abs() < 1e-10);
        // VAF mean of non-ref: (0.05 + 0.10 + 0.20) / 3
        let expected_mean = (0.05 + 0.10 + 0.20) / 3.0;
        assert!((stats.vaf_mean - expected_mean).abs() < 1e-10);
        assert_eq!(stats.vaf_min, 0.05);
        assert_eq!(stats.vaf_max, 0.20);
        assert_eq!(stats.variant_types.get("Substitution"), Some(&1));
        assert_eq!(stats.variant_types.get("Insertion"), Some(&1));
        assert_eq!(stats.variant_types.get("ITD"), Some(&1));
        assert_eq!(stats.variant_types.get("Reference"), None);
    }

    #[test]
    fn test_compute_stats_all_reference() {
        let calls = vec![
            make_call("s1", "TP53", VariantType::Reference, 0.0, 0.0),
            make_call("s1", "BRCA1", VariantType::Reference, 0.0, 0.0),
        ];

        let stats = compute_stats(&calls);
        assert_eq!(stats.total_targets, 2);
        assert_eq!(stats.detected_variants, 0);
        assert_eq!(stats.detection_rate, 0.0);
        assert_eq!(stats.vaf_mean, 0.0);
        assert!(stats.variant_types.is_empty());
    }

    #[test]
    fn test_compute_stats_empty() {
        let stats = compute_stats(&[]);
        assert_eq!(stats.total_targets, 0);
        assert_eq!(stats.detected_variants, 0);
        assert_eq!(stats.detection_rate, 0.0);
    }

    #[test]
    fn test_group_calls_by_sample() {
        let calls = vec![
            make_call("s1", "TP53", VariantType::Substitution, 0.05, 1.0),
            make_call("s2", "TP53", VariantType::Deletion, 0.10, 2.0),
            make_call("s1", "BRCA1", VariantType::Reference, 0.0, 0.0),
        ];

        let groups = group_calls(&calls, "sample").unwrap();
        assert_eq!(groups.len(), 2);
        assert_eq!(groups["s1"].len(), 2);
        assert_eq!(groups["s2"].len(), 1);
    }

    #[test]
    fn test_group_calls_invalid_column() {
        let calls = vec![make_call(
            "s1",
            "TP53",
            VariantType::Substitution,
            0.05,
            1.0,
        )];

        let result = group_calls(&calls, "invalid_col");
        assert!(result.is_err());
    }

    #[test]
    fn test_write_stats_text() {
        let calls = vec![
            make_call("s1", "TP53", VariantType::Substitution, 0.05, 1.0),
            make_call("s1", "BRCA1", VariantType::Reference, 0.0, 0.0),
        ];
        let stats = compute_stats(&calls);

        let mut buf = Vec::new();
        write_stats_text(&stats, &mut buf).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("Total targets:      2"));
        assert!(output.contains("Detected variants:  1"));
        assert!(output.contains("50.0%"));
        assert!(output.contains("Substitution:"));
    }

    #[test]
    fn test_write_stats_tsv() {
        let calls = vec![make_call(
            "s1",
            "TP53",
            VariantType::Substitution,
            0.05,
            1.0,
        )];
        let stats = compute_stats(&calls);

        let mut buf = Vec::new();
        write_stats_tsv(&stats, &mut buf).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.contains("metric\tvalue"));
        assert!(output.contains("total_targets\t1"));
        assert!(output.contains("vaf_mean\t"));
    }
}
