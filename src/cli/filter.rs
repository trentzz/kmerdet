use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result};

use crate::filter::{ExpectedVariant, FilterResult};
use crate::variant::VariantCall;

#[derive(clap::Args, Debug)]
pub struct FilterArgs {
    /// Detection results input file (from detect/merge)
    #[arg(short, long)]
    pub input: PathBuf,

    /// Expected variants file (TSV: CHROM, POS, REF, ALT, TYPE)
    #[arg(long)]
    pub targets: PathBuf,

    /// Match by full ALT_SEQUENCE instead of POS/REF/ALT
    #[arg(long)]
    pub use_alt: bool,

    /// Min k-mer coverage
    #[arg(long, default_value = "3")]
    pub min_coverage: u32,

    /// Min variant allele frequency
    #[arg(long, default_value = "0.0")]
    pub min_vaf: f64,

    /// Min expression value
    #[arg(long, default_value = "0.0")]
    pub min_expression: f64,

    /// Variant types to include [default: all]
    #[arg(long, num_args = 1..)]
    pub types: Vec<String>,

    // ── Composite (adaptive) filter options ───────────────────────────

    /// Minimum QUAL score (Phred-scaled) for composite filtering
    #[arg(long)]
    pub min_qual: Option<f64>,

    /// Maximum p-value for composite filtering
    #[arg(long)]
    pub max_pvalue: Option<f64>,

    /// Enable sequence-context annotations (homopolymer, GC, motif flags)
    #[arg(long)]
    pub context_filter: bool,

    /// Homopolymer length to flag in context annotations [default: 6]
    #[arg(long, default_value = "6")]
    pub homopolymer_threshold: usize,
}

impl FilterArgs {
    /// Apply values from a config file for any args that were not explicitly
    /// set on the command line.
    pub(crate) fn apply_config(
        &mut self,
        cfg: &crate::config::Config,
        matches: &clap::ArgMatches,
    ) {
        let is_set = |id| super::is_subcommand_arg_set(matches, id);

        if !is_set("min_coverage") {
            if let Some(v) = cfg.filter.min_coverage {
                self.min_coverage = v;
            }
        }
        if !is_set("min_vaf") {
            if let Some(v) = cfg.filter.min_vaf {
                self.min_vaf = v;
            }
        }
        if !is_set("min_expression") {
            if let Some(v) = cfg.filter.min_expression {
                self.min_expression = v;
            }
        }
        if !is_set("use_alt") {
            if let Some(v) = cfg.filter.use_alt {
                self.use_alt = v;
            }
        }
        if self.types.is_empty() && !cfg.filter.types.is_empty() {
            self.types = cfg.filter.types.clone();
        }
    }
}

pub fn run(
    mut args: FilterArgs,
    global: &super::GlobalOptions,
    cfg: Option<&crate::config::Config>,
    matches: &clap::ArgMatches,
) -> Result<()> {
    if let Some(cfg) = cfg {
        args.apply_config(cfg, matches);
    }

    // 1. Parse detection results from input TSV file
    let calls = parse_detection_results(&args.input)?;

    // 2. Parse expected variants from --targets TSV file
    let expected = parse_expected_variants(&args.targets)?;

    // 3. Build FilterConfig from CLI args
    let config = crate::filter::FilterConfig {
        min_coverage: args.min_coverage,
        min_vaf: args.min_vaf,
        min_expression: args.min_expression,
        use_alt: args.use_alt,
        types: args.types.clone(),
    };

    // 4. Run filtering
    let mut results = crate::filter::filter_results(&calls, &expected, &config)?;

    // 5. Apply composite (adaptive) filter if any adaptive options are set
    let use_adaptive = args.min_qual.is_some()
        || args.max_pvalue.is_some()
        || args.context_filter;

    if use_adaptive {
        let adaptive_config = crate::filter::adaptive::AdaptiveFilterConfig {
            min_coverage: args.min_coverage as u64,
            min_vaf: args.min_vaf,
            min_expression: args.min_expression,
            allowed_types: args.types.clone(),
            min_qual: args.min_qual,
            max_pvalue: args.max_pvalue,
            context_filter: args.context_filter,
            homopolymer_threshold: args.homopolymer_threshold,
            ..Default::default()
        };

        // For each "Found" result, run the adaptive filter on the matching call
        for result in &mut results {
            if result.found != "Found" {
                continue;
            }

            // Find the matching call to run the adaptive stages on
            let matched_call = calls.iter().find(|c| {
                result.kmer_vaf == Some(c.rvaf)
                    && result.kmer_min_coverage == Some(c.min_coverage)
                    && result.ref_sequence.as_deref() == Some(&c.ref_sequence)
            });

            if let Some(call) = matched_call {
                let stage_result =
                    crate::filter::adaptive::apply_stages(call, &adaptive_config);

                if !stage_result.passed {
                    // Downgrade to "Not Found" and record why
                    result.found = "Not Found".to_string();
                    let reasons: Vec<String> = stage_result
                        .stage_results
                        .iter()
                        .filter(|s| !s.passed)
                        .filter_map(|s| s.reason.clone())
                        .collect();
                    result.filter_notes = reasons.join("; ");
                }

                // Append context flags to notes
                if !stage_result.context_flags.is_empty() {
                    let flag_str = stage_result.context_flags.join(",");
                    if result.filter_notes == "PASS" {
                        result.filter_notes =
                            format!("PASS [context: {}]", flag_str);
                    } else {
                        result.filter_notes =
                            format!("{} [context: {}]", result.filter_notes, flag_str);
                    }
                }
            }
        }
    }

    // 6. Write results as TSV
    let mut output: Box<dyn Write> = match &global.output {
        Some(path) => Box::new(
            File::create(path)
                .with_context(|| format!("creating output file: {}", path.display()))?,
        ),
        None => Box::new(std::io::stdout().lock()),
    };

    write_filter_results(&results, &mut output, global.no_header)?;

    Ok(())
}

/// Parse detection results from a TSV file into VariantCall structs.
fn parse_detection_results(path: &Path) -> Result<Vec<VariantCall>> {
    crate::output::parse_detection_tsv(path)
}

/// Parse expected variants from a TSV file.
///
/// Expected format: CHROM\tPOS\tREF\tALT\tTYPE (tab-separated, with header)
fn parse_expected_variants(path: &Path) -> Result<Vec<ExpectedVariant>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("opening expected variants: {}", path.display()))?;

    let mut variants = Vec::new();
    for record in reader.records() {
        let record = record?;
        variants.push(ExpectedVariant {
            chrom: record.get(0).unwrap_or("").to_string(),
            pos: record
                .get(1)
                .unwrap_or("0")
                .parse()
                .with_context(|| "parsing POS column in expected variants")?,
            ref_allele: record.get(2).unwrap_or("").to_string(),
            alt_allele: record.get(3).unwrap_or("").to_string(),
            variant_type: record.get(4).unwrap_or("").to_string(),
        });
    }
    Ok(variants)
}

/// Column headers for filter output.
const FILTER_HEADERS: &[&str] = &[
    "sample",
    "chrom",
    "pos",
    "ref_allele",
    "alt_allele",
    "variant_type",
    "found",
    "filter_notes",
    "kmer_vaf",
    "kmer_min_coverage",
    "kmer_expression",
    "ref_sequence",
    "variant_sequence",
];

/// Write filter results as TSV.
fn write_filter_results(
    results: &[FilterResult],
    writer: &mut dyn Write,
    no_header: bool,
) -> Result<()> {
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);

    if !no_header {
        wtr.write_record(FILTER_HEADERS)?;
    }

    for r in results {
        wtr.write_record(&[
            &r.sample,
            &r.chrom,
            &r.pos.to_string(),
            &r.ref_allele,
            &r.alt_allele,
            &r.variant_type,
            &r.found,
            &r.filter_notes,
            &r.kmer_vaf
                .map(|v| format!("{:.6}", v))
                .unwrap_or_default(),
            &r.kmer_min_coverage
                .map(|v| v.to_string())
                .unwrap_or_default(),
            &r.kmer_expression
                .map(|v| format!("{:.2}", v))
                .unwrap_or_default(),
            &r.ref_sequence.clone().unwrap_or_default(),
            &r.variant_sequence.clone().unwrap_or_default(),
        ])?;
    }

    wtr.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_write_filter_results_with_header() {
        let results = vec![FilterResult {
            sample: "sample1".to_string(),
            chrom: "chr1".to_string(),
            pos: 12345,
            ref_allele: "A".to_string(),
            alt_allele: "T".to_string(),
            variant_type: "Substitution".to_string(),
            found: "Found".to_string(),
            filter_notes: "PASS".to_string(),
            kmer_vaf: Some(0.05),
            kmer_min_coverage: Some(10),
            kmer_expression: Some(1.5),
            ref_sequence: Some("ACGT".to_string()),
            variant_sequence: Some("TCGT".to_string()),
        }];

        let mut buf = Vec::new();
        write_filter_results(&results, &mut buf, false).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(output.starts_with("sample\tchrom\tpos\t"));
        assert!(output.contains("sample1\tchr1\t12345\t"));
    }

    #[test]
    fn test_write_filter_results_no_header() {
        let results = vec![FilterResult {
            sample: "s1".to_string(),
            chrom: "chr2".to_string(),
            pos: 100,
            ref_allele: "G".to_string(),
            alt_allele: "C".to_string(),
            variant_type: "Substitution".to_string(),
            found: "Not Found".to_string(),
            filter_notes: "coverage 1 < 3".to_string(),
            kmer_vaf: None,
            kmer_min_coverage: None,
            kmer_expression: None,
            ref_sequence: None,
            variant_sequence: None,
        }];

        let mut buf = Vec::new();
        write_filter_results(&results, &mut buf, true).unwrap();
        let output = String::from_utf8(buf).unwrap();
        assert!(!output.starts_with("sample\t"));
        assert!(output.starts_with("s1\t"));
    }
}
