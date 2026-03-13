use std::path::PathBuf;

use clap::{Parser, Subcommand, ValueEnum};

use crate::config::Config;

pub mod benchmark;
pub mod coverage;
pub mod detect;
pub mod filter;
pub mod merge;
pub mod plot;
pub mod pon;
pub mod run;
pub mod stats;

#[derive(Parser, Debug)]
#[command(
    name = "kmerdet",
    version,
    about = "K-mer based variant detection for targeted liquid biopsy monitoring",
    long_about = "kmerdet detects variants by walking k-mer graphs built from jellyfish databases.\n\n\
                  It replaces the km + kmtools + kam ecosystem with a single, fast Rust binary."
)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Command,

    #[command(flatten)]
    pub global: GlobalOptions,
}

#[derive(clap::Args, Debug, Clone)]
pub struct GlobalOptions {
    /// Config file (TOML)
    #[arg(short, long, global = true)]
    pub config: Option<PathBuf>,

    /// Thread count (0 = all available cores)
    #[arg(short = 't', long, global = true, default_value_t = 0)]
    pub threads: usize,

    /// Increase verbosity (-v, -vv, -vvv)
    #[arg(short, long, global = true, action = clap::ArgAction::Count)]
    pub verbose: u8,

    /// Suppress non-error output
    #[arg(short, long, global = true)]
    pub quiet: bool,

    /// Output file [default: stdout]
    #[arg(short, long, global = true)]
    pub output: Option<PathBuf>,

    /// Output format
    #[arg(short, long, global = true, default_value = "tsv")]
    pub format: OutputFormat,

    /// Omit column headers
    #[arg(long, global = true)]
    pub no_header: bool,
}

impl GlobalOptions {
    /// Apply values from a loaded config file. Config values are used only
    /// when the corresponding CLI option was not explicitly provided.
    ///
    /// The `matches` parameter is used to determine whether a CLI flag was
    /// explicitly set by the user (as opposed to being a clap default).
    pub fn apply_config(&mut self, cfg: &Config, matches: &clap::ArgMatches) {
        // threads: only override if the user did not pass --threads explicitly
        if !is_explicitly_set(matches, "threads") {
            if let Some(t) = cfg.runtime.threads {
                self.threads = t;
            }
        }

        // verbose: only override if not explicitly set (count action defaults to 0)
        if !is_explicitly_set(matches, "verbose") {
            if let Some(v) = cfg.runtime.verbose {
                self.verbose = v;
            }
        }

        // output file
        if self.output.is_none() {
            self.output = cfg.output.file.clone();
        }

        // format
        if !is_explicitly_set(matches, "format") {
            if let Some(ref fmt) = cfg.output.format {
                if let Ok(f) = fmt.parse::<OutputFormat>() {
                    self.format = f;
                }
            }
        }

        // no_header
        if !is_explicitly_set(matches, "no_header") {
            if let Some(nh) = cfg.output.no_header {
                self.no_header = nh;
            }
        }
    }
}

/// Load the config file from `--config` if present.
pub fn load_config(global: &GlobalOptions) -> anyhow::Result<Option<Config>> {
    match &global.config {
        Some(path) => {
            let cfg = Config::load(path)?;
            tracing::debug!("loaded config from {}", path.display());
            Ok(Some(cfg))
        }
        None => Ok(None),
    }
}

/// Check whether a given argument was explicitly set on the command line
/// (i.e. not just using the clap default value).
fn is_explicitly_set(matches: &clap::ArgMatches, id: &str) -> bool {
    matches
        .value_source(id)
        .map(|s| s == clap::parser::ValueSource::CommandLine)
        .unwrap_or(false)
}

/// Try to resolve `value_source` for a subcommand arg by walking into the
/// subcommand's matches.
pub fn subcommand_value_source(
    matches: &clap::ArgMatches,
    id: &str,
) -> Option<clap::parser::ValueSource> {
    // First check top-level (only if the id exists there).
    if matches.try_contains_id(id).is_ok() {
        if let Some(source) = matches.value_source(id) {
            return Some(source);
        }
    }
    // Then check within the active subcommand.
    if let Some((_name, sub_matches)) = matches.subcommand() {
        if sub_matches.try_contains_id(id).is_ok() {
            return sub_matches.value_source(id);
        }
    }
    None
}

/// Check if a subcommand argument was explicitly set on the command line.
pub fn is_subcommand_arg_set(matches: &clap::ArgMatches, id: &str) -> bool {
    subcommand_value_source(matches, id)
        .map(|s| s == clap::parser::ValueSource::CommandLine)
        .unwrap_or(false)
}

#[derive(ValueEnum, Debug, Clone, Copy, PartialEq, Eq)]
pub enum OutputFormat {
    Tsv,
    Csv,
    Vcf,
    Json,
    Jsonl,
    Xlsx,
}

impl std::str::FromStr for OutputFormat {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "tsv" => Ok(Self::Tsv),
            "csv" => Ok(Self::Csv),
            "vcf" => Ok(Self::Vcf),
            "json" => Ok(Self::Json),
            "jsonl" => Ok(Self::Jsonl),
            "xlsx" => Ok(Self::Xlsx),
            _ => Err(format!("unknown output format: '{}'", s)),
        }
    }
}

#[derive(Subcommand, Debug)]
pub enum Command {
    /// Run variant detection on targets against a jellyfish database
    Detect(detect::DetectArgs),
    /// Apply tumor-informed filtering to detection results
    Filter(filter::FilterArgs),
    /// Merge results from multiple detection runs
    Merge(merge::MergeArgs),
    /// Generate summary statistics from results
    Stats(stats::StatsArgs),
    /// Generate visualization charts
    Plot(plot::PlotArgs),
    /// Report k-mer coverage statistics for a database
    Coverage(coverage::CoverageArgs),
    /// Full pipeline: detect -> merge -> filter (-> optional stats/plot)
    Run(run::RunArgs),
    /// Compare detection results against ground truth
    Benchmark(benchmark::BenchmarkArgs),
    /// Panel-of-normals: build from normals or filter against a PoN
    Pon(pon::PonArgs),
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::{CommandFactory, FromArgMatches};

    /// Helper: parse CLI args and return both the typed Cli and raw ArgMatches.
    fn parse_cli(args: &[&str]) -> (Cli, clap::ArgMatches) {
        let matches = Cli::command().get_matches_from(args);
        let cli = Cli::from_arg_matches(&matches).unwrap();
        (cli, matches)
    }

    #[test]
    fn test_global_options_apply_config_threads() {
        let (mut cli, matches) = parse_cli(&[
            "kmerdet", "detect", "-d", "test.jf", "-T", "targets/",
        ]);
        let cfg = Config::from_str("[runtime]\nthreads = 8\n").unwrap();
        cli.global.apply_config(&cfg, &matches);
        assert_eq!(cli.global.threads, 8);
    }

    #[test]
    fn test_global_options_cli_overrides_config_threads() {
        let (mut cli, matches) = parse_cli(&[
            "kmerdet", "-t", "4", "detect", "-d", "test.jf", "-T", "targets/",
        ]);
        let cfg = Config::from_str("[runtime]\nthreads = 8\n").unwrap();
        cli.global.apply_config(&cfg, &matches);
        // CLI explicitly set 4, config says 8 -- CLI wins.
        assert_eq!(cli.global.threads, 4);
    }

    #[test]
    fn test_global_options_apply_config_format() {
        let (mut cli, matches) = parse_cli(&[
            "kmerdet", "detect", "-d", "test.jf", "-T", "targets/",
        ]);
        let cfg = Config::from_str("[output]\nformat = \"csv\"\n").unwrap();
        cli.global.apply_config(&cfg, &matches);
        assert_eq!(cli.global.format, OutputFormat::Csv);
    }

    #[test]
    fn test_global_options_cli_overrides_config_format() {
        let (mut cli, matches) = parse_cli(&[
            "kmerdet", "-f", "json", "detect", "-d", "test.jf", "-T", "targets/",
        ]);
        let cfg = Config::from_str("[output]\nformat = \"csv\"\n").unwrap();
        cli.global.apply_config(&cfg, &matches);
        assert_eq!(cli.global.format, OutputFormat::Json);
    }

    #[test]
    fn test_global_options_apply_config_output_file() {
        let (mut cli, matches) = parse_cli(&[
            "kmerdet", "detect", "-d", "test.jf", "-T", "targets/",
        ]);
        let cfg = Config::from_str("[output]\nfile = \"/tmp/out.tsv\"\n").unwrap();
        cli.global.apply_config(&cfg, &matches);
        assert_eq!(
            cli.global.output,
            Some(PathBuf::from("/tmp/out.tsv"))
        );
    }

    #[test]
    fn test_global_options_cli_overrides_config_output_file() {
        let (mut cli, matches) = parse_cli(&[
            "kmerdet", "-o", "/tmp/cli.tsv", "detect", "-d", "test.jf", "-T", "targets/",
        ]);
        let cfg = Config::from_str("[output]\nfile = \"/tmp/cfg.tsv\"\n").unwrap();
        cli.global.apply_config(&cfg, &matches);
        // CLI explicitly set /tmp/cli.tsv, so config is not used.
        assert_eq!(
            cli.global.output,
            Some(PathBuf::from("/tmp/cli.tsv"))
        );
    }

    #[test]
    fn test_global_options_apply_config_no_header() {
        let (mut cli, matches) = parse_cli(&[
            "kmerdet", "detect", "-d", "test.jf", "-T", "targets/",
        ]);
        let cfg = Config::from_str("[output]\nno_header = true\n").unwrap();
        cli.global.apply_config(&cfg, &matches);
        assert!(cli.global.no_header);
    }

    #[test]
    fn test_global_options_apply_config_verbose() {
        let (mut cli, matches) = parse_cli(&[
            "kmerdet", "detect", "-d", "test.jf", "-T", "targets/",
        ]);
        let cfg = Config::from_str("[runtime]\nverbose = 2\n").unwrap();
        cli.global.apply_config(&cfg, &matches);
        assert_eq!(cli.global.verbose, 2);
    }

    #[test]
    fn test_global_options_empty_config_leaves_defaults() {
        let (mut cli, matches) = parse_cli(&[
            "kmerdet", "detect", "-d", "test.jf", "-T", "targets/",
        ]);
        let cfg = Config::from_str("").unwrap();
        cli.global.apply_config(&cfg, &matches);
        assert_eq!(cli.global.threads, 0);
        assert_eq!(cli.global.verbose, 0);
        assert_eq!(cli.global.format, OutputFormat::Tsv);
        assert!(!cli.global.no_header);
        assert!(cli.global.output.is_none());
    }

    #[test]
    fn test_detect_args_apply_config() {
        let (cli, matches) = parse_cli(&[
            "kmerdet", "detect", "-d", "test.jf", "-T", "targets/",
        ]);
        let cfg = Config::from_str(
            r#"
[detect]
count = 5
ratio = 0.10
max_stack = 1000
max_break = 20
max_node = 50000
cluster = true
"#,
        )
        .unwrap();

        if let Command::Detect(mut args) = cli.command {
            args.apply_config(&cfg, &matches);
            assert_eq!(args.count, 5);
            assert!((args.ratio - 0.10).abs() < 1e-10);
            assert_eq!(args.max_stack, 1000);
            assert_eq!(args.max_break, 20);
            assert_eq!(args.max_node, 50000);
            assert!(args.cluster);
        } else {
            panic!("expected Detect command");
        }
    }

    #[test]
    fn test_detect_args_cli_overrides_config() {
        let (cli, matches) = parse_cli(&[
            "kmerdet",
            "detect",
            "-d", "test.jf",
            "-T", "targets/",
            "--count", "10",
            "--ratio", "0.20",
        ]);
        let cfg = Config::from_str(
            r#"
[detect]
count = 5
ratio = 0.10
max_stack = 1000
"#,
        )
        .unwrap();

        if let Command::Detect(mut args) = cli.command {
            args.apply_config(&cfg, &matches);
            // CLI values should win for count and ratio.
            assert_eq!(args.count, 10);
            assert!((args.ratio - 0.20).abs() < 1e-10);
            // max_stack was not set on CLI, so config value should apply.
            assert_eq!(args.max_stack, 1000);
        } else {
            panic!("expected Detect command");
        }
    }

    #[test]
    fn test_filter_args_apply_config() {
        let (cli, matches) = parse_cli(&[
            "kmerdet",
            "filter",
            "--input", "results.tsv",
            "--targets", "expected.tsv",
        ]);
        let cfg = Config::from_str(
            r#"
[filter]
min_coverage = 10
min_vaf = 0.01
min_expression = 1.5
use_alt = true
types = ["Substitution"]
"#,
        )
        .unwrap();

        if let Command::Filter(mut args) = cli.command {
            args.apply_config(&cfg, &matches);
            assert_eq!(args.min_coverage, 10);
            assert!((args.min_vaf - 0.01).abs() < 1e-10);
            assert!((args.min_expression - 1.5).abs() < 1e-10);
            assert!(args.use_alt);
            assert_eq!(args.types, vec!["Substitution"]);
        } else {
            panic!("expected Filter command");
        }
    }

    #[test]
    fn test_filter_args_cli_overrides_config() {
        let (cli, matches) = parse_cli(&[
            "kmerdet",
            "filter",
            "--input", "results.tsv",
            "--targets", "expected.tsv",
            "--min-coverage", "20",
        ]);
        let cfg = Config::from_str(
            r#"
[filter]
min_coverage = 10
min_vaf = 0.01
"#,
        )
        .unwrap();

        if let Command::Filter(mut args) = cli.command {
            args.apply_config(&cfg, &matches);
            // CLI value wins for min_coverage.
            assert_eq!(args.min_coverage, 20);
            // Config value applies for min_vaf (not explicitly set on CLI).
            assert!((args.min_vaf - 0.01).abs() < 1e-10);
        } else {
            panic!("expected Filter command");
        }
    }

    #[test]
    fn test_run_args_apply_config() {
        let (cli, matches) = parse_cli(&[
            "kmerdet",
            "run",
            "-d", "test.jf",
            "--targets-dir", "targets/",
            "--expected", "expected.tsv",
        ]);
        let cfg = Config::from_str(
            r#"
[detect]
count = 5
ratio = 0.10
cluster = true

[filter]
min_coverage = 10
min_vaf = 0.01
use_alt = true
"#,
        )
        .unwrap();

        if let Command::Run(mut args) = cli.command {
            args.apply_config(&cfg, &matches);
            assert_eq!(args.count, 5);
            assert!((args.ratio - 0.10).abs() < 1e-10);
            assert!(args.cluster);
            assert_eq!(args.min_coverage, 10);
            assert!((args.min_vaf - 0.01).abs() < 1e-10);
            assert!(args.use_alt);
        } else {
            panic!("expected Run command");
        }
    }

    #[test]
    fn test_output_format_from_str() {
        assert_eq!("tsv".parse::<OutputFormat>().unwrap(), OutputFormat::Tsv);
        assert_eq!("CSV".parse::<OutputFormat>().unwrap(), OutputFormat::Csv);
        assert_eq!("Json".parse::<OutputFormat>().unwrap(), OutputFormat::Json);
        assert_eq!("JSONL".parse::<OutputFormat>().unwrap(), OutputFormat::Jsonl);
        assert_eq!("xlsx".parse::<OutputFormat>().unwrap(), OutputFormat::Xlsx);
        assert_eq!("vcf".parse::<OutputFormat>().unwrap(), OutputFormat::Vcf);
        assert!("unknown".parse::<OutputFormat>().is_err());
    }

    #[test]
    fn test_is_subcommand_arg_set_explicit() {
        let matches = Cli::command().get_matches_from([
            "kmerdet",
            "detect",
            "-d", "test.jf",
            "-T", "targets/",
            "--count", "10",
        ]);
        assert!(is_subcommand_arg_set(&matches, "count"));
    }

    #[test]
    fn test_is_subcommand_arg_set_default() {
        let matches = Cli::command().get_matches_from([
            "kmerdet",
            "detect",
            "-d", "test.jf",
            "-T", "targets/",
        ]);
        // count has a default of 2, but was not explicitly set.
        assert!(!is_subcommand_arg_set(&matches, "count"));
    }
}
