use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use serde::Deserialize;

/// Top-level TOML configuration file structure.
#[derive(Debug, Default, Deserialize)]
#[serde(default)]
pub struct Config {
    pub detect: DetectConfig,
    pub filter: FilterConfig,
    pub output: OutputConfig,
    pub runtime: RuntimeConfig,
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
pub struct DetectConfig {
    pub db: Option<PathBuf>,
    pub targets: Vec<PathBuf>,
    pub count: Option<u32>,
    pub ratio: Option<f64>,
    pub max_stack: Option<usize>,
    pub max_break: Option<usize>,
    pub max_node: Option<usize>,
    pub cluster: Option<bool>,
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
pub struct FilterConfig {
    pub targets: Option<PathBuf>,
    pub min_coverage: Option<u32>,
    pub min_vaf: Option<f64>,
    pub min_expression: Option<f64>,
    pub use_alt: Option<bool>,
    pub types: Vec<String>,
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
pub struct OutputConfig {
    pub format: Option<String>,
    pub file: Option<PathBuf>,
    pub no_header: Option<bool>,
}

#[derive(Debug, Default, Deserialize)]
#[serde(default)]
pub struct RuntimeConfig {
    pub threads: Option<usize>,
    pub verbose: Option<u8>,
}

impl Config {
    /// Load a TOML config file from the given path.
    pub fn load(path: &Path) -> Result<Self> {
        let contents =
            std::fs::read_to_string(path).with_context(|| format!("reading config: {}", path.display()))?;
        let config: Config =
            toml::from_str(&contents).with_context(|| format!("parsing config: {}", path.display()))?;
        Ok(config)
    }

    /// Parse a TOML string into a Config.
    pub fn from_str(s: &str) -> Result<Self> {
        let config: Config = toml::from_str(s).context("parsing config TOML")?;
        Ok(config)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_load_empty_config() {
        let cfg = Config::from_str("").unwrap();
        assert!(cfg.detect.db.is_none());
        assert!(cfg.detect.targets.is_empty());
        assert!(cfg.detect.count.is_none());
        assert!(cfg.detect.ratio.is_none());
        assert!(cfg.filter.min_coverage.is_none());
        assert!(cfg.output.format.is_none());
        assert!(cfg.runtime.threads.is_none());
    }

    #[test]
    fn test_load_detect_section() {
        let toml = r#"
[detect]
db = "/path/to/db.jf"
targets = ["/path/to/targets"]
count = 5
ratio = 0.10
max_stack = 1000
max_break = 20
max_node = 50000
cluster = true
"#;
        let cfg = Config::from_str(toml).unwrap();
        assert_eq!(cfg.detect.db.as_deref(), Some(Path::new("/path/to/db.jf")));
        assert_eq!(cfg.detect.targets.len(), 1);
        assert_eq!(cfg.detect.count, Some(5));
        assert_eq!(cfg.detect.ratio, Some(0.10));
        assert_eq!(cfg.detect.max_stack, Some(1000));
        assert_eq!(cfg.detect.max_break, Some(20));
        assert_eq!(cfg.detect.max_node, Some(50000));
        assert_eq!(cfg.detect.cluster, Some(true));
    }

    #[test]
    fn test_load_filter_section() {
        let toml = r#"
[filter]
targets = "/path/to/expected.tsv"
min_coverage = 10
min_vaf = 0.01
min_expression = 1.5
use_alt = true
types = ["Substitution", "Insertion"]
"#;
        let cfg = Config::from_str(toml).unwrap();
        assert_eq!(
            cfg.filter.targets.as_deref(),
            Some(Path::new("/path/to/expected.tsv"))
        );
        assert_eq!(cfg.filter.min_coverage, Some(10));
        assert_eq!(cfg.filter.min_vaf, Some(0.01));
        assert_eq!(cfg.filter.min_expression, Some(1.5));
        assert_eq!(cfg.filter.use_alt, Some(true));
        assert_eq!(cfg.filter.types, vec!["Substitution", "Insertion"]);
    }

    #[test]
    fn test_load_output_section() {
        let toml = r#"
[output]
format = "csv"
file = "/tmp/results.csv"
no_header = true
"#;
        let cfg = Config::from_str(toml).unwrap();
        assert_eq!(cfg.output.format.as_deref(), Some("csv"));
        assert_eq!(
            cfg.output.file.as_deref(),
            Some(Path::new("/tmp/results.csv"))
        );
        assert_eq!(cfg.output.no_header, Some(true));
    }

    #[test]
    fn test_load_runtime_section() {
        let toml = r#"
[runtime]
threads = 4
verbose = 2
"#;
        let cfg = Config::from_str(toml).unwrap();
        assert_eq!(cfg.runtime.threads, Some(4));
        assert_eq!(cfg.runtime.verbose, Some(2));
    }

    #[test]
    fn test_load_full_config() {
        let toml = r#"
[detect]
count = 3
ratio = 0.08

[filter]
min_coverage = 5
min_vaf = 0.005

[output]
format = "json"

[runtime]
threads = 8
"#;
        let cfg = Config::from_str(toml).unwrap();
        assert_eq!(cfg.detect.count, Some(3));
        assert_eq!(cfg.detect.ratio, Some(0.08));
        assert_eq!(cfg.filter.min_coverage, Some(5));
        assert_eq!(cfg.filter.min_vaf, Some(0.005));
        assert_eq!(cfg.output.format.as_deref(), Some("json"));
        assert_eq!(cfg.runtime.threads, Some(8));
    }

    #[test]
    fn test_partial_config_leaves_others_default() {
        let toml = r#"
[detect]
count = 10
"#;
        let cfg = Config::from_str(toml).unwrap();
        assert_eq!(cfg.detect.count, Some(10));
        // Everything else should be None/empty/default.
        assert!(cfg.detect.ratio.is_none());
        assert!(cfg.detect.db.is_none());
        assert!(cfg.filter.min_coverage.is_none());
        assert!(cfg.output.format.is_none());
        assert!(cfg.runtime.threads.is_none());
    }

    #[test]
    fn test_load_from_file() {
        let dir = tempfile::tempdir().unwrap();
        let config_path = dir.path().join("test.toml");
        std::fs::write(
            &config_path,
            r#"
[detect]
count = 7
ratio = 0.15
"#,
        )
        .unwrap();

        let cfg = Config::load(&config_path).unwrap();
        assert_eq!(cfg.detect.count, Some(7));
        assert_eq!(cfg.detect.ratio, Some(0.15));
    }

    #[test]
    fn test_load_nonexistent_file() {
        let result = Config::load(Path::new("/nonexistent/config.toml"));
        assert!(result.is_err());
    }

    #[test]
    fn test_load_invalid_toml() {
        let dir = tempfile::tempdir().unwrap();
        let config_path = dir.path().join("bad.toml");
        std::fs::write(&config_path, "this is not [valid toml {{{").unwrap();

        let result = Config::load(&config_path);
        assert!(result.is_err());
    }

    #[test]
    fn test_example_config_is_valid() {
        let example = include_str!("../docs/example_config.toml");
        let cfg = Config::from_str(example).expect("example_config.toml must be valid");
        assert_eq!(cfg.detect.count, Some(2));
        assert_eq!(cfg.detect.ratio, Some(0.05));
        assert_eq!(cfg.filter.min_coverage, Some(3));
        assert_eq!(cfg.output.format.as_deref(), Some("tsv"));
        assert_eq!(cfg.runtime.threads, Some(0));
    }

    #[test]
    fn test_unknown_keys_are_ignored() {
        // serde(default) + deny_unknown_fields is not set, so extra keys
        // should be silently ignored.
        let toml = r#"
[detect]
count = 3
some_future_option = "hello"
"#;
        let cfg = Config::from_str(toml).unwrap();
        assert_eq!(cfg.detect.count, Some(3));
    }
}
