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
    pub sensitivity: SensitivityConfig,
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

/// Sensitivity configuration — controls the sensitivity/specificity tradeoff.
///
/// These parameters can be set individually for fine-grained control, or a named
/// `preset` can be used to set coordinated defaults. Individual values always
/// override the preset.
///
/// # Example TOML
///
/// ```toml
/// [sensitivity]
/// preset = "high"      # sets defaults, individual values below override
/// count = 1            # min absolute k-mer count for walking
/// ratio = 0.0001       # min ratio threshold for walking
/// min_qual = 5.0       # min Phred QUAL score for reporting
/// min_rvaf = 0.00001   # min rVAF for reporting
/// min_coverage = 1     # min k-mer coverage for reporting
/// ```
#[derive(Debug, Default, Deserialize, Clone)]
#[serde(default)]
pub struct SensitivityConfig {
    /// Named preset: "ultra", "high", "standard", "strict".
    /// Sets coordinated defaults; individual fields override.
    pub preset: Option<String>,
    /// Min absolute k-mer count for walk extension.
    pub count: Option<u32>,
    /// Min count ratio threshold for walk extension.
    pub ratio: Option<f64>,
    /// Min Phred QUAL score for reporting a variant call.
    pub min_qual: Option<f64>,
    /// Min rVAF for reporting a variant call.
    pub min_rvaf: Option<f64>,
    /// Min k-mer coverage for reporting a variant call.
    pub min_coverage: Option<u64>,
    /// Enable adaptive thresholds.
    pub adaptive: Option<bool>,
    /// Enable bidirectional walking.
    pub bidirectional: Option<bool>,
    /// Enable bootstrap CIs.
    pub bootstrap: Option<bool>,
}

/// Built-in sensitivity preset values.
#[derive(Debug, Clone, Copy)]
pub struct SensitivityPresetValues {
    pub count: u32,
    pub ratio: f64,
    pub min_qual: f64,
    pub min_rvaf: f64,
    pub min_coverage: u64,
}

impl SensitivityPresetValues {
    /// Look up a named preset. Returns `None` for unknown names.
    pub fn from_name(name: &str) -> Option<Self> {
        match name.to_lowercase().as_str() {
            "ultra" => Some(Self {
                count: 1,
                ratio: 1e-6,
                min_qual: 0.0,
                min_rvaf: 0.0,
                min_coverage: 1,
            }),
            "high" => Some(Self {
                count: 1,
                ratio: 1e-4,
                min_qual: 5.0,
                min_rvaf: 0.00001,
                min_coverage: 1,
            }),
            "standard" => Some(Self {
                count: 2,
                ratio: 0.05,
                min_qual: 10.0,
                min_rvaf: 0.0001,
                min_coverage: 3,
            }),
            "strict" => Some(Self {
                count: 3,
                ratio: 0.05,
                min_qual: 20.0,
                min_rvaf: 0.001,
                min_coverage: 5,
            }),
            _ => None,
        }
    }

    /// Return all known preset names.
    pub fn preset_names() -> &'static [&'static str] {
        &["ultra", "high", "standard", "strict"]
    }
}

/// Resolved sensitivity parameters — all fields are concrete values.
#[derive(Debug, Clone)]
pub struct ResolvedSensitivity {
    pub count: u32,
    pub ratio: f64,
    pub min_qual: f64,
    pub min_rvaf: f64,
    pub min_coverage: u64,
    pub adaptive: bool,
    pub bidirectional: bool,
    pub bootstrap: bool,
}

impl Default for ResolvedSensitivity {
    fn default() -> Self {
        let standard = SensitivityPresetValues::from_name("standard").unwrap();
        Self {
            count: standard.count,
            ratio: standard.ratio,
            min_qual: standard.min_qual,
            min_rvaf: standard.min_rvaf,
            min_coverage: standard.min_coverage,
            adaptive: false,
            bidirectional: false,
            bootstrap: false,
        }
    }
}

impl SensitivityConfig {
    /// Resolve a `SensitivityConfig` into concrete values.
    ///
    /// Resolution order (highest priority first):
    /// 1. Explicit fields in this config
    /// 2. Preset defaults (if `preset` is set)
    /// 3. Hard-coded "standard" defaults
    pub fn resolve(&self) -> Result<ResolvedSensitivity> {
        // Start from the preset or standard defaults.
        let base = if let Some(ref name) = self.preset {
            SensitivityPresetValues::from_name(name).ok_or_else(|| {
                anyhow::anyhow!(
                    "unknown sensitivity preset '{}'; valid presets: {}",
                    name,
                    SensitivityPresetValues::preset_names().join(", ")
                )
            })?
        } else {
            SensitivityPresetValues::from_name("standard").unwrap()
        };

        Ok(ResolvedSensitivity {
            count: self.count.unwrap_or(base.count),
            ratio: self.ratio.unwrap_or(base.ratio),
            min_qual: self.min_qual.unwrap_or(base.min_qual),
            min_rvaf: self.min_rvaf.unwrap_or(base.min_rvaf),
            min_coverage: self.min_coverage.unwrap_or(base.min_coverage),
            adaptive: self.adaptive.unwrap_or(false),
            bidirectional: self.bidirectional.unwrap_or(false),
            bootstrap: self.bootstrap.unwrap_or(false),
        })
    }
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

    #[test]
    fn test_sensitivity_config_preset_only() {
        let toml = r#"
[sensitivity]
preset = "ultra"
"#;
        let cfg = Config::from_str(toml).unwrap();
        let resolved = cfg.sensitivity.resolve().unwrap();
        assert_eq!(resolved.count, 1);
        assert!((resolved.ratio - 1e-6).abs() < 1e-12);
        assert_eq!(resolved.min_qual, 0.0);
        assert_eq!(resolved.min_rvaf, 0.0);
        assert_eq!(resolved.min_coverage, 1);
    }

    #[test]
    fn test_sensitivity_config_preset_with_override() {
        let toml = r#"
[sensitivity]
preset = "high"
count = 2
min_qual = 15.0
"#;
        let cfg = Config::from_str(toml).unwrap();
        let resolved = cfg.sensitivity.resolve().unwrap();
        // count and min_qual are overridden, rest comes from "high" preset
        assert_eq!(resolved.count, 2);
        assert!((resolved.ratio - 1e-4).abs() < 1e-12);
        assert_eq!(resolved.min_qual, 15.0);
        assert!((resolved.min_rvaf - 0.00001).abs() < 1e-12);
    }

    #[test]
    fn test_sensitivity_config_no_preset() {
        let toml = r#"
[sensitivity]
count = 5
ratio = 0.001
min_qual = 30.0
"#;
        let cfg = Config::from_str(toml).unwrap();
        let resolved = cfg.sensitivity.resolve().unwrap();
        assert_eq!(resolved.count, 5);
        assert!((resolved.ratio - 0.001).abs() < 1e-12);
        assert_eq!(resolved.min_qual, 30.0);
        // Unset fields fall back to "standard" defaults
        assert!((resolved.min_rvaf - 0.0001).abs() < 1e-12);
        assert_eq!(resolved.min_coverage, 3);
    }

    #[test]
    fn test_sensitivity_config_unknown_preset() {
        let toml = r#"
[sensitivity]
preset = "nonexistent"
"#;
        let cfg = Config::from_str(toml).unwrap();
        let result = cfg.sensitivity.resolve();
        assert!(result.is_err());
        let err = result.unwrap_err().to_string();
        assert!(err.contains("nonexistent"));
        assert!(err.contains("ultra"));
    }

    #[test]
    fn test_sensitivity_default_is_standard() {
        let cfg = Config::from_str("").unwrap();
        let resolved = cfg.sensitivity.resolve().unwrap();
        assert_eq!(resolved.count, 2);
        assert!((resolved.ratio - 0.05).abs() < 1e-12);
        assert_eq!(resolved.min_qual, 10.0);
    }

    #[test]
    fn test_sensitivity_all_presets() {
        for name in SensitivityPresetValues::preset_names() {
            let values = SensitivityPresetValues::from_name(name);
            assert!(values.is_some(), "preset '{}' should exist", name);
        }
    }

    #[test]
    fn test_sensitivity_strict_values() {
        let toml = r#"
[sensitivity]
preset = "strict"
"#;
        let cfg = Config::from_str(toml).unwrap();
        let resolved = cfg.sensitivity.resolve().unwrap();
        assert_eq!(resolved.count, 3);
        assert!((resolved.ratio - 0.05).abs() < 1e-12);
        assert_eq!(resolved.min_qual, 20.0);
        assert!((resolved.min_rvaf - 0.001).abs() < 1e-12);
        assert_eq!(resolved.min_coverage, 5);
    }
}
