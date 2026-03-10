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
}
