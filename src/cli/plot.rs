use std::path::PathBuf;

use anyhow::Result;

#[derive(clap::Args, Debug)]
pub struct PlotArgs {
    /// Input result file to plot
    #[arg(short, long)]
    pub input: PathBuf,

    /// Chart type to generate
    #[arg(long, default_value = "vaf-histogram")]
    pub chart: ChartType,

    /// Output image file (SVG or PNG)
    #[arg(short, long)]
    pub output: Option<PathBuf>,
}

#[derive(clap::ValueEnum, Debug, Clone, Copy)]
pub enum ChartType {
    /// VAF distribution histogram
    VafHistogram,
    /// Detection rate bar chart
    DetectionBar,
    /// Variant type distribution
    TypeDistribution,
    /// Summary pie chart
    SummaryPie,
}

pub fn run(_args: PlotArgs, _global: &super::GlobalOptions) -> Result<()> {
    anyhow::bail!("plot subcommand not yet implemented (Phase 7)")
}
