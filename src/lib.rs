pub mod benchmark;
pub mod cli;
pub mod confidence;
pub mod config;
pub mod filter;
pub mod graph;
pub mod jellyfish;
pub mod kmer;
pub mod output;
pub mod sequence;
pub mod trace;
pub mod variant;
pub mod walker;

use anyhow::Result;

/// Run kmerdet with the parsed CLI arguments and the raw `ArgMatches`.
///
/// The `matches` are used to determine which CLI flags were explicitly set
/// (as opposed to using their clap default values), so that config-file
/// values can be applied as fallback defaults.
pub fn run(cli: cli::Cli, matches: clap::ArgMatches) -> Result<()> {
    // Load config file early (before logging init) so we can apply runtime settings.
    let cfg = cli::load_config(&cli.global)?;

    let mut global = cli.global.clone();
    if let Some(ref cfg) = cfg {
        global.apply_config(cfg, &matches);
    }

    init_logging(global.verbose, global.quiet);

    if let Some(threads) = std::num::NonZeroUsize::new(global.threads) {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads.get())
            .build_global()
            .ok();
    }

    match cli.command {
        cli::Command::Detect(args) => cli::detect::run(args, &global, cfg.as_ref(), &matches),
        cli::Command::Filter(args) => cli::filter::run(args, &global, cfg.as_ref(), &matches),
        cli::Command::Merge(args) => cli::merge::run(args, &global),
        cli::Command::Stats(args) => cli::stats::run(args, &global),
        cli::Command::Plot(args) => cli::plot::run(args, &global),
        cli::Command::Coverage(args) => cli::coverage::run(args, &global),
        cli::Command::Run(args) => cli::run::run(args, &global, cfg.as_ref(), &matches),
        cli::Command::Benchmark(args) => cli::benchmark::run(args, &global),
        cli::Command::Pon(args) => cli::pon::run(args, &global),
    }
}

fn init_logging(verbose: u8, quiet: bool) {
    use tracing_subscriber::EnvFilter;

    let filter = if quiet {
        EnvFilter::new("error")
    } else {
        match verbose {
            0 => EnvFilter::new("warn"),
            1 => EnvFilter::new("info"),
            2 => EnvFilter::new("debug"),
            _ => EnvFilter::new("trace"),
        }
    };

    tracing_subscriber::fmt()
        .with_env_filter(filter)
        .with_writer(std::io::stderr)
        .init();
}
