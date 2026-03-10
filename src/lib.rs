pub mod cli;
pub mod config;
pub mod filter;
pub mod graph;
pub mod jellyfish;
pub mod kmer;
pub mod output;
pub mod sequence;
pub mod variant;
pub mod walker;

use anyhow::Result;

/// Run kmerdet with the parsed CLI arguments.
pub fn run(cli: cli::Cli) -> Result<()> {
    init_logging(cli.global.verbose, cli.global.quiet);

    if let Some(threads) = std::num::NonZeroUsize::new(cli.global.threads) {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads.get())
            .build_global()
            .ok();
    }

    match cli.command {
        cli::Command::Detect(args) => cli::detect::run(args, &cli.global),
        cli::Command::Filter(args) => cli::filter::run(args, &cli.global),
        cli::Command::Merge(args) => cli::merge::run(args, &cli.global),
        cli::Command::Stats(args) => cli::stats::run(args, &cli.global),
        cli::Command::Plot(args) => cli::plot::run(args, &cli.global),
        cli::Command::Coverage(args) => cli::coverage::run(args, &cli.global),
        cli::Command::Run(args) => cli::run::run(args, &cli.global),
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
