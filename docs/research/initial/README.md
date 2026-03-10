# Initial Research

Research documents for the kmerdet project — a Rust reimplementation of k-mer variant detection from jellyfish databases, inspired by [km](https://github.com/iric-soft/km) and the [kam](https://github.com/trentzz/kam) pipeline.

## Documents

| File | Topic |
|------|-------|
| [01-km-tool-analysis.md](./01-km-tool-analysis.md) | Deep analysis of the km tool: architecture, algorithm, source code walkthrough |
| [02-jellyfish-internals.md](./02-jellyfish-internals.md) | Jellyfish k-mer counter: hash table design, .jf format, Python API, query interface |
| [03-kmer-variant-detection-algorithm.md](./03-kmer-variant-detection-algorithm.md) | The complete algorithm: k-mer walking, graph construction, pathfinding, VAF quantification |
| [04-liquid-biopsy-context.md](./04-liquid-biopsy-context.md) | Liquid biopsy challenges: cfDNA/ctDNA, low VAF detection, parameter tuning |
| [05-kam-tool-analysis.md](./05-kam-tool-analysis.md) | Analysis of the kam pipeline: Docker setup, workflow, tools, design patterns |
| [06-rust-implementation-considerations.md](./06-rust-implementation-considerations.md) | Rust implementation planning: crates, architecture, performance opportunities |
| [07-thesis-use-case-and-findings.md](./07-thesis-use-case-and-findings.md) | Thesis clinical context: 10-patient MRD validation, detection sensitivity, hybrid workflow |
| [08-tool-design-and-ux.md](./08-tool-design-and-ux.md) | CLI design: subcommands, config format, output formats, error handling, progress reporting |
| [09-kmtools-features-to-incorporate.md](./09-kmtools-features-to-incorporate.md) | kmtools reimplementation: chunk, filter, merge, plot, stats — Python to Rust |
| [10-future-improvements.md](./10-future-improvements.md) | Future goals: adaptive filtering, pure Rust jellyfish, UMI-aware counting, ML thresholding |
