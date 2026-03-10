# kam Tool Analysis (trentzz/kam)

## Overview

kam (K-mer Analysis Modules) is a **collection of tools** bundled as a Docker-based pipeline for k-mer-based variant detection. It wraps the km tool and several custom utilities into a cohesive workflow.

- **Repository**: https://github.com/trentzz/kam
- **License**: MIT
- **Languages**: Dockerfile (57.1%), Shell (42.9%)
- **Docker Image**: `trentzz/kam`

## Architecture

kam is primarily a **pipeline orchestrator**, not a standalone analysis tool. It provides:

1. A Docker container with all dependencies pre-installed
2. Example workflows (bash and Nextflow)
3. Custom tools to extend km's functionality

### Docker Container Contents

The Dockerfile installs:

| Component | Purpose |
|-----------|---------|
| Python 3.12 | Runtime for km and Python tools |
| Rust (stable) | Runtime for Rust-based tools |
| Jellyfish 2.2.6 | K-mer counting (built from source with Python bindings) |
| Java 17 (OpenJDK) | Required for Nextflow |
| samtools | BAM/FASTA manipulation |
| HUMID | Read deduplication |
| Nextflow | Workflow engine |

### Tools Included

| Tool | Language | Purpose |
|------|----------|---------|
| **km-walk** (km) | Python | Core k-mer variant detection (from iric-soft) |
| **multiseqex** | Rust | Batch sequence extraction from reference FASTA (multi-core samtools faidx alternative) |
| **kmtools** | Python | km extensions: multithreading, filtering, stats, plotting |
| **refolder** | Rust | Reorganize files into subfolders for parallel processing |
| **vcf2pandas** | Python | Convert VCF files to pandas DataFrames |
| **vcf2xlsx** | Python | Convert VCF files to Excel files |

## Workflow Pipeline

### Example Workflow (example-workflow.sh)

```
Step 1: HUMID deduplication
    humid -d output/humid -m 1 -e R1.fastq.gz R2.fastq.gz
    → Removes PCR duplicates from paired-end reads

Step 2: Jellyfish k-mer counting
    jellyfish count -m 31 -o output.jf -s 100M -t 4 -C -L 2 <(zcat R1) <(zcat R2)
    → Builds k-mer database from deduplicated reads

Step 3: Target sequence extraction
    multiseqex reference.fa --table targets.tsv --flank 35 --output-dir targets/
    → Extracts target regions from reference genome

Step 4: Organize targets for parallel processing
    refolder targets/ --subfolders N --matching "*.fa"
    → Distributes target FASTA files across N subfolders

Step 5: Parallel km execution
    kmtools chunk --threads 4 \
        --km-find-mutation-options "--count 2 --ratio 0.00001" \
        --km-target-directory targets/ \
        --km-jellyfish-file output.jf \
        --output-dir km_output/ \
        --merge --merge-output-file km_results.txt
    → Runs km find_mutation in parallel across target subfolders
    → Merges results into single output file

Step 6: Filter results
    kmtools filter --reference reference-info.tsv \
        --km-output km_results.txt \
        --output km_filtered.tsv
    → Filters km output against reference annotations
```

### Key Parameters Used

- **jellyfish**: `-m 31` (31-mers), `-C` (canonical), `-L 2` (min count 2)
- **multiseqex**: `--flank 35` (35bp flanking on each side of target)
- **km**: `--count 2 --ratio 0.00001` (very sensitive — tuned for low VAF/liquid biopsy)

### Input Files

1. **Paired FASTQ files**: Sequencing reads (R1/R2)
2. **Reference FASTA**: Reference genome
3. **targets.tsv**: Tab-separated file with chromosome and position columns
4. **reference-information.tsv**: Reference annotations for filtering

### targets.tsv Format
```
CHROM, POS
chr1, 100
chr2, 200
```

## Key Design Patterns

1. **Docker-first**: Everything packaged in a single container for reproducibility
2. **Pipeline composition**: Individual tools composed via shell scripts or Nextflow
3. **Parallelism via file organization**: refolder + kmtools chunk enables parallel km execution
4. **Rust for performance-critical tools**: multiseqex (sequence extraction) and refolder are Rust
5. **Python for analysis logic**: km, kmtools, vcf2pandas, vcf2xlsx are Python
6. **Sensitive defaults**: `--count 2 --ratio 0.00001` tuned for liquid biopsy detection

## Relevance to New Rust Tool

The kam pipeline reveals the full workflow context:

1. **Pre-processing**: Deduplication is critical (HUMID)
2. **Target design**: Automated extraction from coordinates (multiseqex)
3. **Parallel execution**: The main bottleneck is km's single-threaded Python execution
4. **Post-processing**: Filtering and format conversion needed

A Rust reimplementation could:
- Natively handle parallelism (no need for refolder + kmtools chunk)
- Directly read jellyfish databases (no Python binding overhead)
- Integrate target extraction (replacing multiseqex as a subcommand)
- Be significantly faster for the core k-mer walking algorithm
- Produce output in multiple formats (TSV, VCF, Excel)
