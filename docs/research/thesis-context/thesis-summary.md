# Thesis Summary: K-mer Variant Detection for Liquid Biopsy ctDNA Monitoring

## Study Overview

The thesis evaluated k-mer-based variant detection (via the kam pipeline) for **minimal residual disease (MRD) monitoring** through liquid biopsy. MRD monitoring tracks residual tumor DNA fragments in patient blood after treatment, providing early warning of relapse without invasive tissue biopsies. The central challenge is that ctDNA in the MRD setting is present at extremely low variant allele frequencies (often <1%), demanding high sensitivity and rigorous artifact filtering.

---

## Study Design

### Validation Cohort

| Parameter | Value |
|-----------|-------|
| Patients | 10 |
| Sequencing method | UMI-based duplex sequencing |
| Panel type | Patient-specific targeted panels (tumor-informed) |
| Sample type | Serial blood draws at multiple timepoints per patient |
| Preprocessing | HUMID UMI deduplication before jellyfish counting |
| Sequencing depth | 2000-5000x typical |

### Tumor-Informed Target Design

For each patient, targets were designed from known somatic mutations identified in the primary tumor biopsy. Flanking reference sequence was extracted via `multiseqex --flank 35`, producing short target FASTA sequences (~150-500 bp) containing both the mutation site and genome-unique anchor k-mers at each end. Both SNVs and INDELs were included where available. This tumor-informed approach maximizes sensitivity by searching for known variants rather than performing de novo discovery.

---

## Pipeline Architecture

The kam pipeline orchestrates six sequential steps from raw FASTQ reads to filtered variant calls:

```
FASTQ (R1/R2)
  --> HUMID dedup (remove PCR duplicates via UMI families)
    --> jellyfish count (build 31-mer database, canonical, min count 2)
      --> multiseqex (extract target sequences from reference, --flank 35)
        --> refolder (organize targets into N subfolders for parallelism)
          --> kmtools chunk (run km find_mutation in parallel, merge results)
            --> kmtools filter (reference-mode filtering against known mutations)
              --> Filtered variant calls (TSV)
```

### Component Roles

| Step | Tool | Language | Function |
|------|------|----------|----------|
| 1. Deduplication | HUMID | C++ | Collapse UMI families to remove PCR duplicates; reduces artifact k-mers |
| 2. K-mer counting | jellyfish | C++ | Build lock-free hash table of all 31-mers and their counts from deduplicated reads |
| 3. Target extraction | multiseqex | Rust | Batch extraction of target regions from reference FASTA (multi-threaded samtools faidx alternative) |
| 4. File organization | refolder | Rust | Distribute target FASTA files across N subfolders for parallel processing |
| 5. Variant detection | km (via kmtools chunk) | Python | Core k-mer walking algorithm: DFS extension, graph construction, pathfinding, variant classification, NNLS quantification |
| 6. Filtering | kmtools filter | Python | Reference-mode filtering: match detected variants against known patient mutations |

### Key Pipeline Parameters

| Parameter | Setting | Rationale |
|-----------|---------|-----------|
| k-mer length (`-m`) | 31 | Balance between specificity (unique k-mers) and INDEL sensitivity; fits in 64-bit word (62 bits) |
| Canonical mode (`-C`) | Enabled | Count both strand orientations together, halving memory |
| Lower count threshold (`-L`) | 2 | Retain low-count k-mers essential for detecting low-VAF variants |
| km `--count` (n_cutoff) | 2 | Absolute minimum k-mer count for extension; very sensitive for liquid biopsy |
| km `--ratio` (cutoff) | 0.00001 | Fraction of sibling counts required for extension; very permissive for low VAF |
| Flank length | 35 bp | Flanking reference on each side of target mutation |
| Threads (kmtools chunk) | 4 | Parallel km execution across target subfolders |

---

## Detection Results by Variant Type

### SNV Performance

| Metric | Value |
|--------|-------|
| Sensitivity | 77% |
| Comparison (alignment-based) | 90-95% |
| Gap | 13-18 percentage points |

SNVs were reliably detected when the VAF exceeded the detection threshold (~0.1%). For a single nucleotide substitution, exactly `k` k-mers (31 k-mers for k=31) differ between the reference and variant paths, making the signal well-defined. False negatives primarily occurred at very low VAF (<0.1%) or in regions with high background noise.

### INDEL Performance

| INDEL category | Approximate sensitivity | Notes |
|----------------|------------------------|-------|
| Overall | 38% | Substantially lower than SNVs |
| Small INDELs (1-7 bp) | Higher (not separately quantified) | More spanning k-mers available |
| Large INDELs (>7 bp) | ~3% | Approaching k/2 = 15.5 bp limit |

### INDEL Failure Root Causes

| Factor | Mechanism | Impact severity |
|--------|-----------|-----------------|
| K-mer length constraint | Insertions >k/2 (>15 bp for k=31) have progressively fewer spanning k-mers, eventually zero | Critical for large INDELs |
| Homopolymer context | Repetitive k-mers cause extension ambiguity; multiple valid paths confound walking | High in repetitive regions |
| Low VAF + INDEL compound effect | Fewer supporting k-mers at the variant junction compounded by already-reduced k-mer coverage for INDELs | High at low VAF |
| PCR stutter | UMI dedup helps but does not fully eliminate stutter artifacts around INDELs | Moderate |
| Classification edge cases | `diff_path_without_overlap` has edge cases for complex indels near target boundaries | Moderate for complex INDELs |

---

## Performance Benchmarks

### Timing Breakdown (~6 Minutes End-to-End)

| Step | Time | Fraction of total |
|------|------|-------------------|
| HUMID deduplication | ~2 min | ~33% |
| jellyfish count | ~1.5 min | ~25% |
| multiseqex + refolder | ~10 sec | ~3% |
| kmtools chunk (km execution) | ~2 min | ~33% |
| kmtools filter | ~15 sec | ~4% |
| **Total** | **~6 min** | **100%** |

Measured on a typical 50-target panel with 4 threads.

### Comparison with Alignment-Based Pipeline

| Step | Alignment pipeline | K-mer pipeline |
|------|-------------------|----------------|
| Read preprocessing / dedup | 10-15 min (Picard MarkDuplicates) | ~2 min (HUMID) |
| Core analysis | 45-90 min (BWA-MEM + GATK Mutect2) | ~3.5 min (jellyfish + km) |
| Post-processing | 10-20 min (annotation + filtering) | ~0.5 min (kmtools filter) |
| **Total** | **1-2 hours** | **~6 minutes** |
| **Speedup** | -- | **10-20x** |

The 10-20x speedup enables same-day clinical reporting, a critical advantage for MRD monitoring.

---

## VAF Detection Threshold

### Established Threshold: ~0.1% VAF

Below 0.1% VAF, the signal-to-noise ratio was insufficient for reliable detection. At 0.1% VAF with typical sequencing depths (2000-5000x), variant k-mers had counts of approximately 2-5, matching the `--count 2` threshold parameter.

### Detection Reliability Formula

```
Detection reliability = f(VAF, sequencing_depth, k, deduplication_rate, target_design)
```

| Factor | Effect on threshold |
|--------|-------------------|
| Higher sequencing depth | Pushes threshold lower (more variant k-mers at same VAF) |
| Better deduplication (more UMI families) | Reduces noise floor |
| Target design quality | Affects anchor k-mer uniqueness and walking success |
| GC content / sequence complexity | Influences jellyfish counting accuracy |

### Parameter Sensitivity Analysis

| Configuration | Characteristics | Use case |
|---------------|----------------|----------|
| ratio=0.00001, count=2 | High sensitivity, moderate false positives | Liquid biopsy screening (thesis default) |
| ratio=0.0001, count=3 | Balanced sensitivity/specificity | Recommended for clinical use |
| ratio=0.001, count=5 | High specificity, misses low-VAF variants | Confirmatory analysis |
| ratio=0.01, count=10 | Too stringent for liquid biopsy | Not suitable for MRD monitoring |

---

## Clinical Validation: Hybrid Strategy

The thesis proposed a hybrid clinical strategy combining k-mer screening with alignment-based confirmation:

### Workflow

1. **Primary screen**: Run kmerdet on targeted panel for rapid results (~6 minutes)
2. **Confirmation**: Run standard pipeline (BWA + GATK) for positive hits (~1-2 hours)
3. **Monitoring**: Use kmerdet for longitudinal tracking of confirmed variants
4. **Escalation**: If kmerdet detects a new variant not in the original panel, trigger full analysis

### Strategy Comparison

| Aspect | K-mer (kmerdet) | Traditional Pipeline | Hybrid |
|--------|-----------------|---------------------|--------|
| Speed | ~6 minutes | 1-2 hours | Fast screening, thorough confirmation |
| SNV sensitivity | 77% | 90-95% | Best of both |
| INDEL detection | Limited (38%) | Better | Traditional fills gaps |
| Resource usage | Low (single node) | High (cluster) | Optimized allocation |
| Turnaround time | Same-day | Next-day | Same-day preliminary |
| Clinical utility | Rapid monitoring | Gold standard | Actionable fast results |

### Clinical Workflow Integration

```
Blood draw
  -> cfDNA extraction
    -> Library prep (UMI-tagged)
      -> Sequencing
        -> HUMID dedup
          -> jellyfish count (.jf)
            -> kmerdet detect (6 min)  -> Preliminary report (same-day)
            -> BWA + GATK (1-2 hours)  -> Confirmatory report
```

---

## Quantification Approach

### NNLS / Linear Regression for rVAF Estimation

The quantification model treats k-mer counts as a linear mixture of contributions from all paths (reference + variant):

```
For each k-mer i:  count[i] = SUM(contrib[i,j] * coef[j])  for all paths j
```

Where:
- `count[i]` = observed jellyfish count for k-mer i
- `contrib[i,j]` = number of times k-mer i appears in path j (0 or 1 for SNVs; can be >1 for ITDs)
- `coef[j]` = expression level (coverage) of path j

**Solution method**: Solve via least squares regression (`numpy.linalg.lstsq`), then refine with gradient descent to eliminate negative coefficients (non-negativity constraint).

**rVAF computation**:
```
rVAF[j] = coef[j] / SUM(all coef)
```

This gives the fraction of reads supporting each path, directly analogous to traditional VAF.

### Cluster Quantification

When multiple mutations overlap in k-mer space, they are grouped into clusters and quantified together. The contribution matrix accounts for shared k-mers between overlapping variant paths, preventing mutual interference between compound heterozygous mutations.

---

## False Positive and False Negative Sources

### False Positives

| Source | Mechanism |
|--------|-----------|
| CHIP mutations | Clonal hematopoiesis variants in cfDNA at low VAF, indistinguishable from tumor ctDNA by sequence alone |
| Sequencing artifacts | Systematic errors at specific sequence contexts (e.g., GGC motifs) |
| Incomplete deduplication | When UMI families are not fully collapsed, artifact k-mers persist |
| Alignment-free limitation | Without mapping context, some germline variants near target regions can be called |

### False Negatives

| Source | Mechanism |
|--------|-----------|
| Below-threshold VAF | Variants at <0.1% VAF not reliably detected |
| INDEL-specific failures | Structural complexity defeats k-mer walking |
| Target design issues | Non-unique anchor k-mers cause walking failures |
| Parameter sensitivity | `ratio` and `count` interaction: aggressive filtering improves specificity but reduces sensitivity |

---

## Summary of Key Numbers

| Metric | Value |
|--------|-------|
| Validation cohort size | 10 patients |
| SNV sensitivity | 77% |
| INDEL sensitivity (overall) | 38% |
| INDEL sensitivity (>7 bp) | ~3% |
| VAF detection threshold | ~0.1% |
| End-to-end runtime | ~6 minutes |
| Speedup vs alignment | 10-20x |
| Typical sequencing depth | 2000-5000x |
| K-mer length | 31 |
| Typical panel size | ~50 targets |

---

## References

- Thesis Chapter 4: Validation Study Design and Methods
- Thesis Chapter 5: Results and Performance Analysis
- Thesis Chapter 6: Discussion, Limitations, and Future Work
- Audoux et al. (2017) -- km: RNA-seq investigation using k-mer decomposition
- Marcais, G., & Kingsford, C. (2011) -- Jellyfish: A fast, lock-free approach for efficient parallel counting of occurrences of k-mers
- Pockrandt et al. -- HUMID: High-performance UMI deduplication
- kam repository: https://github.com/trentzz/kam
