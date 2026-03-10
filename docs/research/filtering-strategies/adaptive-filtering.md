# Adaptive Filtering Strategies

## 1. Current Fixed-Threshold Approach

### How Filtering Works Today

The kmerdet filter engine (`src/filter/mod.rs`) uses fixed thresholds defined in `FilterConfig`:

```rust
pub struct FilterConfig {
    pub min_coverage: u32,    // typically 2
    pub min_vaf: f64,         // typically 0.001
    pub min_expression: f64,  // typically 5.0
    pub use_alt: bool,        // coordinate vs alt-sequence matching
    pub types: Vec<String>,   // variant types to include
}
```

These values are the same for every sample, every target, and every run. The thesis work used `--count 2 --ratio 0.00001` for walking and `--min-vaf 0.001` for filtering, tuned on a 10-patient validation cohort.

### Why Fixed Thresholds Are Problematic

**Coverage variation between samples.** A sample sequenced to 10,000x median coverage has a fundamentally different noise floor than one at 1,000x. At 10,000x, an error k-mer could have count=10 by chance (0.1% error rate). At 1,000x, count=10 would represent 1% of reads -- a strong signal. A fixed count threshold of 2 is appropriate for the low-coverage sample but allows excessive noise in the high-coverage sample.

**Coverage variation between targets.** Within a single sample, different target regions have different effective coverage due to GC content, probe capture efficiency, and sequence composition. A target in a GC-rich region might have 60% of the median coverage, while a GC-neutral target has 120%. Fixed thresholds cannot account for this per-target variation.

**Library quality variation.** Samples with lower library complexity (fewer unique molecules after deduplication) have higher duplicate rates and more PCR artifacts. A fixed threshold does not account for the difference between 5,000 unique molecules and 5,000 reads from 500 molecules.

**Sequencer-specific error profiles.** Different Illumina platforms (NovaSeq, NextSeq, MiSeq) have different error rates and error patterns. NovaSeq two-color chemistry has specific error modes (e.g., poly-G artifacts from lack of signal) that NextSeq does not. Fixed thresholds tuned on one platform may not transfer.

## 2. Estimating Noise from Reference K-mer Coverage

### K-mer Spectrum Analysis

The k-mer frequency spectrum -- the histogram of k-mer counts across all k-mers -- reveals the noise structure of a sample. For a typical targeted sequencing sample:

```
Count:    1   2   3   4   5  ...  50  100  500  1000  5000
Freq:  [high peak]         [valley]    [broad peak at coverage]
       ↑ error k-mers                  ↑ genuine k-mers
```

The low-count peak consists of k-mers arising from sequencing errors. Each error creates k novel k-mers (one for each position where the error base participates in a k-mer). The high-count peak corresponds to genuine genomic k-mers at the sequencing depth.

The valley between these peaks is the natural threshold separating noise from signal. For genome assembly, tools like KmerGenie and Jellyfish's histogram output identify this valley automatically.

### Reference K-mer Distribution

For targeted variant detection, we can focus on the reference k-mers (those derived from the target sequences). These k-mers should be present at approximately the sequencing depth. Their count distribution provides a per-sample coverage estimate:

```rust
fn coverage_stats(ref_kmer_counts: &[u64]) -> CoverageStats {
    let sorted = ref_kmer_counts.to_vec();
    sorted.sort();
    CoverageStats {
        median: sorted[sorted.len() / 2],
        mean: sorted.iter().sum::<u64>() as f64 / sorted.len() as f64,
        p5: sorted[sorted.len() * 5 / 100],
        p95: sorted[sorted.len() * 95 / 100],
        cv: std_dev(&sorted) / mean(&sorted),
    }
}
```

### Error K-mer Distribution

For each reference k-mer, query all four extensions. The three non-reference extensions are predominantly errors. Their count distribution characterizes the noise floor:

```
error_counts = []
for ref_kmer in reference_kmers:
    children = extend_forward(db, ref_kmer)
    for child in children:
        if child.sequence not in reference_kmers:
            error_counts.append(child.count)
```

The distribution of `error_counts` can be modeled as a mixture:
- **Gamma distribution** for the majority (sequencing errors): characterized by shape and rate parameters
- **Outlier component** for genuine low-VAF variants: a few high-count entries

Following the approach of Kelley et al. (Quake) and Guo et al. (statistically-solid k-mers), fit a Gamma distribution to the error component and set the threshold at the 99.9th percentile:

```
threshold = Gamma_quantile(0.999; shape_hat, rate_hat)
```

This automatically adapts to the sample's error rate and coverage depth.

### Practical Implementation

```rust
pub struct AdaptiveFilterConfig {
    /// Target false positive rate for noise k-mers exceeding threshold
    pub noise_fpr: f64,  // default: 0.001 (0.1%)
    /// Minimum threshold regardless of model (safety floor)
    pub min_threshold: u32,  // default: 2
    /// Whether to compute per-target thresholds
    pub per_target: bool,  // default: true
}

pub fn compute_adaptive_threshold(
    ref_kmer_counts: &[u64],
    error_counts: &[u64],
    config: &AdaptiveFilterConfig,
) -> u32 {
    let gamma_params = fit_gamma(error_counts);
    let model_threshold = gamma_quantile(1.0 - config.noise_fpr, gamma_params);
    std::cmp::max(config.min_threshold, model_threshold as u32)
}
```

## 3. Adjusting Thresholds Based on Sequencing Depth

### Depth-Proportional Count Threshold

The most direct adaptation: set the count threshold proportional to coverage and the minimum detectable VAF:

```
adaptive_count = max(2, floor(median_coverage * min_detectable_vaf))
```

| Median Coverage | Min Detectable VAF | Adaptive Count Threshold |
|----------------|-------------------|------------------------|
| 500x | 0.1% | max(2, 0.5) = 2 |
| 1,000x | 0.1% | max(2, 1) = 2 |
| 5,000x | 0.1% | max(2, 5) = 5 |
| 10,000x | 0.1% | max(2, 10) = 10 |
| 50,000x | 0.1% | max(2, 50) = 50 |

At 50,000x (achievable with modern targeted panels), a count threshold of 2 would pass enormous numbers of error k-mers. An adaptive threshold of 50 dramatically reduces false positives while still detecting 0.1% VAF variants.

### Depth-Adjusted Ratio Threshold

Similarly, the ratio threshold can be adjusted. At high coverage, sequencing errors have higher absolute counts but similar relative frequencies. The ratio threshold should reflect the expected error rate rather than being fixed:

```
adaptive_ratio = max(error_rate / 3, base_ratio)
```

Where `error_rate` is estimated from the data (Section 2) or set to a platform default. For Illumina post-dedup: ~0.001 (0.1%), yielding:

```
adaptive_ratio = max(0.001/3, 0.00001) = max(0.00033, 0.00001) = 0.00033
```

This is 33x more stringent than the current liquid biopsy ratio of 0.00001, which is appropriate because the current ratio was chosen for very low coverage situations where false negatives are the primary concern.

### Recommended Depth Tiers

For clinical implementation, define coverage tiers with validated parameter sets. This is the recommendation from standardization efforts for NGS mutation analysis (Jennings et al., 2017):

| Tier | Median Coverage | Count | Ratio | Min VAF | Use Case |
|------|----------------|-------|-------|---------|----------|
| Ultra-low | <500x | 2 | 0.00001 | 0.005 | Exploratory, high sensitivity |
| Standard | 500-2000x | 2 | 0.0001 | 0.002 | Routine MRD monitoring |
| High-depth | 2000-10000x | 5 | 0.0003 | 0.001 | Clinical ctDNA panels |
| Ultra-deep | >10000x | 10 | 0.001 | 0.0005 | Ultra-sensitive MRD |

kmerdet should auto-detect the appropriate tier from median coverage and report which tier was used in the output.

## 4. Per-Target Thresholds Based on Genomic Context

### GC Content Bias

GC content strongly influences sequencing coverage. Extreme GC content (>70% or <30%) typically results in lower coverage, requiring lower thresholds to maintain sensitivity.

From recent work (Zheng et al., 2024; GCfix), GC bias in cfDNA data follows a characteristic pattern where fragment-length-specific GC correction is needed because cfDNA fragments have a narrow size distribution (~167 bp) and GC effects vary by fragment length.

**Implementation:** For each target, compute the GC content and adjust the expected coverage:

```rust
fn gc_coverage_adjustment(gc_fraction: f64, gc_model: &GcBiasModel) -> f64 {
    // GC bias model: fitted from reference k-mer coverage vs GC content
    // Returns a multiplier (1.0 = no adjustment, 0.5 = expect half coverage)
    gc_model.expected_ratio(gc_fraction)
}

fn adjusted_threshold(
    base_threshold: u32,
    gc_fraction: f64,
    gc_model: &GcBiasModel,
) -> u32 {
    let adjustment = gc_coverage_adjustment(gc_fraction, gc_model);
    (base_threshold as f64 * adjustment).round() as u32
}
```

The GC bias model is fitted per-sample from the reference k-mer counts grouped by GC content of their surrounding sequence.

### Repetitiveness

Repetitive regions generate more background k-mer noise because:
1. K-mers from repetitive regions match at multiple genomic locations
2. Sequencing errors in repeats create k-mers that also match at other repeat copies
3. Homopolymer runs cause polymerase slippage artifacts (INDEL noise)

**Quantifying repetitiveness:** For each target, compute the linguistic complexity (ratio of distinct k-mers to total k-mers in the reference sequence). Low complexity indicates repetitiveness.

```rust
fn linguistic_complexity(sequence: &str, k: usize) -> f64 {
    let total_kmers = sequence.len() - k + 1;
    let distinct_kmers: HashSet<&str> = (0..total_kmers)
        .map(|i| &sequence[i..i+k])
        .collect();
    distinct_kmers.len() as f64 / total_kmers as f64
}
```

Targets with complexity < 0.8 should have elevated thresholds (e.g., 2x the base threshold).

### Homopolymer Length

Homopolymer runs are the single largest source of INDEL artifacts in Illumina sequencing. The error rate increases exponentially with homopolymer length:

| Homopolymer Length | Approximate Error Rate (Illumina) |
|-------------------|----------------------------------|
| 1-3 bp | 0.01% (baseline) |
| 4-5 bp | 0.1% |
| 6-7 bp | 1% |
| 8+ bp | 5-10% |

**Implementation:** For each target, find the longest homopolymer run and the total number of homopolymer positions:

```rust
fn max_homopolymer(sequence: &str) -> usize {
    let bytes = sequence.as_bytes();
    let mut max_run = 1;
    let mut current_run = 1;
    for i in 1..bytes.len() {
        if bytes[i] == bytes[i-1] {
            current_run += 1;
            max_run = max_run.max(current_run);
        } else {
            current_run = 1;
        }
    }
    max_run
}
```

For targets with homopolymer runs >= 6 bp:
- Increase INDEL-specific count threshold by 5x
- Flag INDEL calls near homopolymers with a warning annotation
- Apply stricter rVAF thresholds for INDELs (e.g., min_vaf * 5)

### Target Annotation Precomputation

Compute all context features once when loading the target catalog, and store as annotations:

```rust
pub struct TargetAnnotation {
    pub gc_content: f64,
    pub linguistic_complexity: f64,
    pub max_homopolymer: usize,
    pub total_homopolymer_bases: usize,
    pub mappability_score: f64,  // if available from external annotation
    pub adjusted_count_threshold: u32,
    pub adjusted_ratio_threshold: f64,
}
```

These annotations can be serialized to a sidecar file (e.g., `panel_annotations.json`) so they do not need to be recomputed for every sample.

## 5. Panel-of-Normals Approach

### Concept

A Panel of Normals (PoN) is a database of variants found in healthy control samples processed through the same pipeline. Any variant that recurs across multiple normals is likely an artifact (sequencing error, reference genome error, common germline variant, or CHIP) rather than a somatic mutation.

This is standard practice in somatic variant calling. GATK Mutect2's PoN is the most widely used implementation. The Sentieon ctDNA pipeline also incorporates PoN filtering as a core module.

### Building the PoN

1. **Select control samples:** N >= 20 healthy donors, processed with the same library prep, sequencing platform, and target panel as the clinical samples. More donors improve artifact detection sensitivity.

2. **Process each control through kmerdet:**
   ```bash
   for sample in normals/*.jf; do
       kmerdet detect --jf $sample --targets panel/ --output normals_results/
   done
   ```

3. **Aggregate variants across normals:**
   ```rust
   pub struct PonEntry {
       pub target: String,
       pub variant_name: String,
       pub variant_type: String,
       pub ref_allele: String,
       pub alt_allele: String,
       pub frequency: f64,  // fraction of normals with this variant
       pub max_vaf: f64,    // highest VAF seen in any normal
       pub mean_vaf: f64,   // average VAF when present
   }
   ```

4. **Store as a compact database:** Hash map from `(target, normalized_variant)` to `PonEntry`. Serialize to a binary file (e.g., MessagePack or bincode) for fast loading.

### Filtering with the PoN

At detection time, each variant call is checked against the PoN:

```rust
fn pon_filter(call: &VariantCall, pon: &HashMap<PonKey, PonEntry>, max_freq: f64) -> bool {
    let key = PonKey {
        target: call.target.clone(),
        variant_name: normalize_variant(&call.variant_name),
    };
    match pon.get(&key) {
        Some(entry) if entry.frequency > max_freq => false,  // filter
        _ => true,  // pass
    }
}
```

Default `max_freq = 0.05` (variant found in >5% of normals is filtered). This can be tuned:
- Stringent (0.01): filters anything found in even 1 out of 100 normals
- Permissive (0.10): only filters variants found in >10% of normals

### What the PoN Catches

**Sequencer-specific systematic errors.** Certain sequence contexts (e.g., GGC motifs on NovaSeq, poly-G runs on NextSeq) produce recurrent artifactual variants. These appear in every sample and are effectively removed by the PoN.

**Reference genome errors.** Positions where the GRCh38 reference has a rare allele. Every sample will show a "variant" at these positions. The PoN eliminates these without requiring manual curation.

**Common germline variants.** Common SNPs (MAF > 5%) that are present in many normal individuals. While germline filtering is typically done by matching to a buffy coat sample, the PoN provides an additional safety net.

**Clonal hematopoiesis of indeterminate potential (CHIP).** CHIP mutations are age-related somatic mutations in hematopoietic stem cells that are present in cfDNA from blood. As documented by the BLOODPAC consortium's 2024 working group, CHIP variants in genes like DNMT3A, TET2, ASXL1, and TP53 frequently overlap with solid tumor mutations, creating a major confounding factor. A PoN built from age-matched controls captures the most common CHIP variants. Importantly, 50.1% of patients present at least one CHIP mutation among reportable clinical genes (per a 2024 characterization study of 16,812 patients), making this filtering step essential.

### Limitations of the PoN

- **Sample size:** A PoN of 20 samples only captures artifacts present in >= 1/20 = 5% of samples. Rare artifacts are missed. Larger PoN (N=100+) provide better coverage.
- **Platform specificity:** The PoN must be rebuilt when changing sequencing platforms, library prep kits, or target panels.
- **Age matching for CHIP:** CHIP prevalence increases with age. A PoN from young donors may not capture CHIP variants common in elderly patients.
- **Does not catch sample-specific artifacts:** Per-sample noise from library prep stochasticity is not captured by the PoN.

## 6. Machine Learning-Based Filtering

### Feature Engineering

ML-based filtering uses multiple features of each variant call to predict whether it is a true variant or an artifact. Features available from kmerdet output:

| Feature | Source | Rationale |
|---------|--------|-----------|
| rVAF | PathQuant | Primary signal strength indicator |
| min_coverage | Path k-mer counts | Bottleneck support for the variant path |
| expression | PathQuant NNLS | Absolute coverage of the variant |
| ref_expression | PathQuant NNLS | Coverage of the reference allele |
| variant_type | Classification | Different types have different artifact profiles |
| gc_content | Target annotation | GC bias affects coverage and error rates |
| homopolymer_length | Target annotation | Homopolymer proximity indicates INDEL artifact risk |
| linguistic_complexity | Target annotation | Repetitive context increases false positives |
| path_count | Graph | Number of alternative paths (complexity indicator) |
| count_ratio | Extension | Variant vs. reference k-mer count ratio |
| coverage_uniformity | Reference k-mers | CV of coverage across the target |
| pon_frequency | PoN | How often this variant appears in normals |

### Model Selection

A 2025 study published in Scientific Reports demonstrated that Random Forest models can effectively predict high-confidence somatic ctDNA variants in both low and high depth cfDNA NGS data. For kmerdet, gradient-boosted trees (XGBoost/LightGBM-style) offer several advantages:

- **Interpretable:** Feature importance rankings explain which factors drive classification
- **Handles mixed types:** Numeric and categorical features without preprocessing
- **Robust to outliers:** Tree-based methods are not affected by count magnitude
- **Small model size:** A trained model with 100-500 trees is <1 MB, easily shipped with the binary

In Rust, the `linfa` crate provides random forest and gradient boosting implementations. Alternatively, train in Python (scikit-learn or XGBoost) and export the model as a serialized decision tree structure that Rust can evaluate.

### Training Data Challenge

The 10-patient thesis validation set is too small for reliable ML training. Options:

1. **Simulated data:** Generate synthetic k-mer databases with known variants at various VAFs, GC contents, and error rates. Process through kmerdet to create labeled training examples. This provides unlimited data but may not capture all real-world artifact patterns.

2. **Transfer learning from alignment-based callers:** Use GATK Mutect2 or VarScan2 calls as ground truth labels. Variants confirmed by alignment-based callers are labeled "true"; kmerdet-only calls are labeled "candidate artifact." This leverages existing validated tools as an oracle.

3. **Multi-site collaboration:** Pool data from multiple clinical sites running the same panel. Even 100 patients with ~50 targets each provide ~5,000 variant calls for training.

4. **Active learning:** Deploy the rule-based filter initially, manually review borderline cases, and incrementally build a labeled dataset.

### Integration

The ML filter is applied as a post-processing step after detection and before output:

```rust
pub fn ml_filter(
    calls: &[VariantCall],
    model: &TrainedModel,
    threshold: f64,  // e.g., 0.5 for balanced, 0.3 for high sensitivity
) -> Vec<FilteredCall> {
    calls.iter().map(|call| {
        let features = extract_features(call);
        let score = model.predict_proba(&features);
        FilteredCall {
            call: call.clone(),
            ml_score: score,
            ml_pass: score >= threshold,
        }
    }).collect()
}
```

Each call is annotated with a confidence score (0.0-1.0). Users can filter by score threshold, with lower thresholds favoring sensitivity and higher thresholds favoring specificity.

**The ML filter is always optional.** The deterministic rule-based filter remains the default. The ML filter is activated via `--ml-filter` with an optional `--ml-threshold` parameter.

## 7. Composite Filtering Pipeline

### Recommended Filter Order

Apply filters in sequence, with each stage removing progressively more subtle artifacts:

```
Stage 1: Hard thresholds (count >= adaptive_min, rVAF >= adaptive_min_vaf)
    ↓
Stage 2: Type filtering (exclude "Reference" type, apply per-type min_vaf)
    ↓
Stage 3: Context filtering (homopolymer adjustment, GC adjustment)
    ↓
Stage 4: Panel-of-normals (remove recurrent artifacts)
    ↓
Stage 5: ML filter (optional, score-based ranking)
    ↓
Stage 6: Tumor-informed matching (reference mode or alt-sequence mode)
```

Each stage is independently configurable and can be disabled. The output reports which stages each variant passed, enabling retrospective analysis of filter performance.

### Filter Annotations in Output

Add filter columns to the output format:

| Column | Example | Description |
|--------|---------|-------------|
| FILTER | PASS | Overall filter status |
| FILTER_STAGES | 1,2,3,4,6 | Which stages this variant passed |
| ADAPTIVE_THRESHOLD | 5 | Count threshold used for this sample |
| PON_FREQ | 0.00 | Frequency in panel of normals |
| ML_SCORE | 0.92 | ML confidence score (if enabled) |
| CONTEXT_FLAGS | HOMOPOLYMER_5 | Genomic context warnings |

This transparency allows clinicians and bioinformaticians to understand why each variant was called or filtered, which is critical for clinical reporting.

## References

- Kelley et al. (2010) -- Quake: quality-aware detection and correction of sequencing errors (Gamma distribution threshold model)
- Guo et al. (2024) -- Mining statistically-solid k-mers for accurate NGS error correction
- Jennings et al. (2017) -- Standardization of sequencing coverage depth in NGS for cancer diagnostics
- Zheng et al. (2024) -- GCfix: fragment length-specific GC bias correction for cfDNA
- Cibulskis et al. (2013) -- Mutect2: sensitive detection of somatic mutations (Panel of Normals)
- BLOODPAC CH/CHIP Working Group (2024) -- Lexicon for Clonal Hematopoiesis in Liquid Biopsy
- Steele et al. (2025) -- Predicting high confidence ctDNA somatic variants with ensemble machine learning models (Scientific Reports)
- van Allen et al. (2014) -- Clinical interpretation of cancer genomes (artifact filtering strategies)
- panelGC (2024) -- Quantifying and monitoring GC biases in hybridization capture panel sequencing (Briefings in Bioinformatics)
