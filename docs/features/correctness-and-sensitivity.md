# Feature: Correctness Definition and Configurable Sensitivity

## What This Document Defines

This document establishes what "correct" means for kmerdet, how sensitivity and
specificity are measured, and the configurable knobs that let users tune the
sensitivity/specificity tradeoff for their specific clinical or research context.

## Core Principle: kmerdet Is Not a km Clone

kmerdet reimplements the km algorithm but is explicitly designed to **exceed** km's
detection capabilities. km output is NOT ground truth. The project goal is to detect
variants that km misses — particularly at ultra-low VAF and for larger INDELs.

Implications:
- Benchmarking against km measures *concordance*, not *correctness*
- A variant detected by kmerdet but not by km is not automatically a false positive
- Ground truth must come from independent sources (see below)
- Regression testing against km is useful for catching algorithmic bugs, but sensitivity
  improvements that cause kmerdet to diverge from km are expected and desired

## Ground Truth Hierarchy

Correctness is defined relative to ground truth, not relative to km. Ground truth
sources in order of reliability:

### Tier 1: Orthogonal Validation
Variants confirmed by an independent method on the same sample:
- ddPCR (digital droplet PCR) — gold standard for VAF quantification
- Targeted deep sequencing with independent library prep
- Sanger sequencing (for presence/absence, not VAF)
- Primary tumor WGS/WES with high-confidence somatic calls

This is the ultimate arbiter of correctness. When available, this is the only
ground truth that matters.

### Tier 2: Gold-Standard Alignment Workflow
A well-validated alignment-based pipeline (e.g., BWA-MEM2 + GATK HaplotypeCaller +
Mutect2) run on the same sequencing data. This provides an independent computational
ground truth that uses fundamentally different methodology (read alignment vs k-mer
counting).

Limitations: alignment-based callers have their own FP/FN profiles, especially at
ultra-low VAF. Variants found by kmerdet but not by the alignment pipeline are not
necessarily false positives — they may represent kmerdet's sensitivity advantage.

### Tier 3: Simulated Data
Synthetic variants injected into controlled sequences with known VAF, coverage, and
error profiles. Perfect ground truth by construction.

Limitations: May not capture real-world complexity (GC bias, mapping artifacts, UMI
family structure, library-specific error profiles). Results on simulated data are
necessary but not sufficient for validation.

### Using Multiple Tiers Together
The recommended validation approach:
1. Develop and tune on Tier 3 (simulated) data — fast iteration
2. Validate on Tier 2 (alignment workflow) data — computational cross-validation
3. Final validation on Tier 1 (orthogonal) data — clinical-grade evidence

## Correctness Dimensions

### 1. Detection (Binary: Found or Not Found)

The most fundamental question: did kmerdet report a variant that matches a ground
truth entry?

**Matching criteria** (configurable):
- Chromosome and position match (with normalization: chr17 == 17)
- Alleles match (with INDEL normalization: left-alignment equivalence)
- Variant type is compatible (a Substitution call matching a Substitution truth)

**Configurable detection threshold**: A variant is "detected" if:
- It appears in the output at all (most permissive)
- Its rVAF exceeds a minimum threshold (e.g., `--min-rvaf 0.0001`)
- Its QUAL score exceeds a minimum (e.g., `--min-qual 10`)
- Its confidence tier meets a minimum (e.g., `--min-tier MEDIUM`)

This configurability is essential because different clinical contexts have different
requirements:
- **MRD monitoring**: Maximize sensitivity, tolerate higher FP rate (patients are
  already in remission; a FP triggers a repeat test, not immediate treatment)
- **Treatment selection**: Maximize specificity, moderate sensitivity (a FP could
  lead to inappropriate therapy)
- **Screening**: Balance sensitivity and specificity (large cohorts, cost of
  follow-up matters)

### 2. Classification Accuracy

Did kmerdet call the right variant type? This matters because:
- Different variant types have different clinical significance
- Misclassification (e.g., calling an insertion as a complex INDEL) may confuse
  downstream annotation

Measured as: fraction of detected variants with correct type assignment.

### 3. Localization Accuracy

Is the reported position correct? For SNVs this is straightforward. For INDELs,
representation equivalence is a real challenge — the same biological variant can
be reported at different positions depending on left/right alignment.

Measured as: position match after INDEL normalization (left-alignment).

### 4. Quantification Accuracy

How close is the reported rVAF to the true VAF?

Metrics:
- **Bias**: mean(rVAF - true_VAF) across all detected variants
- **RMSE**: sqrt(mean((rVAF - true_VAF)^2))
- **Calibration**: regression slope of rVAF vs true_VAF (ideal = 1.0)
- **Correlation**: Pearson/Spearman correlation of rVAF vs true_VAF

At ultra-low VAF (<0.01%), Poisson noise dominates and point estimates are unreliable.
What matters more is whether the confidence interval contains the true value
(coverage probability).

## Confidence: Beyond rVAF

rVAF is a necessary but insufficient measure of confidence. A variant with rVAF=0.001
could be highly confident (at 100,000x coverage with consistent k-mer support) or
essentially noise (at 1,000x coverage with a single supporting k-mer).

### Confidence Dimensions

**1. Statistical confidence (QUAL / p-value)**
Already implemented: binomial test of variant k-mer counts against error-rate null.
Phred-scaled as QUAL score. This answers: "How unlikely is this observation under
the null hypothesis of no variant?"

**2. rVAF with confidence interval**
Already implemented: bootstrap CI on rVAF. This answers: "What range of VAF values
is consistent with the observed k-mer counts?"

**3. Evidence breadth**
How many independent k-mers support the variant? A variant supported by 1 k-mer
(even at high count) is less reliable than one supported by 15 k-mers. This is
partially captured by `min_coverage` but could be formalized as:
- `n_supporting_kmers`: count of variant-specific k-mers with count > 0
- `kmer_support_fraction`: n_supporting_kmers / expected_n_supporting_kmers

**4. Positional consistency**
Are the supporting k-mer counts uniform across the variant position, or is there
a single outlier driving the call? Already partially captured by positional
uniformity (PU) in the confidence scoring spec. Formalized as:
- Coefficient of variation of variant k-mer counts
- Max/min ratio of variant k-mer counts

**5. Graph topology confidence**
Is the variant path well-connected in the k-mer graph, or is it a fragile single-edge
connection? Metrics:
- Path connectivity: are all edges in the variant path supported by multiple k-mers?
- Branch confidence: at the variant branch point, what is the ratio of variant to
  reference k-mer counts?
- No dead-end paths: does the variant path connect cleanly from source to sink?

**6. Multi-k concordance** (when multi-k mode is used)
Was the variant detected at multiple k values? A variant detected at both k=21 and
k=31 is more reliable than one detected at only one k value. Already captured in
the consensus voting tier system (Tier1 = all k values agree).

### Composite Confidence Score

Combine the above dimensions into a single confidence score. Two approaches:

**Rule-based tiers** (current, transparent):
- HIGH: QUAL >= 30 AND n_supporting_kmers >= 3 AND PU < 1.0
- MEDIUM: QUAL >= 10 AND n_supporting_kmers >= 1
- LOW: everything else

**Learned score** (future, requires validation data):
- Logistic regression or random forest trained on true/false positive labels
- Features: QUAL, rVAF, min_coverage, n_supporting_kmers, PU, graph connectivity
- Output: P(true positive | features)

## Configurable Sensitivity

### The Sensitivity Dial

Users need a single conceptual "dial" that moves between maximum sensitivity
(find everything, accept more FPs) and maximum specificity (only report high-confidence
calls). This manifests as multiple coordinated parameters:

**Detection parameters** (affect what the walker finds):
- `--count`: minimum absolute k-mer count for extension (lower = more sensitive)
- `--ratio`: minimum ratio threshold for extension (lower = more sensitive)
- `--adaptive`: automatically set count/ratio based on sample coverage

**Filtering parameters** (affect what gets reported):
- `--min-rvaf`: minimum rVAF for reporting (lower = more sensitive)
- `--min-qual`: minimum QUAL score for reporting (lower = more sensitive)
- `--min-coverage`: minimum k-mer coverage for reporting

**Preset profiles** for common use cases:
```
kmerdet detect --sensitivity ultra    # count=1, ratio=1e-6, min-qual=0
kmerdet detect --sensitivity high     # count=1, ratio=1e-4, min-qual=5
kmerdet detect --sensitivity standard # count=2, ratio=0.05, min-qual=10  (default)
kmerdet detect --sensitivity strict   # count=3, ratio=0.05, min-qual=20
```

### False Positive Management

False positives in MRD monitoring have real clinical consequences. The FP rate must
be characterizable and controllable.

**FP sources in k-mer variant detection**:
1. Sequencing errors that survive UMI deduplication
2. PCR artifacts (chimeric reads, stutter)
3. Mapping-equivalent regions (k-mers present in multiple genomic loci)
4. Low-complexity / repetitive sequences
5. Systematic errors correlated with sequence context (homopolymers, GC-extreme)

**FP suppression strategies** (layered, each independently configurable):

| Strategy | Parameter | Effect |
|----------|-----------|--------|
| Count threshold | `--count` | Eliminates single-count noise k-mers |
| Ratio threshold | `--ratio` | Requires variant signal above background |
| QUAL filtering | `--min-qual` | Statistical significance threshold |
| Graph pruning | `--no-prune` / prune params | Removes spurious graph branches |
| Panel-of-normals | `--pon` | Removes recurrent technical artifacts |
| Strand bias | (future) | Removes single-strand artifacts |
| Sequence complexity | (future) | Flags low-complexity calls |
| Bootstrap CI | `--bootstrap` | Rejects variants with CI including 0 |

**Recommended approach for tuning**:
1. Run with `--sensitivity ultra` on validation data with known truth
2. Examine FP calls: what QUAL scores do they have? What are their k-mer counts?
3. Set `--min-qual` to a threshold that removes most FPs while retaining TPs
4. Use `--report-dir` to diagnose remaining FPs and identify systematic patterns
5. Apply panel-of-normals if recurrent artifacts are observed

### Reporting Sensitivity Characteristics

Every benchmark run should report:
- LOD50: VAF at which 50% of variants are detected
- LOD80: VAF at which 80% of variants are detected
- LOD95: VAF at which 95% of variants are detected
- FPR at each LOD level
- Sensitivity × Specificity curve (not just sensitivity)

## Acceptance Criteria

- [ ] Ground truth benchmarking uses independent truth, never km output
- [ ] Detection threshold is configurable via `--min-rvaf`, `--min-qual`, `--min-coverage`
- [ ] Sensitivity presets (`--sensitivity ultra/high/standard/strict`) coordinate parameters
- [ ] Confidence score incorporates QUAL, rVAF CI, evidence breadth, and positional consistency
- [ ] Benchmark reports include LOD50/LOD80/LOD95 at specified FPR
- [ ] False positive rate is reported alongside sensitivity at each VAF bin
- [ ] `--report-dir` provides per-variant diagnostic information for FP analysis
- [ ] Documentation clearly states km is not ground truth
