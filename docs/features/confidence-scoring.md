# Feature: Confidence Scoring -- Variant Quality Metrics

## What It Does

Confidence scoring assigns statistical quality metrics to each variant call,
transforming binary detect/no-detect output into calibrated probabilistic
assessments. Each variant receives a Phred-scaled QUAL score, positional
uniformity measure, rVAF confidence interval, and (when strand-aware counting
is available) a strand bias score. A composite quality score integrates all
dimensions into a single ranked confidence value.

## Why It Matters

### No Uncertainty Quantification Today

The current k-mer variant detection pipeline produces three metrics per variant:
rVAF (relative variant allele frequency), expression (absolute path coefficient),
and min_coverage (minimum k-mer count along the variant path). None of these
metrics answers the fundamental clinical question: "How confident should I be
that this variant is real?"

A variant with rVAF=0.2% and min_coverage=3 could be a true low-VAF variant
with strong evidence or a noise-induced false positive. Without a statistical
test against the null hypothesis (no variant present), there is no principled
way to distinguish these cases. Clinical decision-making requires calibrated
confidence -- ideally expressed as a probability or p-value.

### VCF Compatibility

Standard bioinformatics tools expect a QUAL score in Phred scale. Downstream
filtering in bcftools, variant annotation in VEP/Annovar, and visualization in
IGV all use QUAL for ranking and filtering. Without a QUAL score, kmerdet
output cannot participate in standard bioinformatics workflows that combine
calls from multiple tools.

### Strand Bias as Artifact Detection

Research by Guo et al. (2012) demonstrated that strand bias filtering removes
70-90% of false positive variant calls in Illumina data. Sequencing artifacts
(oxidative damage, deamination, polymerase errors) frequently affect only one
strand. A true somatic variant is equally likely to appear on both strands. The
absence of strand information in the current pipeline means this powerful filter
is unavailable.

### MRD Monitoring Requires Calibrated Uncertainty

In MRD (minimal residual disease) monitoring, clinicians track rVAF over time.
A rising rVAF trend may indicate relapse, but without confidence intervals, it
is impossible to distinguish a true trend from measurement noise. Bootstrap
confidence intervals on rVAF enable trend analysis with appropriate uncertainty
quantification.

## Research Backing

### Statistical Frameworks (confidence-metrics.md)

The research document evaluated five statistical frameworks for k-mer enrichment
testing:

**Binomial model**: Given total coverage D and error rate e_k, compute the
probability of observing count c or higher under the null hypothesis (error
only). Simple, well-understood, computationally trivial. Used by SNVer (Wei et
al., 2011) for alignment-based calling.

**Poisson model**: More appropriate for low-count k-mers (count < 20). The
likelihood ratio test between null (lambda = D * e_k) and alternative
(lambda = D * f) can be combined across k-mers using the sum of individual
log-likelihood ratios.

**Negative binomial model**: Extends Poisson to account for overdispersion
(variance > mean), which is consistently observed in sequencing count data
(Love et al., DESeq2, 2014). The overdispersion parameter alpha is estimated
from reference k-mer count variance. Typical values range from alpha ~ 0.01
(well-captured regions) to alpha ~ 0.5-2.0 (GC-biased regions).

**Fisher's exact test**: Constructs a 2x2 contingency table comparing variant
and reference k-mer counts. Also applicable for strand bias (forward vs reverse
counts for variant vs reference). Non-parametric with no distributional
assumptions.

**Bayesian posterior**: Integrates prior knowledge (tumor-informed: variant is
known from biopsy) with k-mer count evidence. Particularly powerful for the
tumor-informed use case where P(variant) = 0.5 is a reasonable prior. The
BATCAVE algorithm (Valeri et al., 2020) and BayVarC (2024) demonstrated this
approach for somatic variant calling and liquid biopsy respectively.

### Alignment-Based Quality Metric Comparison (confidence-metrics.md)

Standard VCF annotations from alignment-based callers include QUAL, GQ, DP, AD,
AF, MQ, FS (Fisher Strand), SOR, ReadPosRankSum, and BaseQRankSum. The k-mer
approach can produce equivalents for QUAL, DP, AD, AF, and FS (with strand-aware
counting). Mapping quality (MQ) and base quality rank sum are inherently
alignment-based and have no k-mer equivalent, but k-mer multiplicity in the
reference genome can serve as a proxy for mapping ambiguity.

### K-merald Approach (Denti et al., 2023)

K-merald constructs k-mer-based sequencing error profiles from the data itself
and genotypes variants based on the empirical error distribution. This achieves
higher genotype quality than edit-distance-based methods, particularly for
indels. The principle -- using the observed error distribution rather than a
fixed rate -- aligns with the adaptive error rate estimation in confidence
scoring.

### Merfin Approach (Formenti et al., 2022)

Merfin evaluates each variant based on expected k-mer multiplicity, checking
whether introducing the variant makes the k-mer spectrum more consistent with
expected copy number. This concept can be adapted: for each candidate variant,
score how well observed k-mer counts match the expected counts under the variant
hypothesis versus the reference hypothesis.

## Design Considerations

### QUAL Score (Phase 1)

The primary quality metric. Computed as a Phred-scaled p-value from a
statistical test against the null hypothesis of no variant.

**Computation per k-mer position**:

For each variant-specific k-mer i with count c_i and local coverage D_i:

```
p_i = 1 - BinomialCDF(c_i - 1, D_i, e_k)
```

Where e_k is the estimated per-k-mer error rate (approximately base_error_rate
/ 3 for a specific substitution).

**Combined p-value across k-mers** (Fisher's method with correlation correction):

```
X = -2 * sum(log(p_i))  for i = 1..m
```

Under independence, X ~ chi-squared(2m). Since adjacent k-mers share reads,
apply correlation correction using effective independent k-mers:

```
m_eff = m * k / read_length  (approximately 0.2m for k=31, read_length=150)
combined_p = 1 - ChiSquaredCDF(X * m_eff / m, df = 2 * m_eff)
```

**Phred scaling**:

```
QUAL = -10 * log10(combined_p)
```

| combined_p | QUAL | Interpretation |
|------------|------|----------------|
| 0.1 | 10 | Low confidence |
| 0.01 | 20 | Marginal |
| 0.001 | 30 | Good (standard filter threshold) |
| 1e-4 | 40 | High confidence |
| 1e-6 | 60 | Very high confidence |
| 1e-10 | 100 | Extremely high confidence |

**Conservative variant**: Rather than Fisher's combined p-value, report the
worst (highest) per-k-mer p-value as the overall variant p-value. This is more
conservative and avoids the independence assumption entirely.

### Strand Bias (Phase 2)

Requires strand-aware k-mer counting (non-canonical mode or separate
forward/reverse databases). When available:

**Fisher's exact test for strand bias**:

```
                Forward    Reverse
Reference:      ref_fwd    ref_rev
Variant:        var_fwd    var_rev
```

The Fisher's exact p-value tests disproportionate strand representation.
Report as Phred-scaled FS (FisherStrand):

```
FS = -10 * log10(fisher_p_value)
```

GATK uses FS > 60 as a strand bias filter. Values > 200 indicate extreme
single-strand artifacts.

**Strand bias ratio** (simpler alternative):

```
SB = |var_fwd - var_rev| / (var_fwd + var_rev)
```

SB close to 0 = balanced (likely true). SB close to 1 = single-strand (likely
artifact).

### Positional Uniformity (Phase 1)

For an SNV, k variant-specific k-mers span the variant position. If the variant
is real, their counts should be approximately uniform (Poisson variation only).
Non-uniform counts suggest positional bias.

**Coefficient of variation**:

```
PU = std(variant_kmer_counts) / mean(variant_kmer_counts)
```

Expected for true variants: PU ~ 0.3-0.5 (Poisson CV = 1/sqrt(mean)).
Values > 1.0 suggest positional bias.

**Chi-squared goodness-of-fit** (more rigorous):

```
chi2 = sum_i (observed_i - expected)^2 / expected
p_positional = 1 - ChiSquaredCDF(chi2, df = k - 1)
```

### rVAF Confidence Interval (Phase 2)

The NNLS decomposition provides a point estimate of rVAF with no uncertainty.
Bootstrap resampling provides confidence intervals:

1. Resample k-mer counts from a Poisson distribution centered on observed counts
2. Solve NNLS on each bootstrap sample
3. Compute rVAF for each replicate
4. Report 2.5th and 97.5th percentiles as 95% CI

For B=1000 bootstrap replicates with k=31, this requires 1000 NNLS solutions.
NNLS for a 31-variable system is computationally trivial (microseconds per
solve), so 1000 replicates add negligible overhead.

The CI width indicates measurement precision. A narrow CI (e.g., 0.08%-0.12%)
means the rVAF estimate is stable. A wide CI (e.g., 0.01%-0.5%) means the
estimate is unreliable, typically at very low coverage or very low VAF.

### Composite Quality Score (Phase 3)

Integrate all quality dimensions into a single score using a trained logistic
regression model:

```
logit(CQS) = w1 * log(p_value)
           + w2 * strand_bias
           + w3 * positional_cv
           + w4 * log(min_coverage)
           + w5 * log(rVAF)
           + w6 * sequence_complexity
           + w7 * homopolymer_proximity
```

Weights are trained from the validation dataset (true positive vs false positive
labels from alignment-based confirmation). The composite score provides a single
ranking metric for clinical prioritization.

### Confidence Tiers

Map QUAL scores to human-readable tiers for clinical reporting:

| Tier | QUAL Range | Label | Action |
|------|-----------|-------|--------|
| HIGH | >= 30 | High confidence | Report without review |
| MEDIUM | 10-30 | Moderate confidence | Report with clinical review |
| LOW | < 10 | Low confidence | Suppress or flag for confirmatory testing |

### VCF INFO/FORMAT Fields

| Field | Type | Description |
|-------|------|-------------|
| QUAL | Float | Phred-scaled variant quality (VCF column 6) |
| FS | Float | Fisher strand bias (Phred-scaled) |
| SB | Float | Strand bias ratio (0-1) |
| PU | Float | Positional uniformity (CV of variant k-mer counts) |
| RVAF_CI | String | 95% CI on rVAF (e.g., "0.08,0.12") |
| CQS | Float | Composite quality score (0-1) |
| TIER | String | Confidence tier (HIGH, MEDIUM, LOW) |

### Error Rate Estimation

The per-k-mer error rate e_k is estimated locally rather than assumed:

**Empirical approach**: Query non-reference extensions at reference k-mer
positions in the target region. The mean count of non-reference extensions
divided by total coverage gives the local error rate.

**Model-based approach**: Use the Illumina error model where the k-mer error
rate is approximately:

```
e_k ~ 1 - (1 - e_base)^k ~ k * e_base
```

For k=31 and Q30 bases (e_base = 0.001), e_k ~ 0.031. The probability of a
specific alternative k-mer is e_base / 3 ~ 0.00033.

The empirical approach is preferred because it captures platform-specific,
sample-specific, and context-specific error rates.

## Acceptance Criteria

### QUAL Score

1. Every variant call includes a QUAL score in Phred scale.

2. QUAL scores are computed from a binomial or Poisson test against the null
   hypothesis of error-only k-mer counts.

3. QUAL scores correlate with true positive rate: variants with QUAL >= 30
   have a higher validation rate than those with QUAL 10-30, which in turn
   are higher than QUAL < 10.

4. QUAL scores are comparable across samples with different coverage depths
   (a variant at 0.5% VAF at 5000x and 0.5% VAF at 2000x should have QUAL
   values reflecting the different statistical power).

### Positional Uniformity

5. The positional uniformity metric (PU) is computed for all variants.

6. PU values for true positive variants are predominantly < 1.0.

7. PU values for false positive variants are enriched above 1.0.

### Strand Bias (when strand-aware counting is available)

8. Fisher strand bias is computed from forward/reverse k-mer counts.

9. Known single-strand artifacts (e.g., oxidative damage G>T) produce
   elevated FS values.

### rVAF Confidence Interval

10. Bootstrap CI is computed with B >= 500 replicates.

11. CI width narrows with increasing coverage (higher precision at higher depth).

12. CI width narrows with increasing rVAF (higher precision at higher signal).

### Composite Score and Tiers

13. The composite quality score ranks true positives higher than false
    positives with AUC > 0.85 on validation data.

14. Confidence tiers (HIGH, MEDIUM, LOW) are assigned based on QUAL thresholds.

### VCF Output

15. All quality metrics are included as INFO fields in VCF output.

16. VCF output with quality fields passes bcftools view validation.
