# Confidence Metrics for K-mer Variant Calls

## Overview

The current k-mer variant detection pipeline produces variant calls with minimal quality annotation: rVAF, expression (path coefficient), and minimum k-mer coverage. These metrics provide basic quantitative information but lack statistical rigor -- there is no p-value, no confidence interval, no strand bias measure, and no composite quality score. This document surveys statistical frameworks and quality metrics that could bring k-mer variant calls to the quality annotation standard expected by clinical genomics.

---

## 1. Current Approach and Limitations

### 1.1 Existing Metrics

#### rVAF (Relative Variant Allele Frequency)

The primary quantitative output of the km algorithm. Computed from NNLS (non-negative least squares) decomposition:

```
For each k-mer i:  count[i] = sum(contrib[i,j] * coef[j])  for all paths j

rVAF[j] = coef[j] / sum(coef)
```

Where `coef[j]` is the expression coefficient for path j, solved by minimizing `||contrib @ coef - counts||^2` subject to `coef >= 0`.

**Strengths**: Directly estimates the fraction of molecules supporting the variant path. Accounts for shared k-mers between overlapping paths. Handles ITDs where k-mers appear multiple times on a path.

**Limitations**: No uncertainty estimate. The NNLS solution is a point estimate with no associated confidence interval. Sensitive to outlier k-mer counts (a single high-count k-mer can dominate the regression). Does not account for k-mer count variance -- treats all counts as equally reliable.

#### Expression (Path Coefficient)

The absolute coefficient `coef[j]` from NNLS, representing the estimated number of molecules supporting path j. This provides the "depth of support" for the variant.

**Limitations**: Dependent on total coverage, making cross-sample comparison difficult. No normalization for target region properties or expected coverage.

#### Min_coverage

The minimum k-mer count along the variant path. Acts as a proxy for the "weakest link" in the variant evidence chain.

**Limitations**: A single low k-mer count can make an otherwise well-supported variant appear weak. Does not account for expected variation in k-mer counts along a path. No statistical calibration -- a min_coverage of 3 at 500x coverage means something very different from 3 at 5000x coverage.

### 1.2 What Is Missing

| Quality Dimension | Available? | Impact |
|---|---|---|
| Statistical significance (p-value) | No | Cannot assess whether observed counts are above noise |
| Confidence interval on rVAF | No | Cannot quantify uncertainty in VAF estimate |
| Strand bias | No | Cannot detect single-strand artifacts |
| Positional bias | No | Cannot detect localized count anomalies |
| Base quality weighting | No | All k-mers treated equally regardless of constituent base qualities |
| Mapping context | N/A | K-mer approach is alignment-free by design |
| Phred-scaled quality score | No | Not VCF-compatible for quality filtering |

---

## 2. Statistical Frameworks for K-mer Enrichment Testing

### 2.1 Binomial Model

The simplest statistical framework: given total coverage D and error rate e, is the observed variant k-mer count c consistent with error alone, or does it require a true variant?

**Null hypothesis (H0)**: The variant k-mers arise from sequencing errors. Expected count under H0: `D * e_k`, where `e_k` is the per-k-mer error probability (probability that a specific base substitution produces this particular k-mer from the reference k-mer).

**Alternative hypothesis (H1)**: A true variant exists at VAF f. Expected count under H1: `D * (f + e_k)`.

**Test statistic**: The one-sided binomial p-value:

```
p = P(X >= c | X ~ Binomial(D, e_k))
  = 1 - BinomialCDF(c - 1, D, e_k)
```

Where:
- `D` = total k-mer coverage at this position (sum of reference + variant counts)
- `c` = observed variant k-mer count
- `e_k` = expected per-k-mer error rate (~0.001 for substitutions, higher for indels)

**Example**: At 3000x coverage with error rate 0.001, a variant k-mer count of 8:
```
p = 1 - BinomialCDF(7, 3000, 0.001)
  = 1 - 0.947
  = 0.053
```
This borderline p-value illustrates the difficulty of calling variants at the noise boundary.

**Advantages**: Simple, well-understood, computationally trivial. Provides a clear threshold for detection.

**Limitations**: Assumes a single error rate for all k-mers (in reality, error rates are sequence-context-dependent). Does not account for overdispersion in k-mer counts. The error rate `e_k` must be estimated, which is itself uncertain.

The binomial model has been used by SNVer (Wei et al., 2011) for alignment-based variant calling and by statistical frameworks for assessing limits of detection in clinical sequencing (Jennings et al., 2017).

### 2.2 Fisher's Exact Test / Chi-squared Test

Construct a 2x2 contingency table comparing observed vs. expected counts for variant and reference k-mers:

```
                  Variant k-mers    Reference k-mers
Observed count:        c_var              c_ref
Expected (error):      D * e_k            D * (1 - e_k)
```

Fisher's exact test computes the probability of observing a table as extreme as the one obtained, given fixed marginal totals. This is directly analogous to the Fisher Strand test used in GATK for strand bias detection.

**For k-mer variant detection**, the contingency table can be constructed per k-mer position along the variant path:

```
                  K-mer i is variant    K-mer i is reference
Position 1:           c_var_1              c_ref_1
Position 2:           c_var_2              c_ref_2
...
```

A chi-squared test across all k-mer positions tests whether the variant signal is consistent across the path (uniform enrichment) or concentrated at specific positions (potential artifact).

**Advantages**: Non-parametric, no distributional assumptions. Well-suited for small sample sizes (Fisher's exact). Can be extended to test multiple hypotheses simultaneously.

**Limitations**: Does not naturally incorporate a prior probability of variant existence. The 2x2 table structure may oversimplify the multi-k-mer evidence structure.

### 2.3 Poisson Model

For low-count k-mers (count < 20), the Poisson distribution provides a natural model:

```
P(count = c | lambda) = (lambda^c * e^(-lambda)) / c!
```

Where `lambda` is the expected count. Under the null hypothesis (error only), `lambda_0 = D * e_k`. Under the alternative (true variant at VAF f), `lambda_1 = D * f`.

**Likelihood ratio test**:

```
LR = P(count = c | lambda_1) / P(count = c | lambda_0)
   = (lambda_1 / lambda_0)^c * exp(lambda_0 - lambda_1)
```

Taking the log:

```
log(LR) = c * log(lambda_1 / lambda_0) + (lambda_0 - lambda_1)
```

**Combined across k-mers**: For k-mers along the variant path with counts c_1, c_2, ..., c_k, the combined log-likelihood ratio is the sum of individual log-LRs (assuming independence):

```
log(LR_combined) = sum_i [c_i * log(lambda_1 / lambda_0) + (lambda_0 - lambda_1)]
```

**Advantages**: Appropriate for low counts where the normal approximation fails. Natural model for count data. The combined LR across k-mers integrates evidence from the entire variant path.

**Limitations**: Assumes counts are independent (violated when k-mers share reads). The Poisson assumption of mean = variance may not hold due to PCR amplification bias and coverage non-uniformity.

### 2.4 Negative Binomial Model

The negative binomial distribution extends the Poisson to allow overdispersion (variance > mean), which is commonly observed in sequencing count data due to PCR amplification bias, GC bias, and capture efficiency variation.

```
P(count = c | mu, alpha) = C(c + 1/alpha - 1, c) * (1/(1 + alpha*mu))^(1/alpha) * (alpha*mu/(1 + alpha*mu))^c
```

Where:
- `mu` = expected count (mean)
- `alpha` = overdispersion parameter (alpha = 0 reduces to Poisson)
- `variance = mu + alpha * mu^2`

**Estimation of overdispersion**: The overdispersion parameter `alpha` can be estimated from the reference k-mer counts along the target region. Reference k-mers should have approximately uniform counts (all representing the reference allele at the same coverage), so the observed variance relative to the mean provides an estimate of `alpha`:

```
alpha_hat = (var(ref_counts) - mean(ref_counts)) / mean(ref_counts)^2
```

**Application to variant detection**: The negative binomial p-value for a variant k-mer count `c` under the null hypothesis:

```
p = 1 - NegBinomCDF(c - 1, mu = D * e_k, alpha = alpha_hat)
```

Research by Love et al. (DESeq2, 2014) and Robinson et al. (edgeR, 2010) established the negative binomial as the standard model for RNA-seq count data, demonstrating that sequence count data consistently exhibit overdispersion relative to the Poisson. While some have argued that the negative binomial is an imperfect fit for sequencing data (Hawinkel et al., 2020), it remains the most widely used model for overdispersed counts in genomics.

**Advantages**: Accounts for overdispersion, providing more conservative (and realistic) p-values than the Poisson model. Well-characterized statistically with established estimation procedures.

**Limitations**: Requires estimating the overdispersion parameter, which itself has uncertainty. With small numbers of reference k-mers per target (~100-300), the estimate of alpha may be imprecise.

### 2.5 Bayesian Posterior Model

A Bayesian framework integrates prior knowledge about variant probability with the k-mer count evidence:

```
P(variant | counts) = P(counts | variant) * P(variant) / P(counts)
```

**Prior probability P(variant)**: For tumor-informed detection (the kam use case), the prior is high -- the variant is already known from the tumor biopsy. A reasonable prior might be P(variant) = 0.5 (reflecting that the variant may or may not be present in this particular blood draw). For de novo discovery, the prior would be much lower (P(variant) ~ 10^-6 per position).

**Likelihood P(counts | variant)**: Computed from the binomial or Poisson model at the maximum likelihood VAF estimate. For k-mer counts c_1, ..., c_k along the variant path:

```
P(counts | variant, VAF = f) = product_i Poisson(c_i | lambda = D * f)
```

The likelihood under H0 (no variant) is:
```
P(counts | no variant) = product_i Poisson(c_i | lambda = D * e_k)
```

**Posterior odds**:
```
posterior_odds = likelihood_ratio * prior_odds
               = [P(counts | variant) / P(counts | no variant)] * [P(variant) / P(no variant)]
```

**Posterior probability**:
```
P(variant | counts) = posterior_odds / (1 + posterior_odds)
```

The BATCAVE algorithm (Valeri et al., 2020) demonstrated the power of tumor-specific priors for somatic variant calling, using high-confidence variants to estimate the mutational prior and mutation rate, then applying this prior to improve calling of low-confidence variants. BayVarC (2024) extended this approach specifically for liquid biopsy, applying Bayesian inference to quantify noise in a locus-specific manner.

**Advantages**: Naturally incorporates prior information from the tumor biopsy. Provides a posterior probability that is directly interpretable. Can be updated as more data (serial blood draws) become available.

**Limitations**: Sensitive to prior specification. Requires choosing between Bayesian model comparison and Bayesian parameter estimation. Computationally more expensive than frequentist tests.

---

## 3. Quality Metrics to Implement

### 3.1 Strand Bias

**Concept**: For each k-mer along the variant path, record whether it came from a forward-strand read or a reverse-strand read. A true variant should have roughly equal support from both strands; a sequencing artifact often shows extreme strand bias.

**Current limitation**: The jellyfish database with `-C` (canonical) counting merges forward and reverse k-mers. Strand information is lost at the counting stage.

**Implementation approach**: During k-mer counting (or as a secondary pass), maintain separate counts for forward-strand and reverse-strand k-mers:

```
strand_bias = |forward_count - reverse_count| / (forward_count + reverse_count)
```

A strand_bias value close to 0 indicates balanced support; close to 1 indicates single-strand support (likely artifact).

**Fisher's exact test for strand bias** (as used in GATK):

```
                Forward    Reverse
Reference:      ref_fwd    ref_rev
Variant:        var_fwd    var_rev
```

The Fisher's exact p-value tests whether the variant allele is disproportionately represented on one strand. GATK reports this as the FisherStrand (FS) annotation, with Phred-scaled values > 60 indicating significant strand bias.

Research by Guo et al. (2012) showed that strand bias filtering removes 70-90% of false positive variant calls in Illumina data, making it one of the most effective single quality filters.

**For kmerdet implementation**: This requires either (a) modifying jellyfish counting to maintain strand-specific databases, or (b) implementing strand-aware k-mer counting in the pure Rust k-mer counter. The latter is the recommended approach for the Rust reimplementation.

### 3.2 Positional Bias (K-mer Uniformity)

**Concept**: For an SNV, exactly k=31 k-mers span the variant position. If the variant is real, all 31 k-mers should have approximately equal counts. If counts are highly non-uniform (e.g., one k-mer has count 20 and the rest have count 2), this suggests the variant evidence is driven by a small number of reads with specific start positions, which is suspicious.

**Metric**: Coefficient of variation of k-mer counts along the variant path:

```
positional_cv = std(variant_kmer_counts) / mean(variant_kmer_counts)
```

Expected value for a true variant: CV ~ 0.3-0.5 (reflecting natural Poisson variation). Values > 1.0 suggest positional bias.

**Alternative metric**: Chi-squared goodness-of-fit test comparing observed k-mer counts to the expected uniform distribution:

```
chi2 = sum_i (observed_i - expected)^2 / expected
```

Where expected = mean(observed counts). With k-1 degrees of freedom, a p-value can be computed.

### 3.3 K-mer Quality Score (Composite)

A composite quality score integrating multiple evidence dimensions:

```
kmer_quality = f(p_value, strand_bias, positional_cv, min_coverage, rVAF)
```

**Proposed composite score** (logistic regression or weighted sum):

```
logit(quality) = w1 * log(p_value)
               + w2 * strand_bias
               + w3 * positional_cv
               + w4 * log(min_coverage)
               + w5 * log(rVAF)
               + w6 * sequence_complexity
               + w7 * homopolymer_proximity
```

Weights can be trained from the validation dataset (10 patients with known true/false calls) or set heuristically based on the relative importance of each feature.

**The Merfin approach** (Formenti et al., 2022) provides precedent for k-mer-based quality evaluation independent of alignment. Merfin evaluates each variant based on the expected k-mer multiplicity in the reads, checking whether introducing the variant makes the k-mer spectrum more consistent with the expected copy number. This concept could be adapted: for each candidate variant, compute how well the observed k-mer counts match the expected counts under the variant hypothesis vs. the reference hypothesis.

### 3.4 Confidence Intervals on rVAF

The NNLS decomposition provides a point estimate of rVAF but no uncertainty. Several approaches can provide confidence intervals:

#### Bootstrap Confidence Intervals

1. Resample the k-mer counts with replacement (or from a parametric distribution fitted to the observed counts)
2. Solve NNLS on each bootstrap sample
3. Compute rVAF for each bootstrap replicate
4. The 2.5th and 97.5th percentiles of the bootstrap rVAF distribution form the 95% CI

**Implementation**: For B=1000 bootstrap replicates with k=31 k-mers, this requires 1000 NNLS solutions per variant. NNLS for a 31-variable system is computationally trivial (microseconds per solve), so 1000 replicates add negligible overhead.

**Considerations for NNLS bootstrap**: Standard bootstrap methods must be adapted to respect the non-negativity constraint. Resampled counts should maintain the constraint that all coefficients >= 0. Research on Bootstrap NNLS (BNNLS) suggests that bootstrap methods that resample from distributions respecting the order/non-negativity restrictions outperform naive resampling, particularly for small sample sizes.

#### Analytical Confidence Intervals

For the linear regression model (before applying non-negativity constraints), the covariance matrix of the coefficients is:

```
Cov(coef) = sigma^2 * (C^T C)^(-1)
```

Where `C` is the contribution matrix and `sigma^2` is estimated from the residual variance:

```
sigma_hat^2 = ||C @ coef - counts||^2 / (n - p)
```

The 95% CI for coefficient j is:

```
coef[j] +/- t_{n-p, 0.025} * sqrt(Cov(coef)[j,j])
```

This CI for `coef[j]` propagates to a CI for rVAF[j] via the delta method:

```
se(rVAF[j]) ~ se(coef[j]) / sum(coef) * sqrt(1 + rVAF[j]^2 * sum(se(coef)^2) / sum(coef)^2)
```

**Limitation**: The analytical CI assumes Gaussian residuals and ignores the non-negativity constraint. When the true coefficient is near zero, the analytical CI may include negative values, which are then clipped to zero, producing asymmetric intervals. The bootstrap approach handles this more naturally.

#### Profile Likelihood Confidence Intervals

Fix rVAF[j] at a grid of values, re-optimize the remaining coefficients subject to the constraint, and find the values where the log-likelihood drops by 1.92 (for 95% CI based on chi-squared approximation). This is the most rigorous approach but computationally expensive for routine use.

### 3.5 Phred-Scaled Quality Score

For VCF compatibility, convert the variant p-value to a Phred-scaled quality score:

```
QUAL = -10 * log10(p_value)
```

Where p_value is from the chosen statistical test (binomial, Poisson, negative binomial, or Bayesian posterior).

| p-value | QUAL | Interpretation |
|---------|------|---------------|
| 0.1 | 10 | Low confidence |
| 0.01 | 20 | Marginal |
| 0.001 | 30 | Good (typical filter threshold) |
| 10^-4 | 40 | High confidence |
| 10^-6 | 60 | Very high confidence |
| 10^-10 | 100 | Extremely high confidence |

**VCF FILTER field**: Based on QUAL and other metrics:

```
PASS:              QUAL >= 30 AND strand_bias < 0.8 AND positional_cv < 1.5
LowQual:           20 <= QUAL < 30
StrandBias:        strand_bias >= 0.8
PositionalBias:    positional_cv >= 1.5
LowCoverage:       min_coverage < 3
BelowThreshold:    QUAL < 20
```

---

## 4. Comparison with Alignment-Based Quality Metrics

### 4.1 Standard VCF Quality Annotations

Alignment-based variant callers (GATK HaplotypeCaller, Mutect2, Strelka2) produce a rich set of quality annotations:

| Metric | Source | Description | K-mer Equivalent? |
|---|---|---|---|
| QUAL | Variant caller | Phred-scaled probability of variant | Computable from k-mer p-value |
| GQ | Genotyper | Phred-scaled genotype quality | Not applicable (no genotyping) |
| DP | Pileup | Read depth at position | Total k-mer coverage at position |
| AD | Pileup | Allelic depths (ref, alt) | Reference vs. variant k-mer counts |
| AF | Caller | Allele frequency estimate | rVAF from NNLS |
| MQ | Aligner | Root mean square mapping quality | Not applicable (alignment-free) |
| MQRankSum | Caller | Rank sum test of MQ for ref vs. alt | Not applicable |
| FS | Caller | Fisher strand bias (Phred-scaled) | Computable with strand-aware counting |
| SOR | Caller | Strand odds ratio | Computable with strand-aware counting |
| ReadPosRankSum | Caller | Rank sum of read positions | Analogous to positional bias metric |
| BaseQRankSum | Caller | Rank sum of base qualities | Not directly available (no per-base quality in k-mer space) |
| TLOD | Mutect2 | Tumor log odds | Computable from Bayesian model |
| NLOD | Mutect2 | Normal log odds | Computable from panel-of-normals |

### 4.2 What Alignment Provides That K-mers Do Not

**Mapping quality (MAPQ)**: Alignment-based methods know whether a read maps uniquely or to multiple locations. Multi-mapped reads reduce confidence in the variant call. K-mer methods have no concept of mapping uniqueness for individual k-mers, though k-mer multiplicity in the reference genome could serve as a proxy.

**Base quality**: Each base in a read has an associated Phred quality score from the sequencer. Alignment-based callers weight variant evidence by base quality, downweighting low-quality bases. K-mer counting discards this information entirely -- a k-mer containing a Q10 base (10% error probability) is counted equally with one containing a Q40 base (0.01% error probability).

**Mate pair information**: Paired-end reads provide geometric constraints (insert size, mate orientation) that can identify artifacts. K-mer approaches lose this paired-end structure.

**Local realignment**: Tools like GATK HaplotypeCaller perform local de Bruijn graph assembly and realignment to resolve complex variants. While this is conceptually similar to k-mer walking, the integrated alignment-reassembly-calling framework allows the caller to consider multiple evidence sources simultaneously.

### 4.3 What K-mers Provide That Alignment Does Not

**No reference bias**: K-mer counts are not influenced by how well reads align to the reference. In regions where the variant allele aligns poorly (large indels, rearrangements), alignment-based callers lose evidence while k-mer counts are unaffected.

**No multi-mapping ambiguity**: K-mer counts are absolute. A k-mer either exists in the database with count N, or it does not. There is no ambiguity about where the reads "came from" -- the k-mer count is the evidence.

**Simultaneous multi-allelic quantification**: The NNLS decomposition naturally handles multiple alleles at the same locus (reference + multiple variant paths), quantifying each simultaneously. Alignment-based callers must explicitly model multi-allelic sites.

**Graph-level evidence**: The connectivity and topology of the k-mer graph provide structural evidence about the variant. A variant supported by a clean, well-connected path with consistent k-mer counts is more credible than one with a fragmented, noisy path -- and this graph structure has no alignment analog.

---

## 5. Recommended Implementation Priority

### Phase 1: Core Statistical Testing

1. **Binomial p-value per variant**: Compute at each k-mer position along the variant path using the local coverage and estimated error rate. Report the worst (highest) p-value as the conservative variant p-value.

2. **Phred-scaled QUAL**: Convert the combined p-value (Fisher's method to combine per-k-mer p-values) to a Phred-scaled score for VCF output.

3. **Coverage-normalized min_coverage**: Report min_coverage as a fraction of median local coverage, providing a metric that is comparable across samples with different depths.

### Phase 2: Extended Quality Metrics

4. **Strand bias**: Implement strand-aware k-mer counting in the Rust k-mer counter. Report Fisher strand bias (FS) in VCF output.

5. **Positional uniformity**: Compute CV of k-mer counts along the variant path. Report as a VCF INFO field.

6. **Bootstrap CI on rVAF**: Implement 1000-replicate bootstrap for NNLS decomposition. Report 95% CI bounds in VCF FORMAT fields.

### Phase 3: Advanced Models

7. **Negative binomial p-value**: Estimate overdispersion from reference k-mer counts. Use NB model for more conservative p-values.

8. **Bayesian posterior**: Implement tumor-informed prior. Report posterior probability alongside frequentist p-value.

9. **Composite quality score**: Train a logistic classifier on the validation dataset using all available quality features.

---

## 6. Mathematical Details for Implementation

### 6.1 Combined P-value Across K-mers (Fisher's Method)

For k-mers with individual p-values p_1, p_2, ..., p_m:

```
X = -2 * sum(log(p_i))
```

Under H0, X follows a chi-squared distribution with 2m degrees of freedom:

```
combined_p = 1 - ChiSquaredCDF(X, df = 2m)
```

**Caveat**: Fisher's method assumes independent tests. K-mers along a path share reads, so they are correlated. A conservative correction: use the effective number of independent k-mers, estimated as `m_eff = m / autocorrelation_factor`, where the autocorrelation factor reflects the read length relative to k.

For 150 bp reads and k=31, each read contributes ~120 overlapping k-mers. K-mers separated by more than 120 positions are approximately independent. For a variant path of length m k-mers, the effective number of independent k-mers is approximately `m_eff = m * k / read_length = m * 31/150 ~ 0.2m`.

### 6.2 Error Rate Estimation

The per-k-mer error rate `e_k` is not constant across the genome. For more accurate testing, estimate `e_k` locally:

**Approach 1 (Empirical)**: In a window around the variant position, count the number of non-reference k-mers at positions where no variant is expected. The average count divided by coverage gives the local error rate.

**Approach 2 (Model-based)**: Use the Illumina error model where error probability depends on base quality, position in read, and local sequence context. The k-mer-level error rate is approximately:

```
e_k ~ 1 - (1 - e_base)^k ~ k * e_base  (for small e_base)
```

For k=31 and average base error rate e_base = 0.001 (Q30), e_k ~ 0.031. But this is the probability that any base in the k-mer is wrong. The probability that the k-mer becomes a specific alternative k-mer (matching the variant) is much lower: approximately e_base / 3 (one specific substitution out of 3 possible).

**Approach 3 (K-merald)**: The K-merald approach (Denti et al., 2023) constructs k-mer-based sequencing error profiles from the data itself, genotyping variants based on the empirical error distribution rather than a fixed rate. This approach achieves higher genotype quality than edit-distance-based methods, particularly for indels.

### 6.3 Negative Binomial Parameter Estimation

Given reference k-mer counts r_1, r_2, ..., r_n along the target region:

**Method of moments**:
```
mu_hat = mean(r_i)
var_hat = var(r_i)
alpha_hat = (var_hat - mu_hat) / mu_hat^2
```

If alpha_hat <= 0, the data is not overdispersed and the Poisson model suffices.

**Maximum likelihood** (more accurate for small n): Solve iteratively using the profile likelihood or Newton-Raphson on the negative binomial log-likelihood.

**Typical values**: For duplex-sequenced cfDNA with UMI deduplication, we observe alpha ~ 0.01-0.1 (mild overdispersion) in well-captured regions and alpha ~ 0.5-2.0 (severe overdispersion) in regions with GC bias or capture non-uniformity.

---

## References

- Wei et al. "[SNVer: a statistical tool for variant calling in analysis of pooled or individual next-generation sequencing data](https://pmc.ncbi.nlm.nih.gov/articles/PMC3201884/)." Nucleic Acids Research, 2011.
- Jennings et al. "[Assessing Limit of Detection in Clinical Sequencing](https://www.jmdjournal.org/article/S1525-1578(21)00006-4/fulltext)." Journal of Molecular Diagnostics, 2017.
- Koboldt et al. "[Statistical modeling for sensitive detection of low-frequency single nucleotide variants](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2905-x)." BMC Genomics, 2016.
- Denti et al. "[Allele detection using k-mer-based sequencing error profiles](https://pmc.ncbi.nlm.nih.gov/articles/PMC10625474/)." Bioinformatics Advances, 2023.
- Formenti et al. "[Merfin: improved variant filtering, assembly evaluation and polishing via k-mer validation](https://pmc.ncbi.nlm.nih.gov/articles/PMC9745813/)." Nature Methods, 2022.
- Valeri et al. "[BATCAVE: calling somatic mutations with a tumor- and site-specific prior](https://pmc.ncbi.nlm.nih.gov/articles/PMC7003682/)." NAR Genomics and Bioinformatics, 2020.
- "[BayVarC: an ultra-sensitive ctDNA variant caller using Bayesian approach](https://www.biorxiv.org/content/10.1101/2024.02.03.578772v1)." bioRxiv, 2024.
- "[Empirical Bayes single nucleotide variant-calling for next-generation sequencing data](https://www.nature.com/articles/s41598-024-51958-z)." Scientific Reports, 2024.
- Guo et al. "[The effect of strand bias in Illumina short-read sequencing data](https://link.springer.com/article/10.1186/1471-2164-13-666)." BMC Genomics, 2012.
- GATK Best Practices. "[FisherStrand annotation](https://gatk.broadinstitute.org/hc/en-us/articles/360040096152-FisherStrand)." Broad Institute.
- GATK Best Practices. "[StrandBiasBySample annotation](https://gatk.broadinstitute.org/hc/en-us/articles/360040096492-StrandBiasBySample)." Broad Institute.
- VarSome. "[Variant calling and quality filters](https://docs.varsome.com/en/variant-calling-and-quality-filters)."
- Love et al. "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome Biology, 2014.
- Robinson et al. "edgeR: a Bioconductor package for differential expression analysis of digital gene expression data." Bioinformatics, 2010.
- Hawinkel et al. "[Sequence count data are poorly fit by the negative binomial distribution](https://pmc.ncbi.nlm.nih.gov/articles/PMC7192467/)." PLOS One, 2020.
- Thesis Chapters 4-6: Validation Study Design, Results, and Discussion.
