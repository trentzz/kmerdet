# Multi-k Detection Strategy

## Motivation

The thesis validation study achieved 77% SNV sensitivity at k=31 but noted significant difficulty with INDELs, particularly insertions >15 bp. The k-mer length tradeoff analysis (see `kmer-length-tradeoffs.md`) shows that no single k value is optimal for all variant types, coverage levels, and sequence contexts. A multi-k approach can exploit the strengths of each k value to improve overall detection.

This document evaluates four concrete approaches and recommends a phased implementation plan.

## Approach 1: Independent Multi-k with Union Merging

### Workflow

Run the complete detection pipeline independently at each k value, then merge the results:

```
FASTQ -> HUMID dedup -> deduplicated FASTQ
  |-> jellyfish count (k=21) -> walk + classify -> results_k21
  |-> jellyfish count (k=31) -> walk + classify -> results_k31
  |-> jellyfish count (k=41) -> walk + classify -> results_k41
       |
       v
  merge(results_k21, results_k31, results_k41) -> final_results
```

### Implementation Details

**Jellyfish counting phase**: Three separate `jellyfish count` invocations from the same deduplicated FASTQ. Each produces an independent .jf database. These can run in parallel if memory permits (~2 GB total for 3 databases from a typical targeted panel).

**Walking phase**: For each target, run the k-mer walking algorithm independently at each k value. Each walk produces its own set of discovered paths and variant calls. Walking at each k is fully independent and parallelizable.

**Merging phase**: Combine results across k values:

1. **Exact variant matching**: Two calls at different k values match if they describe the same genomic variant (same chromosome, position, ref allele, alt allele). The matching must account for left-alignment normalization since different k values may produce slightly different breakpoint coordinates for INDELs.

2. **Union set**: Take the union of all variants detected at any k value. Each variant is annotated with which k values detected it.

3. **Meta-scoring**: Assign a combined confidence score:
   ```
   combined_score = max(score_k21, score_k31, score_k41) * detection_bonus
   where detection_bonus = 1.0 + 0.2 * (n_k_detected - 1)
   ```
   A variant detected at all 3 k values gets a 40% score bonus. A variant detected at only 1 k value retains its original score (no penalty, but no bonus).

### Resource Cost

| Resource | Single-k (k=31) | Triple-k (k=21,31,41) | Overhead |
|----------|-----------------|------------------------|----------|
| Counting time | ~1.5 min | ~4.5 min (serial) or ~1.5 min (parallel) | 1-3x |
| Walking time | ~2 min | ~6 min (serial) or ~2 min (parallel) | 1-3x |
| Memory peak | ~720 MB | ~2.1 GB (all loaded) or ~768 MB (sequential) | 1-3x |
| Disk (JF files) | ~400 MB | ~1.2 GB | 3x |
| Total pipeline time | ~6 min | ~12 min (serial) or ~7 min (parallel) | 1.2-2x |

With parallelization, the overhead is modest: ~1 minute additional for counting (if run in parallel) and walking can be parallelized across k values per target.

### Advantages and Disadvantages

**Advantages**:
- Simplest to implement (reuse existing pipeline 3 times)
- No algorithm changes required
- Full sensitivity at each k value independently
- Easy to debug (each k produces standard output)

**Disadvantages**:
- Highest resource consumption
- Simple union may increase false positive rate if not carefully filtered
- No cross-k information sharing during walking
- Redundant work: most variants are detected at all k values

## Approach 2: Consensus Voting Across k Values

### Voting Framework

Rather than a simple union, use the pattern of detection across k values as a quality signal:

| Detection pattern | Interpretation | Confidence tier |
|-------------------|----------------|-----------------|
| Detected at all 3 k values | High confidence: robust signal | Tier 1 (highest) |
| Detected at 2/3 k values | Moderate confidence | Tier 2 |
| Detected at k=31 only | Standard detection | Tier 3 |
| Detected at k=21 only | Possible low-VAF or large INDEL | Tier 4 (review) |
| Detected at k=41 only | Possible repeat-region artifact | Tier 4 (review) |

### Weighted Voting

Not all k values are equally informative for all variant types. Weight the votes:

**For SNVs:**
```
w_SNV(k=21) = 0.7   (lower specificity, more false positives)
w_SNV(k=31) = 1.0   (standard, well-validated)
w_SNV(k=41) = 0.9   (good specificity, slight coverage loss)

consensus_score_SNV = sum(w_SNV(k) * detected(k)) / sum(w_SNV(k))
```

**For short INDELs (1-10 bp):**
```
w_INDEL_short(k=21) = 1.0   (best junction coverage)
w_INDEL_short(k=31) = 0.9   (good)
w_INDEL_short(k=41) = 0.7   (some junction loss)
```

**For long INDELs (>10 bp):**
```
w_INDEL_long(k=21) = 1.0    (essential for junction spanning)
w_INDEL_long(k=31) = 0.5    (reduced junction coverage)
w_INDEL_long(k=41) = 0.2    (poor junction coverage)
```

### Decision Rules

```
IF consensus_score >= 0.7:
    CALL variant (high confidence)
ELIF consensus_score >= 0.4:
    CALL variant (moderate confidence, flag for review)
ELIF detected at exactly 1 k value AND rVAF > 0.5%:
    CALL variant (low confidence, flag for confirmation)
ELSE:
    NO CALL (insufficient evidence)
```

### Handling Conflicting rVAF Estimates

Different k values may yield different rVAF estimates for the same variant because the k-mer count decomposition operates on different sets of k-mers. Approaches to reconcile:

1. **Weighted average**: `rVAF_combined = sum(w(k) * rVAF(k)) / sum(w(k))` where w(k) is the weight for that variant type.

2. **Minimum variance**: Use the rVAF estimate from the k value that has the lowest variance in k-mer counts across the target region (indicates most uniform coverage).

3. **Coverage-weighted**: Weight by the median k-mer count at each k value, since higher coverage produces more accurate rVAF estimates:
   ```
   rVAF_combined = sum(median_count(k) * rVAF(k)) / sum(median_count(k))
   ```

### Expected Improvement

Based on the thesis data (77% SNV sensitivity at k=31), consensus voting could improve sensitivity by capturing variants missed at k=31 but detected at k=21 (low VAF, short fragments) or k=41 (repetitive contexts):

- **SNV sensitivity**: Estimated 82-87% (5-10% improvement from union detection)
- **Short INDEL sensitivity**: Estimated 15-25% improvement from k=21 contributions
- **Long INDEL sensitivity**: Estimated 20-40% improvement from k=21 contributions
- **Specificity impact**: Consensus voting maintains or slightly improves specificity vs. simple union, because single-k false positives are down-weighted

## Approach 3: Adaptive k Selection Per Target

### Concept

Rather than running all k values for all targets, precompute the optimal k for each target based on its characteristics and run only at that k value. This achieves most of the multi-k benefit at single-k cost.

### Target Analysis Phase (One-Time Preprocessing)

For each target sequence, analyze the reference genome context to determine the best k:

```rust
fn select_k_for_target(target: &Target, reference: &Reference) -> u32 {
    // 1. Check variant type
    if target.variant_type == VariantType::Insertion && target.insertion_length > 15 {
        return 21; // Long insertions need shorter k
    }

    // 2. Check repeat content
    let repeat_fraction = compute_repeat_fraction(target.sequence, reference);
    if repeat_fraction > 0.5 {
        return 41; // High repeat content needs longer k
    }

    // 3. Check anchor uniqueness
    for k in [31, 25, 21, 37, 41] {  // try preferred k first
        let first_kmer = &target.sequence[..k];
        let last_kmer = &target.sequence[target.sequence.len()-k..];
        if is_unique_in_genome(first_kmer, reference)
            && is_unique_in_genome(last_kmer, reference) {
            return k;
        }
    }

    // 4. Default
    return 31;
}
```

### Per-Target k Assignment Rules

| Target characteristic | Assigned k | Rationale |
|----------------------|------------|-----------|
| SNV in non-repetitive region | 31 | Standard, well-validated |
| SNV in low-complexity region | 41 | Need more specificity |
| Short INDEL (1-10 bp) | 31 | Standard k works well |
| Medium INDEL (11-20 bp) | 25 | Better junction coverage |
| Long INDEL (>20 bp) | 21 | Maximum junction spanning |
| ITD (any length) | 21 or 25 | Tandem duplications benefit from shorter k |
| Target in segmental duplication | 43 | Maximum specificity for paralogous regions |
| Target in homopolymer | 31-41 | Longer k reduces ambiguity |

### Implementation in kmerdet

The target catalog (FASTA + metadata) would include a per-target k recommendation:

```toml
# Target catalog metadata
[[targets]]
name = "NPM1_W288fs"
variant_type = "insertion"
insertion_length = 4
recommended_k = 25

[[targets]]
name = "FLT3_ITD"
variant_type = "itd"
recommended_k = 21

[[targets]]
name = "TP53_R175H"
variant_type = "snv"
recommended_k = 31
```

This requires creating jellyfish databases at all k values that any target needs. If the target panel uses k=21, k=25, and k=31, three databases are needed. However, during walking, each target only queries one database, so the walking phase is no slower than single-k.

### Advantages and Disadvantages

**Advantages**:
- Walking cost identical to single-k (each target queries one database)
- Tailored to each variant's characteristics
- Precomputed once per panel design (not per sample)
- Easy to understand and validate (each target has a clear rationale)

**Disadvantages**:
- Still requires multiple jellyfish databases (counting cost scales with distinct k values used)
- No cross-k validation (each variant detected at only one k)
- Requires accurate target metadata (variant type, insertion length)
- New target designs require k selection analysis

## Approach 4: Combined Evidence Scoring (Bayesian Framework)

### Bayesian Model

Formulate variant detection as Bayesian inference with evidence from multiple k values:

```
P(variant | E_k1, E_k2, E_k3) = P(E_k1, E_k2, E_k3 | variant) * P(variant)
                                  / P(E_k1, E_k2, E_k3)
```

Where:
- `E_ki` = evidence from k-mer analysis at k value ki (k-mer counts, rVAF estimate, path structure)
- `P(variant)` = prior probability of the variant (from tumor tissue genotyping)

Assuming conditional independence of evidence across k values given the true variant state (a reasonable approximation since the k-mer databases are derived from the same reads but extract different information):

```
P(variant | E_k1, E_k2, E_k3) proportional to
    P(E_k1 | variant) * P(E_k2 | variant) * P(E_k3 | variant) * P(variant)
```

### Evidence Model at Each k

For each k value, the evidence consists of:

1. **Detection indicator**: d(k) in {0, 1} -- was a variant path found?
2. **rVAF estimate**: v(k) -- the quantified variant allele frequency
3. **Path quality**: q(k) -- number of variant-specific k-mers, coverage uniformity
4. **Walking confidence**: c(k) -- did the walk complete successfully?

The likelihood for each k value can be modeled as:

```
P(E_k | variant present at VAF = f) =
    P(d=1 | f, k) *                    # probability of detection given true VAF
    N(v(k) | f, sigma_k^2) *           # rVAF estimate follows normal around true VAF
    Beta(q(k) | alpha_k, beta_k)       # path quality follows Beta distribution

P(E_k | no variant) =
    P(d=1 | f=0, k) *                  # false positive rate at k
    Uniform(v(k) | 0, noise_floor) *   # rVAF is noise
    Beta(q(k) | alpha_0, beta_0)       # path quality from noise paths
```

### Detection Sensitivity Model

The probability of detecting a variant at VAF = f with k-mer length k can be modeled as a logistic function:

```
P(detect | f, k) = 1 / (1 + exp(-beta_0 - beta_1 * log(f * depth) - beta_2 * k))
```

This captures:
- Higher VAF and higher depth increase detection probability
- Larger k slightly decreases detection probability (via the beta_2 term)
- The logistic function provides a smooth transition from non-detection to detection

Parameters (beta_0, beta_1, beta_2) can be estimated from training data or simulation.

### Handling Conflicting Evidence

The Bayesian framework naturally handles conflicting evidence. If k=21 detects a variant but k=31 does not, the combined posterior probability reflects the relative likelihoods:

- If k=21 detection is at very low rVAF (near noise floor) and k=31 non-detection is expected at that VAF, the posterior is moderate (plausible variant at very low VAF).
- If k=21 detection is at high rVAF but k=31 non-detection is unexpected, the posterior is low (likely k=21 false positive, perhaps from a repeat context).

This nuanced handling of conflicting evidence is the primary advantage over simple voting.

### Practical Considerations

The Bayesian framework requires:
1. **Training data**: Calibrate the detection sensitivity model and noise distributions from known positive and negative samples. The 10-patient thesis dataset provides initial training data but is small.
2. **Variant-type-specific models**: SNVs and INDELs have different detection characteristics at each k; separate models are needed.
3. **Computational cost**: The Bayesian calculation itself is trivial (a few multiplications per variant). The cost is in running the multi-k pipeline to generate the evidence.

### BAYSIC Precedent

The BAYSIC framework (Cantarel et al., 2014) demonstrated that Bayesian combination of variant calls from multiple callers (GATK, FreeBayes, etc.) improves both sensitivity and specificity compared to any single caller. The multi-k approach is analogous: each k value acts as a different "caller" with different strengths, and Bayesian combination extracts maximum information.

## Recommended Implementation Plan

### Phase 1: Adaptive k Per Target (Approach 3)

**Why start here**: Lowest implementation cost, no algorithm changes needed, directly addresses the INDEL sensitivity problem identified in the thesis.

Steps:
1. Add `recommended_k` field to target catalog format
2. Implement reference-based k selection heuristic for new targets
3. Create jellyfish databases at all required k values
4. Route each target to its assigned k database during walking
5. Validate on thesis dataset: expect significant INDEL improvement

### Phase 2: Independent Multi-k with Consensus Voting (Approaches 1+2)

**Why next**: Provides cross-k validation and captures variants missed by the "wrong" k assignment. Build on Phase 1 infrastructure.

Steps:
1. Run walking at all k values for all targets (parallelized)
2. Implement variant matching across k values (normalize coordinates)
3. Implement weighted consensus voting
4. Assign confidence tiers based on detection pattern
5. Validate: expect modest SNV improvement, significant INDEL improvement

### Phase 3: Bayesian Evidence Combination (Approach 4)

**Why last**: Requires training data to calibrate the model. Phase 2 generates the multi-k evidence needed for training.

Steps:
1. Collect multi-k evidence from Phase 2 runs on training samples
2. Fit detection sensitivity and noise models per variant type per k
3. Implement Bayesian posterior calculation
4. Compare to consensus voting on held-out samples
5. If superior, replace voting with Bayesian scoring; otherwise retain voting

### Cost-Benefit Summary

| Approach | Implementation effort | Expected sensitivity gain | Resource overhead |
|----------|----------------------|--------------------------|-------------------|
| Adaptive k (Phase 1) | Low | Moderate (INDELs +15-25%) | Low (counting only) |
| Consensus voting (Phase 2) | Medium | Moderate (SNVs +5-10%, INDELs +20-40%) | Medium (3x walking) |
| Bayesian (Phase 3) | High | Small incremental over voting | Same as Phase 2 |

## References

- Cantarel, B.L. et al. (2014). BAYSIC: a Bayesian method for combining sets of genome variants with improved specificity and sensitivity. *BMC Bioinformatics*, 15, 104.
- Bankevich, A. et al. (2012). SPAdes: A new genome assembly algorithm. *J. Comp. Bio.*, 19(5), 455-477.
- Li, D. et al. (2015). MEGAHIT: An ultra-fast single-node solution for large and complex metagenomics assembly. *Bioinformatics*, 31(10), 1674-1676.
- Audoux, J. et al. (2017). DE-kupl: exhaustive capture of biological variation in RNA-seq data through k-mer decomposition. *Genome Biology*, 18, 243.
- Rizk, G. et al. (2014). Mapping-free variant calling using haplotype reconstruction from k-mer frequencies. *Bioinformatics*, 34(10), 1659-1665.
