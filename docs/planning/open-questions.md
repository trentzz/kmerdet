# Open Questions

Questions we do not yet have answers to, organized by topic. For each question: what it is, why it matters, how to investigate, and which research document is most relevant.

---

## K-mer Counting

### Q1: What is the actual count accuracy of jellyfish at high load factors?

**Why it matters**: The thesis pipeline uses `-s 100M` for hash table sizing. If the number of distinct k-mers exceeds this (which it can for high-depth targeted panels with off-target reads), jellyfish creates overflow files that require merging. `jellyfish-deep-dive.md` notes that split counts across files before merging can cause inaccuracy. If counts are even slightly wrong at the critical low-VAF range (counts of 2-5), variants could be missed or artifacts promoted.

**How to investigate**: Create a synthetic FASTQ with known k-mer composition (specific k-mers at precisely controlled counts). Run jellyfish count with varying `-s` values (undersized, correctly sized, oversized). Query the known k-mers and compare expected vs. observed counts. Measure: (a) absolute accuracy at count=2, 5, 10; (b) whether undersized tables cause count splitting; (c) whether merged databases are bit-for-bit equivalent to correctly-sized single databases.

**Research reference**: `jellyfish-deep-dive.md` Section "Hash Collision Behavior at High Load", specifically the discussion of undersized tables and overflow file merging.

---

### Q2: Does canonical counting meaningfully affect sensitivity compared to non-canonical?

**Why it matters**: Canonical counting (`-C`) merges both strand orientations, effectively doubling the count for most k-mers. This should improve signal-to-noise for low-VAF variants. However, for palindromic or near-palindromic k-mers, canonical counting may introduce artifacts by merging signals from unrelated genomic contexts. The net sensitivity impact has not been measured empirically for ctDNA detection.

**How to investigate**: Run detection at k=31 with both canonical and non-canonical jellyfish databases on the validation cohort. Compare: (a) per-k-mer counts for variant k-mers (should be ~2x higher with canonical); (b) overall sensitivity (canonical should be higher); (c) whether any variants are detected in one mode but not the other; (d) near-palindromic k-mers that might cause problems. Also measure memory: canonical should be ~50% smaller database.

**Research reference**: `jellyfish-deep-dive.md` Section "Canonical Counting Implications" for the tradeoff analysis, noting that canonical counting loses strand information but doubles effective count.

---

### Q3: Is the -L 2 cutoff optimal, or should we use -L 1 (or -L 0)?

**Why it matters**: `-L 2` discards all singleton k-mers from the jellyfish database. At 0.1% VAF with 500x post-dedup coverage, `jellyfish-deep-dive.md` calculates the expected variant k-mer count is ~0.4 -- well below 2, meaning many true variant k-mers are discarded at the counting stage before detection even begins. The sensitivity-landscape analysis estimates that ~30-40% of true variant k-mers at 0.1% VAF have count 1 and are lost to `-L 2`. However, `-L 1` retains enormous numbers of error k-mers (each sequencing error creates a novel singleton k-mer), increasing database size by ~70% and potentially slowing walking.

**How to investigate**: This is part of Research Experiment 4 (Threshold Optimization). Generate .jf databases with -L 1 and -L 2 for the same samples. Compare: (a) database sizes (expect -L 1 to be ~3x larger); (b) number of variant k-mers retained vs. lost; (c) sensitivity at -L 1 vs -L 2 with walking threshold count=2 (which provides the effective floor anyway); (d) walking runtime (more k-mers in database = more extension candidates to evaluate). The key insight from `jellyfish-deep-dive.md` is that filtering should be deferred to the walking phase rather than the counting phase, where contextual information (sibling counts, local coverage) can inform the decision.

**Research reference**: `jellyfish-deep-dive.md` Section "Lower-Count Cutoff (-L Flag)" including the impact table at various VAF/depth combinations, and the alternative approach of post-counting filtering during walking.

---

### Q4: For native Rust counting, will DashMap's overhead be acceptable for production use?

**Why it matters**: `counting-alternatives.md` estimates DashMap at ~3.8 GB for 50M k-mers vs. jellyfish at ~500 MB and a custom hash table at ~750 MB. For workstations with 16 GB RAM holding multiple k-mer databases simultaneously (multi-k at k=21, 31, 41), DashMap could consume >11 GB just for counting, leaving insufficient memory for walking and graph construction.

**How to investigate**: Benchmark DashMap vs. the custom lock-free hash table described in `counting-alternatives.md` on realistic targeted panel data. Measure: (a) memory footprint at 50M k-mers; (b) insertion throughput (k-mers per second per thread); (c) query throughput during walking; (d) scaling with thread count (1, 2, 4, 8 threads). If DashMap overhead is unacceptable, implement the custom hash table described in Phase 5 of the development roadmap.

**Research reference**: `counting-alternatives.md` Section "Memory Requirements" for DashMap vs custom hash table analysis, and the lock-free CAS-based design from CHTKC.

---

## Walking Algorithm

### Q5: What is the optimal max_stack for long INDELs without causing exponential blowup?

**Why it matters**: `walking-improvements.md` documents the tension between `max_stack=500` (too small for >100 bp insertions, which require traversing ~100+ additional nodes) and unlimited stack (which can cause exponential blowup in repetitive regions). FLT3-ITD can be 3-400 bp, requiring 300-400 additional nodes. But a poly-A run of length 50 can generate 4^50 possible paths if unconstrained.

**How to investigate**: Benchmark with increasingly large insertions (10, 20, 50, 100, 200, 400 bp) using synthetic targets:
1. Generate target FASTA with known insertions of each size
2. Generate synthetic .jf databases with controlled coverage
3. Run walking with max_stack = 500, 1000, 2000, 5000
4. Measure: (a) detection success at each stack limit; (b) wall-clock time; (c) maximum observed stack depth for successful detection
5. Separately test on repetitive targets (homopolymer runs, dinucleotide repeats) to find where exponential blowup begins

The dynamic max_stack approach from `walking-improvements.md` Section 4 (`effective_max_stack = config.max_stack + target.expected_indel_size * 2`) should also be tested.

**Research reference**: `walking-improvements.md` Section 4 on long INDEL walking strategies, `sensitivity-landscape.md` Section 1.3 on walking depth limits.

---

### Q6: Does bidirectional walking actually improve sensitivity for deletions?

**Why it matters**: `walking-improvements.md` Section 3 proposes bidirectional walking primarily to handle deletions at target boundaries and variants near anchor k-mers. The theoretical argument is sound (backward walking can bridge gaps that forward walking cannot), but the empirical impact on the validation data is unknown. It is possible that the thesis targets are designed with sufficient flanking (35 bp) to avoid most boundary effects.

**How to investigate**: Implement bidirectional walking per `walking-improvements.md` Section 3 and run on the validation cohort:
1. Compare forward-only vs. bidirectional detection
2. Identify specific variants recovered by bidirectional walking
3. Analyze whether recovered variants are all near target boundaries (confirming the hypothesized mechanism) or in other positions (suggesting additional benefits)
4. Measure the computational overhead (bidirectional is ~2x walking cost)

If bidirectional walking recovers zero or very few additional variants, the complexity may not be justified for the current target design. However, shorter flanking or different panel designs might benefit more.

**Research reference**: `walking-improvements.md` Section 3 on bidirectional walking rationale and implementation, `room-for-improvement.md` Section 10 on target design sensitivity.

---

### Q7: Can we detect when walking has failed vs. when there is truly no variant?

**Why it matters**: Currently, when walking finds no variant path for a target, the output is simply "no variant detected." This is ambiguous: it could mean the variant is absent (true negative at this timepoint) or that the walking algorithm failed to discover a present variant (false negative). In MRD monitoring, this distinction is clinically critical -- a rising VAF trend should not be interrupted by a walking failure masquerading as a true negative.

**How to investigate**: Design a walking quality metric that indicates whether the walk was "healthy" or "troubled":
- Healthy walk: reference path fully traversed source-to-sink, all reference k-mers present with expected coverage, no excessive branching
- Troubled walk: hit max_stack or max_break limits, reference path incomplete, coverage drops in the target region

Implement this metric and correlate with known false negatives from Experiment 2. If troubled walks predict false negatives with high accuracy, the metric can be reported alongside results to flag unreliable no-call results.

**Research reference**: `confidence-metrics.md` for quality metric design, `sensitivity-landscape.md` Section 3.3 on walking as signal loss point.

---

### Q8: How effective is the relaxed threshold during non-reference traversal?

**Why it matters**: `walking-improvements.md` Section 4 proposes reducing the threshold by 50% after 5+ consecutive non-reference extensions, to help the walker traverse long insertion junctions where each k-mer independently must pass threshold. The idea is that the walker is already committed to an INDEL path and should be more permissive about following it through. But this relaxation could also cause the walker to follow error paths deeper than intended.

**How to investigate**: Implement the relaxation and measure on the validation data:
1. Count how many long INDELs (>10 bp) are rescued by relaxation
2. Count how many new false positive paths are created by relaxation
3. Vary the relaxation factor (0.25, 0.5, 0.75) and the engagement threshold (3, 5, 10 consecutive non-ref extensions)
4. Determine the optimal balance between INDEL recovery and false path generation

**Research reference**: `walking-improvements.md` Section 4 for the relaxation proposal, `sensitivity-landscape.md` Section 1.3 on how junction k-mer counts compound with VAF for long INDELs.

---

## Graph Construction

### Q9: Is there a principled way to set edge weights, or will any weighting scheme have failure modes?

**Why it matters**: `graph-improvements.md` Section 2 proposes several alternatives to the binary 0.01/1.0 weighting: coverage-weighted (1/count), log-coverage (-log2(count/total)), and ratio-based (1 - count/expected). Each has theoretical advantages and disadvantages, but none is obviously superior. The recommended hybrid approach (binary for pathfinding, coverage for ranking) sidesteps the question for pathfinding but does not resolve it for cases where Dijkstra's choice between multiple variant paths matters.

**How to investigate**: This is Research Experiment 9 (Graph Weighting Alternatives). Run all four weighting schemes on the validation data and compare path discovery. Specifically look for cases where:
- Coverage weighting finds a variant path that binary weights miss (because binary weights cannot distinguish well-supported from poorly-supported non-reference edges)
- Coverage weighting misses a variant path that binary weights find (because a high-count error k-mer creates a lower-weight path that diverts Dijkstra away from the true variant)
If both cases exist, no single scheme dominates, and the hybrid approach is confirmed as the best strategy.

**Research reference**: `graph-improvements.md` Section 2 for complete analysis of four weighting schemes, `room-for-improvement.md` Section 8 on fixed weight limitations.

---

### Q10: How do we handle targets with >100 paths efficiently?

**Why it matters**: In repetitive regions, the all-shortest-paths enumeration can produce hundreds or thousands of paths, most representing different positions for the same sequencing error within a repeat. This makes NNLS quantification slow (contribution matrix grows exponentially) and output unwieldy. The current `max_node=10000` and `max_break=10` limits prevent runaway computation but may truncate the search before finding the true variant path.

**How to investigate**: Characterize the path count distribution across all targets in the validation data:
1. How many targets produce >10 paths? >50? >100?
2. Are high-path-count targets always repetitive regions?
3. Does graph pruning (Phase 3) reduce path counts below 50 for all targets?
4. Test K-shortest-paths (Yen's, K=10-50) from `graph-improvements.md` Section 5: does K=10 capture all true variants while limiting output?
5. Test if the 2^n cluster limit of 256 paths from `graph-improvements.md` Section 3 is ever reached on real data

If graph pruning and K-shortest-paths together cannot control path count, consider a two-pass approach: fast initial detection with aggressive pruning, followed by targeted re-analysis with relaxed parameters for targets of interest.

**Research reference**: `graph-improvements.md` Section 3 on compound variants (2^n paths), Section 4 on graph pruning strategies, Section 5 on K-shortest-paths algorithms.

---

### Q11: Can A* search reduce pathfinding time for complex graphs?

**Why it matters**: `graph-improvements.md` Section 5 discusses A* search with a heuristic based on remaining reference distance. For the typical kmerdet graph (V < 10,000, E < 40,000), Dijkstra is already fast (milliseconds). A* would provide benefit only for much larger graphs (WGS targets, structural variant detection with extended walking windows). The question is whether Phase 6 SV extensions or very large targets would benefit.

**How to investigate**: Implement A* with the admissible heuristic `h(n) = (ref_length - position(n)) * REF_EDGE_WEIGHT` per `graph-improvements.md`. Benchmark against Dijkstra on:
1. Standard 50-target panels (expect no meaningful difference)
2. Synthetic large targets (1000+ bp) with many branching points
3. SV detection windows (5000+ bp)

If A* provides >2x speedup for SV-scale graphs, include as optional mode. Otherwise, standard Dijkstra suffices.

**Research reference**: `graph-improvements.md` Section 5 on A* with the AStarix heuristic precedent, comparison table of pathfinding algorithms.

---

## Quantification

### Q12: How sensitive is rVAF to the NNLS algorithm vs. true least squares with negative coefficients?

**Why it matters**: `room-for-improvement.md` Section 9 documents that the current approach uses least squares followed by gradient descent to enforce non-negativity, rather than a proper NNLS solver (Lawson-Hanson). The gradient descent correction is a post-hoc fix that may not converge to the true NNLS solution. If the difference between proper NNLS and the current approach is clinically significant (>10% relative error in rVAF), this should be prioritized. If the difference is negligible (rVAF changes by <0.01 percentage points), the current approach is adequate.

**How to investigate**: Implement Lawson-Hanson NNLS in Rust (or use the `nnls` crate if available). Run both the current approach and proper NNLS on all variant calls from the validation data. Compare: (a) rVAF differences (absolute and relative); (b) cases where gradient descent does not converge; (c) cases where negative coefficients persist in the current approach; (d) whether proper NNLS changes any clinical call (variant present/absent at 0.1% threshold).

**Research reference**: `room-for-improvement.md` Section 9, `confidence-metrics.md` Section 1.1 for NNLS limitations, `thesis-summary.md` Section "Quantification Approach" for the current model.

---

### Q13: Should we use maximum likelihood instead of least squares?

**Why it matters**: `room-for-improvement.md` Section 9 notes that the least squares model assumes normally distributed residuals with constant variance (homoscedastic), but k-mer counts follow Poisson or negative binomial distributions where variance scales with the mean. A maximum likelihood estimator (MLE) using the correct distributional assumption could provide better rVAF estimates, particularly at low VAF where the count distribution is most non-Gaussian. However, MLE is computationally more expensive and may not converge reliably for the overdetermined systems (31+ k-mers, 2-5 paths) typical in variant detection.

**How to investigate**: Implement Poisson MLE for the path decomposition model. For each path j with coefficient c_j, the expected count at k-mer i is `lambda_i = sum_j(contrib[i,j] * c_j)`, and the log-likelihood is `sum_i(count_i * log(lambda_i) - lambda_i)`. Maximize using iteratively reweighted least squares (IRLS) or L-BFGS-B with non-negativity constraints. Compare rVAF estimates to NNLS on the validation data. If MLE provides meaningfully better estimates (especially at low VAF), implement as the default; if similar, retain NNLS for simplicity.

**Research reference**: `room-for-improvement.md` Section 9, `confidence-metrics.md` Sections 2.3-2.4 for Poisson and negative binomial models.

---

### Q14: How accurate are rVAF confidence intervals from bootstrap?

**Why it matters**: `confidence-metrics.md` Section 3.4 proposes bootstrap CIs for rVAF, but bootstrap methods for NNLS have known issues. Standard bootstrap can undercover when the true coefficient is near zero (the non-negativity boundary), producing asymmetric intervals that may not have the claimed 95% coverage. For MRD monitoring, where the clinical question is "is rVAF trending up or down," the accuracy of CIs directly affects clinical interpretation.

**How to investigate**: Simulation study:
1. Generate synthetic k-mer count data with known true rVAF values (0.01%, 0.1%, 0.5%, 1%, 5%)
2. Add realistic noise (Poisson counts, overdispersion matching empirical estimates)
3. Run bootstrap NNLS with 1000 replicates
4. Measure CI coverage: what fraction of true rVAF values fall within the 95% CI?
5. Compare naive bootstrap vs. parametric bootstrap (resample from fitted Poisson) vs. percentile bootstrap
6. Test whether coverage improves when using the delta method analytical CI from `confidence-metrics.md` Section 3.4

Target: 90-98% coverage at all VAF levels. If coverage falls below 90% at low VAF, investigate alternative CI methods (profile likelihood from `confidence-metrics.md` Section 3.4).

**Research reference**: `confidence-metrics.md` Section 3.4 for bootstrap methodology, analytical CI via delta method, and profile likelihood alternative.

---

## Clinical

### Q15: What is the minimum sequencing depth for reliable 0.1% VAF detection?

**Why it matters**: The thesis used 2000-5000x depth and achieved a detection threshold of ~0.1% VAF. But `sensitivity-landscape.md` Section 2.2 shows that at 2000x, the probability of observing <2 variant k-mers at 0.1% VAF is ~40%. This means the "0.1% threshold" is actually a stochastic boundary where detection becomes unreliable, not a hard floor. Clinical labs need to know the minimum depth to guarantee detection at 0.1% VAF with, say, 95% probability.

**How to investigate**: Sensitivity simulation:
1. For VAF in {0.01%, 0.05%, 0.1%, 0.2%, 0.5%}:
2. For depth in {500, 1000, 2000, 3000, 5000, 10000}:
3. Compute the probability of variant k-mer count >= threshold (using the binomial model from `sensitivity-landscape.md` Section 1.1):
   `P(count >= threshold) = 1 - BinomCDF(threshold-1, depth, VAF * (1 - (k-1)/read_length))`
4. Find the minimum depth where P(count >= threshold) >= 0.95 for each VAF
5. Validate with actual detection rates on samples with varying coverage (subsample high-coverage data to lower depths)

This produces a depth-vs-VAF detection power table that clinical labs can use to set sequencing targets.

**Research reference**: `sensitivity-landscape.md` Section 1.1 for the probability model and examples, Section 2 for VAF range analysis, `thesis-summary.md` for the 2000-5000x depth context.

---

### Q16: How much does UMI collision affect detection in practice?

**Why it matters**: `humid-analysis.md` Section "UMI Collision Problem in Targeted Panels" calculates that with 8-bp UMIs (65,536 possible) at 5000x unique molecular depth, collisions are almost guaranteed. Estimated 17% sensitivity reduction at 0.1% VAF from collision-mediated variant molecule loss. But the calculation uses birthday paradox approximations -- the actual impact depends on the UMI length used in the specific library prep kit, the true library complexity, and the distribution of molecules across target regions.

**How to investigate**: Empirical measurement on the validation data:
1. Extract UMI sequences from the raw FASTQ (before HUMID)
2. At each target region, count the number of distinct UMIs and compute the expected collision rate
3. For known variant positions, count UMI families supporting the variant vs. reference
4. Compare actual variant UMI count to expected count (from VAF * total families)
5. If UMI collision is reducing variant family count by >10%, test whether longer UMIs (12-16 bp, available in newer library kits) would mitigate the issue

Also relevant for Phase 6: if kmerdet implements UMI-aware counting, it should handle collisions (e.g., by also considering read start position, not just UMI sequence, similar to UMI-tools' position-aware clustering).

**Research reference**: `humid-analysis.md` Section "UMI Collision Problem in Targeted Panels" for the collision probability calculations and impact estimates.

---

### Q17: Could kmerdet work for cfDNA fragmentation analysis?

**Why it matters**: Recent literature shows that cfDNA fragment length distributions carry cancer-specific signals independent of mutation detection. Tumor-derived cfDNA fragments are enriched at shorter sizes (90-150 bp) per `kmer-length-tradeoffs.md` Section "cfDNA Fragment Size Interaction". A k-mer-based approach to fragmentation analysis could complement variant detection: even when no mutations are detected, an abnormal fragment length profile could indicate residual disease.

**How to investigate**: Literature review + feasibility analysis:
1. Review DELFI (Cristiano et al., 2019) and similar fragmentation-based approaches
2. Assess whether k-mer count profiles can serve as a proxy for fragment length distribution (shorter fragments produce fewer k-mers per molecule, but more k-mers per base at the fragment ends)
3. If feasible, design a proof-of-concept: compute k-mer-derived fragment metrics from jellyfish databases and correlate with known MRD status

This is a stretch goal that could significantly expand kmerdet's clinical utility without requiring mutation detection at all.

**Research reference**: `kmer-length-tradeoffs.md` Section "cfDNA Fragment Size Interaction" for fragment size distributions and k-mer yield per fragment.

---

### Q18: What is the false negative rate for CHIP variants in the tumor-informed context?

**Why it matters**: `false-positive-analysis.md` Section 4.3 extensively discusses CHIP as a false positive source, but there is a complementary problem: if a patient's tumor happens to carry a mutation at a known CHIP gene (e.g., TP53 R175H), and the patient also develops CHIP at the same gene (but a different position), the tumor-informed filter would correctly match the tumor variant but miss the CHIP variant. If the CHIP clone's VAF rises above the tumor ctDNA VAF, the total detected VAF could be misinterpreted as rising tumor burden. How common is this scenario?

**How to investigate**: Review the CHIP prevalence data from `false-positive-analysis.md` (10-30% of cancer patients, 50.1% have at least one CHIP mutation in reportable genes). Cross-reference the target panel genes with CHIP-affected genes. For overlapping genes, estimate the probability of coincidental CHIP at a different position than the tracked tumor mutation. This is primarily a clinical workflow question, not a software question, but kmerdet could flag variants in known CHIP genes with a warning annotation.

**Research reference**: `false-positive-analysis.md` Section 4.3 for CHIP prevalence and gene overlap data.

---

## Integration

### Q19: Can we stream results from detect to filter without intermediate files?

**Why it matters**: The current pipeline design uses files to pass data between stages: detect writes TSV/JSON, filter reads it. For a 50-target panel, this is ~50 KB of data -- trivial. But for larger panels or multi-k workflows, intermediate files multiply. More importantly, an in-memory pipeline avoids file I/O latency and simplifies the `run` subcommand.

**How to investigate**: Architecture decision. Two approaches:
1. **Internal streaming**: The `run` subcommand calls detect and filter as library functions, passing data in-memory as `Vec<VariantCall>`. No intermediate files. This is the natural Rust approach.
2. **Unix pipe streaming**: `kmerdet detect ... | kmerdet filter ...` using JSONL format for streaming. Each line is one variant call. This is more flexible but requires spawning separate processes.

The key constraint is whether filter needs all detect results at once (for normalization deduplication, which sorts and deduplicates) or can process them one-at-a-time. If deduplication requires all results, approach 1 is necessary. If per-target processing suffices, approach 2 works.

**Research reference**: No specific research doc, but `10-future-improvements.md` Section on Nextflow integration and `development-roadmap.md` Phase 4 `run` subcommand design.

---

### Q20: How should multi-sample results be combined in a longitudinal monitoring workflow?

**Why it matters**: MRD monitoring involves serial blood draws over months or years. Each draw produces a set of rVAF estimates for tracked mutations. Clinical interpretation requires trend analysis: is the patient's tumor burden stable, declining (treatment response), or rising (relapse)? The question is whether kmerdet should provide this trend analysis or leave it to downstream tools.

**How to investigate**: Clinical workflow design with input from oncologists:
1. What format do clinical reports need? (VAF trend plots with dates, statistical significance of trend, categorical call of "rising/stable/declining")
2. What statistical test is appropriate? (Linear regression on log(rVAF) over time? Change-point detection? Simple comparison of consecutive timepoints?)
3. Should kmerdet produce per-timepoint results that a separate tool combines, or should `kmerdet longitudinal` be a first-class command?
4. What happens when a variant is detected at timepoint 1, not detected at timepoint 2, then detected again at timepoint 3? Is timepoint 2 a true negative (variant absent) or a false negative (below detection limit)?

The walking quality metric from Q7 is relevant here: if timepoint 2 has a "troubled walk" for that target, the non-detection is less trustworthy.

**Research reference**: `10-future-improvements.md` Section "Multi-Sample Joint Calling" for longitudinal mode concept, `thesis-summary.md` Section "Clinical Workflow Integration" for the overall clinical context.

---

### Q21: What is the right balance between km-compatible output and extended output?

**Why it matters**: Users migrating from the km/kam pipeline expect specific output columns (Sample, Target, Variant_name, Type, rVAF, Expression, Min_coverage). Adding confidence metrics (QUAL, strand_bias, positional_cv, bootstrap CI) and adaptive filtering annotations (ADAPTIVE_THRESHOLD, PON_FREQ, ML_SCORE, CONTEXT_FLAGS) enriches the output but changes the schema. Clinical pipelines downstream of km may break if columns change.

**How to investigate**: User feedback + pragmatic compromise:
1. **Default output mode**: TSV format with km-compatible columns only. Existing pipelines work unchanged.
2. **Extended output mode** (`--format extended-tsv`): All km columns plus confidence metrics and filter annotations.
3. **VCF mode** (`--format vcf`): Standard VCF with all metrics in INFO/FORMAT fields per `confidence-metrics.md` Section 4.1.
4. **JSON mode**: All available data in structured format.

The question is which mode should be the default. For clinical deployment, extended TSV or VCF is more appropriate (richer data for clinical decision-making). For backward compatibility during migration, km-compatible TSV is safer.

Recommendation: Default to extended TSV with a `--compat` flag for km-compatible format.

**Research reference**: `confidence-metrics.md` Section 4.1 for VCF field mapping, `adaptive-filtering.md` Section 7 for filter annotation columns.

---

## Algorithm Fundamentals

### Q22: Is the DFS walking approach fundamentally the right strategy, or should kmerdet consider alternatives?

**Why it matters**: The entire detection algorithm is built on iterative DFS walking through the k-mer graph. This approach was designed for RNA-seq fusion detection (km's original use case) and adapted for liquid biopsy. Alternative approaches exist in the literature:
- **2-kupl** (Pelletier et al., 2021): Compares k-mer populations between matched samples rather than walking, enabling detection of variants >100 bp. Referenced in `sensitivity-landscape.md` Section 1.3.
- **GeneBits** (2025): Ultra-sensitive tumor-informed ctDNA monitoring using a different k-mer-based framework. Referenced in `sensitivity-landscape.md`.
- **Direct junction k-mer search**: Instead of walking through the entire variant, specifically search for the junction k-mers at known breakpoints (requires tumor-informed targets). Referenced in `sensitivity-landscape.md` Section 5.2.

**How to investigate**: Implement a prototype junction k-mer search for tumor-informed targets:
1. For each known variant, compute the expected junction k-mers (k-mers spanning the variant breakpoint)
2. Query these k-mers directly from the jellyfish database (O(k) queries per variant, no walking needed)
3. If junction k-mers are present at counts above threshold, call the variant present
4. Compare sensitivity to the walking approach

If junction k-mer search achieves comparable or better sensitivity with much simpler code, it could replace walking for tumor-informed detection entirely. Walking would remain for de novo discovery.

**Research reference**: `sensitivity-landscape.md` Sections 1.3 (2-kupl for long indels) and 5.2 (junction k-mer enrichment).

---

### Q23: How do we handle the error rate estimation chicken-and-egg problem?

**Why it matters**: Adaptive thresholds (`walking-improvements.md` Section 2, `adaptive-filtering.md` Section 2) require an estimate of the per-base sequencing error rate. The proposed approach estimates this from non-reference extensions of reference k-mers. But this estimation itself requires knowing which k-mers are reference (which requires walking or at least target sequence decomposition). Additionally, the error rate is not constant: it varies by sequence context (GGC motifs at 5-10x baseline per `false-positive-analysis.md` Section 2.2), position in read (ends are worse per `false-positive-analysis.md` Section 2.1), and possibly by target region.

**How to investigate**: Implement and compare three error rate estimation approaches:
1. **Global estimate**: Compute error rate across all reference k-mers in all targets. Simple, but ignores context dependence.
2. **Per-target estimate**: Compute error rate from reference k-mers within each target. More accurate but noisier (fewer data points per target).
3. **Fixed platform default**: Use 0.001 for Illumina post-dedup. Simplest, but does not adapt to sample quality.

Measure which approach produces the best-calibrated adaptive threshold (Experiment 5).

**Research reference**: `walking-improvements.md` Section 2 for error rate estimation approach, `adaptive-filtering.md` Section 2 for error k-mer distribution modeling, `false-positive-analysis.md` Sections 2.1-2.2 for context-dependent error rates.

---

### Q24: Is there a way to use the known variant sequence (tumor-informed targets) to improve detection without biasing toward expected results?

**Why it matters**: In tumor-informed detection, we know exactly what variant to look for. The current approach ignores this knowledge during walking and graph construction -- it discovers variants de novo, then checks whether they match the expected variant during filtering. This wastes sensitivity: the walking algorithm might prune the expected variant path due to a single low-count k-mer, even though surrounding k-mers strongly support it.

However, using expected variant information during detection risks confirmation bias: the algorithm might "find" the expected variant even when it is not truly present, by lowering thresholds or guiding the walk specifically toward the expected path.

**How to investigate**: Design a two-tier detection strategy:
1. **Tier 1 (de novo)**: Standard detection with no prior knowledge. This is the current approach.
2. **Tier 2 (guided)**: For targets where Tier 1 produces no result, use the expected variant sequence to:
   - Query specific junction k-mers directly (no walking needed)
   - Relax thresholds for k-mers that match the expected variant path
   - Report as a separate confidence tier ("guided detection")

The guided tier uses prior knowledge but reports it explicitly, so the clinical interpretation can account for the different evidence standard. This is analogous to the "multi-pass detection" suggested in `sensitivity-landscape.md` Section 5.1.

**Research reference**: `walking-improvements.md` Section 4 "Guided Walking Using Partial Sequence Matching" for the guided walking concept, `sensitivity-landscape.md` Section 5.1-5.2 for multi-pass and junction k-mer approaches.

---

## Questions Requiring External Input

### Q25: What are the actual UMI lengths and library prep protocols used in clinical settings?

**Why it matters**: The collision probability calculations in `humid-analysis.md` depend heavily on UMI length (8 bp = 65K possible UMIs, 12 bp = 16.7M possible UMIs). Different clinical labs use different library prep kits with different UMI designs. kmerdet needs to handle the range of UMI configurations in practice.

**How to investigate**: Survey clinical collaborators on their specific protocols. Key parameters: UMI length, dual-index vs single-index, insert size distribution, typical input DNA amount (ng cfDNA), expected library complexity (unique molecules per target region).

---

### Q26: What target panels are in clinical use for MRD monitoring?

**Why it matters**: The thesis used patient-specific panels with ~50 targets. Commercial MRD panels (FoundationOne, Guardant, Signatera) use different designs with different numbers of targets, different variant types, and different flanking strategies. kmerdet's performance may vary across panel designs.

**How to investigate**: Collect panel design specifications from clinical collaborators and public sources. Test kmerdet on panel designs with: (a) different target counts (20, 50, 200, 500); (b) different flanking lengths (20, 35, 50, 100 bp); (c) panels targeting repetitive regions or CHIP-prone genes. This feeds into target design validation (Q7, `room-for-improvement.md` Section 10).

---

### Q27: What turnaround time is clinically acceptable for MRD monitoring?

**Why it matters**: The thesis achieves ~6 minutes end-to-end, which is impressive. But clinical workflows have specific turnaround time targets. If same-day reporting is required (8-hour window), 6 minutes is more than fast enough. If near-real-time monitoring is desired (results within 30 minutes of sequencing completion), pipeline optimization matters more. If batch processing overnight is acceptable, performance optimization has lower priority.

**How to investigate**: Clinical workflow consultation. The answer determines whether Phase 5 performance optimization and Phase 6 UMI integration (which eliminates the HUMID step) are clinically motivated or primarily engineering improvements.

---

## Summary Priority Ranking

| Priority | Questions | Impact | Investigation Effort |
|----------|-----------|--------|---------------------|
| High | Q3 (-L cutoff), Q5 (max_stack), Q12 (NNLS accuracy), Q15 (min depth) | Directly affects sensitivity | Moderate |
| High | Q7 (walk failure detection), Q22 (junction k-mer alternative), Q24 (guided detection) | Could fundamentally improve detection | High |
| Medium | Q1 (JF accuracy), Q6 (bidirectional), Q8 (relaxed threshold), Q9 (edge weights) | Incremental improvements | Moderate |
| Medium | Q13 (MLE), Q14 (bootstrap accuracy), Q16 (UMI collision), Q23 (error estimation) | Better quantification and confidence | Moderate |
| Low | Q2 (canonical), Q4 (DashMap), Q10 (>100 paths), Q11 (A*) | Performance/edge cases | Low-Moderate |
| External | Q17 (fragmentation), Q18 (CHIP FN), Q25-Q27 (clinical context) | Future directions | Varies |

---

## References

All research documents in `docs/research/`:
- `thesis-context/thesis-summary.md`, `thesis-context/room-for-improvement.md`
- `sensitivity-and-confidence/sensitivity-landscape.md`, `confidence-metrics.md`, `false-positive-analysis.md`
- `kmer-length-optimization/kmer-length-tradeoffs.md`, `multi-k-strategy.md`
- `error-correction/humid-analysis.md`, `error-correction-alternatives.md`
- `kmer-counting/jellyfish-deep-dive.md`, `counting-alternatives.md`
- `graph-and-walking/walking-improvements.md`, `graph-improvements.md`
- `filtering-strategies/adaptive-filtering.md`, `indel-filtering.md`
- `initial/10-future-improvements.md`
