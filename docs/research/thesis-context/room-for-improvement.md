# Room for Improvement: Systematic Analysis of Thesis-Identified Limitations

This document catalogs every significant limitation identified in the thesis validation study, with technical root cause analysis, severity assessment, and pointers toward solutions. Each limitation is drawn from specific thesis findings and cross-referenced to other research documents where applicable.

---

## Limitation Overview

| # | Limitation | Severity | Impact Area |
|---|-----------|----------|-------------|
| 1 | SNV sensitivity ceiling (77%) | High | Clinical sensitivity |
| 2 | INDEL sensitivity collapse (38% overall, 3% for >7 bp) | Critical | Clinical sensitivity |
| 3 | VAF detection threshold (~0.1%) | High | MRD monitoring depth |
| 4 | KM's restrictive quality filtering | Moderate | Sample rejection rate |
| 5 | No UMI integration into k-mer counting | High | Noise floor / artifact control |
| 6 | Fixed parameters not adapting to sample characteristics | High | Per-sample optimization |
| 7 | HUMID as bottleneck (60-70% runtime) | Moderate | Throughput |
| 8 | Graph construction limitations | Moderate | Variant path accuracy |
| 9 | Quantification model limitations | Moderate | rVAF accuracy |
| 10 | Target design sensitivity | High | Walking success rate |

---

## 1. SNV Sensitivity Ceiling (77% vs 82-95% Alignment-Based)

### What the Thesis Found

The validation study achieved 77% sensitivity for SNVs across the 10-patient cohort. Alignment-based pipelines (BWA-MEM + GATK Mutect2) achieve 90-95% sensitivity on comparable samples (Thesis Ch5). The 13-18 percentage point gap means roughly 1 in 4-5 true SNVs is missed by the k-mer approach.

False negatives occurred primarily at very low VAF (<0.1%) and in regions with high background noise (Thesis Ch5).

### Why It Happens

**Fixed extension thresholds**: The km algorithm uses a fixed `ratio` (cutoff) and `count` (n_cutoff) for all k-mer extensions across all targets and all samples. The extension rule is:

```
threshold = max(sum_of_sibling_counts * ratio, n_cutoff)
```

With `ratio=0.00001` and `count=2`, a variant k-mer is kept only if its count >= max(total_sibling_counts * 0.00001, 2). At low VAF and moderate coverage, the variant k-mer count may be 1 -- below the count=2 threshold -- causing the walk to terminate before reaching the variant path.

**No quality weighting**: All k-mers are treated equally regardless of the base quality scores of the reads that produced them. A k-mer derived from high-quality bases (Q30+) and one from low-quality bases (Q15) receive the same weight. Alignment-based callers incorporate per-base quality into their probabilistic models, giving them an edge for borderline calls.

**K-mer counting noise**: Sequencing errors create spurious k-mers at a rate of approximately 3% per 31-mer (given ~0.1% per-base error rate, probability of at least one error in 31 bases is ~3%). These error k-mers raise the local background, making it harder to distinguish true low-VAF variant k-mers from noise. After HUMID deduplication, PCR error k-mers are reduced but not eliminated.

**No strand-bias information**: Alignment-based callers can assess whether variant-supporting reads come from both strands (a hallmark of true variants). K-mer counting in canonical mode (`-C`) merges both strands, losing this discriminative signal.

### Severity

**High.** A 77% sensitivity ceiling means the tool cannot be used as a standalone clinical assay without alignment-based confirmation. Every ~4th SNV is missed, which in MRD monitoring could mean delayed relapse detection.

### Approaches to Fix

- **Adaptive thresholds**: Adjust `ratio` and `count` per-sample based on median k-mer coverage and coverage coefficient of variation (see `docs/research/initial/10-future-improvements.md`, Adaptive Filtering Parameters)
- **Multi-pass detection**: Run km with multiple parameter sets and take the union with FDR control (Thesis Ch6.5)
- **Quality-weighted k-mer counting**: Weight k-mer counts by base quality scores during counting or at query time
- **Strand-aware counting**: Maintain separate forward/reverse counts; require signal from both strands for high-confidence calls
- **ML-based adaptive thresholding**: Train a classifier using features beyond raw count (see `docs/research/initial/10-future-improvements.md`, ML section)

---

## 2. INDEL Sensitivity Collapse (38% Overall, ~3% for >7 bp)

### What the Thesis Found

INDEL detection sensitivity was dramatically lower than SNV detection:

| INDEL Size | Approximate Sensitivity | Drop vs SNV |
|------------|------------------------|-------------|
| 1-3 bp | Moderate (not separately quantified) | Significant |
| 4-7 bp | Low | Severe |
| >7 bp | ~3% | Near-total failure |
| Overall | 38% | 39 percentage points below SNV |

Insertions >15 bp (approaching k/2 for k=31) had reduced k-mer coverage across the variant junction. Deletions in homopolymer regions generated ambiguous k-mer paths. Complex indels (simultaneous insertion + deletion) sometimes produced paths that the classification logic misidentified. The graph construction phase occasionally failed to extend through INDEL junctions when coverage was low (Thesis Ch5).

### Why It Happens

**K-mer length constraint (fundamental)**: For an insertion of length L using k-mers of length k, the number of k-mers that span the insertion junction is `k - L`. When L >= k (31 bp), zero k-mers span the junction and the variant is undetectable. Even when L < k, the number of spanning k-mers decreases linearly:

| Insertion length | Spanning k-mers (k=31) | Fraction of full signal |
|-----------------|----------------------|----------------------|
| 1 bp (SNV equivalent) | 31 | 100% |
| 5 bp | 26 | 84% |
| 10 bp | 21 | 68% |
| 15 bp (k/2) | 16 | 52% |
| 20 bp | 11 | 35% |
| 25 bp | 6 | 19% |
| 30 bp | 1 | 3% |
| 31+ bp | 0 | 0% |

At low VAF (0.1-0.5%), the reduction in spanning k-mers compounds with the already-low per-k-mer count, often dropping below the count threshold.

**Homopolymer context**: In homopolymer regions (e.g., AAAAAAA), many k-mers are identical or near-identical, creating ambiguity in the extension. The DFS walker may follow multiple paths without being able to determine the correct INDEL length. This is especially problematic for deletions within or adjacent to homopolymer runs.

**Low VAF + INDEL compound effect**: At 0.1% VAF with 5000x coverage, a variant k-mer has an expected count of ~5. For a 15 bp insertion, only 16 of the 31 possible spanning k-mers exist, and each has count ~5. With counting noise, some of these may fall below the count=2 threshold, breaking the walk.

**PCR stutter around INDELs**: PCR amplification of repetitive regions (especially dinucleotide and trinucleotide repeats) introduces stutter artifacts -- reads with spurious insertions or deletions of repeat units. While HUMID deduplication reduces stutter, it does not fully eliminate it because stutter can occur in the first PCR cycle before UMI attachment.

**Classification edge cases in `diff_path_without_overlap`**: The variant classification logic compares reference and alternative paths by walking from both ends to find the divergence point. For complex indels near target boundaries, the divergence region may extend to the edge of the target sequence, causing the overlap-removal logic to misclassify or fail. Specifically:

- When the INDEL is close to the anchor k-mers, the walk may not have enough reference context to properly identify the start/end of the variant
- Simultaneous insertion + deletion events produce paths where both `end_ref > end_var` and `end_var > end_ref` are partially true, falling into the "Complex Indel" catch-all category

### Severity

**Critical.** The 38% overall INDEL sensitivity (and near-zero for >7 bp) represents the single largest clinical limitation. Many actionable cancer mutations are INDELs (NPM1 4bp insertion, FLT3-ITD, BRCA1/2 frameshifts). A tool that misses 62% of INDELs cannot be relied upon for INDEL detection without alignment-based backup.

### Approaches to Fix

- **Per-target k-mer length**: Use shorter k (e.g., k=21) for targets with known large INDELs while keeping k=31 for SNVs (Thesis Ch6.5, `docs/research/initial/10-future-improvements.md`)
- **Automatic k optimization**: For each target, find the minimum odd k (21, 23, ..., 35) such that anchor k-mers are unique in the reference
- **Multi-k strategy**: Run with multiple k values and merge results
- **Extended walking with higher max_stack**: For known INDEL targets, increase maximum DFS depth to allow longer walks through INDEL junctions
- **Homopolymer-aware extension**: Modify the extension threshold in homopolymer regions to be more permissive
- **Improved classification logic**: Handle edge cases in `diff_path_without_overlap` for complex indels near target boundaries

---

## 3. VAF Detection Threshold (~0.1%)

### What the Thesis Found

The validation study established a practical detection threshold of **VAF > 0.1%**. Below this, the signal-to-noise ratio was insufficient for reliable detection. At 0.1% VAF with typical sequencing depths (2000-5000x), variant k-mers had counts of approximately 2-5. Above 0.1% VAF, detection was consistent and rVAF correlated well with expected values (Thesis Ch5).

### Why It Happens

**Sequencing depth limitation**: At 0.1% VAF with 3000x effective coverage (after dedup), a variant produces ~3 supporting reads. Each read generates up to 31 variant k-mers (for an SNV), so the expected variant k-mer count is ~3. With Poisson-distributed sampling noise, some k-mers may have count 1 or 0, breaking the walk.

**K-mer counting noise floor**: Sequencing errors create low-count k-mers at a rate far exceeding true ultra-low-VAF variants. At 0.05% VAF with 3000x coverage, a true variant k-mer has expected count ~1.5 -- indistinguishable from error-derived k-mers without additional information (base quality, strand bias, molecular origin).

**Fixed count thresholds**: The `--count 2` parameter is a hard floor. A variant k-mer with count 1 is discarded regardless of context. Lowering to count=1 would dramatically increase false positives by following error k-mers. There is no way to express "this count-1 k-mer is credible because its neighbors all have count 2-3 and it sits in a known variant context."

**HUMID deduplication ceiling**: HUMID collapses UMI families, but the number of unique molecules is limited by library complexity. If only 3 unique molecules carry the variant at 0.1% VAF, no amount of sequencing depth can increase the molecule count beyond 3.

### Severity

**High.** MRD monitoring ideally requires detection down to 0.01% VAF for early relapse detection. The 0.1% threshold means the tool can only detect ctDNA when the tumor burden is already somewhat elevated, potentially missing the earliest signs of relapse. Competing technologies (duplex sequencing with error correction, CAPP-Seq) claim detection limits of 0.01-0.02% VAF.

### Approaches to Fix

- **UMI-aware k-mer counting**: Count unique molecules per k-mer rather than raw read counts; a k-mer supported by 2 unique molecules is far more credible than one supported by 20 reads from a single PCR family (see `docs/research/initial/10-future-improvements.md`)
- **Adaptive per-sample thresholds**: Set `count` based on achieved sequencing depth and deduplication rate (Thesis Ch6.5)
- **Context-aware thresholding**: Allow lower counts in known variant regions while maintaining strict thresholds elsewhere
- **Panel-of-normals**: Build a database of k-mer paths from healthy controls; paths seen in normals are artifacts, allowing lower detection thresholds for paths not seen in normals (see `docs/research/initial/10-future-improvements.md`)
- **Statistical model**: Replace hard count threshold with a probabilistic model incorporating coverage, error rate, and expected VAF

---

## 4. KM's Restrictive Quality Filtering Rejecting Valid Samples

### What the Thesis Found

The km tool's fixed filtering parameters (particularly `n_cutoff=500` in the default configuration for RNA-seq, but even the liquid biopsy parameters `count=2, ratio=0.00001`) sometimes rejected valid samples or valid targets within samples. Targets with legitimate but unusual coverage profiles (very high or very low coverage relative to median) could produce no results or false results. The thesis noted that km's behavior at parameter boundaries was not always predictable, with some samples requiring manual parameter adjustment to produce results (Thesis Ch6).

### Why It Happens

**Hard-coded defaults designed for RNA-seq**: The km tool was originally designed for RNA-seq analysis of leukemia mutations, where expression levels are typically high (hundreds to thousands of reads). The default `n_cutoff=500` is far too high for liquid biopsy, where variant k-mer counts may be 2-10. While the kam pipeline overrides this to `count=2`, the override is global -- the same value is used for every target in every sample.

**No per-target or per-sample adaptation**: A 50-target panel may have targets with coverage ranging from 500x to 10,000x depending on capture efficiency, GC content, and library complexity. A single `count` threshold cannot be optimal for all targets. High-coverage targets could use stricter thresholds (reducing false positives), while low-coverage targets need more permissive thresholds.

**Coverage coefficient of variation**: The thesis found that coverage variability within a sample (CV) was significant. Samples with high CV had targets where the k-mer counts were far below the sample median, causing those targets to fall below the threshold even when the variant was present.

### Severity

**Moderate.** This primarily affects samples with unusual coverage profiles rather than the majority. However, in a clinical setting, rejecting a valid sample or missing a variant due to parameter rigidity is a significant failure mode -- it introduces a systematic bias against certain targets and sample types.

### Approaches to Fix

- **Adaptive thresholds per target**: Compute coverage statistics per target and set `count` proportional to local median coverage (e.g., count = max(2, median_coverage * 0.001)) (see `docs/research/initial/10-future-improvements.md`, Adaptive Filtering Parameters)
- **Coverage-aware ratio**: High-coverage targets can use a higher ratio (more specific); low-coverage targets use a lower ratio (more sensitive)
- **Automatic parameter grid search**: For each target, try multiple parameter combinations and select the one that produces the most biologically plausible result (requires a plausibility scoring function)
- **Fallback mechanism**: If primary parameters produce no result for a target, automatically retry with more permissive parameters and flag the result as low-confidence

---

## 5. No UMI Integration into K-mer Counting

### What the Thesis Found

The kam pipeline uses HUMID as a separate preprocessing step: UMI-based deduplication collapses read families, then jellyfish counts k-mers from the deduplicated output. The UMI information is discarded after deduplication and is not available during k-mer counting or walking. The thesis identified this as a significant missed opportunity for improved artifact detection (Thesis Ch6.5).

### Why It Happens

**Architectural separation**: HUMID produces deduplicated FASTQ files; jellyfish counts k-mers from those files. Neither tool is aware of the other's internals. The UMI group membership is lost at the FASTQ output boundary.

**Jellyfish's design**: Jellyfish is a general-purpose k-mer counter. It has no concept of read groups, UMI families, or molecular identity. It counts raw k-mer occurrences with no metadata.

**Historical precedent**: The km tool was designed for RNA-seq without UMIs. UMI support was added later at the pipeline level (via HUMID) rather than at the algorithmic level.

### Why It Matters

A k-mer at a true variant site supported by 3 unique molecules (3 different UMIs) is far more credible than one supported by 20 reads all from a single PCR duplicate family. With standard counting, both scenarios produce very different raw counts but the same molecular evidence.

Conversely, a PCR error that occurs in the first amplification cycle is propagated to many reads, creating a k-mer with high count but only 1 supporting molecule. Standard counting treats this as a strong signal; molecule-level counting would correctly identify it as a single-molecule observation.

| Scenario | Read count | Molecule count | True signal? |
|----------|-----------|---------------|-------------|
| True variant, 3 UMI families | 20 | 3 | Yes |
| PCR error, 1 UMI family | 20 | 1 | No |
| True variant, low coverage | 3 | 3 | Yes |
| Sequencing error | 1 | 1 | Probably not |

### Severity

**High.** Without molecule-level counting, the tool cannot distinguish PCR artifacts from true variants at the k-mer level. This directly contributes to both false positives (high-count artifacts) and the inability to lower the VAF detection threshold (must use conservative thresholds to avoid artifacts).

### Approaches to Fix

- **UMI-aware k-mer counting**: Implement a custom k-mer counter that tracks `molecule_count[K] = number of distinct UMI families containing k-mer K` alongside `read_count[K]` (see `docs/research/initial/10-future-improvements.md`)
- **Approximate molecule counting**: Use HyperLogLog or Bloom filters for memory-efficient approximate unique molecule counting per k-mer
- **Streaming integration**: Process UMI grouping and k-mer counting in a single pass, eliminating the intermediate FASTQ files
- **Dual-threshold model**: Require both `read_count >= R` and `molecule_count >= M` for k-mer extension; e.g., `read_count >= 2 AND molecule_count >= 2`

---

## 6. Fixed Parameters Not Adapting to Sample Characteristics

### What the Thesis Found

The kam pipeline uses the same `ratio=0.00001` and `count=2` for all samples regardless of coverage depth, library complexity, deduplication rate, or GC bias profile. The thesis suggested that at >5000x median coverage, `count >= 5` would be appropriate, while at <500x, `count <= 2` is necessary. High coverage coefficient of variation (>0.5) should halve the effective ratio (Thesis Ch6.5).

### Why It Happens

**One-size-fits-all design**: The km algorithm takes parameters as command-line arguments and applies them uniformly to every k-mer extension decision. There is no mechanism to inspect the k-mer database characteristics before setting thresholds.

**No pre-analysis step**: The pipeline does not compute sample-level statistics (median coverage, CV, dedup rate) before running km. Without this information, adaptive thresholding is not possible.

### Severity

**High.** Fixed parameters create a fundamental tradeoff: parameters tuned for low-coverage samples produce excess false positives on high-coverage samples, while parameters tuned for high-coverage samples miss true variants on low-coverage samples. In a clinical lab processing samples with variable coverage (due to cfDNA input variation, extraction efficiency, etc.), this leads to inconsistent performance across samples.

### Approaches to Fix

- **Pre-analysis coverage estimation**: Before running km, query the jellyfish database for reference k-mer counts to estimate median coverage and CV
- **Per-target thresholds**: Compute local coverage at each target's reference k-mers and set count/ratio proportionally
- **Recommended parameter table** (from thesis):

| Median coverage | Recommended count | Recommended ratio | Rationale |
|----------------|-------------------|-------------------|-----------|
| <500x | 2 | 0.00001 | Maximize sensitivity at low depth |
| 500-2000x | 2-3 | 0.00005 | Balanced |
| 2000-5000x | 3-5 | 0.0001 | Reduce noise at moderate depth |
| >5000x | 5-10 | 0.001 | Can afford strict thresholds |

- **Coverage CV adjustment**: If CV > 0.5, halve the effective ratio to avoid losing low-coverage targets
- **Deduplication rate as complexity proxy**: High dedup rate (>80%) indicates low library complexity; use more conservative thresholds since effective molecule count is low

---

## 7. HUMID as Bottleneck (60-70% of Preprocesssing Runtime)

### What the Thesis Found

HUMID deduplication consumed approximately 2 minutes of the ~6-minute pipeline, representing ~33% of total wall-clock time. However, for the preprocessing-only portion (HUMID + jellyfish = 3.5 minutes), HUMID represents ~57% of preprocessing time. The thesis noted that HUMID is single-threaded and writes intermediate FASTQ files to disk, creating an I/O bottleneck (Thesis Ch5, Ch6.5). The thesis estimated ~30% wall-clock improvement from overlapping deduplication and counting.

### Why It Happens

**Single-threaded execution**: HUMID processes reads sequentially, not utilizing multiple CPU cores. Modern servers have 16-64+ cores, but HUMID uses only one.

**Disk I/O overhead**: HUMID reads FASTQ input and writes deduplicated FASTQ output. Jellyfish then re-reads the deduplicated FASTQ. For a typical liquid biopsy sample, the intermediate deduplicated FASTQ files can be >10 GB, creating substantial disk I/O.

**Separate process**: HUMID runs as a completely separate binary before jellyfish. There is no pipelining or streaming between the two.

### Severity

**Moderate.** The absolute time (2 minutes) is small enough that it does not prevent same-day reporting. However, it represents a significant fraction of total runtime and is easily improvable -- making it a high-value optimization target relative to effort.

### Approaches to Fix

- **Streaming integration**: Pipe HUMID output directly into jellyfish count via Unix named pipes (FIFOs), eliminating disk I/O for intermediate files (see `docs/research/initial/10-future-improvements.md`)
- **In-process deduplication**: Implement UMI-aware dedup in Rust within kmerdet, feeding deduplicated reads directly to the k-mer counter in memory
- **Parallel deduplication**: Implement multi-threaded UMI grouping using a concurrent hash table (similar to jellyfish's approach)
- **Skip deduplication when UMI-aware counting is used**: If k-mer counting is molecule-aware (Limitation #5 fix), separate deduplication becomes unnecessary

---

## 8. Graph Construction Limitations

### What the Thesis Found

The graph construction phase uses fixed edge weights: reference edges get weight 0.01, non-reference edges get weight 1.0. Dijkstra's algorithm then finds shortest paths from source ("BigBang") to sink ("BigCrunch"). This reference-biased weighting causes the pathfinder to prefer the reference path even when variant evidence is strong, and the fixed weights do not account for actual k-mer count differences between paths (Thesis Ch6).

The graph construction occasionally failed to extend through INDEL junctions when coverage was low. The maximum DFS depth (`max_stack=500`) and maximum branching points (`max_break=10`) can be reached before completing the walk in complex regions.

### Why It Happens

**Fixed edge weights**: The 0.01 vs 1.0 weighting is a heuristic, not derived from the data. A non-reference edge with 10,000 supporting k-mer counts gets the same weight (1.0) as one with 2 counts. This throws away the coverage information during pathfinding, even though it was used during the walking phase.

**Reference-biased pathfinding**: By giving reference edges 100x lower weight, the algorithm strongly prefers the reference path. This is appropriate for finding the reference as a baseline, but it means the pathfinding step is not well-calibrated for discovering alternative paths -- those are found only by removing reference edges entirely and re-running, a somewhat crude approach.

**DFS depth limits**: For complex variants (large INDELs, compound mutations), the k-mer walk may need to traverse many nodes. The `max_stack=500` limit was chosen to prevent runaway computation in repetitive regions, but it also prevents legitimate discovery of complex variants.

### Severity

**Moderate.** The fixed-weight scheme works adequately for simple SNVs (which make up the majority of clinical targets) but contributes to INDEL failures and can produce suboptimal path selection in complex regions. It also prevents the system from using coverage information to weight path confidence.

### Approaches to Fix

- **Coverage-weighted edges**: Set edge weight inversely proportional to k-mer count: `weight = 1 / log(count + 1)`. Higher-count k-mers get lower-weight edges, guiding pathfinding toward well-supported paths
- **Bidirectional shortest paths**: Instead of Dijkstra from source and sink separately, use a bidirectional search that meets in the middle, improving efficiency for long targets
- **Adaptive DFS limits**: Increase `max_stack` for known INDEL targets; decrease for simple SNV targets
- **All-paths enumeration**: Instead of shortest paths only, enumerate all paths up to a coverage threshold and score them by total supporting evidence
- **Remove reference bias from initial pathfinding**: Use uniform weights for initial graph exploration, then classify paths relative to the reference afterward

---

## 9. Quantification Model Limitations

### What the Thesis Found

The NNLS/linear regression model for rVAF estimation uses least squares to decompose observed k-mer counts into per-path contributions. The thesis noted that negative coefficients can occur (biologically impossible -- a path cannot have negative expression) and must be iteratively removed via gradient descent refinement. The model assumes uniform counting noise, which does not hold in practice due to GC bias, PCR amplification bias, and sequence-context-dependent error rates (Thesis Ch5, Ch6.5).

### Why It Happens

**Least squares assumptions**: The model minimizes `||contrib @ coef - counts||^2`, which assumes:
1. Normally distributed residuals (counting noise is actually Poisson or negative binomial)
2. Homoscedastic errors (variance is actually proportional to count -- higher-count k-mers have higher absolute noise)
3. Independent observations (adjacent k-mers share reads and are highly correlated)

**Negative coefficient handling**: When least squares produces negative coefficients, the current approach uses gradient descent to push them to zero. This is a post-hoc fix rather than a principled approach. Proper non-negative least squares (NNLS) algorithms exist (e.g., Lawson-Hanson) but are not used.

**GC and amplification bias**: K-mers with different GC content are amplified and sequenced at different rates. A path through a GC-rich region may have systematically lower counts than expected, biasing the coefficient estimate downward.

### Severity

**Moderate.** For clinical purposes, exact rVAF quantification is less important than detection (present/absent). However, inaccurate rVAF can affect MRD monitoring over time -- if longitudinal rVAF estimates are noisy, it becomes harder to distinguish a rising trend (relapse) from measurement noise. The thesis reported good correlation between rVAF and expected VAF above 0.1%, suggesting the model works adequately in its operating range despite theoretical limitations.

### Approaches to Fix

- **Maximum likelihood estimation**: Replace least squares with a Poisson or negative binomial likelihood model that accounts for count-proportional variance (Thesis Ch6.5)
- **Proper NNLS solver**: Use Lawson-Hanson NNLS directly instead of least squares followed by gradient descent correction
- **GC bias correction**: Estimate and correct for GC-dependent count bias before quantification
- **Weighted least squares**: Weight each k-mer's contribution inversely by its expected variance: `weight[i] = 1 / count[i]`
- **Bootstrap confidence intervals**: Resample k-mer counts to estimate uncertainty in rVAF, providing clinically useful confidence intervals

---

## 10. Target Design Sensitivity

### What the Thesis Found

Target design quality significantly affected walking success. Poorly designed targets with non-unique anchor k-mers (first and last k-mers of the target sequence) caused walking failures because the DFS could not establish a clear source-to-sink path. Flank length (35 bp on each side) affected the number of reference k-mers available for context, and insufficient flanking could cause the variant to be too close to the target boundary for proper classification (Thesis Ch5, Ch6).

### Why It Happens

**Non-unique anchor k-mers**: The first and last k-mers of the target serve as the source ("BigBang") and sink ("BigCrunch") of the graph. If these k-mers appear elsewhere in the genome, the walk may follow paths through unrelated genomic regions, producing spurious variant calls or failing to converge. This is particularly problematic for:
- Targets in repetitive regions (SINEs, LINEs, segmental duplications)
- Targets near pseudogenes with high sequence similarity
- Short targets where flank length is insufficient to reach unique k-mers

**Fixed flank length**: The `--flank 35` parameter adds 35 bp of reference context on each side. For some genomic regions, 35 bp is insufficient to reach a unique k-mer. The optimal flank length depends on the local repetitiveness of the genome.

**No validation of target quality**: The pipeline does not check whether anchor k-mers are unique before running km. A non-unique anchor silently produces a bad result rather than raising an error.

### Severity

**High.** Target design failures are silent -- the pipeline produces output, but the output may be wrong or missing. In a clinical setting, this means certain patient mutations may be systematically undetectable due to their genomic context, without any warning to the operator. This creates a hidden bias in the assay's sensitivity profile.

### Approaches to Fix

- **Anchor uniqueness validation**: Before running km, check that the first and last k-mers of each target are unique in the reference genome; flag or reject targets with non-unique anchors
- **Adaptive flank extension**: Automatically extend flanking regions until unique anchor k-mers are found (up to a maximum, e.g., 200 bp)
- **Per-target k optimization**: For targets where k=31 does not produce unique anchors, try larger k values (33, 35, etc.) until uniqueness is achieved (see `docs/research/initial/10-future-improvements.md`)
- **Target quality score**: Compute a quality score for each target based on anchor uniqueness, GC content, homopolymer frequency, and k-mer complexity; report this score alongside results
- **Multi-anchor strategy**: Use multiple anchor points (not just first/last k-mer) to establish the graph source and sink, providing redundancy against non-unique anchors

---

## Cross-Limitation Interactions

Several limitations interact and compound each other:

| Interaction | Effect |
|-------------|--------|
| Fixed thresholds (#6) + Low VAF (#3) | Fixed count=2 cannot be lowered for high-coverage samples where count=1 might be valid |
| No UMI integration (#5) + VAF threshold (#3) | Without molecule-level counting, cannot distinguish count-2 true variant from count-2 artifact, keeping threshold high |
| INDEL collapse (#2) + Fixed k (#6) | A per-target k would directly address large INDEL sensitivity |
| HUMID bottleneck (#7) + No UMI integration (#5) | Solving #5 (UMI-aware counting) eliminates #7 entirely |
| Graph weights (#8) + INDEL collapse (#2) | Fixed weights contribute to INDEL walk failures at low coverage |
| Target design (#10) + SNV ceiling (#1) | Some missed SNVs are actually target design failures, not algorithm failures |
| Quantification (#9) + VAF threshold (#3) | Noisy rVAF estimates at low VAF make it hard to set a reliable detection threshold |

---

## Priority Ranking for kmerdet Implementation

Based on clinical impact, implementation effort, and interaction effects:

| Priority | Limitation | Rationale |
|----------|-----------|-----------|
| P1 | #6 Adaptive parameters | Low effort, high impact, unlocks improvements for #1, #3, #4 |
| P1 | #10 Target validation | Low effort, eliminates silent failures |
| P1 | #2 INDEL sensitivity (per-target k) | Medium effort, addresses the most critical clinical gap |
| P2 | #5 UMI-aware counting | High effort, but solves #5, #7, and partially #3 simultaneously |
| P2 | #3 VAF threshold (statistical model) | Medium effort, directly improves MRD monitoring depth |
| P2 | #1 SNV sensitivity (multi-pass, quality weighting) | Medium effort, closes gap with alignment-based tools |
| P3 | #8 Graph construction (coverage-weighted edges) | Medium effort, moderate impact |
| P3 | #9 Quantification model (MLE) | Medium effort, moderate impact |
| P3 | #7 HUMID bottleneck | Resolved automatically if #5 is implemented |
| P3 | #4 Restrictive filtering | Resolved automatically if #6 is implemented |

---

## References

- Thesis Chapter 4: Validation Study Design and Methods
- Thesis Chapter 5: Results and Performance Analysis
- Thesis Chapter 6: Discussion, Limitations, and Future Work
- Thesis Chapter 6.5: Recommendations for Future Work
- `docs/research/initial/01-km-tool-analysis.md` -- km algorithm parameters and defaults
- `docs/research/initial/03-kmer-variant-detection-algorithm.md` -- Algorithm details and limitations
- `docs/research/initial/04-liquid-biopsy-context.md` -- Clinical context and challenges
- `docs/research/initial/05-kam-tool-analysis.md` -- Pipeline architecture
- `docs/research/initial/10-future-improvements.md` -- Future improvements and priority ranking
