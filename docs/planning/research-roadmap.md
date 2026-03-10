# Research Roadmap: Experiments and Analyses

This document defines a sequence of experiments to systematically improve kmerdet's detection performance, informed by findings from all research documents. Each experiment answers specific questions, has clear success criteria, and builds on prior results.

---

## Dependency Graph

```
Experiment 1 (Baseline)
    |
    +--> Experiment 2 (Bottleneck Analysis)
    |        |
    |        +--> Experiment 5 (Adaptive Thresholds)
    |        +--> Experiment 8 (Error Correction)
    |
    +--> Experiment 3 (K-mer Length Sweep)
    |        |
    |        +--> Experiment 4 (Threshold Optimization)
    |
    +--> Experiment 7 (INDEL Normalization)
    |
    +--> Experiment 6 (Confidence Calibration)
    |
    +--> Experiment 9 (Graph Weighting)
    |
    +--> Experiment 10 (Panel-of-Normals)
```

Experiments 1 is the prerequisite for all others. Experiments 2, 3, 6, 7, 9, and 10 can proceed in parallel once baseline is established. Experiments 4, 5, and 8 depend on earlier results as shown.

---

## Timeline

| Phase | Experiments | Prerequisites | Estimated Duration |
|-------|------------|---------------|-------------------|
| Foundation | 1 (Baseline) | Working detect + filter subcommands | 2-3 weeks |
| Parallel Wave 1 | 2, 3, 7 | Experiment 1 complete | 3-4 weeks |
| Parallel Wave 2 | 4, 6, 9 | Experiments 1, 2, 3 complete | 3-4 weeks |
| Parallel Wave 3 | 5, 8, 10 | Experiments 2, 4 complete | 4-6 weeks |
| Total | All 10 | | 12-17 weeks |

---

## Experiment 1: Baseline Reproduction

### Questions Answered
- Does kmerdet match km output on the same inputs?
- What is the exact performance baseline for kmerdet before any improvements?
- Are there any bugs in the reimplementation that cause divergence from km?

### Method
1. Obtain the 10-patient validation dataset used in the thesis (HUMID-deduplicated FASTQ, jellyfish .jf databases, target FASTA files, reference mutation files).
2. Run kmerdet `detect` on each patient's .jf database with the thesis parameters: `--count 2 --ratio 0.00001 --k 31`.
3. Run kmerdet `filter` in reference mode against each patient's known mutations.
4. Compare output against the original km/kmtools output field-by-field: variant name, variant type, rVAF, expression, min_coverage.

### Reference Documents
- `thesis-summary.md`: Expected results (77% SNV sensitivity, 38% INDEL sensitivity, ~0.1% VAF threshold)
- `room-for-improvement.md`: Known limitations to verify kmerdet reproduces

### Requirements
- Working `detect` and `filter` subcommands (Phase 1-2 of development)
- Test data access (10-patient .jf files + targets + reference mutations)
- km/kmtools installed for generating ground truth comparison output

### Success Criteria
- rVAF values match within 1% absolute tolerance (e.g., km reports 0.15%, kmerdet reports 0.14-0.16%)
- Same variants detected and classified with same variant type labels
- Same variants pass/fail filter in reference mode
- SNV sensitivity matches thesis at 77% (+/- 2 percentage points)
- INDEL sensitivity matches thesis at 38% (+/- 3 percentage points)
- Wall-clock time competitive with km (~2 min for walking phase on 50-target panel)

### Outputs
- Per-patient concordance table (matching/mismatching calls)
- Summary statistics table (sensitivity, specificity, rVAF correlation)
- List of any discrepancies with root cause analysis
- Baseline performance numbers for comparison in subsequent experiments

---

## Experiment 2: Sensitivity Bottleneck Analysis

### Questions Answered
- Which pipeline stage loses the most true variants?
- Is the primary bottleneck in walking (threshold too aggressive), counting (k-mers dropped by -L 2), graph construction (paths not connected), or filtering (representation mismatch)?
- What is the theoretical maximum sensitivity achievable with perfect parameters at each stage?

### Method
1. For each patient, identify all expected variants from the reference mutation file.
2. Instrument each pipeline stage with detailed logging:
   - **Counting stage**: For each expected variant, query all k variant k-mers from the .jf database. Record how many have count >= 2, how many have count 1 (lost to -L 2), how many have count 0. Reference: `jellyfish-deep-dive.md` Section on -L flag impact at low VAF.
   - **Walking stage**: Log every extension decision. For failed extensions, record the k-mer, its count, the threshold applied, and the sibling counts. Reference: `walking-improvements.md` Section 1 on threshold analysis and `sensitivity-landscape.md` Section 3.3 on walking as signal loss point.
   - **Graph stage**: After graph construction, check path connectivity for each expected variant path. Record whether the variant path is fully connected source-to-sink. Reference: `graph-improvements.md` Section 1 on binary weight limitations and `sensitivity-landscape.md` Section 3.4 on disconnected paths.
   - **Classification stage**: Check whether discovered paths are correctly classified. Record misclassifications (e.g., INDEL classified as Complex). Reference: `sensitivity-landscape.md` Section 3.5 on classification edge cases.
   - **Filter stage**: Check whether correctly detected variants match in the filter. Record INDEL representation mismatches. Reference: `indel-filtering.md` on the representation problem.
3. Aggregate across all patients to produce per-stage loss counts.

### Reference Documents
- `sensitivity-landscape.md`: Stage-by-stage analysis of where variants are lost, including probability calculations for k-mer survival at various VAF levels
- `walking-improvements.md`: Threshold analysis showing how fixed ratio/count parameters filter variants at different coverage levels
- `graph-improvements.md`: Path connectivity requirements and how single missing k-mers break paths
- `indel-filtering.md`: INDEL representation mismatch as a filter-stage loss source

### Requirements
- Experiment 1 complete (baseline established)
- Instrumented build of kmerdet with verbose per-stage logging
- Known ground truth variants for all 10 patients

### Success Criteria
- Every false negative is attributed to a specific pipeline stage
- Loss percentages sum to the total false negative rate within rounding error
- Identify the single largest loss source (expected: walking threshold based on `sensitivity-landscape.md` Section 3.3)
- Produce a Sankey-style flow diagram showing variants entering and surviving each stage

### Outputs
- Per-stage variant survival table (variants entering / surviving / lost at each stage)
- Root cause classification for each lost variant (threshold, connectivity, misclassification, representation mismatch, below VAF floor)
- Ranked list of loss sources by impact
- Recommendations for which stage to optimize first

---

## Experiment 3: K-mer Length Sweep

### Questions Answered
- What is the optimal k for SNV detection? For INDEL detection?
- How does per-k-mer coverage change across k values at 0.1% VAF?
- Is the sensitivity improvement from shorter k for INDELs as large as theory predicts?
- Does k=21 introduce specificity problems (non-unique anchor k-mers)?

### Method
1. Generate jellyfish databases at k=21, 25, 31, 37, 43 from the same HUMID-deduplicated FASTQ for each patient:
   ```bash
   for k in 21 25 31 37 43; do
       jellyfish count -m $k -s 100M -t 4 -C -L 2 dedup.fastq -o counts_k${k}.jf
   done
   ```
   Reference: `kmer-length-tradeoffs.md` Section "Jellyfish Multi-k Considerations" for resource costs.

2. Run kmerdet `detect` at each k value for all 10 patients.

3. For each k, compute:
   - SNV sensitivity and specificity
   - INDEL sensitivity by size category (1-3 bp, 4-7 bp, 8-15 bp, 16-30 bp, >30 bp)
   - Per-k-mer coverage at variant sites (expected vs. observed)
   - Number of targets where anchor k-mers are non-unique

4. Compute the multi-k union (detect at any k) and intersection (detect at all k) sensitivity.

### Reference Documents
- `kmer-length-tradeoffs.md`: Theoretical framework for k-mer uniqueness vs. sensitivity, expected counts per k at 0.1% VAF, fragment length interaction (cfDNA at 167 bp loses 15% of k-mers going from k=21 to k=43)
- `multi-k-strategy.md`: Four approaches to multi-k detection (independent union, consensus voting, adaptive per-target, Bayesian), expected improvement estimates (SNVs +5-10%, short INDELs +15-25%, long INDELs +20-40%)
- `room-for-improvement.md` Section 2: INDEL sensitivity collapse root causes, k-mer length constraint table showing spanning k-mers decrease linearly with insertion length

### Requirements
- Experiment 1 complete (baseline at k=31)
- Disk space for 5 .jf databases per patient (~2 GB per patient)
- Compute time for 50 detection runs (10 patients x 5 k values)

### Success Criteria
- Identify optimal k for SNVs (expected: k=31, possibly k=25 at very low VAF per `kmer-length-tradeoffs.md` Section on 0.05% VAF)
- Demonstrate >20% INDEL sensitivity improvement at k=21 vs k=31 for INDELs >7 bp (from ~3% to >25% per `multi-k-strategy.md` expected improvements)
- Quantify the specificity cost of k=21 (expected: minimal for targeted panels per `kmer-length-tradeoffs.md` uniqueness table showing k=21 is "essentially unique")
- Determine whether multi-k union improves overall sensitivity enough to justify 3x resource cost

### Outputs
- Sensitivity x k-value matrix by variant type
- Specificity x k-value analysis
- Per-target anchor uniqueness report at each k
- Recommendation for default k value and multi-k strategy
- Resource cost measurements (time, memory, disk) per k

---

## Experiment 4: Threshold Optimization

### Questions Answered
- What is the sensitivity/specificity tradeoff curve for ratio and count parameters?
- Are the thesis defaults (ratio=0.00001, count=2) optimal for the validation data?
- Is there a single parameter set that dominates, or is the optimal set coverage-dependent?
- What is the sensitivity at 0.05% VAF with optimized thresholds?

### Method
1. Using k=31 (or the optimal k from Experiment 3), run kmerdet `detect` across a grid of parameters:
   - ratio in {0.000001, 0.00001, 0.0001, 0.001, 0.01}
   - count in {1, 2, 3, 5, 10}
   - Total: 25 parameter combinations per patient

2. For each combination, compute:
   - SNV sensitivity and false positive count
   - INDEL sensitivity and false positive count
   - True positive rate at each VAF bin (>1%, 0.5-1%, 0.1-0.5%, <0.1%)

3. Construct ROC-like curves: sensitivity vs. false positive rate as parameters vary.

4. Identify the Pareto-optimal parameter sets (no other set achieves both higher sensitivity and lower FP).

### Reference Documents
- `adaptive-filtering.md`: Depth-proportional count thresholds, recommended depth tiers table (ultra-low <500x through ultra-deep >10000x), GC content and homopolymer adjustments
- `sensitivity-landscape.md` Section 1.1: Analysis of how fixed thresholds interact with coverage and VAF, probability calculations showing ~40% of 0.1% VAF variants have k-mer count <2 at 2000x coverage
- `room-for-improvement.md` Section 6: Recommended parameter table by coverage tier, thesis suggestion that count=5 is appropriate above 5000x coverage
- `walking-improvements.md` Section 1: How ratio threshold at low VAF is dominated by the absolute minimum count, making ratio nearly irrelevant for ctDNA

### Requirements
- Experiment 1 and 3 complete
- 250 detection runs (10 patients x 25 parameter combinations)
- Ground truth variant calls with VAF estimates

### Success Criteria
- Identify parameter set(s) achieving >85% SNV sensitivity with zero false positives on the validation data (current: 77%)
- Demonstrate that count=1 with appropriate ratio can recover variants lost at count=2 without unacceptable FP increase
- Show that coverage-stratified parameters outperform any single fixed parameter set
- Produce a recommended parameter table indexed by coverage tier

### Outputs
- 25-cell grid: sensitivity at each (ratio, count) combination
- ROC curves for SNV and INDEL detection
- Pareto frontier of optimal parameter sets
- Coverage-stratified parameter recommendations
- False positive analysis at each parameter set

---

## Experiment 5: Adaptive Threshold Validation

### Questions Answered
- Does coverage-adaptive thresholding improve sensitivity without increasing false positives?
- Is the Poisson-based threshold model (`walking-improvements.md`) better than fixed thresholds?
- Does per-target threshold adaptation (using local coverage) outperform per-sample adaptation?

### Method
1. Implement the adaptive threshold from `walking-improvements.md` Section 2:
   - Estimate local coverage from reference k-mer counts in a sliding window
   - Estimate error rate from non-reference extension counts
   - Set threshold at the Poisson 99th percentile of error k-mer distribution
   - Floor at min_threshold=2

2. Implement per-target threshold adaptation from `adaptive-filtering.md` Section 4:
   - Compute GC content, linguistic complexity, and max homopolymer length per target
   - Adjust threshold based on target annotations
   - Apply the depth-tier table from `adaptive-filtering.md` Section 3

3. Run three configurations on all 10 patients:
   - A: Fixed thresholds (thesis defaults)
   - B: Per-sample adaptive thresholds (based on median coverage)
   - C: Per-target adaptive thresholds (based on local coverage and context)

4. Compare sensitivity and false positive rates.

### Reference Documents
- `adaptive-filtering.md`: Full framework for adaptive thresholding, including error k-mer distribution modeling (Gamma distribution from Quake), GC bias correction, and homopolymer-specific adjustments
- `walking-improvements.md` Section 2: Poisson-based statistical model for setting thresholds, formula for computing threshold from coverage + error rate + false extension rate
- `room-for-improvement.md` Section 6: Thesis finding that fixed parameters create systematic bias against low-coverage targets
- `confidence-metrics.md` Section 2.1-2.4: Statistical models (binomial, Poisson, negative binomial) for k-mer count significance testing

### Requirements
- Experiments 1, 2, and 4 complete
- Adaptive threshold implementation in kmerdet walker
- Per-target annotation computation

### Success Criteria
- Per-target adaptive thresholds achieve sensitivity >= fixed thresholds on every patient (no patient is worse)
- Overall sensitivity improves by >= 5 percentage points for SNVs (target: 82%+)
- False positive rate remains at zero on the validation data
- Demonstrate that the Poisson model correctly separates error k-mers from true variant k-mers (validated by Experiment 2 instrumentation data)

### Outputs
- Three-way comparison table: fixed vs. per-sample vs. per-target sensitivity
- Per-target threshold values and their derivation
- Calibration plot: predicted vs. observed error k-mer counts
- Analysis of which targets benefit most from adaptive thresholds

---

## Experiment 6: Confidence Metric Calibration

### Questions Answered
- Do Phred-scaled QUAL scores predict true vs. false positive variant calls?
- What is the calibration of the p-value: does QUAL=30 actually correspond to 1:1000 odds?
- Which statistical model (binomial, Poisson, negative binomial) best fits the k-mer count data?
- Do bootstrap confidence intervals on rVAF have correct coverage?

### Method
1. Implement the binomial p-value per variant as described in `confidence-metrics.md` Section 2.1:
   - For each variant k-mer along the variant path, compute the binomial p-value given local coverage and estimated error rate
   - Combine per-k-mer p-values using Fisher's method with the autocorrelation correction (effective m_eff = m * k / read_length)
   - Convert combined p-value to Phred-scaled QUAL

2. Implement the negative binomial model from `confidence-metrics.md` Section 2.4:
   - Estimate overdispersion alpha from reference k-mer count variance
   - Compute NB p-values as an alternative to binomial

3. For each true positive and false positive call from Experiment 1:
   - Record QUAL score from both models
   - Record strand bias (if available), positional uniformity (CV of k-mer counts), and min_coverage
   - Record rVAF and its bootstrap 95% CI (1000 replicates)

4. Construct ROC curves using QUAL as a classifier for true/false positives.
5. Construct calibration plots (predicted probability vs. observed frequency of true positives).
6. Evaluate bootstrap CI coverage: what fraction of true VAFs fall within the 95% CI?

### Reference Documents
- `confidence-metrics.md`: Complete framework for statistical testing (binomial, Fisher's exact, Poisson, NB, Bayesian), strand bias computation, positional uniformity metric, composite quality score design, bootstrap CI methodology
- `false-positive-analysis.md` Section 2: Error k-mer count model at various coverages, showing that error k-mers at GGC positions can have counts of 25+ at 5000x coverage
- `sensitivity-landscape.md` Section 2.2: Signal-to-noise overlap at 0.1-1% VAF, showing where statistical models must discriminate

### Requirements
- Experiment 1 complete (true/false positive labels available)
- Statistical test implementations in kmerdet
- Bootstrap NNLS implementation

### Success Criteria
- QUAL score separates true positives from false positives with AUROC >= 0.90
- Calibration is within 2x of nominal (QUAL=30 corresponds to at most 2:1000 odds, not worse)
- Negative binomial model provides better calibration than binomial in high-coverage regions (overdispersion present per `confidence-metrics.md` Section 2.4 typical values alpha=0.01-0.1)
- Bootstrap 95% CI achieves 90-98% coverage of true VAF values
- Identify a QUAL threshold that achieves zero false positives while maximizing sensitivity

### Outputs
- ROC curves for QUAL as a TP/FP classifier (binomial vs. NB model)
- Calibration plots for both models
- Recommended QUAL threshold for clinical use
- Bootstrap CI coverage analysis
- Overdispersion parameter estimates by target

---

## Experiment 7: INDEL Normalization Impact

### Questions Answered
- How many true INDELs are missed by the filter due to representation mismatch (not detection failure)?
- Does left-alignment normalization resolve all representation mismatches, or do some require ALT sequence matching?
- What fraction of INDEL false negatives are filter-stage losses vs. detection-stage losses?

### Method
1. Using Experiment 2 instrumentation data, identify all INDELs that were:
   - Correctly detected by walking but filtered as "Not Found" (representation mismatch)
   - Correctly detected with a different position/allele representation than the reference file

2. Implement INDEL left-alignment normalization in kmerdet filter as described in `indel-filtering.md` Section 2:
   - Normalize both detected variants and expected variants before comparison
   - Use the Tan et al. (2015) algorithm: trim suffix, trim prefix, left-align, re-trim

3. Implement ALT sequence comparison as fallback from `indel-filtering.md` Section 3.

4. Run three filter configurations on all INDEL calls:
   - A: Raw coordinate matching (current behavior)
   - B: Normalized coordinate matching
   - C: Normalized + ALT sequence fallback

5. Validate normalization correctness against bcftools norm output.

### Reference Documents
- `indel-filtering.md`: Complete analysis of the INDEL representation problem, normalization algorithm with Rust implementation, ALT sequence matching approach, test cases for validation
- `room-for-improvement.md` Section 2: INDEL sensitivity collapse (38% overall) with classification edge cases as a contributing factor
- `sensitivity-landscape.md` Section 1.2-1.3: INDEL-specific failure modes including homopolymer ambiguity and classification edge cases

### Requirements
- Experiment 1 complete (baseline INDEL filter results)
- Normalization implementation in kmerdet
- Reference genome sequence for targets (available from target FASTA with flanking)
- bcftools installed for validation

### Success Criteria
- Identify >= 5 INDELs across the 10-patient cohort where normalization rescues a missed match (based on thesis finding that representation mismatch contributed to INDEL false negatives)
- Normalized coordinate matching recovers >= 80% of representation-mismatch false negatives
- Normalized + ALT fallback recovers >= 95%
- Normalization output matches bcftools norm exactly on all test cases
- Zero false matches introduced by normalization (no incorrect variant pairings)

### Outputs
- Count of INDELs rescued by normalization vs. ALT matching
- Per-patient impact on INDEL filter sensitivity
- Validation report against bcftools norm
- Recommended default filter strategy (normalize-then-alt per `indel-filtering.md` Section 4)

---

## Experiment 8: Error Correction Pre-Processing

### Questions Answered
- Does k-mer error correction (Lighter) before jellyfish counting reduce the noise floor?
- Does error correction remove true variant signal at 0.1% VAF?
- Is context-aware correction (error-profile-based) safer than spectrum-based correction?
- Does walking-phase tip removal and bulge smoothing reduce false positive paths without affecting true variants?

### Method
1. **Lighter preprocessing** (from `error-correction-alternatives.md` Sections on Lighter):
   - Run Lighter on HUMID-deduplicated FASTQ with conservative parameters (alpha adjusted for 5000x depth)
   - Generate jellyfish database from corrected reads
   - Compare k-mer counts at variant sites: before vs. after correction
   - Run kmerdet detect and compare sensitivity

2. **Walking-phase graph cleaning** (from `error-correction-alternatives.md` Section "Lightweight Error Correction During Walking"):
   - Implement tip removal: remove dead-end branches shorter than k/2 with coverage below 10% of main path
   - Implement bubble smoothing: collapse parallel paths shorter than k with >100:1 coverage ratio
   - Run kmerdet detect with and without graph cleaning
   - Count paths before and after cleaning, measure impact on false positive paths

3. **Context-aware correction** (from `error-correction-alternatives.md` Section "Strategy: HUMID + Context-Aware Correction"):
   - For each low-count k-mer, check 1-Hamming-distance neighbors
   - Flag as probable error only if: single high-count neighbor, base change matches known error profile (G>T for oxidative damage), count ratio >100:1
   - Apply this conservative correction and measure impact

### Reference Documents
- `error-correction-alternatives.md`: Complete evaluation of Lighter, BFC, Musket, CARE 2.0 for targeted panel data; risk assessment showing these tools may "correct" true variant k-mers; UMI-aware counting as alternative; walking-phase tip removal and bulge smoothing with Rust implementation sketches
- `false-positive-analysis.md` Section 2: Systematic error k-mers at GGC motifs reaching counts of 25+ at 5000x coverage
- `humid-analysis.md`: HUMID limitations including no consensus calling within UMI families, single-read representative selection preserving errors

### Requirements
- Experiment 2 complete (to know which stage to target)
- Lighter installed and parameterized for targeted panel depth
- Tip removal and bubble smoothing implementations in kmerdet

### Success Criteria
- Walking-phase graph cleaning reduces false positive paths by >= 30% without removing any true variant path
- Lighter preprocessing does NOT significantly reduce sensitivity (< 2% loss) while reducing noise k-mer counts by >= 20%
- If Lighter removes true variant k-mers (expected at 0.1% VAF per `error-correction-alternatives.md` critical insight), document the VAF threshold below which correction is unsafe
- Context-aware correction removes >= 50% of known systematic error k-mers (GGC motif errors) without affecting true variants

### Outputs
- Before/after k-mer count distributions at variant sites for each correction method
- Sensitivity impact table (3 methods x 10 patients)
- False positive path count reduction from graph cleaning
- Recommendation on which correction method(s) to integrate into the default pipeline
- Safety analysis: minimum VAF at which each method preserves variant signal

---

## Experiment 9: Graph Weighting Alternatives

### Questions Answered
- Does coverage-weighted edge weighting improve path finding compared to binary 0.01/1.0 weights?
- Does log-coverage weighting provide better discrimination between variant paths and error paths?
- Does the recommended hybrid approach (binary pathfinding + coverage-based ranking) outperform alternatives?
- For K-shortest-paths (Yen's or Eppstein's), what value of K captures all true variants while limiting noise paths?

### Method
1. Implement three weighting schemes from `graph-improvements.md` Section 2:
   - **Binary** (current): ref=0.01, alt=1.0
   - **Coverage-weighted**: weight(A->B) = 1 / count(B)
   - **Log-coverage**: weight(A->B) = -log2(count(B) / total_coverage)
   - **Ratio-based**: weight(A->B) = 1 - min(1, count(B) / expected_count)

2. For each scheme, run pathfinding on all 10 patients and measure:
   - Number of paths found per target
   - Are all true variant paths among the discovered paths?
   - Are true variant paths ranked higher (lower weight) than noise paths?
   - Runtime for pathfinding

3. Implement K-shortest-paths (Yen's algorithm, K=10) from `graph-improvements.md` Section 5 as alternative to all-shortest-paths enumeration:
   - Compare: does K=10 capture all true variants that all-shortest-paths finds?
   - Does limiting to K paths reduce noise path output?

4. Test the hybrid approach: binary weights for pathfinding, coverage-based scoring for ranking.

### Reference Documents
- `graph-improvements.md`: Analysis of binary weight limitations, four alternative weighting schemes with formulas and examples, K-shortest-paths algorithms (Yen's O(K*V*(V log V + E)), Eppstein's O(E + V log V + K)), graph pruning strategies
- `room-for-improvement.md` Section 8: Fixed edge weights throwing away coverage information during pathfinding
- `walking-improvements.md` Section 5: Branch pruning as complementary strategy

### Requirements
- Experiment 1 complete (baseline pathfinding results)
- Alternative weighting and K-shortest-paths implementations in kmerdet

### Success Criteria
- At least one alternative weighting scheme ranks true variant paths higher than binary weights in >= 80% of cases
- K=10 captures all true variant paths found by all-shortest-paths enumeration
- K-shortest-paths reduces noise path output by >= 50% in repetitive targets
- The hybrid approach (binary pathfinding + coverage ranking) achieves equivalent sensitivity to the best alternative weighting with simpler implementation
- Runtime for all schemes is within 2x of current binary weights

### Outputs
- Path count and ranking comparison table across weighting schemes
- True variant path rank distribution for each scheme
- K-shortest-paths analysis: recall at K=5, 10, 20, 50
- Runtime benchmarks
- Recommendation for default weighting strategy

---

## Experiment 10: Panel-of-Normals Construction

### Questions Answered
- How many normal samples are needed for effective artifact filtering?
- What is the false positive reduction from PoN filtering?
- Does the PoN catch systematic errors (GGC motifs, reference errors) as predicted?
- What is the sensitivity cost of PoN filtering (are any true somatic variants incorrectly filtered)?

### Method
1. Obtain N >= 20 healthy control samples processed through the same pipeline (same library prep, sequencing platform, target panel). If fewer are available, test with 5, 10, and 20 to measure marginal improvement.

2. Process each control through kmerdet detect with the optimized parameters from Experiments 4-5:
   ```bash
   for normal in normals/*.jf; do
       kmerdet detect --jf $normal --targets panel/ --output pon_results/
   done
   ```

3. Build the PoN database following `adaptive-filtering.md` Section 5:
   - For each target, record all non-reference variant paths found in normals
   - Compute frequency (fraction of normals with each variant)
   - Store as compact hash map: (target, normalized_variant) -> PonEntry

4. Test incremental PoN sizes (5, 10, 20 normals) to measure marginal FP reduction.

5. Apply PoN filtering to the 10-patient validation data:
   - Filter variants found in > 5% of normals (default threshold)
   - Test thresholds: 1%, 5%, 10%
   - Check whether any true somatic variants are incorrectly filtered

6. Characterize what the PoN catches: categorize filtered variants as sequencing error, reference error, CHIP, or common germline.

### Reference Documents
- `adaptive-filtering.md` Section 5: PoN construction methodology, filtering criteria, what the PoN catches vs. misses
- `false-positive-analysis.md`: Comprehensive catalog of false positive sources including systematic sequencing errors (GGC motifs at 0.5% error rate), reference genome errors (368K positions), CHIP variants (10-30% of cancer patients, 50.1% of patients have at least one CHIP mutation in reportable genes)
- `room-for-improvement.md` Section 3: VAF detection threshold and how PoN enables lower thresholds for paths not seen in normals
- `10-future-improvements.md`: PoN as P1 priority feature

### Requirements
- >= 5 healthy control samples (ideally >= 20) with same target panel and library prep
- Optimized detection parameters from Experiments 4-5
- PoN database storage implementation

### Success Criteria
- PoN with 20 normals reduces false positives by >= 50% when running with more permissive detection thresholds (count=1 instead of count=2)
- Zero true somatic variants incorrectly filtered at the 5% frequency threshold
- Incremental analysis shows diminishing returns: 5->10 normals provides large improvement, 10->20 provides modest improvement, suggesting 20 is sufficient for clinical use
- PoN correctly captures >= 90% of GGC-motif systematic errors identified in Experiment 2
- Storage requirements are minimal (< 100 KB per target panel per `false-positive-analysis.md` estimate of ~50 bytes per target-variant pair)

### Outputs
- PoN database for the target panel
- FP reduction curve as a function of PoN size (5, 10, 20 normals)
- Categorized list of PoN-filtered variants (error type breakdown)
- Sensitivity impact analysis (any true variants lost?)
- Recommended PoN size and frequency threshold for clinical deployment

---

## Cross-Experiment Analysis Plan

After completing all experiments, synthesize results into:

1. **Optimized default configuration**: The parameter set that maximizes sensitivity while maintaining zero false positives on the validation data. This becomes the kmerdet default for liquid biopsy ctDNA monitoring.

2. **Sensitivity improvement summary**: Waterfall chart showing cumulative sensitivity gain from each improvement (adaptive thresholds, multi-k, INDEL normalization, graph improvements, PoN filtering).

3. **Remaining gap analysis**: After all improvements, what is the residual sensitivity gap vs. alignment-based callers? What are the irreducible limitations of the k-mer approach?

4. **Clinical recommendation**: Updated hybrid strategy (k-mer screening + alignment-based confirmation) based on the achieved sensitivity profile.

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
