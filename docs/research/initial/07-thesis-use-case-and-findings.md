# Thesis Use Case and Findings

## Clinical Context: Liquid Biopsy MRD Monitoring

The thesis work applied k-mer-based variant detection to **minimal residual disease (MRD) monitoring** via liquid biopsy. MRD monitoring tracks residual tumor DNA fragments in a patient's blood after treatment, providing early warning of relapse without invasive tissue biopsies.

The key challenge: ctDNA in the MRD setting is present at extremely low variant allele frequencies (often <1%), requiring both high-sensitivity detection and rigorous artifact filtering.

## 10-Patient Validation Study

### Study Design

A validation cohort of **10 patients** was used to evaluate the kam pipeline for ctDNA monitoring using UMI-based duplex sequencing:

- **Sequencing method**: UMI-based duplex sequencing (provides error correction via molecular barcodes)
- **Panel**: Targeted panel covering patient-specific somatic mutations identified from tumor tissue
- **Samples**: Serial blood draws at multiple timepoints per patient
- **Preprocessing**: HUMID deduplication to collapse UMI families before jellyfish counting
- **Pipeline**: Full kam workflow (HUMID -> jellyfish count -> multiseqex -> refolder -> kmtools chunk -> kmtools filter)

### Patient-Specific Target Design

For each patient, targets were designed from:

1. Known somatic mutations identified in the primary tumor biopsy
2. Flanking reference sequence extracted via `multiseqex` with `--flank 35`
3. Both SNVs and INDELs included where available

This tumor-informed approach maximizes sensitivity by searching for known variants rather than performing de novo discovery.

## Detection Performance

### SNV Detection

- **Sensitivity**: 77% for SNVs across the validation cohort
- SNVs were reliably detected when the VAF exceeded the detection threshold
- False negatives primarily occurred at very low VAF (<0.1%) or in regions with high background noise
- The k-mer walking algorithm performed well for single-base substitutions, as exactly `k` k-mers differ between reference and variant paths

### INDEL Detection Challenges

INDEL detection proved significantly more difficult:

- **Lower sensitivity** compared to SNVs
- Insertions >15 bp (approaching k/2 for k=31) had reduced k-mer coverage across the variant junction
- Deletions in homopolymer regions generated ambiguous k-mer paths
- Complex indels (simultaneous insertion + deletion) sometimes produced paths that the classification logic misidentified
- The graph construction phase occasionally failed to extend through INDEL junctions when coverage was low

### Root Causes of INDEL Difficulty

| Factor | Impact |
|--------|--------|
| K-mer length constraint | Insertions >k/2 have fewer spanning k-mers |
| Homopolymer context | Repetitive k-mers cause extension ambiguity |
| Low VAF + INDEL | Compound effect — fewer supporting k-mers at variant junction |
| PCR stutter | UMI deduplication helps but does not fully eliminate stutter artifacts around indels |
| Classification logic | `diff_path_without_overlap` has edge cases for complex indels near target boundaries |

## VAF Detection Threshold

### VAF > 0.1% Finding

The validation study established a practical detection threshold of **VAF > 0.1%**:

- Below 0.1% VAF, the signal-to-noise ratio was insufficient for reliable detection
- At 0.1% VAF with typical sequencing depths (2000-5000x), variant k-mers had counts of ~2-5
- The km parameters used (`--count 2 --ratio 0.00001`) were tuned to detect at this level
- Above 0.1% VAF, detection was consistent and quantification (rVAF) correlated well with expected values

### Factors Affecting the Threshold

```
Detection reliability = f(VAF, sequencing_depth, k, deduplication_rate, target_design)
```

- Higher sequencing depth pushes the threshold lower (more variant k-mers at same VAF)
- Better deduplication (more UMI families) reduces noise floor
- Target design quality affects anchor k-mer uniqueness and walking success
- GC content and sequence complexity influence jellyfish counting accuracy

## Hybrid Clinical Strategy

### Combining K-mer Detection with Standard Approaches

The thesis proposed a **hybrid strategy** that uses k-mer detection alongside traditional variant calling:

1. **Primary screen**: Run kmerdet on targeted panel for rapid results (~6 minutes)
2. **Confirmation**: Run standard pipeline (alignment + variant calling) for positive hits
3. **Monitoring**: Use kmerdet for longitudinal tracking of confirmed variants
4. **Escalation**: If kmerdet detects a new variant not in the original panel, trigger full analysis

### Advantages of the Hybrid Approach

| Aspect | K-mer (kmerdet) | Traditional Pipeline | Hybrid |
|--------|-----------------|---------------------|--------|
| Speed | ~6 minutes | Hours | Fast screening, thorough confirmation |
| Sensitivity (SNV) | 77% | 90-95% | Best of both |
| INDEL detection | Limited | Better | Traditional fills gaps |
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
            -> kmerdet detect (6 min) -> Preliminary report
            -> BWA + GATK (hours)     -> Confirmatory report
```

## 6-Minute Workflow Performance

### Timing Breakdown

The kam pipeline achieved end-to-end k-mer analysis in approximately **6 minutes** for a typical liquid biopsy panel:

| Step | Time | Notes |
|------|------|-------|
| HUMID deduplication | ~2 min | Depends on read count and UMI complexity |
| jellyfish count | ~1.5 min | 31-mers, canonical, from deduplicated reads |
| multiseqex + refolder | ~10 sec | Target extraction and organization |
| kmtools chunk (km execution) | ~2 min | Parallel across targets (4 threads) |
| kmtools filter | ~15 sec | Reference-based filtering |
| **Total** | **~6 min** | Typical 50-target panel, 4 threads |

### Comparison with Traditional Pipelines

Traditional alignment-based pipelines for the same data:
- **BWA-MEM alignment**: 15-30 minutes
- **Duplicate marking**: 10-15 minutes
- **GATK HaplotypeCaller / Mutect2**: 30-60 minutes
- **Annotation and filtering**: 10-20 minutes
- **Total**: 1-2 hours minimum

The 10-20x speedup from k-mer analysis enables same-day clinical reporting.

## The kam Pipeline's Role

The kam pipeline was essential for enabling this workflow:

1. **HUMID integration**: UMI-aware deduplication before k-mer counting removes PCR artifacts that would otherwise create false positive k-mer paths
2. **Sensitive parameters**: `--count 2 --ratio 0.00001` tuned specifically for liquid biopsy VAF ranges
3. **Parallelization**: `kmtools chunk` with `refolder` enabled processing 50+ targets in parallel
4. **Tumor-informed filtering**: `kmtools filter` with reference mode matched results against known patient mutations
5. **Reproducible environment**: Docker container ensured consistent results across runs

## Key Findings: False Positives and Negatives

### False Positives

Sources of false positive calls:

- **CHIP mutations**: Clonal hematopoiesis variants present in cfDNA at low VAF, indistinguishable from tumor ctDNA by sequence alone
- **Sequencing artifacts**: Systematic errors at specific sequence contexts (e.g., GGC motifs with certain sequencers)
- **Incomplete deduplication**: When UMI families are not fully collapsed, artifact k-mers persist
- **Alignment-free limitation**: Without mapping context, some germline variants near target regions can be called

### False Negatives

Sources of false negative calls:

- **Below threshold VAF**: Variants at <0.1% VAF not reliably detected
- **INDEL-specific failures**: As described above, structural complexity defeats k-mer walking
- **Target design issues**: Poorly designed targets with non-unique anchor k-mers cause walking failures
- **Parameter sensitivity**: The `ratio` and `count` parameters interact — aggressive filtering improves specificity but reduces sensitivity

### Parameter Sensitivity Analysis

```
ratio=0.00001, count=2:  High sensitivity, moderate false positives
ratio=0.0001,  count=3:  Balanced — recommended for clinical use
ratio=0.001,   count=5:  High specificity, misses low-VAF variants
ratio=0.01,    count=10: Too stringent for liquid biopsy (misses most ctDNA)
```

The optimal parameters depend on:
- Expected VAF range for the clinical application
- Sequencing depth achieved
- Whether the assay is screening (favor sensitivity) or confirmatory (favor specificity)

## Thesis Chapter 6 Recommendations

### Immediate Improvements

1. **Adaptive filtering parameters**: Automatically adjust `ratio` and `count` based on median k-mer coverage in each sample, rather than using fixed values across all samples
2. **INDEL-specific walking strategy**: Use shorter k-mers (k=21) for INDEL targets while keeping k=31 for SNVs
3. **Multi-pass detection**: Run km with multiple parameter sets and take the union (with appropriate FDR control)

### Medium-Term Goals

4. **Parallel HUMID integration**: Integrate deduplication directly into the k-mer counting step rather than as a separate preprocessing stage
5. **Improved quantification**: Replace linear regression with a maximum-likelihood model that accounts for k-mer count variance
6. **Panel-of-normals**: Build a database of k-mer counts from healthy controls to filter recurrent artifacts

### Long-Term Vision

7. **Rust reimplementation**: Port the core algorithm to Rust for performance (this project — kmerdet)
8. **UMI-aware k-mer counting**: Incorporate UMI information directly into the k-mer database, counting unique molecules rather than reads
9. **Machine learning integration**: Train a classifier to distinguish true variants from artifacts using features beyond just k-mer count (e.g., strand bias, positional bias, sequence context)
10. **Clinical validation**: Expanded multi-site validation study with standardized reference materials

## Summary

The thesis demonstrated that k-mer-based variant detection is viable for clinical ctDNA monitoring, with particular strength in SNV detection (77% sensitivity) and rapid turnaround (~6 minutes). The approach works best as part of a hybrid strategy alongside traditional pipelines, with the k-mer method providing fast preliminary results. Key limitations center on INDEL detection and the VAF >0.1% threshold, both of which are targets for improvement in the Rust reimplementation.

## References

- Thesis Chapter 4: Validation Study Design and Methods
- Thesis Chapter 5: Results and Performance Analysis
- Thesis Chapter 6: Discussion, Limitations, and Future Work
- HUMID: High-performance UMI deduplication (Pockrandt et al.)
- kam repository: https://github.com/trentzz/kam
