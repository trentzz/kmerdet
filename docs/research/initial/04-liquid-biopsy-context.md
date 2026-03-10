# Liquid Biopsy Context for K-mer Variant Detection

## What is Liquid Biopsy?

Liquid biopsy analyzes cell-free DNA (cfDNA) from blood to detect tumor-derived mutations. The tumor-derived fraction is called circulating tumor DNA (ctDNA).

### Key Characteristics of cfDNA/ctDNA

1. **Fragment size**: cfDNA is highly fragmented (~167 bp modal size, corresponding to nucleosome-wrapped DNA)
2. **Low tumor fraction**: ctDNA typically represents 0.01%-10% of total cfDNA
   - Early-stage cancers: often <1%
   - Advanced cancers: can be 10-50%
   - Post-treatment monitoring: can be <0.1%
3. **Low variant allele frequency (VAF)**: Mutations may be present at 0.1-5% VAF
4. **Low input DNA**: Total cfDNA yield from a blood draw is typically 5-30 ng
5. **High background noise**: Clonal hematopoiesis of indeterminate potential (CHIP) mutations are common in cfDNA

## Challenges for Variant Detection

### 1. Ultra-Low VAF Detection
- Standard variant callers are designed for VAF >5-10%
- Liquid biopsy requires detection down to 0.1-0.5% VAF
- At 0.1% VAF with 1000x coverage, only ~1 supporting read
- Sequencing error rate (~0.1-1%) overlaps with true variant signal

### 2. PCR and Sequencing Artifacts
- PCR amplification introduces duplicate molecules and errors
- Sequencing errors can mimic real variants at low VAF
- Solutions: UMI (unique molecular identifiers), error correction
- **HUMID** (used in the kam pipeline) deduplicates reads before k-mer counting

### 3. Short Fragment Lengths
- cfDNA fragments (~167 bp) are shorter than typical sequencing reads
- This means fewer k-mers per fragment
- Overlapping paired-end reads can be beneficial (error correction at overlap)
- Fragment length itself can be informative (tumor cfDNA tends to be shorter)

### 4. Coverage Requirements
- Need very high sequencing depth (1000-10,000x) for low VAF detection
- Higher coverage = more k-mers = better sensitivity for k-mer approach
- But also more noise/errors to filter

## Why K-mer Approach May Be Advantageous for Liquid Biopsy

### Advantages

1. **No alignment bias**: cfDNA fragments are short; alignment can be problematic
   - Short fragments may map to multiple locations
   - K-mer counting doesn't suffer from mapping ambiguity

2. **Direct count-based quantification**: K-mer counts directly reflect molecular abundance
   - No need for complex probabilistic models of alignment quality
   - Linear relationship between k-mer count and allele frequency

3. **Speed**: Pre-computed k-mer database enables rapid targeted queries
   - Important for clinical turnaround time
   - Can query many targets simultaneously

4. **Sensitivity**: K-mer thresholding can be tuned for low-VAF detection
   - The `n_cutoff` parameter in km directly controls minimum detectable signal
   - At 1000x coverage, a 0.5% VAF variant would have ~5 supporting k-mers

5. **Handles complex variants well**: ITDs, insertions, fusions
   - These are notoriously difficult with alignment at low VAF
   - K-mer walking naturally discovers these paths

### Challenges for K-mer Approach

1. **k-mer length vs fragment size**: With 31-mers and ~167 bp fragments:
   - Each fragment produces ~137 k-mers
   - A single variant is covered by up to 31 k-mers per fragment
   - Fragment length is sufficient for k-mer analysis

2. **Error k-mers**: Sequencing errors create spurious k-mers
   - At 0.1% error rate, each k-mer position has ~0.1% chance of error
   - For a 31-mer, probability of at least one error ≈ 3%
   - These error k-mers have low counts and are filtered by threshold
   - Deduplication (HUMID) before counting reduces error k-mers

3. **Coverage depth requirements**: Need enough k-mers to distinguish signal from noise
   - For 0.5% VAF at 1000x: expected variant k-mer count ≈ 5
   - For 0.1% VAF at 5000x: expected variant k-mer count ≈ 5
   - The `n_cutoff` parameter must be set appropriately for the expected coverage

## Parameter Tuning for Liquid Biopsy

Based on the km algorithm, key parameters to adjust:

| Parameter | Standard RNA-seq | Liquid Biopsy Recommendation |
|-----------|-----------------|------------------------------|
| k (k-mer length) | 31 | 31 (fine for 167bp fragments) |
| ratio (cutoff) | 0.30 | 0.001-0.01 (much lower for low VAF) |
| count (n_cutoff) | 500 | 2-10 (much lower for low VAF) |
| jellyfish -L | 2 | 2 (keep low-count k-mers) |
| jellyfish -C | yes | yes (canonical counting) |

The kam workflow already uses `--count 2 --ratio 0.00001`, showing these parameters have been tuned for sensitive detection.

## The kam Pipeline for Liquid Biopsy

From the kam example workflow:

```
1. HUMID deduplication    → Remove PCR duplicates
2. jellyfish count        → Build k-mer database (k=31, -C canonical, -L 2)
3. multiseqex             → Extract target sequences from reference
4. refolder               → Organize targets into subfolders for parallelism
5. kmtools chunk          → Run km find_mutation in parallel across targets
6. kmtools filter         → Filter results against reference information
```

This pipeline adds:
- **Deduplication** (HUMID) to reduce PCR artifacts
- **Parallelization** (kmtools chunk + refolder) for many targets
- **Filtering** (kmtools filter) to reduce false positives
- **Flexible target design** (multiseqex) from reference coordinates

## Clinical Applications

### Monitoring (MRD - Minimal Residual Disease)
- Track known mutations over time
- Very low VAF (0.01-0.1%)
- K-mer approach: pre-define targets for patient-specific mutations

### Screening
- Detect hotspot mutations across a panel
- Moderate VAF (0.1-5%)
- K-mer approach: use catalog of common cancer mutations as targets

### Treatment Response
- Track VAF changes during therapy
- K-mer approach: quantitative rVAF tracking over serial samples

### Resistance Monitoring
- Detect emergence of resistance mutations
- K-mer approach: include known resistance mutation targets

## References

- eVIDENCE: Variant filtering for low-frequency detection in cfDNA (Scientific Reports, 2019)
- ARTEMIS: K-mer analysis of repeat elements in cfDNA (Science Translational Medicine)
- Impact of cfDNA Reference Materials on Liquid Biopsy NGS Assays (PMC, 2023)
