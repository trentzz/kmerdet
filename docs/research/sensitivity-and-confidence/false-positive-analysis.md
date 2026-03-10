# False Positive Analysis: Sources and Suppression Strategies

## Overview

The thesis validation achieved zero false positives across the 10-patient cohort using the full kam pipeline (HUMID deduplication + kmtools filter). This is a remarkable result, but it was achieved on a small cohort with tumor-informed targets and conservative filtering. Understanding the sources of false positives -- and the conditions under which they would emerge -- is essential for building a robust production system. This document catalogs the sources of false positive k-mer variant calls, analyzes why the thesis achieved zero FP, and identifies scenarios where this would break down.

---

## 1. PCR Artifacts Surviving Deduplication

### 1.1 How PCR Errors Propagate Through Amplification

PCR amplification during library preparation introduces base substitution errors at a rate of approximately 10^-5 to 10^-6 per base per cycle for high-fidelity polymerases (Q5, KAPA HiFi) and 10^-4 to 10^-5 for standard Taq polymerase. Over 10-15 cycles of PCR, these errors accumulate:

**Early-cycle errors**: An error introduced in cycle 1 is amplified in all subsequent cycles, producing a large family of reads carrying the error. With 14 cycles of amplification, an early error produces up to 2^13 = 8192 copies. These error copies are indistinguishable from true variant reads in terms of sequence content.

**Late-cycle errors**: An error in cycle 14 produces only 1 copy. These are easily handled by UMI deduplication (they are singletons within their UMI family).

**The critical window**: Errors in cycles 1-5 produce large families that can survive deduplication and create k-mers with counts well above the detection threshold. At 3000x post-dedup coverage, a cycle-1 error could contribute ~100 reads (3000 * 2^1 / 2^14 ~ 0.4 reads per cycle-1 error origin, but considering the target represents a small fraction of the genome, multiple cycle-1 errors at the same position from different molecules are extremely unlikely). The real concern is when multiple independent PCR errors converge on the same base change.

**Error spectrum bias**: PCR errors are not uniformly distributed across all possible substitutions. Oxidative damage during library preparation preferentially creates 8-oxo-guanine, leading to G>T (and complementary C>A) artifacts. This oxidation artifact is a well-documented source of false positive variant calls in cfDNA sequencing, particularly in FFPE samples but also in plasma cfDNA.

### 1.2 Why HUMID Helps but Does Not Eliminate All Artifacts

HUMID performs UMI-based deduplication by clustering reads with matching UMI sequences and similar start positions. Within each UMI family, a consensus sequence is generated (or a representative read is selected). This eliminates PCR duplicates -- reads that are exact copies of a single original molecule.

**What HUMID catches**:
- Late-cycle PCR duplicates (identical UMI, identical sequence)
- Most PCR errors (present in only a subset of reads within a UMI family, removed by consensus)
- Optical duplicates (adjacent clusters on the flowcell)

**What HUMID misses**:

1. **Pre-PCR errors**: Errors introduced during end-repair, A-tailing, or adapter ligation occur before UMI tagging. These errors are present in the original tagged molecule and are indistinguishable from true variants. Every read in the UMI family carries the error, so consensus cannot correct it.

2. **UMI errors**: Sequencing errors in the UMI sequence itself can cause reads from one family to be assigned to a different family, or cause a single family to be split into two. Split families may each be too small for reliable consensus, while merged families combine molecules with different sequences.

3. **Cycle-1 PCR errors with single-molecule UMI families**: When a UMI family contains only 1 read (common at low input amounts), no consensus correction is possible. If that single read carries a PCR error, it is reported as the consensus.

4. **Chimeric molecules**: During PCR, incomplete extension products can anneal to different templates in subsequent cycles, creating chimeric reads that combine sequence from two different molecules. These chimeras have valid UMIs and valid sequence structure but contain fusion artifacts that can mimic structural variants.

### 1.3 PCR Stutter at Microsatellites

Microsatellites (short tandem repeats, STRs) are particularly susceptible to PCR artifacts due to polymerase slippage. During extension through a repeat region, the nascent strand can dissociate and re-anneal at a different repeat unit, adding or removing repeat copies. This creates "stutter" products that differ from the template by one or more repeat units.

**Impact on k-mer detection**: Stutter products generate k-mers spanning the modified repeat boundary. These k-mers are indistinguishable from true indel variant k-mers at the repeat locus. At a dinucleotide repeat (CA)_n, the stutter rate is typically 1-5% per allele per cycle, cumulating to 10-50% over 10 cycles. This means stutter k-mers can have counts comparable to true variant k-mers at 1-5% VAF.

Research by Woerner et al. (2017) demonstrated that UMIs reduce stutter-associated genotyping errors by ~75%, but residual stutter persists at rates of 0.5-2% even after UMI deduplication. For cfDNA applications targeting microsatellite instability (MSI), this residual stutter is a significant source of false positives.

**Low-temperature PCR**: Rasmussen et al. (2019) showed that isothermal amplification at lower temperatures "drastically reduces stutter artifact formation," suggesting that library preparation chemistry choices directly impact false positive rates for k-mer-based detection at repeat loci.

### 1.4 Chimeric Reads from PCR

Chimeric PCR artifacts arise from several mechanisms:

**Template switching**: During extension, the polymerase switches from one template to another at a region of microhomology. The resulting chimeric read contains sequence from two different molecules. If one molecule carries a true variant and the other is reference, the chimeric read may create novel junction k-mers that are interpreted as structural variants.

**Incomplete extension and re-annealing**: A partially extended product from one cycle acts as a primer for a different template in the next cycle. This is more common at high PCR cycle numbers and with damaged template DNA (both common in cfDNA protocols).

**Impact on k-mer walking**: Chimeric reads create k-mers that bridge sequences from different genomic locations. If the chimeric junction happens to fall within a target region, it creates a novel k-mer path that the walking algorithm may interpret as a variant. However, chimeric k-mers are typically at low count (they arise from rare events), so the `--count 2` threshold provides some protection.

---

## 2. Sequencing Error K-mers

### 2.1 Error Rate by Position in Read

Illumina sequencing error rates are not uniform across the read:

**First 5-10 bases**: Elevated error rate (~0.5-1%) due to phasing effects during cluster generation and incomplete strand synthesis initialization.

**Middle of read (positions 10-100)**: Lowest error rate (~0.1-0.3% for modern instruments like NovaSeq 6000).

**End of read (last 20-50 bases)**: Progressively increasing error rate (0.5-2%), caused by signal decay from phasing accumulation, incomplete extension, and dephasing. For 150 bp reads, positions 130-150 can have error rates 5-10x higher than the middle of the read.

**Impact on k-mer detection**: K-mers derived from the ends of reads are more likely to contain errors and thus more likely to be error k-mers. A 31-mer from the last 31 bases of a 150 bp read has a probability of containing at least one error of approximately `1 - (1 - 0.01)^31 ~ 27%`, compared to `1 - (1 - 0.002)^31 ~ 6%` for a k-mer from the middle of the read.

Since k-mer counting does not weight by base quality or read position, these error k-mers contribute equally to the database. At very high coverage (10,000x), a 1% error rate at a specific position generates ~100 error k-mers, well above any reasonable count threshold.

### 2.2 Systematic Errors by Sequence Context

Illumina platforms exhibit well-characterized sequence-context-dependent error patterns (Schirmer et al., NAR 2015; Manley et al., NAR Genomics 2021):

**GGC motif**: The most consistently problematic motif across Illumina platforms. Positions preceded by GG show elevated substitution error rates, particularly G>T transitions. The mechanism involves altered enzyme kinetics on the lagging strand when encountering the GGC pattern. Not every GGC instance triggers an error, but across a large number of reads, the elevated rate creates a consistent population of error k-mers at GGC-containing positions.

**Inverted repeats**: Short inverted repeats can form hairpin structures in single-stranded DNA on the flow cell, inhibiting proper incorporation and causing elevated error rates at the bases flanking the repeat.

**GC-extreme regions**: Both very high (>70%) and very low (<30%) GC content regions show elevated error rates. High-GC regions have reduced polymerase efficiency and increased phasing; low-GC regions have reduced cluster density and signal intensity.

**Homopolymer-adjacent positions**: The base immediately following a homopolymer run of length >= 4 has an elevated error rate due to signal bleed-through from the preceding identical bases.

**Impact on k-mer variant detection**: These systematic errors create error k-mers that are **reproducible across reads**. Unlike random errors (which produce diverse error k-mers at low count each), systematic errors produce the **same** error k-mer from many independent reads. At high coverage, this creates error k-mers with counts of 10-50, indistinguishable from true variant k-mers at 0.3-1.5% VAF.

**Quantitative example**: At a GGC position with a systematic G>T error rate of 0.5%, 5000x coverage produces ~25 error k-mers. With `--count 2`, these easily pass the threshold. The resulting false positive would have an apparent rVAF of ~0.5% -- right in the clinically relevant range for ctDNA monitoring.

### 2.3 How Error K-mers Pass Counting Thresholds at High Coverage

The relationship between coverage, error rate, and expected error k-mer count is:

```
expected_error_count = coverage * error_rate_at_position
```

| Coverage | Error Rate | Expected Error Count | Passes --count 2? | Passes --count 5? |
|----------|-----------|---------------------|-------|-------|
| 1000x | 0.1% | 1 | No | No |
| 2000x | 0.1% | 2 | Borderline | No |
| 5000x | 0.1% | 5 | Yes | Borderline |
| 5000x | 0.5% (GGC) | 25 | Yes | Yes |
| 10000x | 0.1% | 10 | Yes | Yes |
| 10000x | 0.5% (GGC) | 50 | Yes | Yes |

At the high coverages typical of cfDNA targeted sequencing (5000-10,000x), **error k-mers routinely exceed counting thresholds** at positions with systematic error patterns. This is a fundamental challenge: increasing coverage improves sensitivity for true variants but also amplifies systematic errors proportionally.

### 2.4 K-mer Error Correction as Mitigation

K-mer error correction tools from the genome assembly field could be adapted for variant detection preprocessing:

**Spectrum-based correction** (Musket, Liu et al. 2013): Counts all k-mers and identifies low-frequency k-mers as likely errors. Corrects reads by replacing error k-mers with the nearest high-frequency k-mer. For variant detection, this approach would remove true variant k-mers (which are low-frequency by definition), making it unsuitable in its standard form.

**Bayesian correction** (BayesHammer, Nikolenko et al. 2013): Uses Hamming graph clustering to separate error k-mers from rare true k-mers. BayesHammer's subclustering approach can distinguish k-mers that come from different instances of a near-repeat, which is conceptually similar to distinguishing variant k-mers from error k-mers. However, at 0.1% VAF, the variant k-mer count profile is indistinguishable from an error k-mer profile.

**Bloom filter correction** (Lighter, Song et al. 2014): Uses sampling-based correction without explicit counting. Very memory-efficient but not designed to preserve rare true variants.

**Recommended approach for kmerdet**: Rather than applying generic error correction (which risks removing true variants), implement **error-aware thresholding**: estimate the expected error k-mer count at each position based on the local sequence context and coverage, and set the detection threshold adaptively. This is effectively a position-specific version of the `--count` parameter.

---

## 3. Repetitive Region Ambiguity

### 3.1 K-mers Mapping to Multiple Genomic Locations

With k=31, approximately 85% of all possible 31-mers in the human genome are unique (appear exactly once). The remaining 15% map to multiple locations, with some k-mers appearing thousands of times (in satellite DNA, ribosomal DNA, etc.).

When a target region overlaps a non-unique k-mer:
- The k-mer count in the jellyfish database reflects contributions from ALL genomic locations, not just the target
- A variant k-mer that happens to match a common k-mer elsewhere in the genome will have an inflated count
- Reference k-mers from repetitive regions have elevated counts, potentially masking variant signal

**For targeted cfDNA panels**: Targets are typically designed in non-repetitive regions (exonic sequences of cancer genes), so this is less of a concern. However, some clinically important loci are near repeats:
- **FLT3-ITD**: The internal tandem duplication region in FLT3 is inherently repetitive (it creates a tandem repeat)
- **TP53 exon 7**: Contains a GC-rich region with sequence similarity to processed pseudogenes
- **NPM1 exon 12**: The 4-bp TCTG insertion site is in a relatively unique context, but the flanking sequence contains short repeats

### 3.2 Homopolymer Runs Creating Ambiguous Paths

In a homopolymer of length L >= k+1, all k-mers within the run are identical. This creates a self-loop in the k-mer graph: the k-mer AAAA...A (31 A's) has extension A, which produces the same k-mer. The walking algorithm cannot determine how many times to traverse this self-loop, making the path length through the homopolymer ambiguous.

**Consequences for variant detection**:
- A 1-bp insertion in a long poly-A run produces no new k-mers (the same AAAA...A k-mer spans both the reference and variant)
- A 1-bp deletion similarly produces no distinguishing k-mers
- Only the junction k-mers at the boundaries of the homopolymer differ, and these may have ambiguous extensions due to the repeat

**Prevalence**: The human genome contains ~1 million homopolymer runs of length >= 6 and ~100,000 of length >= 10. In coding regions targeted by cfDNA panels, homopolymer runs are less common but still present (e.g., poly-A tracts in BRCA1 exon 11, poly-T in APC).

### 3.3 Low-Complexity Sequences

Beyond homopolymers, other low-complexity sequences cause k-mer ambiguity:
- **Dinucleotide repeats** (CA)_n: Very common in the human genome (~200,000 instances). Each 31-mer within the repeat is one of approximately 2 distinct k-mers (depending on phase), creating a simple cycle in the k-mer graph.
- **Trinucleotide repeats** (CAG)_n: Associated with repeat expansion diseases. Similar k-mer graph structure to dinucleotide repeats.
- **AT-rich regions**: Low sequence complexity means fewer distinct k-mers per region, increasing the chance of k-mer collisions between different genomic locations.

**K-mer complexity score**: A useful metric for assessing the reliability of a variant call is the linguistic complexity of the k-mers along the variant path. This can be computed as the ratio of distinct k-mers to total k-mers in the path. Low values (<0.7) indicate repetitive sequence context and reduced confidence.

---

## 4. Reference Artifacts

### 4.1 Reference Genome Errors Creating False Variant Calls

The human reference genome (GRCh38) contains known errors that can produce false variant calls:

**Falsely duplicated regions**: GRCh38 contains approximately 1.2 Mbp of falsely duplicated sequence and 8.04 Mbp of collapsed regions, affecting 33 protein-coding genes including 12 with medical relevance (Aganezov et al., Science 2022; Formenti et al., 2022). When a target falls in a falsely duplicated region, reads map to both copies but k-mers are counted once (from the canonical sequence), potentially creating apparent coverage anomalies.

**Minor alleles in the reference**: At some positions, the reference genome carries the minor allele rather than the major allele in the global population. When a patient is homozygous for the major allele, the k-mer analysis reports a "variant" that is actually the normal state. This is particularly problematic for non-European populations whose genetic background is less represented in GRCh38.

**Assembly gaps and errors**: GRCh38 still contains 819 gaps. Regions adjacent to gaps may have reduced sequence accuracy, with bases that are incorrect or ambiguous. K-mers spanning these regions may not match the patient's actual sequence, creating spurious variant paths.

**The T2T-CHM13 reference** (Nurk et al., 2022) resolved many GRCh38 errors, identifying 368,574 heterozygous SNVs in the autosomes when CHM13 reads were aligned to GRCh38 -- indicating positions where GRCh38 has the wrong base. The FixItFelix tool (Rhie et al., 2023) specifically addresses fixing reference errors that impact variant calling.

**Impact on k-mer detection**: If a target sequence in kmerdet is designed from GRCh38 and the patient has the correct (non-GRCh38) base, the k-mer walking algorithm will detect a "variant" that is actually the reference genome error. For tumor-informed detection, this is partially mitigated by comparing variant calls against the patient's germline, but for de novo panels, reference errors can create systematic false positives.

### 4.2 Common Polymorphisms Not in the Reference

Positions where common polymorphisms (MAF > 1%) exist but the reference carries one allele are a source of germline false positives. The patient may be heterozygous for a common polymorphism that is misinterpreted as a somatic variant.

**Mitigation in the kam pipeline**: The kmtools filter uses reference mode to compare detected variants against a database of known polymorphisms. However, this requires a comprehensive and up-to-date polymorphism database (dbSNP, gnomAD).

**Residual risk**: Novel or population-specific polymorphisms not yet in databases can still cause false positives. This is particularly relevant for patients from underrepresented populations.

### 4.3 CHIP (Clonal Hematopoiesis of Indeterminate Potential)

CHIP variants represent the most challenging source of false positives for cfDNA analysis because they are real somatic mutations present in the blood but not derived from the tumor.

**Prevalence and impact**: CHIP mutations are found in the cfDNA of 10-30% of cancer patients (Razavi et al., Nature Medicine 2019; Hu et al., JCO Precision Oncology 2020). In late-stage NSCLC patients, 72% had at least one CHIP mutation detectable by targeted panel NGS (Hu et al., 2020). CHIP mutations are age-associated, with prevalence increasing from <5% in individuals under 50 to >20% in those over 70.

**Commonly affected genes**: CHIP most frequently affects DNA methylation genes (DNMT3A, TET2, ASXL1), which are also common cancer driver genes. Other affected genes include TP53, JAK2, KRAS, PPM1D, and SF3B1. A recent characterization of 16,812 advanced cancer patients (2025) found that CHIP variants in DNA repair genes like ATM and CHEK2 can be misinterpreted as tumor-derived actionable mutations, potentially leading to inappropriate treatment decisions.

**Why CHIP is uniquely problematic for k-mer detection**:
- CHIP variants are present in hematopoietic cells, which are the primary source of cfDNA (>90% of cfDNA comes from hematopoietic cell turnover)
- CHIP VAFs are typically 0.5-10%, well within the detection range of the k-mer pipeline
- The variants are real somatic mutations with proper allelic balance and strand support, so they pass all quality filters designed to catch technical artifacts
- Without a matched white blood cell control, CHIP variants are indistinguishable from ctDNA variants

**Mitigation strategies**:
1. **Paired blood control**: Sequence the patient's buffy coat (white blood cells) alongside plasma cfDNA. Variants present in both are CHIP, not tumor-derived. This is the gold standard approach adopted by Foundation Medicine, Guardant Health, and other commercial platforms.
2. **CHIP variant databases**: Filter against known CHIP-associated genes and hotspot positions. The BLOODPAC Consortium's CH/CHIP Working Group (2024) has established standardized terminology and curation approaches.
3. **Panel-of-normals**: CHIP variants that are recurrent across healthy controls can be filtered as likely CHIP rather than tumor-derived.
4. **Tumor-informed approach**: The kam pipeline uses tumor-informed targets -- variants known from the primary tumor. CHIP variants are unlikely to match the patient's specific tumor mutation profile, providing inherent protection against CHIP false positives.

---

## 5. The Thesis's Zero False Positives Achievement

### 5.1 How It Was Achieved

The zero false positive rate in the 10-patient validation was the result of multiple layered safeguards:

**Layer 1 -- UMI deduplication (HUMID)**: Collapsed PCR duplicates into consensus sequences, removing the vast majority of PCR-introduced errors. With duplex UMI consensus, the effective error rate drops from ~10^-3 (raw reads) to ~10^-4 (simplex consensus) or ~10^-7 (duplex consensus).

**Layer 2 -- Jellyfish -L 2 threshold**: Removed all singleton k-mers from the database, eliminating most random sequencing errors (which produce unique k-mers that appear exactly once).

**Layer 3 -- Walking thresholds**: The ratio threshold (0.30 sibling sum or the adjusted ratio/count parameters for low-VAF mode) prevented the walking algorithm from extending into low-confidence branches that might represent noise.

**Layer 4 -- Graph pathfinding**: The shortest-path algorithm favored reference paths (weight 0.01) over variant paths (weight 1.0), requiring strong evidence to call a variant. Fragmented or poorly-supported variant paths were not returned.

**Layer 5 -- Tumor-informed targets**: Only variants already known from the tumor biopsy were searched for. This dramatically reduces the search space and eliminates most false positive sources -- the pipeline does not discover new variants, only confirms/quantifies known ones.

**Layer 6 -- kmtools filter (reference mode)**: Compared detected variants against a reference database, filtering known polymorphisms, common artifacts, and variants inconsistent with the expected tumor mutation profile.

**Layer 7 -- Conservative parameter choice**: The parameter set `ratio=0.0001, count=3` (balanced recommendation from thesis Chapter 5) was designed to prioritize specificity over sensitivity in the clinical validation context.

**Combined effect**: Each layer independently reduces the false positive rate by 1-2 orders of magnitude:

```
Raw FP rate:                    ~10^-2  (1 per 100 positions)
After HUMID:                    ~10^-4  (removes PCR artifacts)
After -L 2:                     ~10^-5  (removes singleton errors)
After walking threshold:        ~10^-6  (removes low-confidence branches)
After graph pathfinding:        ~10^-7  (requires connected path)
After tumor-informed filter:    ~10^-8  (only known variants searched)
After kmtools filter:           ~10^-9  (explicit artifact removal)
```

With a 50-target panel and ~200 bp per target, the total search space is ~10,000 positions. An FP rate of ~10^-9 per position yields an expected 10^-5 false positives per sample -- effectively zero.

### 5.2 What Scenarios Could Break This

The zero FP result is robust within the conditions of the validation study but would not hold under several realistic scenarios:

#### Scenario 1: De novo variant discovery

If the pipeline searches for unknown variants (not tumor-informed), the search space explodes from ~50 known positions to ~10,000 positions across the target panel. Without the tumor-informed filter (Layer 5), the FP rate per position must be <10^-4 to maintain zero FP across the panel. This is at the boundary of what the remaining layers achieve.

#### Scenario 2: Lower-quality sequencing data

The validation used high-quality UMI-tagged duplex sequencing. With standard library preparation (no UMIs), Layer 1 is dramatically weakened. PCR duplicate removal based on position alone (e.g., Picard MarkDuplicates) is far less effective than UMI-based deduplication. The FP rate after deduplication could be 10-100x higher.

#### Scenario 3: Higher sequencing depth

As discussed in Section 2.3, higher coverage amplifies systematic error k-mers proportionally. At 10,000x coverage (vs. the 2000-5000x in the validation), systematic error k-mers at GGC positions could have counts of 50+, creating confident-looking false positive calls that pass all threshold filters.

#### Scenario 4: Expanded target panel

With a 500-gene panel (vs. the 50-target panel in the validation), the probability of encountering a GGC error hotspot, a homopolymer artifact, or a reference genome error within the target set increases proportionally. A 10x larger panel requires a 10x lower per-position FP rate to maintain the same per-sample FP rate.

#### Scenario 5: Different patient populations

The 10-patient validation cohort may not represent the full diversity of genomic backgrounds. Patients from populations poorly represented in GRCh38 (African, South Asian, East Asian) may have more positions where the reference allele differs from the patient's germline, creating more false positive opportunities.

#### Scenario 6: CHIP accumulation with age

Older patients (>70) have significantly higher CHIP burden. In the validation cohort, if patients were younger (e.g., 40-60), CHIP false positives would be naturally lower. In an older cohort, CHIP variants at TP53, DNMT3A, and JAK2 could produce false positives that pass all technical filters.

#### Scenario 7: Cross-contamination

Sample cross-contamination during library preparation or sequencing introduces k-mers from other samples. If another sample on the same sequencing run carries a variant at one of the target positions, cross-contamination at rates as low as 0.01-0.1% can produce detectable k-mer counts.

### 5.3 Panel-of-Normals as Additional Safeguard

A panel-of-normals (PoN) is a database of variant calls made in a cohort of healthy control samples processed through the identical pipeline. Any variant recurrently observed in normals is flagged as a likely artifact (technical or biological) rather than a true tumor variant.

**Construction**: Process N >= 20 healthy control samples through the full pipeline (HUMID + jellyfish + km + kmtools). For each target position, record:
- K-mer paths found in normals
- Variant frequencies observed in normals
- K-mer count distributions at each position

**Filtering criteria**: Flag a variant if:
- It appears in > 5% of normals (configurable via `--pon-max-freq`)
- Its k-mer count distribution overlaps with the normal distribution (indicating noise rather than signal)
- The variant matches a known CHIP hotspot and appears at a CHIP-typical VAF (0.5-5%)

**What the PoN catches**:

| False Positive Source | Caught by PoN? | Explanation |
|---|---|---|
| Systematic sequencing errors (GGC) | Yes | Same errors appear in every sample |
| Reference genome errors | Yes | Same "variant" appears in all normals |
| Common CHIP variants | Partially | Recurrent CHIP variants at high frequency |
| Rare CHIP variants | No | Individual-specific, not recurrent |
| PCR stutter at microsatellites | Yes | Stutter pattern is reproducible |
| Cross-contamination | No | Sample-specific |
| PCR chimeras | Partially | Some chimeric junctions are reproducible |

**Implementation in kmerdet**: Store the PoN as a compact database mapping `(target_id, variant_path_hash) -> frequency_in_normals`. At detection time, look up each candidate variant path and filter if frequency exceeds the threshold. Storage requirement: ~50 bytes per target-variant pair, or ~2.5 KB for a 50-target panel with ~1 non-reference path per target in normals. Even for a 500-target panel, the PoN database is trivially small.

---

## 6. Comprehensive False Positive Mitigation Strategy

### 6.1 Pre-analytical Mitigations

| Strategy | Mechanism | FP Reduction |
|---|---|---|
| UMI-tagged library prep | Enables molecular deduplication | 10-100x |
| Duplex UMI consensus | Error correction from both strands | 1000x |
| Low-cycle PCR (8-10 cycles) | Fewer PCR-introduced errors | 2-5x |
| PCR-free protocols (where possible) | Eliminates PCR artifacts entirely | 10x |
| Oxidation damage repair (e.g., 8-oxoG cleanup) | Reduces G>T artifacts | 2-5x |
| Paired buffy coat sequencing | Enables CHIP filtering | Eliminates CHIP FP |

### 6.2 Analytical Mitigations (in kmerdet)

| Strategy | Implementation | Phase |
|---|---|---|
| Adaptive error-aware thresholds | Position-specific count thresholds based on local error rate | Phase 1-2 |
| Strand bias filter | Require balanced strand support (FS < 60) | Phase 2-3 |
| Positional uniformity check | CV of k-mer counts along variant path < 1.5 | Phase 2-3 |
| Statistical significance test | Binomial/NB p-value < 0.001 (QUAL >= 30) | Phase 2-3 |
| Panel-of-normals filter | Filter recurrent artifacts seen in healthy controls | Phase 3-4 |
| CHIP annotation | Annotate variants in known CHIP genes | Phase 3-4 |
| Reference error database | Flag positions with known GRCh38 errors | Phase 4 |
| Sequence context annotation | Flag GGC, homopolymer, and low-complexity positions | Phase 2 |
| Cross-sample contamination check | Compare k-mer profiles across samples on same run | Phase 5+ |

### 6.3 Post-analytical Mitigations

| Strategy | Mechanism |
|---|---|
| Tumor-informed interpretation | Only report variants known from tumor biopsy |
| Longitudinal consistency | Flag variants that appear/disappear inconsistently across timepoints |
| Population frequency filter | Filter variants with gnomAD MAF > 0.1% (likely germline) |
| Paired buffy coat CHIP filter | Remove variants present in matched white blood cells |
| Expert review | Manual review of borderline calls with visualization |

---

## 7. Expected False Positive Rates by Configuration

| Configuration | Expected FP Rate (per sample) | Notes |
|---|---|---|
| Full pipeline (UMI + tumor-informed + filter + PoN) | ~0 | As demonstrated in thesis |
| UMI + tumor-informed + filter (no PoN) | ~0 | Tumor-informed targets provide strong protection |
| UMI + de novo discovery + filter | 0.1-1 per sample | No tumor-informed filter; depends on panel size |
| No UMI + tumor-informed + filter | 0.01-0.1 per sample | PCR artifacts survive deduplication |
| No UMI + de novo + no filter | 5-50 per sample | Minimal protection; not recommended |
| High coverage (>10,000x) + any config | 2-10x higher | Systematic errors amplified |

---

## References

- Schirmer et al. "[Insight into biases and sequencing errors for amplicon sequencing with the Illumina MiSeq platform](https://academic.oup.com/nar/article/43/6/e37/2453415)." Nucleic Acids Research, 2015.
- Manley et al. "[Sequencing error profiles of Illumina sequencing instruments](https://academic.oup.com/nargab/article/3/1/lqab019/6193612)." NAR Genomics and Bioinformatics, 2021.
- Pfeiffer et al. "[Systematic evaluation of error rates and causes in short samples in next-generation sequencing](https://www.nature.com/articles/s41598-018-29325-6)." Scientific Reports, 2018.
- Meacham et al. "[Identification and correction of systematic error in high-throughput sequence data](https://pmc.ncbi.nlm.nih.gov/articles/PMC3295828/)." BMC Bioinformatics, 2011.
- Stoler and Nekrutenko. "[Sequencing error profiles of Illumina sequencing instruments](https://pmc.ncbi.nlm.nih.gov/articles/PMC3141275/)." Nucleic Acids Research, 2011.
- Woerner et al. "[Reducing noise and stutter in short tandem repeat loci with unique molecular identifiers](https://pubmed.ncbi.nlm.nih.gov/33429137/)." Forensic Science International: Genetics, 2021.
- Rasmussen et al. "[Low temperature isothermal amplification of microsatellites drastically reduces stutter artifact formation](https://academic.oup.com/nar/article/47/21/e141/5570702)." Nucleic Acids Research, 2019.
- Razavi et al. "High-intensity sequencing reveals the sources of plasma circulating cell-free DNA variants." Nature Medicine, 2019.
- Hu et al. "[Clonal Hematopoiesis in Late-Stage Non-Small-Cell Lung Cancer](https://ascopubs.org/doi/10.1200/PO.20.00046)." JCO Precision Oncology, 2020.
- Spoor et al. "[Liquid biopsy in esophageal cancer: a case report of false-positive circulating tumor DNA detection due to clonal hematopoiesis](https://pmc.ncbi.nlm.nih.gov/articles/PMC8421960/)." Annals of Translational Medicine, 2021.
- "[Accurate and reliable detection of clonal hematopoiesis in plasma cell-free DNA](https://ashpublications.org/bloodadvances/article/10/2/413/547996/)." Blood Advances, 2025.
- "[Characterization of Plasma Cell-Free DNA Variants as of Tumor or Clonal Hematopoiesis Origin](https://aacrjournals.org/clincancerres/article/31/13/2710/763079/)." Clinical Cancer Research, 2025.
- "[Lexicon for Clonal Hematopoiesis in Liquid Biopsy](https://pmc.ncbi.nlm.nih.gov/articles/PMC12741886/)." PMC, 2024.
- Aganezov et al. "[A complete reference genome improves analysis of human genetic variation](https://www.science.org/doi/10.1126/science.abl3533)." Science, 2022.
- Rhie et al. "[FixItFelix: improving genomic analysis by fixing reference errors](https://link.springer.com/article/10.1186/s13059-023-02863-7)." Genome Biology, 2023.
- Liu et al. "[Musket: a multistage k-mer spectrum-based error corrector for Illumina sequence data](https://academic.oup.com/bioinformatics/article/29/3/308/257257)." Bioinformatics, 2013.
- Nikolenko et al. "[BayesHammer: Bayesian clustering for error correction in single-cell sequencing](https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-14-S1-S7)." BMC Genomics, 2013.
- Song et al. "[Lighter: fast and memory-efficient sequencing error correction without counting](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0509-9)." Genome Biology, 2014.
- Narzisi et al. "[Indel variant analysis of short-read sequencing data with Scalpel](https://pmc.ncbi.nlm.nih.gov/articles/PMC5507611/)." Current Protocols in Bioinformatics, 2018.
- "[Sensitivity, specificity, and accuracy of a liquid biopsy approach utilizing molecular amplification pools](https://www.nature.com/articles/s41598-021-89592-8)." Scientific Reports, 2021.
- "[Dual-molecular barcode sequencing detects rare variants in tumor and cell free DNA](https://www.nature.com/articles/s41598-020-60361-3)." Scientific Reports, 2020.
- Illumina. "[Unique Molecular Identifiers (UMIs)](https://www.illumina.com/techniques/sequencing/ngs-library-prep/multiplexing/unique-molecular-identifiers.html)."
- Thesis Chapters 4-6: Validation Study Design, Results, and Discussion.
