# Feature: VCF 4.3 Output -- Bioinformatics Integration

## What It Does

VCF output writes variant calls in VCF 4.3 format (Variant Call Format),
the standard interchange format for genomic variant data. This enables kmerdet
output to be consumed by the entire ecosystem of bioinformatics tools:
bcftools for filtering and manipulation, IGV for visualization, GATK for
downstream analysis, VEP/Annovar for annotation, and clinical reporting
pipelines that expect VCF input.

## Why It Matters

### The Interoperability Gap

The current output format (TSV with km-compatible columns) is readable only by
custom scripts and the kamiltool filtering pipeline. It cannot be loaded into
any standard bioinformatics tool without format conversion. This isolation means:

- **No IGV visualization**: Clinicians cannot view kmerdet calls alongside BAM
  alignment data in IGV, losing the ability to manually inspect questionable
  calls.

- **No bcftools filtering**: The rich filtering language of `bcftools filter`
  (e.g., `QUAL>30 && INFO/DP>100 && INFO/FS<60`) is unavailable.

- **No annotation**: VEP, Annovar, and SnpEff require VCF input for functional
  annotation (gene name, amino acid change, clinical significance).

- **No multi-caller integration**: Clinical pipelines often combine calls from
  multiple callers and require VCF format for intersection and union operations.

- **No ACMG/AMP classification**: Downstream variant classification tools
  (InterVar, Franklin) expect VCF input.

### Clinical Adoption Barrier

VCF is the lingua franca of clinical genomics. Any variant calling tool that
does not produce VCF output is effectively excluded from standard clinical
workflows. For kmerdet to be adopted in clinical laboratories, VCF output is
not optional -- it is a prerequisite.

The CAP (College of American Pathologists) and CLIA laboratory standards
increasingly reference VCF as the expected format for clinical sequencing
results. Accredited labs that adopt kmerdet will need VCF output for their
reporting pipelines.

## Research Backing

### VCF 4.3 Specification

The VCF 4.3 specification (SAMtools/hts-specs) defines the standard fields:
CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT, and per-sample columns.
kmerdet must map its internal variant representation to these fields.

### Tool Design Document (08-tool-design-and-ux.md)

The UX design document specifies VCF 4.3 as one of six supported output
formats, with INFO tags RVAF, EXPR, MINCOV, and KMER. The `noodles-vcf` crate
or manual formatting are identified as implementation options.

### Confidence Metrics Document (confidence-metrics.md)

The confidence metrics research maps kmerdet fields to standard VCF quality
annotations:

| VCF Field | K-mer Equivalent | Source |
|-----------|-----------------|--------|
| QUAL | Phred-scaled p-value | Binomial/Poisson test |
| DP | Total k-mer coverage | Reference k-mer counts |
| AD | Ref/alt k-mer counts | Reference vs variant path counts |
| AF | rVAF | NNLS decomposition |
| FS | Fisher strand bias | Strand-aware counting |

Fields with no k-mer equivalent (MQ, BaseQRankSum) are omitted. Fields unique
to kmerdet (MINCOV, EXPR, KMER_K, VTYPE) are added as custom INFO fields.

## Design Considerations

### VCF Header

The header must include:

**File format declaration**:
```
##fileformat=VCFv4.3
```

**Source and version**:
```
##source=kmerdet v0.1.0
##kmerdetCommand=detect -d sample.jf -t targets/ --count 2
```

**Reference genome**:
```
##reference=GRCh38
```

**Contig definitions** (derived from target catalog or reference info file):
```
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr7,length=159345973>
...
```

**INFO field definitions**:
```
##INFO=<ID=RVAF,Number=1,Type=Float,Description="Relative variant allele frequency from NNLS decomposition">
##INFO=<ID=EXPR,Number=1,Type=Float,Description="Variant path expression coefficient from NNLS">
##INFO=<ID=MINCOV,Number=1,Type=Integer,Description="Minimum k-mer count along variant path">
##INFO=<ID=KMER_K,Number=1,Type=Integer,Description="K-mer length used for detection">
##INFO=<ID=VTYPE,Number=1,Type=String,Description="kmerdet variant type (SNV,Insertion,Deletion,ITD,Complex)">
##INFO=<ID=TARGET,Number=1,Type=String,Description="Target name from panel catalog">
##INFO=<ID=PU,Number=1,Type=Float,Description="Positional uniformity (CV of variant k-mer counts)">
##INFO=<ID=FS,Number=1,Type=Float,Description="Fisher strand bias (Phred-scaled)">
##INFO=<ID=SB,Number=1,Type=Float,Description="Strand bias ratio (0-1)">
##INFO=<ID=RVAF_CI,Number=2,Type=Float,Description="95% confidence interval on rVAF (lower,upper)">
##INFO=<ID=CQS,Number=1,Type=Float,Description="Composite quality score">
##INFO=<ID=TIER,Number=1,Type=String,Description="Confidence tier (HIGH,MEDIUM,LOW)">
```

**FILTER definitions**:
```
##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=LowQual,Description="QUAL below threshold">
##FILTER=<ID=LowCoverage,Description="Minimum k-mer coverage below threshold">
##FILTER=<ID=StrandBias,Description="Significant strand bias (FS > 60)">
##FILTER=<ID=PositionalBias,Description="Non-uniform variant k-mer distribution (PU > 1.5)">
##FILTER=<ID=LowVAF,Description="rVAF below minimum detectable VAF">
##FILTER=<ID=PanelOfNormals,Description="Variant found in panel of normals">
```

**FORMAT field definitions** (for multi-sample VCF):
```
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype (0/1 for heterozygous variant detected)">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for ref and alt alleles">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (from k-mer coverage)">
##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele frequency (rVAF)">
```

### Field Mapping

**CHROM**: Chromosome from target catalog or reference info. Extracted from the
target's genomic coordinates. Must match contig IDs in the header.

**POS**: 1-based position of the variant on the chromosome. For SNVs, the
position of the changed base. For INDELs, the position of the base preceding
the event (VCF convention). Derived from the variant classification's
Start_offset combined with the target's genomic anchor position.

**ID**: Set to "." (no dbSNP ID) unless kmerdet is configured with an rsID
annotation database. In future phases, could be populated from ClinVar or dbSNP
lookups.

**REF**: Reference allele. For SNVs, a single base. For deletions, the deleted
sequence plus the preceding base. For insertions, the preceding base.

**ALT**: Alternate allele. For SNVs, the substituted base. For deletions, the
preceding base only. For insertions, the preceding base followed by the inserted
sequence.

**QUAL**: Phred-scaled variant quality from confidence scoring. If confidence
scoring is not yet implemented, set to "." (missing).

**FILTER**: PASS if the variant passes all configured filters, or a
semicolon-separated list of filter names that failed.

**INFO**: Key-value pairs for kmerdet-specific and standard annotations.

### INDEL Normalization

VCF convention requires INDELs to be left-aligned and parsimoniously represented.
kmerdet's variant classification produces coordinates relative to the target
sequence, which may not be left-aligned in genomic coordinates.

Left-alignment algorithm:

```rust
fn left_align(chrom: &str, pos: u64, ref_allele: &str, alt_allele: &str,
              reference: &Reference) -> (u64, String, String) {
    let mut pos = pos;
    let mut ref_seq = ref_allele.to_string();
    let mut alt_seq = alt_allele.to_string();

    // Trim matching suffix
    while ref_seq.len() > 1 && alt_seq.len() > 1
          && ref_seq.ends_with(alt_seq.chars().last().unwrap()) {
        ref_seq.pop();
        alt_seq.pop();
    }

    // Left-shift while bases match
    while pos > 1 {
        let prev_base = reference.base_at(chrom, pos - 1);
        if ref_seq.starts_with(prev_base) && alt_seq.starts_with(prev_base) {
            // Already at leftmost position
            break;
        }
        if ref_seq.chars().last() == alt_seq.chars().last() {
            ref_seq.pop();
            alt_seq.pop();
            ref_seq.insert(0, prev_base);
            alt_seq.insert(0, prev_base);
            pos -= 1;
        } else {
            break;
        }
    }

    (pos, ref_seq, alt_seq)
}
```

This normalization ensures that kmerdet's VCF output is compatible with bcftools
norm and that variants can be correctly matched across tools.

### Multi-Sample VCF

When merging results from multiple samples (via the `merge` subcommand), kmerdet
produces a multi-sample VCF with one column per sample:

```
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT  sample1  sample2  sample3
chr7   ... .  A   T   45   PASS   ...  GT:AD:DP:AF 0/1:4500,5:4505:0.001  0/0:.:.:. 0/1:3200,8:3208:0.002
```

Samples where the variant was not detected receive `0/0:.:.:..` (homozygous
reference, missing counts). The INFO field contains aggregate statistics
(e.g., number of samples with detection).

### Genotype Representation

kmerdet detects somatic variants in cfDNA, where traditional diploid genotypes
(0/0, 0/1, 1/1) do not apply. However, VCF requires a GT field. The convention:

- **0/1**: Variant detected (heterozygous is the closest representation for a
  somatic variant present at sub-clonal frequency)
- **0/0**: Variant not detected at this locus in this sample
- **1/1**: Not used (would imply 100% VAF, unlikely for somatic ctDNA variants)

The AF (allele frequency) field in FORMAT provides the actual rVAF, which is
the clinically meaningful quantity.

### Coordinate System

kmerdet internally tracks variant positions relative to target sequences (0-based
offsets). VCF uses 1-based genomic coordinates. The conversion requires:

1. Target-to-genome coordinate mapping (from the target catalog metadata or
   reference info file)
2. Adding the target's genomic start position to the internal offset
3. Converting from 0-based internal to 1-based VCF convention

The target catalog must include genomic coordinates (chromosome, start, end)
for each target to enable VCF output. If genomic coordinates are unavailable,
VCF output is not possible (an error is raised with a clear message).

### Implementation Approach

Two implementation options:

**Option A: noodles-vcf crate**: The `noodles` library provides a Rust-native
VCF writer with type-safe header construction and record building. Advantages:
correct VCF formatting guaranteed, header validation, support for VCF 4.3
features. Disadvantage: additional dependency.

**Option B: Manual formatting**: Write VCF lines directly using Rust's
`write!` macro. Simpler dependency graph but requires careful attention to
VCF formatting rules (tab separation, semicolons in INFO, colons in FORMAT).

The recommendation is Option A (noodles-vcf) for correctness, with a fallback
to manual formatting if the noodles dependency is problematic.

## Acceptance Criteria

### Format Compliance

1. Output passes `bcftools view` without errors. Specifically:
   `bcftools view output.vcf > /dev/null` exits with code 0.

2. Output passes `vcf-validator` (EBI GA4GH tool) without errors.

3. VCF header includes all required meta-information lines: fileformat,
   source, reference, contig definitions, INFO definitions, FILTER definitions.

4. All INFO fields used in data lines are defined in the header.

### Field Correctness

5. CHROM, POS, REF, ALT are correctly derived from kmerdet's internal variant
   representation with proper coordinate transformation.

6. INDELs are left-aligned before output (VCF convention).

7. QUAL contains the Phred-scaled quality score when confidence scoring is
   available, or "." when not computed.

8. FILTER contains "PASS" for variants passing all filters, or the list of
   failed filter names.

### Tool Compatibility

9. Output is loadable in IGV (Integrative Genomics Viewer) and variants
   display at correct genomic positions.

10. `bcftools filter -e 'QUAL<30'` correctly filters low-quality variants.

11. `bcftools query -f '%CHROM\t%POS\t%INFO/RVAF\n'` correctly extracts
    kmerdet-specific fields.

12. VEP (Variant Effect Predictor) can annotate the VCF output with gene
    and consequence information.

### Multi-Sample

13. Multi-sample VCF produced by the `merge` subcommand contains correct
    per-sample FORMAT fields.

14. Samples with no detection at a locus have GT=0/0 and missing data fields.

### Content Accuracy

15. For all variants in the VCF, the REF allele matches the reference genome
    at the specified CHROM:POS.

16. AF (allele frequency) in FORMAT matches RVAF in INFO for each sample.

17. AD (allelic depths) are consistent with DP (AD_ref + AD_alt = DP, within
    rounding).
