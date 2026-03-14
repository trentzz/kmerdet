# Real Data Benchmarking Guide

Instructions for running km vs kmerdet comparisons on real clinical sequencing data.

## Prerequisites

- **HUMID** (UMI deduplication): https://github.com/jfjlaros/HUMID
- **jellyfish** >= 2.0: https://github.com/gmarcais/Jellyfish
- **samtools** >= 1.15
- **bedtools** >= 2.30
- **km** (Python): `pip install km-walk`
- **kmerdet** (Rust): built from this repository

## Pipeline Overview

```
Raw FASTQ → HUMID dedup → jellyfish count → .jf database
                                                ↓
Reference genome + target BED → target FASTA extraction
                                                ↓
Tumor VCF → ground truth TSV
                                                ↓
                    km find_mutation + kmerdet detect → comparison
```

## Step 1: UMI Deduplication with HUMID

HUMID performs UMI-based deduplication on paired-end reads. This removes PCR duplicates while preserving unique molecules, critical for accurate k-mer counting at low VAF.

```bash
# Standard duplex UMI dedup
humid \
    --umi-length 8 \
    --consensus \
    --min-family-size 2 \
    -o deduped_R1.fastq.gz \
    -p deduped_R2.fastq.gz \
    raw_R1.fastq.gz raw_R2.fastq.gz

# For single-strand UMI (simpler protocols)
humid \
    --umi-length 8 \
    --min-family-size 1 \
    -o deduped_R1.fastq.gz \
    -p deduped_R2.fastq.gz \
    raw_R1.fastq.gz raw_R2.fastq.gz
```

**Important**: Duplex consensus sequences dramatically reduce sequencing errors, which directly improves k-mer-based variant detection at low VAF.

## Step 2: Jellyfish K-mer Counting

Create a jellyfish database from the deduplicated reads. Both km and kmerdet read the same `.jf` file.

```bash
# Convert FASTQ to FASTA (jellyfish reads FASTA)
zcat deduped_R1.fastq.gz deduped_R2.fastq.gz \
    | awk 'NR%4==1{print ">"substr($0,2)} NR%4==2{print}' \
    > deduped_reads.fa

# Count k-mers with jellyfish
#   -m 31  : k-mer length (must match detection parameters)
#   -C     : canonical mode (REQUIRED — both tools expect this)
#   -s 100M: initial hash size (increase for large samples)
#   -L 2   : minimum count threshold (filters singletons)
#   -t 4   : threads
jellyfish count \
    -m 31 \
    -C \
    -s 100M \
    -L 2 \
    -t 4 \
    -o sample.jf \
    deduped_reads.fa

# Verify the database
jellyfish stats sample.jf
# Expected output:
#   Unique:    ~10-50M k-mers
#   Distinct:  ~10-50M k-mers
#   Total:     ~100-500M k-mer occurrences
#   Max_count: typically 5000-50000
```

**Notes**:
- Always use `-C` (canonical mode). Both km and kmerdet expect canonical k-mers.
- `-L 2` is standard. It filters singleton error k-mers and reduces database size.
- For ultra-low VAF work (< 0.1%), consider `-L 1` and let the detection tools handle filtering.
- If jellyfish produces multiple files (`sample_0.jf`, `sample_1.jf`), merge them:
  ```bash
  jellyfish merge -o sample.jf sample_*.jf
  ```

## Step 3: Target FASTA Extraction

Create target sequences from the reference genome. Each target is a short (~150-200bp) sequence centered on a variant of interest, flanked by enough context for k-mer walking.

### From a BED file of target regions

```bash
# targets.bed should have flanking regions already
# Format: chrom  start  end  name
# Example: chr7  55241607  55241807  EGFR_T790M

bedtools getfasta \
    -fi reference.fa \
    -bed targets.bed \
    -name \
    -fo targets/all_targets.fa

# Split into individual files (one per target, as km expects)
mkdir -p targets/
awk '/^>/{name=$1; gsub(">","",name); file="targets/"name".fa"} {print > file}' \
    targets/all_targets.fa
```

### From a VCF of known variants

Extract flanking regions around each variant:

```bash
# For each variant, extract ±75bp flanking sequence
# (total ~150bp target, appropriate for k=31 with ~35bp walking room)
FLANK=75

while IFS=$'\t' read -r chrom pos id ref alt rest; do
    [[ "$chrom" == "#"* ]] && continue
    start=$((pos - FLANK - 1))  # BED is 0-based
    end=$((pos + ${#ref} + FLANK))
    name="${chrom}_${pos}_${ref}_${alt}"
    echo -e "${chrom}\t${start}\t${end}\t${name}"
done < variants.vcf > targets.bed

bedtools getfasta -fi reference.fa -bed targets.bed -name -fo targets/all_targets.fa
```

### Target sizing guidelines

| k-mer length | Min target size | Recommended target size |
|-------------|----------------|------------------------|
| 21          | 80 bp          | 120-150 bp             |
| 25          | 100 bp         | 140-170 bp             |
| 31          | 120 bp         | 150-200 bp             |

The target must have at least `k` bases of unique sequence on each side of the variant to anchor the walking algorithm.

## Step 4: Ground Truth Creation

Create a ground truth TSV from tumor-informed variant calls.

### From a tumor VCF

```bash
# Extract variants and format as ground truth TSV
echo -e "chrom\tpos\tref\talt\ttype\ttrue_vaf\tcategory" > ground_truth.tsv

bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/TYPE\t%INFO/AF\t%INFO/TYPE\n' \
    tumor_variants.vcf \
    >> ground_truth.tsv
```

### Manual ground truth TSV format

```tsv
chrom	pos	ref	alt	type	true_vaf	category
chr7	55241707	T	A	Substitution	0.05	Substitution
chr7	55249071	C	CGGAAT	Insertion	0.03	Insertion
chr17	7577548	GCAC	G	Deletion	0.08	Deletion
FLT3	1788	A	AGATATATA	ITD	0.15	ITD
```

**Column definitions**:
- `chrom`: Chromosome or gene name
- `pos`: 1-based position (matching the target FASTA coordinate)
- `ref`: Reference allele (VCF-style, with anchor base for INDELs)
- `alt`: Alternate allele
- `type`: Variant type — must be one of: `Substitution`, `Insertion`, `Deletion`, `ITD`, `Complex`
- `true_vaf`: True variant allele frequency as a decimal (0.05 = 5%)
- `category`: Same as type (used for grouping in analysis)

**Including absent variants**: For specificity testing, include variants known to be absent with `true_vaf = 0.0`:

```tsv
chr12	25398284	C	T	Substitution	0.0	Substitution
```

## Step 5: Run the Comparison

```bash
# Simple comparison with default parameters
bash docs/benchmarking/framework/compare_km_kmerdet.sh \
    --db sample.jf \
    --targets targets/ \
    --truth ground_truth.tsv \
    --output-dir results/patient01

# Full parameter sweep
bash docs/benchmarking/framework/compare_km_kmerdet.sh \
    --db sample.jf \
    --targets targets/ \
    --truth ground_truth.tsv \
    --output-dir results/patient01 \
    --param-sweep
```

## Expected Directory Structure

```
data/
├── real/
│   ├── patient01/
│   │   ├── sample.jf              # Jellyfish database
│   │   ├── targets/               # Target FASTA files
│   │   │   ├── EGFR_T790M.fa
│   │   │   ├── KRAS_G12D.fa
│   │   │   └── ...
│   │   └── ground_truth.tsv       # Truth file
│   ├── patient02/
│   │   └── ...
│   └── cohort_summary/
│       └── all_patients_truth.tsv  # Aggregated truth
├── results/
│   ├── patient01/
│   │   ├── km/                    # km outputs
│   │   ├── kmerdet/               # kmerdet outputs
│   │   └── comparison/            # Side-by-side analysis
│   └── aggregated/
│       └── comparison_report.md   # Multi-patient report
```

## Multi-Patient Aggregation

To compare km vs kmerdet across a cohort:

```bash
# Run comparison for each patient
for patient_dir in data/real/patient*/; do
    patient=$(basename "$patient_dir")
    echo "=== Processing $patient ==="

    bash docs/benchmarking/framework/compare_km_kmerdet.sh \
        --db "$patient_dir/sample.jf" \
        --targets "$patient_dir/targets/" \
        --truth "$patient_dir/ground_truth.tsv" \
        --output-dir "results/$patient" \
        --param-sweep
done

# Aggregate results across patients
python3 - <<'PYEOF'
import json
import pathlib
import pandas as pd

results_root = pathlib.Path("results")
all_summaries = []

for patient_dir in sorted(results_root.glob("patient*")):
    summary_file = patient_dir / "comparison" / "comparison_summary.tsv"
    if summary_file.exists():
        df = pd.read_csv(summary_file, sep="\t")
        df["patient"] = patient_dir.name
        all_summaries.append(df)

if all_summaries:
    agg = pd.concat(all_summaries, ignore_index=True)
    agg.to_csv(results_root / "aggregated" / "all_patients_comparison.tsv",
               sep="\t", index=False)

    # Print overall summary
    overall = agg[agg["category"] == "overall"]
    for metric in ["sensitivity", "f1"]:
        sub = overall[overall["metric"] == metric]
        print(f"\n{metric}:")
        print(f"  km mean:      {sub['km'].mean():.3f} (std {sub['km'].std():.3f})")
        print(f"  kmerdet mean: {sub['kmerdet'].mean():.3f} (std {sub['kmerdet'].std():.3f})")
PYEOF
```

## Troubleshooting

### km produces no output
- Verify jellyfish database: `jellyfish query sample.jf ACGTACGTACGTACGTACGTACGTACGTACG` (should return a count)
- Check that target FASTA files have sequences long enough (>= 2k for k-mer length)
- km requires targets as a directory of individual FASTA files

### kmerdet can't read the .jf file
- Ensure the jellyfish-reader crate version matches your jellyfish version
- kmerdet supports jellyfish 2.x binary format via the `jellyfish-reader` crate
- Check: `kmerdet detect --jf sample.jf --targets targets/ --count 2 --ratio 0.05 2>&1`

### Mismatched variant positions
- Ensure ground truth positions match the coordinate system of your target FASTA files
- For simulated data: positions are relative to the target sequence (1-based)
- For real data: positions should be genomic coordinates matching the reference genome

### Low sensitivity on real data
- Check UMI deduplication quality: HUMID should report consensus family statistics
- Verify coverage: `jellyfish stats sample.jf` — total k-mer count should be proportional to coverage
- Try relaxing parameters: `--count 1 --ratio 0.00001` for maximum sensitivity
