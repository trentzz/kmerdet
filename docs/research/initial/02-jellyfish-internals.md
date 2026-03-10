# Jellyfish K-mer Counter Internals

## Overview

Jellyfish is a fast, multi-threaded k-mer counter for DNA sequences. It uses a lock-free hash table to count k-mer occurrences with high parallelism and memory efficiency.

- **Author**: Guillaume Marcais & Carl Kingsford (University of Maryland)
- **Paper**: "A fast, lock-free approach for efficient parallel counting of occurrences of k-mers" (Bioinformatics, 2011)
- **Repository**: https://github.com/gmarcais/Jellyfish
- **Current Version**: 2.2.x

## Core Architecture

### Lock-Free Hash Table

Jellyfish's key innovation is a **lock-free concurrent hash table** that allows multiple threads to count k-mers simultaneously without traditional mutex locks:

1. **Compare-and-Swap (CAS)**: Uses the CPU's atomic CAS instruction to update counters
   - Thread reads current value
   - Computes new value (current + 1)
   - Atomically swaps only if value hasn't changed since read
   - If another thread modified it, retry
2. **No locks**: Eliminates mutex contention, enabling near-linear scaling with threads
3. **Open addressing**: Uses probing (not chaining) for collision resolution

### K-mer Encoding

- **2-bit encoding**: Each nucleotide stored as 2 bits (A=00, C=01, G=10, T=11)
- A 31-mer requires only 62 bits (fits in a 64-bit word)
- **Canonical form**: By default, stores the lexicographically smaller of a k-mer and its reverse complement
  - This halves memory requirements
  - Handles both strand orientations automatically
  - Controlled by `-C` flag during counting

### Variable-Length Counters

- K-mers with low counts use small counters
- K-mers with high counts automatically expand to use multiple hash entries
- This optimizes memory usage since most k-mers have low counts (Zipf-like distribution)

## The .jf Database Format

### Binary Format

The `.jf` file is a binary file containing:

1. **Header**: JSON metadata including:
   - K-mer length
   - Hash table size
   - Whether canonical mode was used
   - Counter length
   - Other configuration parameters
2. **Hash table**: The actual k-mer count data in binary form

### Header Parsing (from km's Jellyfish.py)

```python
with open(filename, mode='rb') as f:
    header = f.readline().decode("ascii", errors='ignore')
header = "{" + header.split("{", 1)[1]
# Parse nested JSON to extract 'canonical' flag and other metadata
header = json.loads(header)
self.canonical = header["canonical"]
```

### Querying

```python
# Python API via jellyfish bindings
jf = jellyfish.QueryMerFile(filename)
k = jellyfish.MerDNA.k()  # Get k-mer length from DB

# Query a specific k-mer
kmer = jellyfish.MerDNA(seq)
if canonical:
    kmer.canonicalize()
count = jf[kmer]
```

## Key Commands

### jellyfish count
```bash
jellyfish count -m 31 -s 100M -t 4 -C -L 2 reads.fastq
```
- `-m 31`: k-mer length of 31
- `-s 100M`: initial hash table size (100 million entries)
- `-t 4`: number of threads
- `-C`: canonical mode (count both strands together)
- `-L 2`: lower count threshold (discard k-mers with count < 2, reduces noise)
- Output: binary `.jf` file

### jellyfish query
```bash
jellyfish query database.jf ACGTACGTACGTACGTACGTACGTACGTACG
```
- Returns the count of a specific k-mer

### jellyfish dump
```bash
jellyfish dump database.jf > counts.fa
```
- Converts binary DB to human-readable FASTA-like format
- Each entry: k-mer sequence + count

### jellyfish merge
```bash
jellyfish merge db1.jf db2.jf -o merged.jf
```
- Combines multiple databases

### jellyfish histo / stats
- Generate count histograms and summary statistics

## How km Uses Jellyfish

### Python Bindings

km uses jellyfish's Python bindings (`dna_jellyfish` or `jellyfish` module):

```python
import dna_jellyfish as jellyfish

# Open DB
jf = jellyfish.QueryMerFile("sample.jf")

# Query individual k-mers
kmer = jellyfish.MerDNA("ACGT...31bp...")
kmer.canonicalize()
count = jf[kmer]
```

### K-mer Walking via Jellyfish

The critical operation for variant detection is **k-mer extension**:

```python
def get_child(self, seq, forward=True):
    child = []
    sum = 0
    for c in ['A', 'C', 'G', 'T']:
        if forward:
            c_seq = seq[1:] + c      # Shift window right, add new base
        else:
            c_seq = c + seq[0:-1]    # Shift window left, add new base
        c_count = self.query(c_seq)
        child.append((c_seq, c_count))
        sum += c_count
    # Keep only children above threshold
    threshold = max(sum * cutoff, n_cutoff)
    return [x[0] for x in child if x[1] >= threshold]
```

This queries 4 candidate extensions (one per nucleotide), filters by count threshold, and returns viable extensions. This is how the algorithm "walks" through the k-mer graph built from sequencing data.

## Performance Characteristics

- **Counting speed**: ~1 billion k-mers per minute (multi-threaded)
- **Memory**: Proportional to number of distinct k-mers, not total k-mers
- **Scaling**: Near-linear with CPU threads for counting
- **Query speed**: O(1) per k-mer lookup from the hash table

## Relevance to Rust Implementation

For a Rust implementation, key considerations:

1. **Reading .jf files**: Need to parse the jellyfish binary format (header JSON + hash table)
   - Alternative: Use jellyfish C API via FFI bindings
   - Alternative: Use `rust-jellyfish` crate if available
   - Alternative: Shell out to `jellyfish query` (slow for many queries)
2. **Canonical k-mer handling**: Must match jellyfish's canonical representation
3. **2-bit encoding**: Rust can efficiently implement this with bitwise operations
4. **Batch querying**: The km algorithm queries many k-mers; batch or memory-mapped access would be ideal

## References

- Marcais, G., & Kingsford, C. (2011). A fast, lock-free approach for efficient parallel counting of occurrences of k-mers. Bioinformatics, 27(6), 764-770.
- https://github.com/gmarcais/Jellyfish
- https://www.cs.cmu.edu/~ckingsf/software/jellyfish/
