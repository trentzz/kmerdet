#!/usr/bin/env bash
# run_benchmarks.sh — Main entry point for the kmerdet benchmarking suite.
#
# Runs accuracy benchmarks (kmerdet detect + benchmark subcommand) and
# performance benchmarks (timing, memory) for one or more datasets.
#
# Usage:
#   bash run_benchmarks.sh [OPTIONS]
#
# Options:
#   --config FILE          Path to benchmark_config.toml
#                          [default: $(dirname $0)/benchmark_config.toml]
#   --dataset NAME         Run only the named dataset (repeatable)
#   --data-dir DIR         Override data directory for simulated datasets
#   --output-dir DIR       Where to write results [default: results/latest]
#   --kmerdet PATH         Path to kmerdet binary [default: kmerdet]
#   --km PATH              Path to km binary for comparison [default: km]
#   --threads N            Number of threads [default: 4]
#   --quick                Run a quick subset (1 rep, fewer parameter combos)
#   --accuracy-only        Skip performance benchmarks
#   --perf-only            Skip accuracy benchmarks
#   --no-comparison        Skip km vs kmerdet comparison
#   --vaf-bins LIST        Comma-separated VAF bin boundaries
#   --sweep-vaf LIST       Comma-separated VAF thresholds for ROC sweep
#   -h, --help             Show this message and exit
#
# Examples:
#   # Quick smoke test on simulated data
#   bash run_benchmarks.sh --quick --dataset simulated_snv_indel
#
#   # Full benchmark suite with custom output directory
#   bash run_benchmarks.sh --output-dir results/$(date +%Y%m%d)-v0.2.0
#
#   # Performance scaling test only
#   bash run_benchmarks.sh --perf-only --dataset simulated_snv_indel

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="${SCRIPT_DIR}/benchmark_config.toml"
OUTPUT_DIR="results/latest"
KMERDET="kmerdet"
KM="km"
THREADS=4
TIMING_REPS=3
VAF_BINS="0.0,0.001,0.01,0.05,0.1,0.5,1.0"
SWEEP_VAF="0.0,0.0001,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5"
QUICK=false
ACCURACY_ONLY=false
PERF_ONLY=false
NO_COMPARISON=false
DATA_DIR=""
SELECTED_DATASETS=()

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)        CONFIG_FILE="$2"; shift 2 ;;
        --dataset)       SELECTED_DATASETS+=("$2"); shift 2 ;;
        --data-dir)      DATA_DIR="$2"; shift 2 ;;
        --output-dir)    OUTPUT_DIR="$2"; shift 2 ;;
        --kmerdet)       KMERDET="$2"; shift 2 ;;
        --km)            KM="$2"; shift 2 ;;
        --threads)       THREADS="$2"; shift 2 ;;
        --quick)         QUICK=true; TIMING_REPS=1; shift ;;
        --accuracy-only) ACCURACY_ONLY=true; shift ;;
        --perf-only)     PERF_ONLY=true; shift ;;
        --no-comparison) NO_COMPARISON=true; shift ;;
        --vaf-bins)      VAF_BINS="$2"; shift 2 ;;
        --sweep-vaf)     SWEEP_VAF="$2"; shift 2 ;;
        -h|--help)
            sed -n '/^# Usage:/,/^[^#]/p' "$0" | grep '^#' | sed 's/^# \?//'
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            exit 1
            ;;
    esac
done

# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------
log()  { echo "[$(date '+%H:%M:%S')] $*"; }
info() { echo "[$(date '+%H:%M:%S')] INFO  $*"; }
warn() { echo "[$(date '+%H:%M:%S')] WARN  $*" >&2; }
die()  { echo "[$(date '+%H:%M:%S')] ERROR $*" >&2; exit 1; }

require_bin() {
    local bin="$1"
    if ! command -v "$bin" &>/dev/null; then
        warn "Required binary not found: $bin"
        return 1
    fi
    return 0
}

# Run a command under /usr/bin/time -v or hyperfine, recording wall-clock
# time and peak RSS. Outputs a JSON fragment to the given file.
#
# time_cmd <output_json> <label> [command...]
time_cmd() {
    local out_json="$1"
    local label="$2"
    shift 2
    local cmd=("$@")

    local wall_sec peak_kb

    if [[ "$QUICK" == "false" ]] && command -v hyperfine &>/dev/null; then
        # Use hyperfine for accurate timing
        local hf_out
        hf_out=$(mktemp --suffix=.json)
        hyperfine \
            --runs "$TIMING_REPS" \
            --warmup 1 \
            --export-json "$hf_out" \
            -- "${cmd[*]}" 2>/dev/null || true
        wall_sec=$(python3 -c "
import json, sys
d = json.load(open('${hf_out}'))
r = d['results'][0]
print(r['median'])
" 2>/dev/null || echo "null")
        peak_kb="null"
        rm -f "$hf_out"
    elif command -v /usr/bin/time &>/dev/null; then
        # Fallback to /usr/bin/time -v
        local time_out
        time_out=$(mktemp)
        local start_ns
        start_ns=$(date +%s%N)
        /usr/bin/time -v -- "${cmd[@]}" 2>"$time_out" || true
        local end_ns
        end_ns=$(date +%s%N)
        wall_sec=$(echo "scale=3; ($end_ns - $start_ns) / 1000000000" | bc)
        peak_kb=$(grep "Maximum resident" "$time_out" | grep -oE '[0-9]+' | tail -1 || echo "null")
        rm -f "$time_out"
    else
        # Bare timing
        local start_ns
        start_ns=$(date +%s%N)
        "${cmd[@]}" 2>/dev/null || true
        local end_ns
        end_ns=$(date +%s%N)
        wall_sec=$(echo "scale=3; ($end_ns - $start_ns) / 1000000000" | bc)
        peak_kb="null"
    fi

    # Write JSON fragment
    python3 - <<PYEOF
import json
with open('${out_json}', 'w') as f:
    json.dump({
        "label": "${label}",
        "wall_sec": ${wall_sec},
        "peak_rss_kb": ${peak_kb},
    }, f, indent=2)
PYEOF
}

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
log "kmerdet benchmarking suite"
log "Output directory: $OUTPUT_DIR"
log "Threads: $THREADS"
log "Quick mode: $QUICK"

mkdir -p "$OUTPUT_DIR"/{accuracy,performance,logs}

# Resolve kmerdet binary
if [[ "$KMERDET" == "kmerdet" ]] && ! command -v kmerdet &>/dev/null; then
    # Try to find it in the project's target/release directory
    if [[ -x "target/release/kmerdet" ]]; then
        KMERDET="$(pwd)/target/release/kmerdet"
        info "Using kmerdet at: $KMERDET"
    else
        die "kmerdet binary not found. Run 'cargo build --release' first, or pass --kmerdet /path/to/kmerdet"
    fi
fi

require_bin "$KMERDET" || die "kmerdet not available"

KMERDET_VERSION=$("$KMERDET" --version 2>&1 | head -1 || echo "unknown")
info "kmerdet version: $KMERDET_VERSION"

# Record run metadata
python3 - <<PYEOF
import json, datetime, subprocess, os

git_hash = "unknown"
try:
    git_hash = subprocess.check_output(
        ["git", "rev-parse", "--short", "HEAD"],
        stderr=subprocess.DEVNULL
    ).decode().strip()
except Exception:
    pass

meta = {
    "run_timestamp": datetime.datetime.utcnow().isoformat() + "Z",
    "kmerdet_version": "${KMERDET_VERSION}",
    "git_hash": git_hash,
    "threads": ${THREADS},
    "quick_mode": ${QUICK},
    "vaf_bins": "${VAF_BINS}",
    "sweep_vaf": "${SWEEP_VAF}",
    "hostname": os.uname().nodename,
}
with open("${OUTPUT_DIR}/run_metadata.json", "w") as f:
    json.dump(meta, f, indent=2)
print(json.dumps(meta, indent=2))
PYEOF

# ---------------------------------------------------------------------------
# Discover datasets
# ---------------------------------------------------------------------------
# For this script, dataset definitions come either from the config TOML or
# are specified inline below for the built-in "simulated" type.
#
# A dataset entry must provide:
#   DATASET_NAME       - short identifier
#   DATASET_DATA_DIR   - directory with mock .jf database and targets
#   DATASET_TRUTH_TSV  - ground truth TSV file
#   DATASET_JF_DB      - path to jellyfish database file
#   DATASET_TARGETS    - path to targets FASTA or directory of FASTA files

run_accuracy_benchmark() {
    local name="$1"
    local jf_db="$2"
    local targets="$3"
    local truth_tsv="$4"
    local out_prefix="${OUTPUT_DIR}/accuracy/${name}"
    local count="${5:-2}"
    local ratio="${6:-0.00001}"

    info "--- Accuracy benchmark: $name (count=$count, ratio=$ratio) ---"

    local detect_out="${out_prefix}_count${count}_ratio${ratio}_detected.tsv"
    local bench_out="${out_prefix}_count${count}_ratio${ratio}_benchmark.json"

    # Run detection
    log "Running: kmerdet detect --jf $jf_db --targets $targets ..."
    "$KMERDET" detect \
        --jf "$jf_db" \
        --targets "$targets" \
        --count "$count" \
        --ratio "$ratio" \
        --threads "$THREADS" \
        --format tsv \
        --output "$detect_out" \
        2>"${out_prefix}_detect.log" \
        || { warn "kmerdet detect failed for $name — see ${out_prefix}_detect.log"; return 1; }

    log "Detected $(wc -l < "$detect_out") lines in $detect_out"

    # Run benchmark comparison
    log "Running: kmerdet benchmark --results $detect_out --truth $truth_tsv ..."
    "$KMERDET" benchmark \
        --results "$detect_out" \
        --truth "$truth_tsv" \
        --vaf-bins "$VAF_BINS" \
        --sweep-vaf "$SWEEP_VAF" \
        --format json \
        --output "$bench_out" \
        2>"${out_prefix}_benchmark.log" \
        || { warn "kmerdet benchmark failed for $name — see ${out_prefix}_benchmark.log"; return 1; }

    log "Benchmark report written to: $bench_out"

    # Print summary
    if command -v python3 &>/dev/null; then
        python3 - <<PYEOF
import json, sys
with open('${bench_out}') as f:
    r = json.load(f)
s = r.get('summary', {})
c = r.get('confusion', {})
print(f"  Overall: sensitivity={s.get('sensitivity', float('nan')):.3f}  "
      f"precision={s.get('precision', float('nan')):.3f}  "
      f"F1={s.get('f1', float('nan')):.3f}  "
      f"(TP={c.get('tp',0)}, FP={c.get('fp',0)}, FN={c.get('fn_count',0)})")
print("  Per-type:")
for t in r.get('per_type', []):
    print(f"    {t['variant_type']:20s}  sensitivity={t['sensitivity']:.3f}  "
          f"present={t['present']}  TP={t['tp']}  FN={t['fn_count']}")
print("  VAF bins:")
for b in r.get('per_vaf_bin', []):
    if b['present'] > 0:
        print(f"    {b['bin_label']:20s}  sensitivity={b['sensitivity']:.3f}  "
              f"present={b['present']}  TP={b['tp']}  FN={b['fn_count']}")
PYEOF
    fi
}

run_performance_benchmark() {
    local name="$1"
    local jf_db="$2"
    local targets="$3"
    local out_prefix="${OUTPUT_DIR}/performance/${name}"

    info "--- Performance benchmark: $name ---"

    # Full-dataset timing
    local timing_json="${out_prefix}_timing.json"
    log "Timing kmerdet detect on full dataset..."
    local detect_out_perf="${out_prefix}_perf_output.tsv"

    time_cmd "$timing_json" "detect_${name}" \
        "$KMERDET" detect \
            --jf "$jf_db" \
            --targets "$targets" \
            --count 2 \
            --ratio 0.00001 \
            --threads "$THREADS" \
            --format tsv \
            --output "$detect_out_perf"

    python3 - <<PYEOF
import json
with open('${timing_json}') as f:
    t = json.load(f)
print(f"  detect wall_sec={t['wall_sec']:.2f}  peak_rss_kb={t.get('peak_rss_kb', 'N/A')}")
PYEOF

    if [[ "$QUICK" == "true" ]]; then
        info "Quick mode: skipping scaling benchmarks"
        return
    fi

    # Thread scaling: vary --threads
    log "Thread scaling benchmark..."
    local thread_results=()
    for nthreads in 1 2 4 8; do
        if [[ "$nthreads" -gt "$(nproc 2>/dev/null || echo 4)" ]]; then
            continue
        fi
        local thread_json="${out_prefix}_threads${nthreads}.json"
        local thread_out="${out_prefix}_threads${nthreads}_output.tsv"
        time_cmd "$thread_json" "detect_${name}_t${nthreads}" \
            "$KMERDET" detect \
                --jf "$jf_db" \
                --targets "$targets" \
                --count 2 \
                --ratio 0.00001 \
                --threads "$nthreads" \
                --format tsv \
                --output "$thread_out"
        python3 - <<PYEOF
import json
with open('${thread_json}') as f:
    t = json.load(f)
print(f"  threads={${nthreads}}  wall_sec={t['wall_sec']:.2f}")
PYEOF
    done
}

# ---------------------------------------------------------------------------
# Main benchmark loop
# ---------------------------------------------------------------------------
# For each dataset, we check whether it is simulated (mock db) or real.
# Simulated datasets are expected to have been generated by generate_simulated.py
# and contain a mock_db.tsv (or .jf) and a targets/ subdirectory.

DATASETS_RAN=0

run_simulated_dataset() {
    local name="$1"
    local data_dir="$2"
    local truth_tsv="$3"

    # Validate that the simulated dataset exists
    if [[ ! -d "$data_dir" ]]; then
        warn "Simulated dataset directory not found: $data_dir"
        warn "Run: python3 datasets/generate_simulated.py --output-dir $data_dir"
        return 1
    fi

    if [[ ! -f "$truth_tsv" ]]; then
        warn "Ground truth file not found: $truth_tsv"
        return 1
    fi

    # Find the jellyfish database (mock or real)
    local jf_db=""
    if [[ -f "${data_dir}/mock_db.tsv" ]]; then
        # Mock database: kmerdet will use its built-in mock mode
        jf_db="${data_dir}/mock_db.tsv"
    elif [[ -f "${data_dir}/database.jf" ]]; then
        jf_db="${data_dir}/database.jf"
    else
        # Try to find any .jf file
        jf_db=$(find "$data_dir" -name "*.jf" -type f | head -1 || true)
        if [[ -z "$jf_db" ]]; then
            warn "No jellyfish database found in $data_dir"
            return 1
        fi
    fi

    # Find targets directory
    local targets=""
    if [[ -d "${data_dir}/targets" ]]; then
        targets="${data_dir}/targets"
    elif [[ -f "${data_dir}/targets.fa" ]]; then
        targets="${data_dir}/targets.fa"
    else
        warn "No targets found in $data_dir"
        return 1
    fi

    if [[ "$PERF_ONLY" != "true" ]]; then
        run_accuracy_benchmark "$name" "$jf_db" "$targets" "$truth_tsv"

        if [[ "$QUICK" != "true" ]]; then
            # Additional accuracy runs: vary count and ratio
            for count in 2 5; do
                for ratio in 0.00001 0.0001; do
                    run_accuracy_benchmark \
                        "${name}_c${count}_r${ratio/./p}" \
                        "$jf_db" "$targets" "$truth_tsv" \
                        "$count" "$ratio"
                done
            done
        fi
    fi

    if [[ "$ACCURACY_ONLY" != "true" ]]; then
        run_performance_benchmark "$name" "$jf_db" "$targets"
    fi

    DATASETS_RAN=$((DATASETS_RAN + 1))
}

# ---------------------------------------------------------------------------
# Handle --data-dir override for simulated datasets
# ---------------------------------------------------------------------------
run_dataset() {
    local name="$1"

    # Check if this is a selected dataset
    if [[ ${#SELECTED_DATASETS[@]} -gt 0 ]]; then
        local found=false
        for sel in "${SELECTED_DATASETS[@]}"; do
            if [[ "$sel" == "$name" ]]; then
                found=true
                break
            fi
        done
        if [[ "$found" != "true" ]]; then
            return
        fi
    fi

    case "$name" in
        simulated_snv_indel)
            local dd="${DATA_DIR:-data/simulated/snv_indel}"
            run_simulated_dataset "$name" "$dd" "${dd}/ground_truth.tsv"
            ;;
        simulated_large_indels)
            local dd="${DATA_DIR:-data/simulated/large_indels}"
            run_simulated_dataset "$name" "$dd" "${dd}/ground_truth.tsv"
            ;;
        simulated_ultra_low_vaf)
            local dd="${DATA_DIR:-data/simulated/ultra_low_vaf}"
            run_simulated_dataset "$name" "$dd" "${dd}/ground_truth.tsv"
            ;;
        *)
            warn "Unknown dataset: $name (add a case to run_dataset())"
            ;;
    esac
}

# ---------------------------------------------------------------------------
# Run benchmarks
# ---------------------------------------------------------------------------
info "Starting benchmark runs"

if [[ ${#SELECTED_DATASETS[@]} -gt 0 ]]; then
    for ds in "${SELECTED_DATASETS[@]}"; do
        run_dataset "$ds"
    done
else
    # Default: run all built-in simulated datasets
    for name in simulated_snv_indel simulated_large_indels simulated_ultra_low_vaf; do
        run_dataset "$name"
    done
fi

# ---------------------------------------------------------------------------
# Aggregate results
# ---------------------------------------------------------------------------
if [[ "$DATASETS_RAN" -gt 0 ]]; then
    info "Aggregating results..."
    python3 - <<PYEOF
import json
import pathlib
import sys

results_dir = pathlib.Path("${OUTPUT_DIR}/accuracy")
aggregate = []

for bench_file in sorted(results_dir.glob("*.json")):
    try:
        with open(bench_file) as f:
            r = json.load(f)
        s = r.get("summary", {})
        aggregate.append({
            "dataset": bench_file.stem,
            "sensitivity": s.get("sensitivity"),
            "precision": s.get("precision"),
            "f1": s.get("f1"),
            "tp": r.get("confusion", {}).get("tp"),
            "fp": r.get("confusion", {}).get("fp"),
            "fn": r.get("confusion", {}).get("fn_count"),
        })
    except Exception as e:
        print(f"  Warning: could not parse {bench_file}: {e}", file=sys.stderr)

if aggregate:
    agg_file = pathlib.Path("${OUTPUT_DIR}/aggregate_accuracy.json")
    with open(agg_file, "w") as f:
        json.dump(aggregate, f, indent=2)

    # Print summary table
    print("\n=== Accuracy Summary ===")
    hdr = f"{'Dataset':<45} {'Sensitivity':>11} {'Precision':>9} {'F1':>7} {'TP':>5} {'FP':>5} {'FN':>5}"
    print(hdr)
    print("-" * len(hdr))
    for row in aggregate:
        sens = f"{row['sensitivity']:.3f}" if row['sensitivity'] is not None else "  N/A"
        prec = f"{row['precision']:.3f}"   if row['precision']    is not None else "  N/A"
        f1   = f"{row['f1']:.3f}"          if row['f1']           is not None else "  N/A"
        print(f"{row['dataset']:<45} {sens:>11} {prec:>9} {f1:>7} "
              f"{row['tp']:>5} {row['fp']:>5} {row['fn']:>5}")
    print(f"\nFull results: ${OUTPUT_DIR}/aggregate_accuracy.json")
PYEOF
fi

# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
log "Benchmark suite complete."
log "Results in: $OUTPUT_DIR"
log ""
log "Next steps:"
log "  Analyze:  python3 docs/benchmarking/framework/analyze_results.py --results-dir $OUTPUT_DIR"
log "  Plot:     python3 docs/benchmarking/framework/plot_results.py --results-dir $OUTPUT_DIR"
