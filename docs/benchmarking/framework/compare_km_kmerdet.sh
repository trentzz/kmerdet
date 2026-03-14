#!/usr/bin/env bash
# compare_km_kmerdet.sh — Head-to-head comparison of km vs kmerdet.
#
# Runs both tools on the same jellyfish database and targets with identical
# parameters, then compares accuracy and performance side-by-side.
#
# Usage:
#   bash compare_km_kmerdet.sh [OPTIONS]
#
# Required:
#   --db PATH              Path to jellyfish .jf database
#   --targets PATH         Path to targets directory or FASTA file
#   --truth PATH           Path to ground truth TSV
#
# Options:
#   --output-dir DIR       Where to write results [default: comparison_results]
#   --kmerdet PATH         Path to kmerdet binary [default: kmerdet]
#   --km PATH              Path to km binary [default: km]
#   --threads N            Number of threads [default: 4]
#   --params SPEC          Parameter spec: "count=N,ratio=F" [default: count=2,ratio=0.05]
#   --param-sweep          Run all parameter combinations (see below)
#   --vaf-bins LIST        Comma-separated VAF bin boundaries
#   -h, --help             Show this message and exit
#
# Parameter sweep (with --param-sweep):
#   (count=2, ratio=0.05)     — km liquid biopsy default
#   (count=2, ratio=0.00001)  — kmerdet sensitive default
#   (count=5, ratio=0.05)     — conservative
#
# Examples:
#   # Quick comparison with default parameters
#   bash compare_km_kmerdet.sh --db data.jf --targets targets/ --truth truth.tsv
#
#   # Full parameter sweep
#   bash compare_km_kmerdet.sh --db data.jf --targets targets/ --truth truth.tsv \
#       --param-sweep --output-dir /tmp/comparison

set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="comparison_results"
KMERDET="kmerdet"
KM="km"
THREADS=4
PARAM_SWEEP=false
VAF_BINS="0.0,0.001,0.01,0.05,0.1,0.5,1.0"
SWEEP_VAF="0.0,0.0001,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1,0.2,0.5"

# Default parameters
DEFAULT_PARAMS=("count=2,ratio=0.05")

# Sweep parameter sets
SWEEP_PARAMS=(
    "count=2,ratio=0.05"
    "count=2,ratio=0.00001"
    "count=5,ratio=0.05"
)

# Required arguments
JF_DB=""
TARGETS=""
TRUTH=""

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
    case "$1" in
        --db)           JF_DB="$2"; shift 2 ;;
        --targets)      TARGETS="$2"; shift 2 ;;
        --truth)        TRUTH="$2"; shift 2 ;;
        --output-dir)   OUTPUT_DIR="$2"; shift 2 ;;
        --kmerdet)      KMERDET="$2"; shift 2 ;;
        --km)           KM="$2"; shift 2 ;;
        --threads)      THREADS="$2"; shift 2 ;;
        --params)       DEFAULT_PARAMS=("$2"); shift 2 ;;
        --param-sweep)  PARAM_SWEEP=true; shift ;;
        --vaf-bins)     VAF_BINS="$2"; shift 2 ;;
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

parse_params() {
    # Parse "count=N,ratio=F" → sets COUNT and RATIO
    local spec="$1"
    COUNT=$(echo "$spec" | sed 's/.*count=\([^,]*\).*/\1/')
    RATIO=$(echo "$spec" | sed 's/.*ratio=\([^,]*\).*/\1/')
}

time_and_run() {
    # Run a command, capture wall time and peak RSS
    local out_json="$1"
    local label="$2"
    shift 2
    local cmd=("$@")

    local wall_sec peak_kb
    local time_out
    time_out=$(mktemp)
    local start_ns
    start_ns=$(date +%s%N)

    if command -v /usr/bin/time &>/dev/null; then
        /usr/bin/time -v -- "${cmd[@]}" 2>"$time_out" || true
        local end_ns
        end_ns=$(date +%s%N)
        wall_sec=$(echo "scale=3; ($end_ns - $start_ns) / 1000000000" | bc)
        peak_kb=$(grep "Maximum resident" "$time_out" | grep -oE '[0-9]+' | tail -1 || echo "null")
    else
        "${cmd[@]}" 2>"$time_out" || true
        local end_ns
        end_ns=$(date +%s%N)
        wall_sec=$(echo "scale=3; ($end_ns - $start_ns) / 1000000000" | bc)
        peak_kb="null"
    fi
    rm -f "$time_out"

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
# Validation
# ---------------------------------------------------------------------------
[[ -n "$JF_DB" ]]  || die "--db is required"
[[ -n "$TARGETS" ]] || die "--targets is required"
[[ -n "$TRUTH" ]]  || die "--truth is required"
[[ -f "$JF_DB" ]]  || die "Database file not found: $JF_DB"
[[ -e "$TARGETS" ]] || die "Targets not found: $TARGETS"
[[ -f "$TRUTH" ]]  || die "Truth file not found: $TRUTH"

# Check for binaries
HAVE_KM=true
HAVE_KMERDET=true

if ! command -v "$KM" &>/dev/null; then
    warn "km binary not found: $KM"
    HAVE_KM=false
fi

if ! command -v "$KMERDET" &>/dev/null; then
    # Try project build directory
    if [[ -x "target/release/kmerdet" ]]; then
        KMERDET="$(pwd)/target/release/kmerdet"
    else
        warn "kmerdet binary not found: $KMERDET"
        HAVE_KMERDET=false
    fi
fi

[[ "$HAVE_KM" == "true" || "$HAVE_KMERDET" == "true" ]] || \
    die "Neither km nor kmerdet found. Install at least one tool."

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
mkdir -p "$OUTPUT_DIR"/{km,kmerdet,comparison,logs}

log "=== km vs kmerdet Head-to-Head Comparison ==="
log "Database:   $JF_DB"
log "Targets:    $TARGETS"
log "Truth:      $TRUTH"
log "Output:     $OUTPUT_DIR"
log "Threads:    $THREADS"
log "km found:   $HAVE_KM"
log "kmerdet:    $HAVE_KMERDET"

# Determine parameter sets to run
if [[ "$PARAM_SWEEP" == "true" ]]; then
    PARAM_SETS=("${SWEEP_PARAMS[@]}")
    log "Parameter sweep mode: ${#PARAM_SETS[@]} parameter sets"
else
    PARAM_SETS=("${DEFAULT_PARAMS[@]}")
fi

# Record run metadata
python3 - <<PYEOF
import json, datetime, os
meta = {
    "run_timestamp": datetime.datetime.utcnow().isoformat() + "Z",
    "jf_db": "${JF_DB}",
    "targets": "${TARGETS}",
    "truth": "${TRUTH}",
    "threads": ${THREADS},
    "param_sweep": ${PARAM_SWEEP},
    "have_km": ${HAVE_KM},
    "have_kmerdet": ${HAVE_KMERDET},
    "hostname": os.uname().nodename,
}
with open("${OUTPUT_DIR}/run_metadata.json", "w") as f:
    json.dump(meta, f, indent=2)
PYEOF

# ---------------------------------------------------------------------------
# Run benchmarks for each parameter set
# ---------------------------------------------------------------------------
run_km() {
    local count="$1"
    local ratio="$2"
    local suffix="c${count}_r${ratio}"
    local out_raw="${OUTPUT_DIR}/km/km_raw_${suffix}.tsv"
    local out_norm="${OUTPUT_DIR}/km/km_normalized_${suffix}.tsv"
    local out_timing="${OUTPUT_DIR}/km/km_timing_${suffix}.json"
    local out_bench="${OUTPUT_DIR}/km/km_benchmark_${suffix}.json"

    info "Running km find_mutation (count=$count, ratio=$ratio)..."

    # Resolve targets: km needs individual FASTA files or a directory
    local target_arg=""
    if [[ -d "$TARGETS" ]]; then
        # km's -t flag expects a directory
        target_arg="-t $TARGETS"
    else
        target_arg="-t $TARGETS"
    fi

    # Time km find_mutation
    time_and_run "$out_timing" "km_detect_${suffix}" \
        "$KM" find_mutation \
            $target_arg \
            -c "$count" \
            -r "$ratio" \
            "$JF_DB" \
        > "$out_raw" \
        2>"${OUTPUT_DIR}/logs/km_${suffix}.log"

    local n_lines
    n_lines=$(wc -l < "$out_raw" 2>/dev/null || echo 0)
    log "  km produced $n_lines lines"

    # Normalize km output to kmerdet format
    python3 "${SCRIPT_DIR}/parse_km_output.py" \
        --input "$out_raw" \
        --output "$out_norm" \
        2>"${OUTPUT_DIR}/logs/km_parse_${suffix}.log"

    # Run kmerdet benchmark on normalized km output
    if [[ "$HAVE_KMERDET" == "true" ]]; then
        "$KMERDET" benchmark \
            --results "$out_norm" \
            --truth "$TRUTH" \
            --vaf-bins "$VAF_BINS" \
            --sweep-vaf "$SWEEP_VAF" \
            --format json \
            --output "$out_bench" \
            2>"${OUTPUT_DIR}/logs/km_benchmark_${suffix}.log" \
            || warn "kmerdet benchmark failed on km output (${suffix})"
    fi

    # Print km timing summary
    python3 - <<PYEOF
import json
with open('${out_timing}') as f:
    t = json.load(f)
print(f"  km timing: wall_sec={t['wall_sec']:.2f}  peak_rss_kb={t.get('peak_rss_kb', 'N/A')}")
PYEOF
}

run_kmerdet() {
    local count="$1"
    local ratio="$2"
    local suffix="c${count}_r${ratio}"
    local out_detect="${OUTPUT_DIR}/kmerdet/kmerdet_detected_${suffix}.tsv"
    local out_timing="${OUTPUT_DIR}/kmerdet/kmerdet_timing_${suffix}.json"
    local out_bench="${OUTPUT_DIR}/kmerdet/kmerdet_benchmark_${suffix}.json"

    info "Running kmerdet detect (count=$count, ratio=$ratio)..."

    # Time kmerdet detect
    time_and_run "$out_timing" "kmerdet_detect_${suffix}" \
        "$KMERDET" detect \
            --jf "$JF_DB" \
            --targets "$TARGETS" \
            --count "$count" \
            --ratio "$ratio" \
            --threads "$THREADS" \
            --format tsv \
            --output "$out_detect"

    local n_lines
    n_lines=$(wc -l < "$out_detect" 2>/dev/null || echo 0)
    log "  kmerdet produced $n_lines lines"

    # Run kmerdet benchmark
    "$KMERDET" benchmark \
        --results "$out_detect" \
        --truth "$TRUTH" \
        --vaf-bins "$VAF_BINS" \
        --sweep-vaf "$SWEEP_VAF" \
        --format json \
        --output "$out_bench" \
        2>"${OUTPUT_DIR}/logs/kmerdet_benchmark_${suffix}.log" \
        || warn "kmerdet benchmark failed (${suffix})"

    # Print kmerdet timing summary
    python3 - <<PYEOF
import json
with open('${out_timing}') as f:
    t = json.load(f)
print(f"  kmerdet timing: wall_sec={t['wall_sec']:.2f}  peak_rss_kb={t.get('peak_rss_kb', 'N/A')}")
PYEOF
}

# ---------------------------------------------------------------------------
# Main comparison loop
# ---------------------------------------------------------------------------
for param_spec in "${PARAM_SETS[@]}"; do
    log ""
    log "============================================"
    log "Parameters: $param_spec"
    log "============================================"

    parse_params "$param_spec"

    # Run km
    if [[ "$HAVE_KM" == "true" ]]; then
        run_km "$COUNT" "$RATIO" || warn "km run failed for $param_spec"
    else
        warn "Skipping km (not available)"
    fi

    # Run kmerdet
    if [[ "$HAVE_KMERDET" == "true" ]]; then
        run_kmerdet "$COUNT" "$RATIO" || warn "kmerdet run failed for $param_spec"
    else
        warn "Skipping kmerdet (not available)"
    fi
done

# ---------------------------------------------------------------------------
# Run comparison analysis
# ---------------------------------------------------------------------------
log ""
log "Running comparison analysis..."

python3 "${SCRIPT_DIR}/compare_results.py" \
    --km-dir "${OUTPUT_DIR}/km" \
    --kmerdet-dir "${OUTPUT_DIR}/kmerdet" \
    --truth "$TRUTH" \
    --output-dir "${OUTPUT_DIR}/comparison" \
    2>"${OUTPUT_DIR}/logs/compare_results.log" \
    || warn "Comparison analysis failed — see ${OUTPUT_DIR}/logs/compare_results.log"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
log ""
log "=== Comparison Complete ==="
log "Results:      $OUTPUT_DIR"
log "  km output:      $OUTPUT_DIR/km/"
log "  kmerdet output:  $OUTPUT_DIR/kmerdet/"
log "  Comparison:     $OUTPUT_DIR/comparison/"
log ""
log "Next steps:"
log "  View comparison: cat $OUTPUT_DIR/comparison/comparison_summary.tsv"
log "  Detailed report: cat $OUTPUT_DIR/comparison/comparison_report.md"
log "  Plot comparison: python3 ${SCRIPT_DIR}/plot_results.py \\"
log "      --results-dir $OUTPUT_DIR/comparison --comparison $OUTPUT_DIR/comparison/comparison_summary.tsv"
