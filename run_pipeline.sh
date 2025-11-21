#!/bin/bash
#
# VIROME PIPELINE RUNNER v2.1
# Production-ready automation script
#


# Load java as needed for nextflow 
module load Java/21.0.2

# Gives tower access to seqera so can connect to the website and can be tracked
export TOWER_ACCESS_TOKEN="eyJ0aWQiOiAxMzAxOX0uMmM3NjYyYWJjYmE1NWU4NzkwZmE5OWYwNTQzMWVhZjA3M2Q4ZmRjOQ=="

set -euo pipefail

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Default parameters
INPUT_DIR="data"
OUTPUT_DIR="results_$(date +%Y%m%d_%H%M%S)"
PROFILE="slurm"
RESUME=""
DRY_RUN=false
HELP=false
USE_CONDA=false

# Array to collect extra Nextflow-specific flags (like -with-tower)
NF_EXTRA_ARGS=()

print_header() {
    echo -e "${BLUE}=========================================="
    echo "  VIROME ANALYSIS PIPELINE v2.1"
    echo "==========================================${NC}"
}

log_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }

print_usage() {
    cat << EOF
USAGE: $0 [OPTIONS] [NEXTFLOW_FLAGS]

This script is a wrapper for the virome Nextflow pipeline.

OPTIONS:
    -i, --input DIR   Input directory (default: data)
    -o, --output DIR  Output directory (default: results_TIMESTAMP)
    -p, --profile PROF Nextflow profile (default: slurm)
    -c, --conda     Enable conda for Nextflow
    -r, --resume    Resume previous run
    -n, --dry-run    Show command but do not run
    -h, --help     Show help

NEXTFLOW_FLAGS:
    Any flags not listed above (e.g., -with-tower, -name "my_run") are passed 
    directly to the 'nextflow run' command.

EXAMPLES:
    $0 -i my_data -o results -with-tower
    $0 -p test -name "Test_Run_1" -process.max_cpus 4
EOF
}

# ---------------- Argument Parsing ----------------

while [[ "$#" -gt 0 ]]; do
    case "$1" in
        -i|--input)
            INPUT_DIR="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -p|--profile)
            PROFILE="$2"
            shift 2
            ;;
        -c|--conda)
            USE_CONDA=true
            shift 1
            ;;
        -r|--resume)
            RESUME="-resume"
            shift 1
            ;;
        -n|--dry-run)
            DRY_RUN=true
            shift 1
            ;;
        -h|--help)
            HELP=true
            shift 1
            ;;
        *)
            # Instead of erroring out, collect the argument and pass it to Nextflow
            log_warn "Argument not recognized by wrapper, passing to Nextflow: $1"
            NF_EXTRA_ARGS+=("$1")
            shift 1
            ;;
    esac
done

# ---------------- Help ----------------
if [ "$HELP" = true ]; then
    print_usage
    exit 0
fi

# ---------------- Prepare Nextflow Environment Variables ----------------
# Check if TOWER_ACCESS_TOKEN is set for -with-tower flag
if [[ "${NF_EXTRA_ARGS[*]}" =~ "-with-tower" ]]; then
    # More robust check - verify token exists AND has content
    if [ -z "${TOWER_ACCESS_TOKEN:-}" ] || [ "${TOWER_ACCESS_TOKEN}" == "" ]; then
        log_error "ERROR: The '-with-tower' flag was provided, but the TOWER_ACCESS_TOKEN environment variable is not set or is empty."
        log_error "Please run 'export TOWER_ACCESS_TOKEN=\"YOUR_TOKEN\"' before executing the pipeline."
        log_error ""
        log_error "To verify your token is set, run: echo \$TOWER_ACCESS_TOKEN"
        exit 1
    fi
    
    # Debug: Show token length (not the token itself for security)
    log_info "TOWER_ACCESS_TOKEN found (length: ${#TOWER_ACCESS_TOKEN} characters)"
    
    # Ensure the token is exported so Nextflow can access it
    export TOWER_ACCESS_TOKEN
    log_info "TOWER_ACCESS_TOKEN successfully exported for Nextflow."
fi


# ---------------- Display run info ----------------
print_header
log_info "Starting pipeline with:"
echo -e " Input:     ${YELLOW}$INPUT_DIR${NC}"
echo -e " Output:    ${YELLOW}$OUTPUT_DIR${NC}"
echo -e " Profile:   ${YELLOW}$PROFILE${NC}"
echo -e " Conda:    ${YELLOW}$USE_CONDA${NC}"
echo -e " Resume:    ${YELLOW}${RESUME:-false}${NC}"
echo -e " Nextflow Args: ${YELLOW}${NF_EXTRA_ARGS[*]:-none}${NC}"
echo -e " Dry Run:   ${YELLOW}$DRY_RUN${NC}"
echo "------------------------------------------"

# ---------------- Build NF command ----------------
NF_CMD_ARRAY=(
    "run"
    "main.nf"
    "-profile" "$PROFILE"
)

# Input/output rules
if [ "$PROFILE" != "test" ]; then
    NF_CMD_ARRAY+=( "--input_dir" "$INPUT_DIR" )
    NF_CMD_ARRAY+=( "--outdir" "$OUTPUT_DIR" )
else
    log_warn "Profile 'test': using internal test paths."
fi

# Conda support
if [ "$USE_CONDA" = true ]; then
    NF_CMD_ARRAY+=( "-with-conda" )
fi

# Resume
if [ -n "$RESUME" ]; then
    NF_CMD_ARRAY+=( "$RESUME" )
fi

# Append all unhandled arguments
NF_CMD_ARRAY+=("${NF_EXTRA_ARGS[@]}")

# ---------------- Dry Run ----------------
if [ "$DRY_RUN" = true ]; then
    log_warn "DRY RUN - Not executing."
    echo -e "${BLUE}nextflow ${NF_CMD_ARRAY[*]}${NC}"
    exit 0
fi

# ---------------- Execute ----------------
log_info "Executing:"
echo -e "${BLUE}nextflow ${NF_CMD_ARRAY[*]}${NC}"
echo "------------------------------------------"

nextflow "${NF_CMD_ARRAY[@]}"

status=$?

if [ $status -eq 0 ]; then
    log_info "Pipeline finished successfully."
    log_info "Results saved to: ${YELLOW}$OUTPUT_DIR${NC}"
else
    log_error "Pipeline FAILED (exit code $status)"
    exit $status
fi