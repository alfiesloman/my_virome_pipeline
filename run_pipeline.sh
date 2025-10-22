#!/bin/bash
#
# VIROME PIPELINE RUNNER v2.0
# Production-ready automation script
#

set -euo pipefail

# Color codes
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

# Default parameters
# Note: OUTPUT_DIR uses a timestamp for unique run results
INPUT_DIR="data"
OUTPUT_DIR="results_$(date +%Y%m%d_%H%M%S)"
PROFILE="slurm"
RESUME=""
DRY_RUN=false
HELP=false

print_header() {
    echo -e "${BLUE}=========================================="
    echo "   VIROME ANALYSIS PIPELINE v2.0"
    echo "==========================================${NC}"
}

log_info() { echo -e "${GREEN}[INFO]${NC} $1"; }
log_warn() { echo -e "${YELLOW}[WARN]${NC} $1"; }
log_error() { echo -e "${RED}[ERROR]${NC} $1"; }

print_usage() {
    # Fix: Added EOF delimiter and closing brace
    cat << EOF
USAGE: $0 [OPTIONS]

This script is a wrapper for the virome Nextflow pipeline.
It translates these simple options into the required Nextflow parameters.

OPTIONS:
    -i, --input DIR     Input directory containing FASTA files (default: data)
    -o, --output DIR    Output directory (default: results_TIMESTAMP)   
    -p, --profile PROF  Nextflow profile (default: slurm)
    -r, --resume        Resume the previous Nextflow run
    -n, --dry-run       Show the Nextflow command but do not run it
    -h, --help          Show this help message

EXAMPLES:
    $0                       # Run with defaults (slurm profile, 'data' input)
    $0 -i my_data -o results   # Custom paths
    $0 -p test               # Run using the 'test' profile
    $0 -r                      # Resume the last run
EOF
}

# --- Argument Parsing ---
# Parse command-line arguments
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
        -r|--resume)
            RESUME="-resume" # This is the Nextflow flag
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
            log_error "Unknown option: $1"
            print_usage
            exit 1
            ;;
    esac
done

# --- Main Execution ---

# Show help and exit if -h is used
if [ "$HELP" = true ]; then
    print_usage
    exit 0
fi

# Print the header
print_header

# Log the parameters
log_info "Starting pipeline with the following settings:"
echo -e "  Input (wrapper):     ${YELLOW}$INPUT_DIR${NC}"
echo -e "  Output (wrapper):    ${YELLOW}$OUTPUT_DIR${NC}"
echo -e "  Profile:             ${YELLOW}$PROFILE${NC}"
echo -e "  Resume:              ${YELLOW}${RESUME:-false}${NC}"
echo -e "  Dry Run:             ${YELLOW}$DRY_RUN${NC}"
echo "------------------------------------------"

# Construct the Nextflow command in an array for safety
NF_CMD_ARRAY=(
    "nextflow"
    "run"
    "main.nf"
    "-profile" "$PROFILE"
)

# ONLY add input/output args if NOT using the test profile
# The test profile provides its own --input_dir and --outdir
if [ "$PROFILE" != "test" ]; then
    NF_CMD_ARRAY+=("--input_dir" "$INPUT_DIR")
    NF_CMD_ARRAY+=("--outdir" "$OUTPUT_DIR")
else
    log_warn "Running with 'test' profile. Using 'test_data' input and 'test_results' output."
fi

# Add the -resume flag if it was set
if [ -n "$RESUME" ]; then
    NF_CMD_ARRAY+=("$RESUME")
fi

# Handle the dry-run flag
if [ "$DRY_RUN" = true ]; then
    log_warn "DRY RUN enabled. The command will not be executed."
    log_info "Generated Nextflow command:"
    # Print the command array elements
    echo -e "${BLUE}${NF_CMD_ARRAY[*]}${NC}"
    exit 0
fi

# Execute the pipeline
log_info "Executing Nextflow command:"
echo -e "${BLUE}${NF_CMD_ARRAY[*]}${NC}"
echo "------------------------------------------"

# Run the command
"${NF_CMD_ARRAY[@]}"

# Check exit status
if [ $? -eq 0 ]; then
    log_info "Pipeline finished successfully."
    log_info "Results are in: ${YELLOW}$OUTPUT_DIR${NC}"
else
    log_error "Pipeline failed. Check Nextflow logs for details."
    exit 1
fi