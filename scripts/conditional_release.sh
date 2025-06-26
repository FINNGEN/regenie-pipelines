#!/bin/bash

# simplified_gcp_sync.sh
# Copies files from specified local list files to Google Cloud Storage.
# Supports an optional --test mode, which limits processing to the first 10 items in each list.
# GCS destination paths remain the same regardless of test mode.

# Exit on error, unset variables, and pipeline failures.
set -euo pipefail

# --- Configuration ---
ALL_CHAINS_LIST_FILE="conditional_analysis.all_chains"
PHENO_CHAINS_LIST_FILE="conditional_analysis.pheno_chains"
ALL_OUTPUTS_LIST_FILE="conditional_analysis.all_outputs"

ALL_CHAINS_GCS_DEST="gs://r13-data/conditional_analysis/release/all_chains/"
PHENO_CHAINS_GCS_DEST="gs://r13-data/conditional_analysis/release/pheno_chains/"
REGENIE_OUTPUTS_GCS_DEST="gs://r13-data/conditional_analysis/release/regenie_outputs/"

# --- Functions ---

# Displays script usage.
usage() {
    echo "Usage: $0 [--test]"
    echo "  --test: Limits processing to first 10 items for all lists."
    exit 1
}

# Copies files from a list file to GCS.
copy_files_from_list() {
    local list_file="$1"
    local gcs_dest="$2"
    local input_cmd="$3" # "head -n10" or "cat"

    if [ -f "$list_file" ]; then
        echo "INFO: Copying from '$list_file' to $gcs_dest"
        eval "$input_cmd" "$list_file" | while read -r f; do
            [ -n "$f" ] && gsutil cp "$f" "$gcs_dest" || echo "WARNING: Failed to upload '$f'. Continuing."
        done
    else
        echo "WARNING: List file '$list_file' not found. Skipping."
    fi
}

# --- Initial Checks & Setup ---

# Check for 'gsutil'
command -v gsutil &> /dev/null || { echo "ERROR: 'gsutil' not installed. Please install Google Cloud SDK."; exit 1; }

TEST_MODE=false
# If first argument is --test, activate test mode. Otherwise, it's normal mode.
if [[ "${1:-}" == "--test" ]]; then
    TEST_MODE=true
    echo "INFO: Running in test mode (first 10 items)."
    # Shift arguments so that the script doesn't try to interpret --test as a file path later if more args were expected
    shift
fi

# Determine the command to limit inputs based on TEST_MODE
INPUT_LIMIT_CMD=$( [ "$TEST_MODE" = true ] && echo "head -n10" || echo "cat" )

# --- Main Logic: Copy Files to GCS ---
echo ""
echo "--- Copying Files to GCS ---"

# Copy files from provided list files
echo ""
copy_files_from_list "$ALL_CHAINS_LIST_FILE" "$ALL_CHAINS_GCS_DEST" "$INPUT_LIMIT_CMD"
echo ""
copy_files_from_list "$PHENO_CHAINS_LIST_FILE" "$PHENO_CHAINS_GCS_DEST" "$INPUT_LIMIT_CMD"
echo ""
copy_files_from_list "$ALL_OUTPUTS_LIST_FILE" "$REGENIE_OUTPUTS_GCS_DEST" "$INPUT_LIMIT_CMD"

echo ""
echo "--- Script finished ---"
