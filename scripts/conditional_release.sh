#!/bin/bash

# combined_cromwell_gcp_sync.sh
# This script extracts outputs from a Cromwell JSON, creating local files.
# It then copies files from specified lists to GCS.
# It supports an optional --test mode, which limits all list file processing
# to the first 10 items. Cromwell outputs are always fully extracted locally.
# GCS destination paths remain the same regardless of test mode.

# Exit on error, unset variables, and pipeline failures.
set -euo pipefail

# --- Configuration ---
# CROMWELL_EXTRACT_GCS_DEST is removed as its functionality is no longer needed for uploads.

ALL_CHAINS_LIST_FILE="conditional_analysis.all_chains"
PHENO_CHAINS_LIST_FILE="conditional_analysis.pheno_chains"
ALL_OUTPUTS_LIST_FILE="conditional_analysis.all_outputs"

ALL_CHAINS_GCS_DEST="gs://r13-data/conditional_analysis/release/all_chains/"
PHENO_CHAINS_GCS_DEST="gs://r13-data/conditional_analysis/release/pheno_chains/"
REGENIE_OUTPUTS_GCS_DEST="gs://r13-data/conditional_analysis/release/regenie_outputs/"

# --- Functions ---

# Displays script usage.
usage() {
    echo "Usage: $0 <path_to_cromwell_output.json> [--test]"
    echo "  <path_to_cromwell_output.json>: Path to Cromwell workflow output JSON."
    echo "  --test: Limits processing to first 10 items for list files only."
    echo "          Cromwell outputs are extracted locally but NOT uploaded to GCS."
    exit 1
}

# Copies files from a list file to GCS.
# This function now only handles 'list_file' type copying.
copy_paths_to_gcs() {
    local list_file="$1"
    local gcs_dest="$2"
    local is_test_mode="$3" # "true" or "false"

    echo "INFO: Copying to $gcs_dest"

    if [ -f "$list_file" ]; then
        if [ "$is_test_mode" = true ]; then
            head -n10 "$list_file" | while read -r f; do
                [ -n "$f" ] && gsutil cp "$f" "$gcs_dest" || echo "WARNING: Failed to upload '$f'. Continuing."
            done
        else
            cat "$list_file" | while read -r f; do
                [ -n "$f" ] && gsutil cp "$f" "$gcs_dest" || echo "WARNING: Failed to upload '$f'. Continuing."
            done
        fi
    else
        echo "WARNING: List file '$list_file' not found. Skipping."
    fi
}

# --- Initial Checks & Setup ---

# Check for 'gsutil' (jq is checked before its use)
command -v gsutil &> /dev/null || { echo "ERROR: 'gsutil' not installed. Please install Google Cloud SDK."; exit 1; }

# Validate arguments
[ -z "${1:-}" ] && usage

CROMWELL_JSON_FILE="$1"
TEST_MODE=false
# If second argument is --test, activate test mode. Otherwise, it's normal mode.
if [[ "${2:-}" == "--test" ]]; then
    TEST_MODE=true
    echo "INFO: Running in test mode (first 10 items for list files)."
    # Shift arguments to remove both the cromwell file and --test if present
    shift 2
else
    # Shift once to remove the cromwell file if only one arg is present
    shift 1
fi


# --- Stage 1: Extract Cromwell Outputs to Local Files ---
echo ""
echo "--- Stage 1: Extracting Cromwell Outputs Locally ---"
echo "INFO: Processing '$CROMWELL_JSON_FILE' to create local output files."


# This stage remains, creating local files from Cromwell outputs.
# It does NOT apply the test mode limit, as user requested outputs are always created locally.
if ! jq -e '.outputs | type == "object"' "$CROMWELL_JSON_FILE" &> /dev/null; then
    echo "WARNING: No valid 'outputs' object found in '$CROMWELL_JSON_FILE'. Skipping local extraction."
else
    jq_output=$(jq -r '.outputs | to_entries[] | "\(.key)\t\(.value | tojson)"' "$CROMWELL_JSON_FILE" | tr -d '\0')
    jq_status=${PIPESTATUS[0]}

    if [ "$jq_status" -ne 0 ]; then
        echo "ERROR: 'jq' command failed (exit code $jq_status) during local Cromwell output extraction."
    elif [ -z "$jq_output" ]; then
        echo "INFO: 'outputs' section is empty. No local files to create."
    else
        echo "$jq_output" | while IFS=$'\t' read -r key json_value_str; do
            FILENAME="$key"

            if echo "$json_value_str" | jq -e '.[0]' &> /dev/null; then
                echo "$json_value_str" | jq -r '.[]' > "$FILENAME"
            else
                echo "$json_value_str" | jq -r 'fromjson' > "$FILENAME"
            fi

            [ $? -eq 0 ] && echo "INFO: Created local file: $FILENAME" || {
                echo "ERROR: Failed to create local file: $FILENAME. Skipping."
            }
        done
        echo "INFO: Local Cromwell output files processed."
    fi
fi

# --- Stage 2: Copy Files from Lists to GCS ---
echo ""
echo "--- Stage 2: Copying Files from Lists to GCS ---"

# Removed section 2.1 for copying extracted Cromwell outputs to GCS.
echo "INFO: Local Cromwell output files are created but not uploaded to GCS by this script."

# Copy files from provided list files
echo ""
copy_paths_to_gcs "$ALL_CHAINS_LIST_FILE" "$ALL_CHAINS_GCS_DEST" "$TEST_MODE"
echo ""
copy_paths_to_gcs "$PHENO_CHAINS_LIST_FILE" "$PHENO_CHAINS_GCS_DEST" "$TEST_MODE"
echo ""
copy_paths_to_gcs "$ALL_OUTPUTS_LIST_FILE" "$REGENIE_OUTPUTS_GCS_DEST" "$TEST_MODE"

echo ""
echo "--- Script finished ---"
