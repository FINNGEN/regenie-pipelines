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
    echo "          Cromwell outputs are extracted locally but NOT uploaded to GCS in this mode."
    exit 1
}

# Copies files from a list file to GCS using gsutil's batch input (-I) and parallel (-m) options.
copy_paths_to_gcs() {
    local list_file="$1"
    local gcs_dest="$2"
    local is_test_mode="$3" # "true" or "false"

    echo "INFO: Processing list file: $list_file" # Added minimal logging

    if [ ! -f "$list_file" ]; then return 0; fi
    local files_to_process_cmd; [ "$is_test_mode" = true ] && files_to_process_cmd="head -n10 \"$list_file\"" || files_to_process_cmd="cat \"$list_file\""

    # Create a temporary file to store the paths, ensuring only non-empty lines are included.
    local temp_paths_file=$(mktemp)
    eval "$files_to_process_cmd" | grep -v '^[[:space:]]*$' > "$temp_paths_file"

    local total_files=$(wc -l < "$temp_paths_file" | tr -d '[:space:]')

    if [ "$total_files" -eq 0 ]; then
        rm -f "$temp_paths_file"
        echo "INFO: No files found in $list_file to copy." # Added minimal logging
        return 0
    fi

    # Execute gsutil -m cp -I, piping the list of files from the temporary file.
    # Redirect stdout and stderr to /dev/null to keep the output clean.
    if gsutil -m cp -I "$gcs_dest" < "$temp_paths_file" >/dev/null 2>&1; then
        echo "INFO: Successfully copied files from $list_file." # Added minimal logging
    else
        echo "WARNING: Failed to copy some files from $list_file. Check gsutil logs." >&2 # Added minimal logging
    fi

    rm -f "$temp_paths_file" # Clean up the temporary file
}

# --- Initial Checks & Setup ---

# Validate arguments
if [ -z "${1:-}" ]; then
    usage
fi

CROMWELL_JSON_FILE="$1"
TEST_MODE=false; [[ "${2:-}" == "--test" ]] && { TEST_MODE=true; shift 2; } || shift 1

# --- Stage 1: Extract Cromwell Outputs to Local Files ---
echo ""
echo "--- Stage 1: Extracting Cromwell Outputs Locally ---"

# This stage creates local files from Cromwell outputs.
# It does NOT apply the test mode limit, as user requested outputs are always created locally.
if jq -e '.outputs | type == "object"' "$CROMWELL_JSON_FILE" &> /dev/null; then
    # Process Cromwell JSON outputs to create local files.
    # The 'jq' expression handles both single string and array values, outputting each path on a new line.
    jq -r '.outputs | to_entries[] | "\(.key)\t\(.value | tojson)"' "$CROMWELL_JSON_FILE" | tr -d '\0' | while IFS=$'\t' read -r key json_value_str; do
        FILENAME="$key"
        # The 'fromjson' filter was removed as 'json_value_str' already contains parsed JSON value.
        # This prevents the 'array (...) only strings can be parsed' error.
        echo "$json_value_str" | jq -r 'if type == "array" then .[] else . end' > "$FILENAME" || true
    done
fi

# --- Stage 2: Copy Files from Lists to GCS ---
echo ""
echo "--- Stage 2: Copying Files from Lists to GCS ---"
echo "INFO: Local Cromwell output files were created but are not uploaded to GCS by this script."

copy_paths_to_gcs "$ALL_CHAINS_LIST_FILE" "$ALL_CHAINS_GCS_DEST" "$TEST_MODE"
copy_paths_to_gcs "$PHENO_CHAINS_LIST_FILE" "$PHENO_CHAINS_GCS_DEST" "$TEST_MODE"
copy_paths_to_gcs "$ALL_OUTPUTS_LIST_FILE" "$REGENIE_OUTPUTS_GCS_DEST" "$TEST_MODE"

echo ""
echo "--- Script finished ---"
