#!/bin/bash

# Default values
PYTHON_SCRIPT="compute_symmetric_channel_binder_rmsd.py"
TARGET_SEQ="IPDAFWWAVVTMTTVGYGD"
OUTPUT_FILE="analysis_results.csv"

# Automatically detect number of CPU cores
NUM_PROCESSES=$(python3 -c "import multiprocessing as mp; print(mp.cpu_count())")
echo "Detected $NUM_PROCESSES CPU cores"

# Function to print usage
usage() {
    echo "Usage: $0 [options]"
    echo "Required options:"
    echo "  -r, --ref-dir DIR        Directory containing reference files"
    echo "  -2, --af2-dir DIR        Directory containing AF2 model files"
    echo "  -3, --af3-dir DIR        Directory containing AF3 model files"
    echo "Optional:"
    echo "  -p, --python-script PATH  Path to the Python script (default: $PYTHON_SCRIPT)"
    echo "  -s, --sequence SEQ       Target sequence (default: $TARGET_SEQ)"
    echo "  -o, --output FILE        Output CSV file (default: $OUTPUT_FILE)"
    echo "  -n, --num-processes N    Number of parallel processes (default: auto-detected = $NUM_PROCESSES)"
    echo "  -h, --help              Show this help message"
    exit 1
}

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -r|--ref-dir)
            REF_DIR="$2"
            shift 2
            ;;
        -2|--af2-dir)
            AF2_DIR="$2"
            shift 2
            ;;
        -3|--af3-dir)
            AF3_DIR="$2"
            shift 2
            ;;
        -p|--python-script)
            PYTHON_SCRIPT="$2"
            shift 2
            ;;
        -s|--sequence)
            TARGET_SEQ="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        -n|--num-processes)
            NUM_PROCESSES="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Check if required directories are provided and exist
for dir_var in REF_DIR AF2_DIR AF3_DIR; do
    if [ -z "${!dir_var}" ]; then
        echo "Error: ${dir_var% _DIR} directory not specified"
        usage
    fi
    if [ ! -d "${!dir_var}" ]; then
        echo "Error: ${dir_var% _DIR} directory does not exist: ${!dir_var}"
        exit 1
    fi
done

# Check if Python script exists
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "Error: Python script not found: $PYTHON_SCRIPT"
    exit 1
fi

# Create temporary directory for parallel processing
TEMP_DIR=$(mktemp -d)
trap 'rm -rf "$TEMP_DIR"' EXIT

# Create output CSV header
echo "Model,AF_Version,Binder_RMSD,Channel_RMSD,Binder_Length_Ref,Binder_Length_Model,Channel_Length,Reference_Path,Model_Path" > "$OUTPUT_FILE"

# Create a lock file for atomic writes to the output file
exec {LOCKFD}>"$TEMP_DIR/output.lock"

# Function to find reference file
find_reference() {
    local model_num="$1"
    local sample_num="$2"
    local ref_pattern="*${model_num}*sample_${sample_num}*/interface.pdb"
    find "$REF_DIR" -type f -path "*/$ref_pattern" -print -quit
}

# Function to safely write to output file
safe_write() {
    flock $LOCKFD echo "$1" >> "$OUTPUT_FILE"
}

# Function to process a single model
process_model() {
    local model_path="$1"
    local af_version="$2"
    local model_name=$(basename "$model_path")
    local process_num="$3"

    # Extract model and sample numbers based on AF version
    if [[ "$af_version" == "AF2" ]]; then
        model_num=$(echo "$model_name" | grep -o '_[0-9]\+_sample_[0-9]\+_' | cut -d'_' -f2)
        sample_num=$(echo "$model_name" | grep -o '_sample_[0-9]\+_' | cut -d'_' -f3)
    else
        model_num=$(echo "$model_name" | grep -o '_[0-9]\+_[0-9]\+_' | cut -d'_' -f2)
        sample_num=$(echo "$model_name" | grep -o '_[0-9]\+_[0-9]\+_' | cut -d'_' -f3)
    fi

    if [ -z "$model_num" ] || [ -z "$sample_num" ]; then
        echo "Process $process_num: Warning - Could not extract model/sample numbers from $model_name"
        return
    fi

    reference_path=$(find_reference "$model_num" "$sample_num")

    if [ -z "$reference_path" ]; then
        echo "Process $process_num: Warning - No reference found for model $model_name"
        return
    fi

    echo "Process $process_num: Processing ${af_version} model ${model_num}_${sample_num}..."

    # Create a temporary output file for this model
    local temp_output="$TEMP_DIR/output_${process_num}_${model_num}_${sample_num}"

    # Run analysis
    python3 "$PYTHON_SCRIPT" \
        -refpdb "$reference_path" \
        -testpdb "$model_path" \
        -targetseq "$TARGET_SEQ" > "$temp_output" 2>&1

    # Extract results and write to CSV if successful
    if [ $? -eq 0 ]; then
        binder_rmsd=$(tail -n 6 "$temp_output" | grep "Binder RMSD:" | awk '{print $3}')
        channel_rmsd=$(tail -n 5 "$temp_output" | grep "Channel RMSD:" | awk '{print $3}')
        binder_length_ref=$(tail -n 4 "$temp_output" | grep "Binder length (reference):" | awk '{print $4}')
        binder_length_model=$(tail -n 3 "$temp_output" | grep "Binder length (test):" | awk '{print $4}')
        channel_length=$(tail -n 2 "$temp_output" | grep "Channel length:" | awk '{print $3}')

        safe_write "${model_num}_${sample_num},$af_version,$binder_rmsd,$channel_rmsd,$binder_length_ref,$binder_length_model,$channel_length,$reference_path,$model_path"
    else
        echo "Process $process_num: Warning - Analysis failed for model $model_name"
    fi

    rm -f "$temp_output"
}

export -f process_model
export -f find_reference
export -f safe_write
export LOCKFD
export TEMP_DIR
export OUTPUT_FILE
export PYTHON_SCRIPT
export TARGET_SEQ
export REF_DIR

# Create list of all models
echo "Creating model list..."
{
    find "$AF2_DIR" -type f -name "*.pdb" -printf '%p,AF2\n'
    find "$AF3_DIR" -type f -name "*.cif" -printf '%p,AF3\n'
} > "$TEMP_DIR/all_models.txt"

total_models=$(wc -l < "$TEMP_DIR/all_models.txt")
echo "Total models to process: $total_models"
echo "Using $NUM_PROCESSES parallel processes"

# Process all models in parallel using GNU Parallel
if command -v parallel >/dev/null 2>&1; then
    # Use GNU Parallel if available
    parallel --progress --jobs "$NUM_PROCESSES" \
        process_model {1} {2} {#} :::: "$TEMP_DIR/all_models.txt" --colsep ','
else
    # Fallback to xargs if GNU Parallel is not available
    export MAX_PROCS=$NUM_PROCESSES
    cat "$TEMP_DIR/all_models.txt" | xargs -P "$NUM_PROCESSES" -I{} bash -c \
        'IFS=","; set -- {}; process_model "$1" "$2" "$$"'
fi

echo "Analysis complete. Results saved to $OUTPUT_FILE"

# Print summary statistics
echo -e "\nSummary Statistics:"
echo "===================="
echo "Number of models processed:"
echo "AF2: $(grep "AF2" "$OUTPUT_FILE" | wc -l)"
echo "AF3: $(grep "AF3" "$OUTPUT_FILE" | wc -l)"
echo -e "\nAverage RMSDs:"
echo "Binder RMSD: $(awk -F',' 'NR>1 {sum+=$3; count++} END {printf "%.2f\n", sum/count}' "$OUTPUT_FILE")"
echo "Channel RMSD: $(awk -F',' 'NR>1 {sum+=$4; count++} END {printf "%.2f\n", sum/count}' "$OUTPUT_FILE")"
