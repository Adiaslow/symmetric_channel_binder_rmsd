#!/bin/bash

# Default values
PYTHON_SCRIPT="symmetric_channel_binder_rmsd.py"
TARGET_SEQ="IPDAFWWAVVTMTTVGYGD"
OUTPUT_FILE="analysis_results.csv"

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

# Create output CSV header
echo "Model,AF_Version,Binder_RMSD,Channel_RMSD,Binder_Length_Ref,Binder_Length_Model,Channel_Length,Reference_Path,Model_Path" > "$OUTPUT_FILE"

# Function to find reference file
find_reference() {
    local model_num="$1"
    local sample_num="$2"
    local ref_pattern="*${model_num}*sample_${sample_num}*/interface.pdb"
    find "$REF_DIR" -type f -path "*/$ref_pattern" -print -quit
}

# Function to extract model and sample numbers
extract_numbers() {
    local filename="$1"
    local af_version="$2"
    local model_num=""
    local sample_num=""

    if [[ $af_version == "AF2" ]]; then
        model_num=$(echo "$filename" | grep -o '_[0-9]\+_sample_[0-9]\+_' | cut -d'_' -f2)
        sample_num=$(echo "$filename" | grep -o '_sample_[0-9]\+_' | cut -d'_' -f3)
    else  # AF3
        # Updated pattern for AF3 CIF files
        model_num=$(echo "$filename" | grep -o '_[0-9]\+_[0-9]\+_' | cut -d'_' -f2)
        sample_num=$(echo "$filename" | grep -o '_[0-9]\+_[0-9]\+_' | cut -d'_' -f3)
    fi

    echo "$model_num $sample_num"
}

# Process model function
process_model() {
    local model_path="$1"
    local af_version="$2"
    local model_name=$(basename "$model_path")

    # Extract model and sample numbers
    read -r model_num sample_num <<< "$(extract_numbers "$model_name" "$af_version")"

    if [ -z "$model_num" ] || [ -z "$sample_num" ]; then
        echo "Warning: Could not extract model/sample numbers from $model_name"
        return
    fi

    # Find corresponding reference
    reference_path=$(find_reference "$model_num" "$sample_num")

    if [ -z "$reference_path" ]; then
        echo "Warning: No reference found for model $model_name"
        return
    fi

    echo "Processing $af_version model ${model_num}_${sample_num}..."

    # Run the Python script and capture output
    output=$(python3 "$PYTHON_SCRIPT" \
        -refpdb "$reference_path" \
        -refchanchain B \
        -refbindchain A \
        -testpdb "$model_path" \
        -testchanchain A \
        -testbindchain B \
        -targetseq "$TARGET_SEQ")

    # Check if the Python script executed successfully
    if [ $? -ne 0 ]; then
        echo "Warning: Analysis failed for model $model_name"
        return
    fi

    # Extract values from output
    binder_rmsd=$(echo "$output" | grep "Binder RMSD:" | awk '{print $3}')
    channel_rmsd=$(echo "$output" | grep "Channel RMSD:" | awk '{print $3}')
    binder_length_ref=$(echo "$output" | grep "Binder length (reference):" | awk '{print $4}')
    binder_length_model=$(echo "$output" | grep "Binder length (test):" | awk '{print $4}')
    channel_length=$(echo "$output" | grep "Channel length:" | awk '{print $3}')

    # Add to CSV
    echo "${model_num}_${sample_num},$af_version,$binder_rmsd,$channel_rmsd,$binder_length_ref,$binder_length_model,$channel_length,$reference_path,$model_path" >> "$OUTPUT_FILE"
}

# Process AF2 models
echo "Processing AF2 models..."
find "$AF2_DIR" -type f -name "*.pdb" | while read -r model_path; do
    process_model "$model_path" "AF2"
done

# Process AF3 models (updated to look for CIF files)
echo "Processing AF3 models..."
find "$AF3_DIR" -type f -name "*.cif" | while read -r model_path; do
    process_model "$model_path" "AF3"
done

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
