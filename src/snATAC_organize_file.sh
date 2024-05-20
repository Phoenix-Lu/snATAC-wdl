#!/bin/bash

# Function to count occurrences of '1' and '2' in file names
count_occurrences() {
    local directory="$1"
    local max_1=0
    local max_2=0
    local file_1=""
    local file_2=""

    # Loop through files in directory
    while IFS= read -r file; do
        # Count occurrences of '1' and '2' in file name
        count_1=$(basename "$file" | tr -cd '1' | wc -c)
        count_2=$(basename "$file" | tr -cd '2' | wc -c)

        # Update max counts and corresponding file names
        if [ "$count_1" -gt "$max_1" ]; then
            max_1="$count_1"
            file_1="$file"
        fi

        if [ "$count_2" -gt "$max_2" ]; then
            max_2="$count_2"
            file_2="$file"
        fi
    done < <(find "$directory" -type f)

    echo "$file_1" "$file_2"
}

# Main script
main() {
    local directory="$1"
    declare -A result

    # Initialize arrays
    file1_array=()
    file2_array=()

    # Loop through subdirectories
    for subdir in "$directory"/*/; do
        subdir_name=$(basename "$subdir")
        result["sample"]+="\"$subdir_name\","  # Collect sample names

        # Count occurrences of '1' and '2' in file names
        files=$(count_occurrences "$subdir")

        # Get absolute paths for files
        file_1_abs=$(realpath "${files%% *}" 2>/dev/null)
        file_2_abs=$(realpath "${files##* }" 2>/dev/null)

        file1_array+=("\"$file_1_abs\"")  # Collect file1 paths
        file2_array+=("\"$file_2_abs\"")  # Collect file2 paths
    done

    # Remove trailing comma and assemble JSON output
    result["sample"]="${result["sample"]%,}"
    result["File1"]="${file1_array[@]}"
    result["File2"]="${file2_array[@]}"

    # Convert result to JSON format
    printf "{\n"
    printf '  "Read2": ['
    for ((i=0; i<${#file2_array[@]}; i++)); do
        if [ $i -eq $((${#file2_array[@]} - 1)) ]; then
            printf '%s' "${file2_array[i]}"
        else
            printf '%s, ' "${file2_array[i]}"
        fi
    done
    printf '],\n'
    printf '  "Read1": ['
    for ((i=0; i<${#file1_array[@]}; i++)); do
        if [ $i -eq $((${#file1_array[@]} - 1)) ]; then
            printf '%s' "${file1_array[i]}"
        else
            printf '%s, ' "${file1_array[i]}"
        fi
    done
    printf '],\n'
    printf '  "sample": [%s]\n' "${result["sample"]}"
    printf "}\n"
}

# Check for directory argument
if [ $# -ne 1 ]; then
    echo "Usage: $0 directory"
    exit 1
fi

# Run main script
main "$1"
