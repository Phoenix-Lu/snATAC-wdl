#!/bin/bash

usage() {
    echo "Usage: $0 <fq1> <fq2> <portion>"
    echo "  fq1: Fastq file 1 (with .fq.gz extension)"
    echo "  fq2: Fastq file 2 (with .fq.gz extension)"
    echo "  portion: Portion to downsample (e.g., 0.1 for 10%, 0.5 for 50%)"
}

if [ "$#" -ne 3 ]; then 
    echo "Error: Invalid number of arguments."
    usage 
    exit 1
fi

fq1=$1
fq2=$2
portion=$3

fq1_basename=$(basename "$fq1" '.fq.gz')
fq2_basename=$(basename "$fq2" '.fq.gz')
fq1_downsampled="${fq1_basename}_${portion}_sampled.fq"
fq2_downsampled="${fq2_basename}_${portion}_sampled.fq"

# Check if input files exist
if [ ! -f "$fq1" ] || [ ! -f "$fq2" ]; then
    echo "Error: Input files not found."
    usage
    exit 1
fi

seqtk sample -s 11 "$fq1" "$portion" > "$fq1_downsampled"
seqtk sample -s 11 "$fq2" "$portion" > "$fq2_downsampled"

echo "Downsampled files created:"
echo "$fq1_downsampled"
echo "$fq2_downsampled"