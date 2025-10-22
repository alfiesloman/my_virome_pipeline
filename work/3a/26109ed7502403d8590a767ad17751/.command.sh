#!/bin/bash -ue
# Create output directory
mkdir -p scaffolds_checkv

# Check if input file is not empty
if [ -s "scaffolds_virus.fna" ]; then
    # Run CheckV
    checkv end_to_end \
        scaffolds_virus.fna \
        scaffolds_checkv \
        -t 8

    # Rename output for consistency
    if [ -f "scaffolds_checkv/proviruses.fna" ] && [ -f "scaffolds_checkv/viruses.fna" ]; then
        cat scaffolds_checkv/proviruses.fna scaffolds_checkv/viruses.fna > scaffolds_checkv/viruses_combined.fna
        mv scaffolds_checkv/viruses_combined.fna scaffolds_checkv/viruses.fna
    elif [ -f "scaffolds_checkv/proviruses.fna" ]; then
        mv scaffolds_checkv/proviruses.fna scaffolds_checkv/viruses.fna
    elif [ ! -f "scaffolds_checkv/viruses.fna" ]; then
        touch scaffolds_checkv/viruses.fna
    fi
else
    echo "Warning: Empty virus FASTA file for scaffolds"
    mkdir -p scaffolds_checkv
    touch scaffolds_checkv/viruses.fna
fi
