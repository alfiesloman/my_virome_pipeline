#!/bin/bash -ue
export CHECKVDB="/mnt/scratch/projects/biol-soilv-2024/checkv-db-v1.5"
# Create output directory
mkdir -p 003_checkv

# Check if input file is not empty
if [ -s "003_virus.fna" ]; then
    # Run CheckV
    checkv end_to_end \
        003_virus.fna \
        003_checkv \
        -t 4

    # Rename output for consistency
    if [ -f "003_checkv/proviruses.fna" ] && [ -f "003_checkv/viruses.fna" ]; then
        cat 003_checkv/proviruses.fna 003_checkv/viruses.fna > 003_checkv/viruses_combined.fna
        mv 003_checkv/viruses_combined.fna 003_checkv/viruses.fna
    elif [ -f "003_checkv/proviruses.fna" ]; then
        mv 003_checkv/proviruses.fna 003_checkv/viruses.fna
    elif [ ! -f "003_checkv/viruses.fna" ]; then
        touch 003_checkv/viruses.fna
    fi
else
    echo "Warning: Empty virus FASTA file for 003"
    mkdir -p 003_checkv
    touch 003_checkv/viruses.fna
fi
