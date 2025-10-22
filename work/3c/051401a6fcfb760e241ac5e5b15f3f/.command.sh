#!/bin/bash -ue
export CUDA_VISIBLE_DEVICES=""  

# Create output directory
mkdir -p 003_genomad

# Run GeNomad
genomad end-to-end \
    --cleanup \
    --splits 8 \
    003.fasta \
    003_genomad \
    /mnt/scratch/projects/biol-soilv-2024/databases/genomad_db

# Validate output
if [ ! -f "003_genomad/003_summary/003_virus.fna" ]; then
    echo "Warning: No viruses found by GeNomad for 003"
    touch "003_genomad/003_summary/003_virus.fna"
fi
