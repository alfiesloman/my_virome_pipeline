#!/bin/bash -ue
# Create output directory
mkdir -p scaffolds_genomad

# Run GeNomad
genomad end-to-end \
    --cleanup \
    --splits 8 \
    scaffolds.fasta \
    scaffolds_genomad \
    /mnt/scratch/projects/biol-soilv-2024/databases/genomad_db

# Validate output
if [ ! -f "scaffolds_genomad/scaffolds_summary/scaffolds_virus.fna" ]; then
    echo "Warning: No viruses found by GeNomad for scaffolds"
    touch "scaffolds_genomad/scaffolds_summary/scaffolds_virus.fna"
fi
