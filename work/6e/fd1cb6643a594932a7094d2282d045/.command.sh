#!/bin/bash -ue
# Check if input file has content
if [ -s "viruses.fna" ]; then
    # Run BACPHLIP
    bacphlip \
        --multi_fasta viruses.fna \
        --local_hmmsearch /mnt/scratch/users/as3243/miniconda/envs/bacphlip/bin/hmmsearch \
        --threads 16 \
        --output_dir bacphlip_temp

    # Move results to expected location
    if [ -f "bacphlip_temp/bacphlip_results.tsv" ]; then
        mv bacphlip_temp/bacphlip_results.tsv 003_bacphlip_results.tsv
    else
        echo -e "sequence_id\tprediction\tconfidence" > 003_bacphlip_results.tsv
        echo "Warning: BACPHLIP produced no results for 003" >> 003_bacphlip_results.tsv
    fi
else
    echo "Warning: Empty CheckV FASTA file for 003"
    echo -e "sequence_id\tprediction\tconfidence" > 003_bacphlip_results.tsv
    echo "No predictions - empty input file" >> 003_bacphlip_results.tsv
fi
