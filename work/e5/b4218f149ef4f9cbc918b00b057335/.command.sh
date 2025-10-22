#!/bin/bash -ue
# Check if input file has content
if [ -s "viruses.fna" ]; then
    # Create working directory
    mkdir -p 003_vcontact_temp

    # Run Prodigal to predict genes
    prodigal \
        -i viruses.fna \
        -a 003_vcontact_temp/003_proteins.faa \
        -o 003_vcontact_temp/003_genes.gff \
        -p meta

    # Create gene2genome mapping
    vcontact2_gene2genome \
        -p 003_vcontact_temp/003_proteins.faa \
        -o 003_vcontact_temp/gene2genome.csv \
        -s 'Prodigal-FAA'

    # Run vConTACT2
    vcontact2 \
        --raw-proteins 003_vcontact_temp/003_proteins.faa \
        --rel-mode 'Diamond' \
        --proteins-fp 003_vcontact_temp/gene2genome.csv \
        --db 'ProkaryoticViralRefSeq211-Merged' \
        --pcs-mode MCL \
        --vcs-mode ClusterONE \
        --c1-bin /mnt/scratch/users/as3243/miniconda/envs/vcontact2/lib/cluster_one-v1.0.jar \
        --output-dir 003_vcontact_out \
        --threads 16

    # Clean up temporary files
    rm -rf 003_vcontact_temp
else
    echo "Warning: Empty CheckV FASTA file for 003"
    mkdir -p 003_vcontact_out
    echo "No analysis - empty input file" > 003_vcontact_out/analysis_log.txt
fi
