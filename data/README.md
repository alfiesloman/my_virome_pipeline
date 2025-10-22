# Input Data Directory

This directory should contain your input FASTA files for virome analysis.

## Expected File Format:
- Files should be in FASTA format with `.fasta` extension
- Each file represents a sample
- File naming convention: `sample_name.fasta`

## Example:
```
data/
├── sample1.fasta
├── sample2.fasta
└── sample3.fasta
```

## File Requirements:
- Files must contain assembled contigs/sequences
- Headers should follow standard FASTA format
- Files should be readable and not corrupted

## To add your data:
1. Copy your FASTA files to this directory
2. Ensure they have `.fasta` extension
3. Run the pipeline: `./run_pipeline.sh`
