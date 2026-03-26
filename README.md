# Virome Analysis Pipeline (Nextflow)

Production-ready virome analysis workflow implemented in Nextflow DSL2.

## What it runs

Given per-sample contig FASTA files (`*.fasta`), the pipeline runs the following stages (from `main.nf`):

1. `GeNomad` (`RUN_GENOMAD`) - virus identification and annotation
2. `CheckV` (`RUN_CHECKV`) - virus quality assessment
3. `iPHoP` (`RUN_IPHOP`) - host prediction for viruses
4. `BACPHLIP` (`RUN_BACPHLIP`) - phage lifestyle prediction
5. `vConTACT2` chain
   - `Prodigal` (`RUN_VCONTACT_PRODIGAL`) - gene/protein prediction
   - `vcontact2_gene2genome` (`RUN_VCONTACT_GENE2GENOME`) - gene2genome mapping
   - `vcontact2` main clustering/taxonomy (`RUN_VCONTACT_MAIN`)

Tools included in this repo:
- GeNomad
- CheckV
- iPHoP
- BACPHLIP
- vConTACT2

## Inputs

Place input FASTA files in your input directory (default is `data/`).

Expected format (see `data/README.md`):
- FASTA files with the `.fasta` extension
- Each file corresponds to one sample
- File naming convention: `sample_name.fasta`

The Nextflow workflow reads input as:
- `${params.input_dir}/*.fasta`
- `sample_id = file.baseName`

## Databases and configured paths

Database/tool paths are defined in `nextflow.config` under `params` (and validated at runtime by `main.nf` for these parameters):
- `--genomad_db`
- `--iphop_db`
- `--bacphlip_hmmsearch`
- `--cluster_one_jar`

Other configured paths (used by the workflow):
- `--checkv_db` (used by CheckV; defined in `nextflow.config`)

If you need to override these paths at run time, `main.nf` exposes the following options in its help message:
- `--genomad_db`
- `--iphop_db`
- `--bacphlip_hmmsearch`
- `--cluster_one_jar`

## Outputs

Outputs are published under your `--outdir` (configured in `nextflow.config` as `params.outdir`), and per-sample directories are created for each tool stage.

Important paths:
- `--outdir/<sample_id>/...` (tool outputs)
- `--outdir/pipeline_info/` execution reports and traces (from `nextflow.config`):
  - `execution_report.html`
  - `execution_timeline.html`
  - `execution_trace.txt`
  - `pipeline_dag.html`

Notes on typical filenames produced by stages:
- GeNomad: `${sample_id}_genomad/.../${contig}_virus.fna` and summary TSVs
- CheckV: `viruses.fna` (and the CheckV output directory)
- iPHoP: host prediction CSV (e.g. `Host_prediction_to_genus_m90.csv` if inputs are non-empty)
- BACPHLIP: `${sample_id}_bacphlip_results.tsv` (and a `${sample_id}_bacphlip.log`)
- vConTACT2 chain:
  - `${sample_id}_proteins.faa` and `${sample_id}_genes.gff` (Prodigal)
  - `gene2genome.csv` (gene2genome step)
  - `${sample_id}_vcontact_out/` (vConTACT2 main outputs)

Exact on-disk layout depends on Nextflow `publishDir` behavior, but all results will be under `--outdir`.

## Running the pipeline

This repo includes a wrapper script: `run_pipeline.sh`.

### Wrapper usage

By default, the wrapper runs:
- profile: `slurm`
- input: `data/`
- output: `results_YYYYmmdd_HHMMSS` (timestamped)

Run:
```bash
./run_pipeline.sh
```

Common options:
```bash
./run_pipeline.sh -i data -o results -p slurm
./run_pipeline.sh -i my_fastas -o results -p test
./run_pipeline.sh -i data -o results -p conda -c
./run_pipeline.sh -i data -o results -p slurm -r
./run_pipeline.sh -i data -o results -p slurm -n   # dry-run
```

Wrapper options (from `run_pipeline.sh`):
- `-i, --input DIR` (default: `data`)
- `-o, --output DIR` (default: `results_TIMESTAMP`)
- `-p, --profile PROF` (default: `slurm`)
- `-c, --conda` enable Nextflow `-with-conda`
- `-r, --resume` pass Nextflow `-resume`
- `-n, --dry-run` print the Nextflow command but do not execute
- `-h, --help` show wrapper help

Any unrecognized flags are passed directly to `nextflow run main.nf ...`.
For example: `-process.max_cpus 4` or `-with-tower`.

### Profiles

`nextflow.config` defines profiles including:
- `slurm` (cluster execution; `process.executor = 'slurm'`)
- `test` (uses internal `test_data/` and reduced resource limits, but still validates database paths at startup)
- `conda` / `mamba` (enables conda environments; Singularity/Docker disabled)
- `viking2` (includes `conf/viking2.config`; you will need that file in your worktree for this profile to work)

## Resumability, tracing, and work directory

From `nextflow.config`:
- `resume = true` (Nextflow will attempt to resume)
- `cleanup = true`
- `workDir` is set to:
  - `/mnt/scratch/users/as3243/nextflow_work`
- Tracing is enabled and written into `--outdir/pipeline_info/`

## Important note (debug selector)

`nextflow.config` currently includes:
- `process.selector = 'RUN_BACPHLIP'` (marked as “Delete after trouble shooting”)

If left enabled, Nextflow may restrict execution to only the `RUN_BACPHLIP` process.
Remove/comment this line to run the full pipeline stages.

