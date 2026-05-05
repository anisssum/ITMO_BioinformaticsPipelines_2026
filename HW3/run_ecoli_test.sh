#!/bin/bash
# ============================================================
#  run_ecoli_test.sh — example launch commands for all profiles
# ============================================================

mkdir -p test_data

# --- LOCAL profile (uses existing conda env) -----------------
nextflow run main.nf \
    -profile local \
    --sra_id SRR33727007 \
    --assembly true \
    --outdir ./results_ecoli_local \
    --threads 4 \
    --conda_env_path /home/user/miniconda3/envs/bioenv \
    -resume

# --- CLUSTER profile (SLURM + conda from environment.yml) ----
# nextflow run main.nf \
#     -profile cluster \
#     --sra_id SRR33727007 \
#     --assembly true \
#     --outdir ./results_ecoli_cluster \
#     --threads 8 \
#     -resume

# --- CONTAINER profile (Docker) ------------------------------
# nextflow run main.nf \
#     -profile container \
#     --sra_id SRR33727007 \
#     --assembly true \
#     --outdir ./results_ecoli_docker \
#     --threads 4 \
#     -resume

# --- CONTAINER_SINGULARITY profile (HPC + Singularity) -------
# nextflow run main.nf \
#     -profile container_singularity \
#     --sra_id SRR33727007 \
#     --assembly true \
#     --outdir ./results_ecoli_singularity \
#     --threads 4 \
#     -resume
