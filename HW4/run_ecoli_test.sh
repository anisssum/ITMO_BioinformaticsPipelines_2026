#!/bin/bash
# ============================================================
#  Example launch commands for legacy and multi-sample modes
# ============================================================

# --- MULTI-SAMPLE CSV design check (no tools required) -------
nextflow run main.nf \
    -profile stub \
    -stub-run \
    --samplesheet samplesheet.csv \
    --outdir ./results_stub \
    -resume

# --- LOCAL legacy single-sample profile ----------------------
# nextflow run main.nf \
#     -profile local \
#     --sra_id SRR33727007 \
#     --assembly true \
#     --outdir ./results_ecoli_local \
#     --threads 4 \
#     --conda_env_path /home/user/miniconda3/envs/bioenv \
#     -resume

# --- LOCAL multi-sample profile ------------------------------
# nextflow run main.nf \
#     -profile local \
#     --samplesheet samplesheet.csv \
#     --outdir ./results_multisample_local \
#     --conda_env_path /home/user/miniconda3/envs/bioenv \
#     -resume

# --- CLUSTER profile (SLURM + bundled env) -------------------
# nextflow run main.nf \
#     -profile cluster \
#     --samplesheet samplesheet.csv \
#     --outdir ./results_multisample_cluster \
#     -resume

# --- CONTAINER profile (Docker) ------------------------------
# nextflow run main.nf \
#     -profile container \
#     --samplesheet samplesheet.csv \
#     --outdir ./results_multisample_docker \
#     -resume

# --- CONTAINER_SINGULARITY profile ---------------------------
# nextflow run main.nf \
#     -profile container_singularity \
#     --samplesheet samplesheet.csv \
#     --outdir ./results_multisample_singularity \
#     -resume
