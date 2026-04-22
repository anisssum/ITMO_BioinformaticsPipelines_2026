#!/bin/bash

mkdir -p test_data

nextflow run main.nf \
    --sra_id SRR33727007 \
    --assembly true \
    --outdir ./results_ecoli \
    --threads 4 \
    -resume