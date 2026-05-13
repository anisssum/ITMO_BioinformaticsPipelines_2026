#!/bin/bash

nextflow run main.nf \
    -profile local \
    --input_reads_folder ./results_local/raw_reads/ \
    --assembly true \
    --outdir results_local
	
nextflow run main.nf \
    -profile cluster \
    --input_reads_folder ./results_cluster/raw_reads/ \
    --assembly true \
    --outdir results_cluster

nextflow run main.nf \
    -profile container \
    --input_reads_folder ./results_container/raw_reads/ \
    --assembly true \
    --outdir results_container