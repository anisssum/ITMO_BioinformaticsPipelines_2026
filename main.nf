#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Pipeline parameters
params.input_reads_folder = ''
params.sra_id = ''
params.reference = ''
params.assembly = false
params.outdir = './results'
params.threads = 2

// Process 1: Download SRA data
process download_sra {
    publishDir "${params.outdir}/raw_reads", mode: 'copy'

    input:
        val sra_id
    
    output:
        path "downloaded_reads.fastq"
    
    script:
    """
    echo "Downloading ${sra_id}..."
    prefetch ${sra_id} --max-size 500M
    fasterq-dump ${sra_id} --threads ${task.cpus}
    
    if [ -f ${sra_id}.fastq ]; then
        cp ${sra_id}.fastq downloaded_reads.fastq
    elif [ -f ${sra_id}_1.fastq ]; then
        cat ${sra_id}_1.fastq ${sra_id}_2.fastq > downloaded_reads.fastq
    else
        echo "ERROR: No FASTQ files found"
        exit 1
    fi
    
    echo "Download complete"
    """
}

// Process 2: Run QC
process run_qc {
    publishDir "${params.outdir}/${reads_type}_qc", mode: 'copy'

    input:
        val reads_type
        path reads
    
    output:
        path "${reads_type}_qc_report/"
    
    script:
    def cpus = task.cpus ?: 2
    """
    mkdir -p ${reads_type}_qc_report
    fastqc -o ${reads_type}_qc_report/ -t ${cpus} ${reads}
    """
}

// Process 3: Trim reads with fastp
process trim_reads {
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

    input:
        path reads
    
    output:
        path "trimmed.fastq"
    
    script:
    def cpus = task.cpus ?: 2
    """
    # Less stringent trimming: lower quality threshold, shorter min length
    fastp -i ${reads} -o trimmed.fastq \
        -q 15 \
        -l 30 \
        --thread ${cpus} \
        --length_required 30
    
    echo "Trimmed reads size:"
    wc -l trimmed.fastq
    """
}

// Process 4: Assemble reads
process assemble_reads {
    publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
        path reads

    output:
        path "assembly_contigs.fasta"

    script:
    def cpus = task.cpus ?: 2
    """
    if [ ! -s ${reads} ]; then
        echo "ERROR: Reads file is empty, cannot perform assembly"
        exit 1
    fi

    # Run SPAdes assembly
    spades.py -s ${reads} -o assembly \
        -t ${cpus} \
        -m 8 \
        --only-assembler \
        -k 21,33,55 \
        --cov-cutoff auto

    if [ -f assembly/contigs.fasta ] && [ -s assembly/contigs.fasta ]; then
        cp assembly/contigs.fasta assembly_contigs.fasta
        echo "Assembly successful. Number of contigs:"
        grep -c "^>" assembly_contigs.fasta
    else
        echo "ERROR: Assembly failed, please check input reads and SPAdes logs"
        exit 1
    fi
    """
}

// Process 5: Map reads to reference
process map_reads {
    publishDir "${params.outdir}/mapping", mode: 'copy'

    input:
        path reads
        path reference
    
    output:
        path "mapped.bam"
    
    script:
    def cpus = task.cpus ?: 2
    """
    # Validate reference file
    if [ ! -s ${reference} ]; then
        echo "ERROR: Reference file is empty or does not exist: ${reference}"
        exit 1
    fi
    
    if ! grep -q "^>" ${reference}; then
        echo "ERROR: Reference file does not appear to be in FASTA format"
        exit 1
    fi
    
    bwa index ${reference}
    bwa mem -t ${cpus} ${reference} ${reads} > aligned.sam
    samtools view -bS aligned.sam > unsorted.bam
    samtools sort -@ ${cpus} unsorted.bam -o mapped.bam
    samtools index mapped.bam
    
    echo "Mapping statistics:"
    samtools flagstat mapped.bam
    """
}

// Process 6: Plot coverage with R
process plot_coverage {
    publishDir "${params.outdir}/coverage", mode: 'copy'

    input:
        path bam
    
    output:
        path "coverage_plot.png"
        path "coverage_stats.txt"
    
    script:
    """
    echo "Calculating coverage..."
    samtools depth ${bam} > depth.txt
    
    # Create R script using a heredoc with 'EOF' (prevents Groovy interpolation)
    cat > plot.R << 'EOF'
    library(ggplot2)
    d <- read.table("depth.txt", header=FALSE)
    if(nrow(d) == 0) {
        d <- data.frame(V1="dummy", V2=1:100, V3=0)
    }
    names(d) <- c("contig", "pos", "depth")
    if(nrow(d) > 50000) {
        set.seed(123)
        d <- d[sample(1:nrow(d), 50000),]
    }
    p <- ggplot(d, aes(x=pos, y=depth)) + 
        geom_line(color="steelblue") + 
        theme_bw()
    ggsave("coverage_plot.png", p, width=10, height=6)
    stats <- data.frame(
        mean_depth = mean(d[["depth"]]),
        median_depth = median(d[["depth"]]),
        max_depth = max(d[["depth"]]),
        min_depth = min(d[["depth"]])
    )
    write.table(stats, "coverage_stats.txt", row.names=FALSE, quote=FALSE, sep="\t")
    cat("Coverage plot saved\\n")
    EOF
    
    Rscript plot.R
    """
}

// Workflow for initial QC
workflow qc_initial {
    take:
        reads
    
    main:
        run_qc('initial', reads)
    
    emit:
        qc = run_qc.out
}

// Workflow for trimmed QC
workflow qc_trimmed {
    take:
        reads
    
    main:
        run_qc('trimmed', reads)
    
    emit:
        qc = run_qc.out
}

// Workflow for assembly
workflow assemble {
    take:
        reads
    
    main:
        assemble_reads(reads)
    
    emit:
        reference = assemble_reads.out
}

// Main workflow
workflow {
    // Step 1: Get input reads
    if (params.sra_id) {
        raw_reads = download_sra(params.sra_id)
    } else if (params.input_reads_folder) {
        raw_reads = Channel.fromPath("${params.input_reads_folder}/*.fastq")
    } else {
        error "Either --sra_id or --input_reads_folder must be provided"
    }
    
    // Step 2: QC on raw reads
    qc_initial(raw_reads)
    
    // Step 3: Trim reads
    trimmed_reads = trim_reads(raw_reads)
    
    // Step 4: QC on trimmed reads
    qc_trimmed(trimmed_reads)
    
    // Step 5: Get reference
    if (params.reference) {
        reference_ch = Channel.fromPath(params.reference)
    } else if (params.assembly) {
        reference_ch = assemble(trimmed_reads).reference
    } else {
        error "Either --reference or --assembly must be true"
    }
    
    // Step 6: Map trimmed reads to reference
    mapping_result = map_reads(trimmed_reads, reference_ch)
    
    // Step 7: Plot coverage
    plot_coverage(mapping_result)
}