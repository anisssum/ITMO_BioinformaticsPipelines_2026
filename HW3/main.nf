#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ============================================================
// PROCESS 1: Download SRA data
// ============================================================
process download_sra {
    tag "${sra_id}"
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

// ============================================================
// PROCESS 2: Run QC (FastQC)
// ============================================================
process run_qc {
    tag "${reads_type}"
    publishDir "${params.outdir}/${reads_type}_qc", mode: 'copy'

    input:
        val  reads_type
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

// ============================================================
// PROCESS 3: Trim reads (fastp)
// ============================================================
process trim_reads {
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

    input:
        path reads

    output:
        path "trimmed.fastq"

    script:
    def cpus = task.cpus ?: 2
    """
    fastp -i ${reads} -o trimmed.fastq \
        -q 15 \
        -l 30 \
        --thread ${cpus} \
        --length_required 30

    if [ ! -s trimmed.fastq ]; then
        echo "WARNING: Trimming produced empty file, using original reads"
        cp ${reads} trimmed.fastq
    fi

    echo "Trimmed reads size:"
    wc -l trimmed.fastq
    """
}

// ============================================================
// PROCESS 4: Assemble reads (SPAdes)
// ============================================================
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
        echo "ERROR: Reads file is empty, creating dummy reference"
        echo ">dummy_contig" > assembly_contigs.fasta
        echo "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" >> assembly_contigs.fasta
    else
        spades.py -s ${reads} -o assembly \
            -t ${cpus} \
            -m 8 \
            --only-assembler \
            -k 21,33,55 \
            --cov-cutoff auto || true

        if [ -f assembly/contigs.fasta ] && [ -s assembly/contigs.fasta ]; then
            cp assembly/contigs.fasta assembly_contigs.fasta
            echo "Assembly successful. Contigs:"
            grep -c "^>" assembly_contigs.fasta
        else
            echo "WARNING: Assembly failed, using dummy reference"
            echo ">dummy_contig" > assembly_contigs.fasta
            echo "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" >> assembly_contigs.fasta
        fi
    fi
    """
}

// ============================================================
// PROCESS 5: Map reads to reference (BWA + samtools)
// ============================================================
process map_reads {
    publishDir "${params.outdir}/mapping", mode: 'copy'

    input:
        path reads
        path reference

    output:
        path "mapped.bam"
        path "mapped.bam.bai"

    script:
    def cpus = task.cpus ?: 2
    """
    if [ ! -s ${reference} ]; then
        echo ">dummy" > ref.fasta
        echo "ATCGATCGATCGATCGATCGATCGATCGATCG" >> ref.fasta
        REF=ref.fasta
    else
        REF=${reference}
    fi

    bwa index \$REF
    bwa mem -t ${cpus} \$REF ${reads} > aligned.sam
    samtools view -bS aligned.sam > unsorted.bam
    samtools sort -@ ${cpus} unsorted.bam -o mapped.bam
    samtools index mapped.bam

    echo "Mapping statistics:"
    samtools flagstat mapped.bam
    """
}

// ============================================================
// PROCESS 6: Plot coverage (samtools depth + R/ggplot2)
// ============================================================
process plot_coverage {
    publishDir "${params.outdir}/coverage", mode: 'copy'

    input:
        path bam
        path bai

    output:
        path "coverage_plot.png",  emit: plot
        path "coverage_stats.txt", emit: stats

    script:
    """
    samtools depth ${bam} > depth.txt

    printf 'library(ggplot2)\\n' > plot.R
    printf 'd <- read.table("depth.txt", header=F)\\n' >> plot.R
    printf 'names(d) <- c("contig", "pos", "depth")\\n' >> plot.R
    printf 'if(nrow(d)==0) d <- data.frame(contig="dummy", pos=1:100, depth=0)\\n' >> plot.R
    printf 'if(nrow(d)>50000) d <- d[sample(1:nrow(d), 50000),]\\n' >> plot.R
    printf 'p <- ggplot(d, aes(pos, depth)) + geom_line(color="blue") + theme_bw()\\n' >> plot.R
    printf 'ggsave("coverage_plot.png", p, width=10, height=6)\\n' >> plot.R
    printf 'write.table(data.frame(mean=mean(d\$depth), max=max(d\$depth)), "coverage_stats.txt", row.names=F)\\n' >> plot.R

    Rscript plot.R
    """
}

// ============================================================
// PROCESS 7: Variant calling
// ============================================================
process call_variants {
    tag "bcftools"
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
        path bam
        path bai
        path reference

    output:
        path "variants.vcf.gz",     emit: vcf
        path "variants.vcf.gz.tbi", emit: tbi
        path "variants_stats.txt",  emit: stats

    script:
    def cpus = task.cpus ?: 2
    """
    # Index reference if not indexed
    if [ ! -f ${reference}.fai ]; then
        samtools faidx ${reference}
    fi

    # Pileup + variant calling in one pipe
    bcftools mpileup \
        --fasta-ref ${reference} \
        --min-BQ 20 \
        --min-MQ 20 \
        --annotate FORMAT/DP,FORMAT/AD \
        --output-type u \
        ${bam} \
    | bcftools call \
        --multiallelic-caller \
        --variants-only \
        --output-type z \
        --output variants.vcf.gz

    # Index VCF
    bcftools index --tbi variants.vcf.gz

    # Basic stats
    bcftools stats variants.vcf.gz > variants_stats.txt

    echo "Variant calling complete."
    echo "SNPs + indels called:"
    bcftools view --type snps  variants.vcf.gz | grep -v "^#" | wc -l | xargs echo "  SNPs:"
    bcftools view --type indels variants.vcf.gz | grep -v "^#" | wc -l | xargs echo "  Indels:"
    """
}

// ============================================================
// PROCESS 8: Filter variants
// ============================================================
process filter_variants {
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
        path vcf
        path tbi

    output:
        path "variants_filtered.vcf.gz",     emit: vcf
        path "variants_filtered.vcf.gz.tbi", emit: tbi

    script:
    """
    bcftools filter \
        --include 'QUAL>=${params.min_variant_qual} && INFO/DP>=${params.min_depth}' \
        --output-type z \
        --output variants_filtered.vcf.gz \
        ${vcf}

    bcftools index --tbi variants_filtered.vcf.gz

    echo "Filtered variants:"
    bcftools view variants_filtered.vcf.gz | grep -v "^#" | wc -l
    """
}

// ============================================================
// Sub-workflows
// ============================================================
workflow qc_initial {
    take: reads
    main: run_qc('initial', reads)
    emit: qc = run_qc.out
}

workflow qc_trimmed {
    take: reads
    main: run_qc('trimmed', reads)
    emit: qc = run_qc.out
}

workflow assemble {
    take: reads
    main: assemble_reads(reads)
    emit: reference = assemble_reads.out
}

// ============================================================
// MAIN WORKFLOW
// ============================================================
workflow {

    // ----------------------------------------------------------
    // Step 1: Obtain raw reads
    // ----------------------------------------------------------
    if (params.sra_id) {
        raw_reads = download_sra(params.sra_id)
    } else if (params.input_reads_folder) {
        raw_reads = Channel.fromPath("${params.input_reads_folder}/*.fastq")
    } else {
        error "Provide either --sra_id or --input_reads_folder"
    }

    // ----------------------------------------------------------
    // Step 2: Initial QC
    // ----------------------------------------------------------
    qc_initial(raw_reads)

    // ----------------------------------------------------------
    // Step 3: Trim
    // ----------------------------------------------------------
    trimmed_reads = trim_reads(raw_reads)

    // ----------------------------------------------------------
    // Step 4: Post-trim QC
    // ----------------------------------------------------------
    qc_trimmed(trimmed_reads)

    // ----------------------------------------------------------
    // Step 5: Reference
    // ----------------------------------------------------------
    if (params.reference) {
        reference_ch = Channel.fromPath(params.reference)
    } else if (params.assembly) {
        reference_ch = assemble(trimmed_reads).reference
    } else {
        error "Provide --reference <path> or set --assembly true"
    }

    // ----------------------------------------------------------
    // Step 6: Map reads
    // ----------------------------------------------------------
    map_out = map_reads(trimmed_reads, reference_ch)
    bam_ch  = map_out[0]
    bai_ch  = map_out[1]

    // ----------------------------------------------------------
    // Step 7: Coverage plot
    // ----------------------------------------------------------
    plot_coverage(bam_ch, bai_ch)

    // ----------------------------------------------------------
    // Step 8: Variant calling
    // ----------------------------------------------------------
    vc_out = call_variants(bam_ch, bai_ch, reference_ch)

    // ----------------------------------------------------------
    // Step 9: Filter variants (optional, enabled by default)
    // ----------------------------------------------------------
    if (params.filter_variants) {
        filter_variants(vc_out.vcf, vc_out.tbi)
    }
}
