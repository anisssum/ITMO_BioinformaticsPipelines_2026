#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.input_reads_folder = ''
params.sra_id = ''
params.reference = ''
params.assembly = false
params.outdir = './results'
params.threads = 2

process download_sra {

    publishDir "${params.outdir}/raw_reads", mode: 'copy'

    input:
    val sra_id

    output:
    path "*.fastq"

    script:
    """
    /home/anya/sratoolkit.3.4.1-ubuntu64/bin/fasterq-dump ${sra_id}
    """
}

process qc_raw {

    publishDir "${params.outdir}/qc_raw", mode: 'copy'

    input:
    path reads

    output:
    path "*.html"

    script:
    """
    fastqc ${reads}
    """
}

process trim_reads {

    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    path reads

    output:
    path "trimmed.fastq"

    script:
    """
    fastp -i ${reads} -o trimmed.fastq
    """
}

process assemble_reads {

    publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
    path reads

    output:
    path "contigs.fasta"

    script:
    """
    spades.py \
        -s ${reads} \
        -o spades_out

    cp spades_out/contigs.fasta contigs.fasta
    """
}

process map_reads {

    publishDir "${params.outdir}/mapping", mode: 'copy'

    input:
    path reads
    path reference

    output:
    path "mapped.bam"

    script:
    """
    bwa index ${reference}

    bwa mem ${reference} ${reads} > aln.sam

    samtools view -bS aln.sam > aln.bam

    samtools sort aln.bam -o mapped.bam

    samtools index mapped.bam
    """
}

process call_variants {

    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    path bam
    path reference

    output:
    path "variants.vcf"

    script:
    """
    samtools faidx ${reference}

    bcftools mpileup \
        -Ou \
        -f ${reference} \
        ${bam} | \
    bcftools call \
        -mv \
        -Ov \
        -o variants.vcf
    """
}

workflow {

    if (params.sra_id) {
        reads_ch = download_sra(params.sra_id)
    }
    else {
        reads_ch = Channel.fromPath("${params.input_reads_folder}/*.fastq")
    }

    qc_raw(reads_ch)

    trimmed_ch = trim_reads(reads_ch)

    if (params.reference) {
        ref_ch = Channel.fromPath(params.reference)
    }
    else if (params.assembly) {
        ref_ch = assemble_reads(trimmed_ch)
    }

    bam_ch = map_reads(trimmed_ch, ref_ch)

    call_variants(bam_ch, ref_ch)
}