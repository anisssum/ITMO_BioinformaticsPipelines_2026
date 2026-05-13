#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.sra_id = ''
params.input_reads_folder = ''
params.reference = ''
params.assembly = false
params.outdir = './results'
params.threads = 2

include { download_sra } from './modules/processes/download_sra'
include { qc_raw  ; qc_trimmed} from './modules/processes/qc_raw'
include { trim_reads } from './modules/processes/trim_reads'
include { assemble_reads } from './modules/processes/assemble_reads'
include { map_reads } from './modules/processes/map_reads'
include { plot_coverage } from './modules/processes/plot_coverage'
include { BCFTOOLS_CALL } from './modules/nf-core/bcftools/call/main'
include { BCFTOOLS_MPILEUP } from './modules/nf-core/bcftools/mpileup/main'

workflow {

    if (params.sra_id) {
        reads_ch = download_sra(params.sra_id)
    }
    else {
        reads_ch = Channel.fromPath("${params.input_reads_folder}/*.fastq")
    }

    split_reads_ch = reads_ch.multiMap { reads ->
        qc_in: reads
        trim_in: reads
    }

    qc_raw(split_reads_ch.qc_in)

    trimmed_ch = trim_reads(split_reads_ch.trim_in)

    qc_trimmed(trimmed_ch)

    if (params.reference) {
        ref_ch = Channel.fromPath(params.reference)
    }
    else if (params.assembly) {
        ref_ch = assemble_reads(trimmed_ch)
    }

    bam_ch = map_reads(trimmed_ch, ref_ch)

    plot_coverage(bam_ch)

    bam_ref_ch = bam_ch.join(ref_ch)

    def split_channels = bam_ref_ch.multiMap { sample_id, bam, reference ->
        bam_input: tuple( [id: sample_id, single_end: false], bam, [], [] )
        fasta_input: tuple( [id: 'reference'], reference, [] )
    }

    mpileup_out = BCFTOOLS_MPILEUP(
        split_channels.bam_input,
        split_channels.fasta_input.first(),
        false)

    vcf_ch = mpileup_out.vcf.map { meta, vcf ->
        tuple( [id: meta.id], vcf )}

        BCFTOOLS_CALL(
        vcf_ch,                        
        split_channels.fasta_input.first(),  
        [],                            
        [])
}