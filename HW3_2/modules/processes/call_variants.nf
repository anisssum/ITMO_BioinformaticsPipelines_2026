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
