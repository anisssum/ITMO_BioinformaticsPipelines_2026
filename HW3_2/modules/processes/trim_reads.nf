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