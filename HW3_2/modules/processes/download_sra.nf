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