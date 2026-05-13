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

workflow qc_trimmed {
    take: reads
    main: qc_raw(reads)
    emit: qc = qc_raw.out
}