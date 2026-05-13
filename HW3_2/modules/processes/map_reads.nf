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