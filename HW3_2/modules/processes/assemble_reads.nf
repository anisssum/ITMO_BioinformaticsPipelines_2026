
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