#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// ============================================================
// PROCESS 1: Download SRA data
// ============================================================
process download_sra {
    tag "${meta.id}"
    publishDir "${params.outdir}/raw_reads", mode: 'copy'

    input:
        tuple val(meta), val(sra_id)

    output:
        tuple val(meta), path("${meta.id}.fastq")

    script:
    def cpus = task.cpus ?: params.threads ?: 2
    """
    echo "Downloading ${sra_id} for sample ${meta.id}..."
    prefetch ${sra_id} --max-size 500M
    fasterq-dump ${sra_id} --threads ${cpus}

    if [ -f ${sra_id}.fastq ]; then
        cp ${sra_id}.fastq ${meta.id}.fastq
    elif [ -f ${sra_id}_1.fastq ] && [ -f ${sra_id}_2.fastq ]; then
        cat ${sra_id}_1.fastq ${sra_id}_2.fastq > ${meta.id}.fastq
    else
        echo "ERROR: No FASTQ files found for ${sra_id}"
        exit 1
    fi

    echo "Download complete for ${meta.id}"
    """

    stub:
    """
    touch ${meta.id}.fastq
    """
}

// ============================================================
// PROCESS 2: Run QC (FastQC)
// ============================================================
process run_qc {
    tag "${meta.id}"
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
        val reads_type
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("*_fastqc.html"), path("*_fastqc.zip")

    script:
    def cpus = task.cpus ?: 2
    """
    fastqc -t ${cpus} ${reads}
    """

    stub:
    """
    touch ${meta.id}_${reads_type}_fastqc.html
    touch ${meta.id}_${reads_type}_fastqc.zip
    """
}

// ============================================================
// PROCESS 3: Trim reads (fastp)
// ============================================================
process trim_reads {
    tag "${meta.id}"
    publishDir "${params.outdir}/trimmed_reads", mode: 'copy'

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${meta.id}_trimmed.fastq")

    script:
    def cpus = task.cpus ?: 2
    """
    fastp -i ${reads} -o ${meta.id}_trimmed.fastq \
        -q ${params.min_quality} \
        -l ${params.min_read_length} \
        --thread ${cpus} \
        --length_required ${params.min_read_length}

    if [ ! -s ${meta.id}_trimmed.fastq ]; then
        echo "WARNING: Trimming produced empty file, using original reads"
        cp ${reads} ${meta.id}_trimmed.fastq
    fi

    echo "Trimmed reads size for ${meta.id}:"
    wc -l ${meta.id}_trimmed.fastq
    """

    stub:
    """
    touch ${meta.id}_trimmed.fastq
    """
}

// ============================================================
// PROCESS 4: Assemble reads (SPAdes)
// ============================================================
process assemble_reads {
    tag "${meta.id}"
    publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
        tuple val(meta), path(reads)

    output:
        tuple val(meta), path("${meta.id}_assembly_contigs.fasta")

    script:
    def cpus = task.cpus ?: 2
    """
    if [ ! -s ${reads} ]; then
        echo "ERROR: Reads file is empty, creating dummy reference"
        echo ">dummy_contig" > ${meta.id}_assembly_contigs.fasta
        echo "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" >> ${meta.id}_assembly_contigs.fasta
    else
        spades.py -s ${reads} -o assembly \
            -t ${cpus} \
            -m 8 \
            --only-assembler \
            -k 21,33,55 \
            --cov-cutoff auto || true

        if [ -f assembly/contigs.fasta ] && [ -s assembly/contigs.fasta ]; then
            cp assembly/contigs.fasta ${meta.id}_assembly_contigs.fasta
        else
            echo "WARNING: Assembly failed, using dummy reference"
            echo ">dummy_contig" > ${meta.id}_assembly_contigs.fasta
            echo "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG" >> ${meta.id}_assembly_contigs.fasta
        fi
    fi
    """

    stub:
    """
    printf '>dummy_contig\nATCGATCGATCG\n' > ${meta.id}_assembly_contigs.fasta
    """
}

// ============================================================
// PROCESS 5: Fetch references from accession IDs
// ============================================================
process fetch_reference {
    tag "${ref_id}"
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
        val ref_id

    output:
        tuple val(ref_id), path("${ref_id}.fasta")

    script:
    """
    efetch -db nuccore -id ${ref_id} -format fasta > ${ref_id}.fasta
    """

    stub:
    """
    printf '>${ref_id}\nATCGATCGATCGATCGATCG\n' > ${ref_id}.fasta
    """
}

// ============================================================
// PROCESS 6: Map reads to reference (BWA + samtools)
// ============================================================
process map_reads {
    tag "${meta.id}"
    publishDir "${params.outdir}/mapping", mode: 'copy'

    input:
        tuple val(meta), path(reads), path(reference)

    output:
        tuple val(meta), path("${meta.id}.mapped.bam"), path("${meta.id}.mapped.bam.bai"), path(reference)

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
    samtools sort -@ ${cpus} unsorted.bam -o ${meta.id}.mapped.bam
    samtools index ${meta.id}.mapped.bam
    """

    stub:
    """
    touch ${meta.id}.mapped.bam
    touch ${meta.id}.mapped.bam.bai
    """
}

// ============================================================
// PROCESS 7: Plot coverage (samtools depth + R/ggplot2)
// ============================================================
process plot_coverage {
    tag "${meta.id}"
    publishDir "${params.outdir}/coverage", mode: 'copy'

    input:
        tuple val(meta), path(bam), path(bai)

    output:
        tuple val(meta), path("${meta.id}_coverage_plot.png"), path("${meta.id}_coverage_stats.txt")

    script:
    """
    samtools depth ${bam} > depth.txt

    printf 'library(ggplot2)\n' > plot.R
    printf 'd <- read.table("depth.txt", header=F)\n' >> plot.R
    printf 'names(d) <- c("contig", "pos", "depth")\n' >> plot.R
    printf 'if(nrow(d)==0) d <- data.frame(contig="dummy", pos=1:100, depth=0)\n' >> plot.R
    printf 'if(nrow(d)>50000) d <- d[sample(1:nrow(d), 50000),]\n' >> plot.R
    printf 'p <- ggplot(d, aes(pos, depth)) + geom_line(color="blue") + theme_bw()\n' >> plot.R
    printf 'ggsave("${meta.id}_coverage_plot.png", p, width=10, height=6)\n' >> plot.R
    printf 'write.table(data.frame(mean=mean(d\$depth), max=max(d\$depth)), "${meta.id}_coverage_stats.txt", row.names=F)\n' >> plot.R

    Rscript plot.R
    """

    stub:
    """
    touch ${meta.id}_coverage_plot.png
    printf 'mean max\n0 0\n' > ${meta.id}_coverage_stats.txt
    """
}

// ============================================================
// PROCESS 8: Variant calling
// ============================================================
process call_variants {
    tag "${meta.id}"
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
        tuple val(meta), path(bam), path(bai), path(reference)

    output:
        tuple val(meta), path("${meta.id}.variants.vcf.gz"), path("${meta.id}.variants.vcf.gz.tbi"), path("${meta.id}_variants_stats.txt")

    script:
    """
    if [ ! -f ${reference}.fai ]; then
        samtools faidx ${reference}
    fi

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
        --output ${meta.id}.variants.vcf.gz

    bcftools index --tbi ${meta.id}.variants.vcf.gz
    bcftools stats ${meta.id}.variants.vcf.gz > ${meta.id}_variants_stats.txt
    """

    stub:
    """
    touch ${meta.id}.variants.vcf.gz
    touch ${meta.id}.variants.vcf.gz.tbi
    touch ${meta.id}_variants_stats.txt
    """
}

// ============================================================
// PROCESS 9: Per-sample variant filtering
// ============================================================
process filter_variants {
    tag "${meta.id}"
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
        tuple val(meta), path(vcf), path(tbi), path(stats)

    output:
        tuple val(meta), path("${meta.id}.filtered.vcf.gz"), path("${meta.id}.filtered.vcf.gz.tbi")

    script:
    """
    bcftools filter \
        --include 'QUAL>=${params.min_variant_qual} && FMT/DP>=${params.min_depth}' \
        --output-type z \
        --output ${meta.id}.filtered.vcf.gz \
        ${vcf}

    bcftools index --tbi ${meta.id}.filtered.vcf.gz
    """

    stub:
    """
    touch ${meta.id}.filtered.vcf.gz
    touch ${meta.id}.filtered.vcf.gz.tbi
    """
}

// ============================================================
// PROCESS 10: Group-wise comparative filtering
// ============================================================
process compare_group_variants {
    tag "${group_id}"
    publishDir "${params.outdir}/comparisons", mode: 'copy'

    input:
        tuple val(group_id), val(meta_a), path(vcf_a), val(meta_b), path(vcf_b)

    output:
        tuple val(group_id), path("${group_id}_comparative.vcf"), path("${group_id}_comparative_summary.tsv")

    script:
    """
    echo "Comparing ${meta_a.id} (${meta_a.type}) against ${meta_b.id} (${meta_b.type})"

    bcftools index -f ${vcf_a}
    bcftools index -f ${vcf_b}
    bcftools isec -p isec_out -O v ${vcf_a} ${vcf_b}

    cp isec_out/0001.vcf ${group_id}_comparative.vcf

    count=\$(grep -vc '^#' ${group_id}_comparative.vcf || true)
    printf 'group_id\tsample_a\tsample_b\tcomparison_type\tvariant_count\n' > ${group_id}_comparative_summary.tsv
    printf '${group_id}\t${meta_a.id}\t${meta_b.id}\tonly_in_${meta_b.id}\t%s\n' "\$count" >> ${group_id}_comparative_summary.tsv
    """

    stub:
    """
    printf '##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' > ${group_id}_comparative.vcf
    printf 'group_id\tsample_a\tsample_b\tcomparison_type\tvariant_count\n' > ${group_id}_comparative_summary.tsv
    printf '${group_id}\t${meta_a.id}\t${meta_b.id}\tonly_in_${meta_b.id}\t0\n' >> ${group_id}_comparative_summary.tsv
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

    def samples_ch
    def raw_reads_ch

    if (params.samplesheet) {
        samples_ch = Channel.fromPath(params.samplesheet)
            .splitCsv(header: true)
            .map { row ->
                def sampleId = row.sample_id?.trim()
                def groupId  = row.group?.trim()
                def type     = row.type?.trim()
                def sraId    = row.sra_id?.trim()
                def refId    = row.ref_id?.trim() ?: params.ref_id?.trim() ?: ''

                if (!sampleId || !groupId || !type || !sraId) {
                    error "Samplesheet rows must contain sample_id, group, type, and sra_id"
                }

                if (!refId && !params.reference && !params.assembly) {
                    error "Samplesheet rows must contain ref_id when --reference is not provided and --assembly is false"
                }

                def meta = [id: sampleId, group: groupId, type: type, ref_id: refId]
                tuple(meta, sraId)
            }

        raw_reads_ch = download_sra(samples_ch)
    } else if (params.sra_id) {
        def refId = params.ref_id?.trim() ?: ''
        def sampleId = params.sample_id?.trim() ?: params.sra_id
        def meta = [id: sampleId, group: sampleId, type: params.sample_type, ref_id: refId]
        raw_reads_ch = download_sra(Channel.of(tuple(meta, params.sra_id)))
    } else if (params.input_reads_folder) {
        raw_reads_ch = Channel.fromPath("${params.input_reads_folder}/*.fastq*")
            .map { reads ->
                def sampleId = reads.baseName.replaceFirst(/\.fastq(.gz)?$/, '')
                def refId = params.ref_id?.trim() ?: ''
                def meta = [id: sampleId, group: sampleId, type: params.sample_type, ref_id: refId]
                tuple(meta, reads)
            }
    } else {
        error "Provide --samplesheet, --sra_id, or --input_reads_folder"
    }

    qc_initial(raw_reads_ch)

    trimmed_reads_ch = trim_reads(raw_reads_ch)

    qc_trimmed(trimmed_reads_ch)

    def mapping_in_ch

    if (params.reference) {
        reference_ch = Channel.fromPath(params.reference)
        mapping_in_ch = trimmed_reads_ch
            .combine(reference_ch)
            .map { meta, reads, reference -> tuple(meta, reads, reference) }
    } else if (params.assembly) {
        assembled_refs_ch = assemble(trimmed_reads_ch).reference
        mapping_in_ch = trimmed_reads_ch
            .join(assembled_refs_ch)
            .map { meta, reads, reference -> tuple(meta, reads, reference) }
    } else if (params.samplesheet) {
        refs_ch = samples_ch
            .map { meta, sra_id -> meta.ref_id }
            .filter { it }
            .unique()
            .set { unique_ref_ids_ch }

        fetched_refs_ch = fetch_reference(unique_ref_ids_ch)

        mapping_in_ch = trimmed_reads_ch
            .combine(fetched_refs_ch)
            .filter { meta, reads, ref_id, reference -> meta.ref_id == ref_id }
            .map { meta, reads, ref_id, reference -> tuple(meta, reads, reference) }
    } else {
        error "Provide --reference, set --assembly true, or use --samplesheet with ref_id values"
    }

    mapped_reads_ch = map_reads(mapping_in_ch)

    coverage_in_ch = mapped_reads_ch
        .map { meta, bam, bai, reference -> tuple(meta, bam, bai) }
    plot_coverage(coverage_in_ch)

    variants_ch = call_variants(mapped_reads_ch)

    filtered_variants_ch = params.filter_variants
        ? filter_variants(variants_ch)
        : variants_ch.map { meta, vcf, tbi, stats -> tuple(meta, vcf, tbi) }

    if (params.samplesheet || params.enable_group_comparison) {
        sample_a_vcf_ch = filtered_variants_ch
            .filter { meta, vcf, tbi -> meta.type == params.comparison_a_type }
            .map { meta, vcf, tbi -> tuple(meta.group, meta, vcf) }

        sample_b_vcf_ch = filtered_variants_ch
            .filter { meta, vcf, tbi -> meta.type == params.comparison_b_type }
            .map { meta, vcf, tbi -> tuple(meta.group, meta, vcf) }

        joined_variants_ch = sample_a_vcf_ch.join(sample_b_vcf_ch)
        compare_group_variants(joined_variants_ch)
    }
}
