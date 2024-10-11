process SAMTOOLS_MERGE_SORT_INDEX {
    tag "samtools merge, sort and index"
    conda 'bioconda::samtools'
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path bam_files  // channel of BAM files
    val prefix

    output:
    path "${prefix}.phased.sorted.bam", emit: bam_file
    path "${prefix}.phased.sorted.bam.bai", emit: bam_index

    script:
    """
    samtools merge -@ ${params.threads} ${prefix}.phased.bam ${bam_files.join(' ')}
    samtools sort -@ ${params.threads} ${prefix}.phased.bam -o ${prefix}.phased.sorted.bam
    samtools index -@ ${params.threads} ${prefix}.phased.sorted.bam
    """
}