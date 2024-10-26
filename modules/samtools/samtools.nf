process SAMTOOLS_MERGE_SORT_INDEX {
    tag "samtools merge, sort and index"
    conda 'bioconda::samtools==1.17'
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

process SAMTOOLS_MERGE_SORT_INDEX_EXON {
    tag "samtools merge, sort and index"
    conda 'bioconda::samtools==1.17'
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

process SAMTOOLS_MERGE_SORT_INDEX_ONLY {
    tag "samtools merge, sort and index"
    conda 'bioconda::samtools==1.17'
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

process SAMTOOLS_INDEX {
    tag "samtools index"
    conda 'bioconda::samtools==1.17'

    input:
    path bam_file
    val prev_ch

    output:
    path "${bam_file}.bai", emit: bam_index

    script:
    """
    samtools index -@ ${params.threads} ${bam_file}
    """
}

process SAMTOOLS_FAIDX {
    tag "samtools faidx"
    conda 'bioconda::samtools==1.17'

    input:
    path fasta_file
    val prev_ch

    output:
    path "${fasta_file}.fai", emit: fasta_index

    script:
    """
    samtools faidx ${fasta_file}
    """
}
