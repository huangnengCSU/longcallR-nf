process BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR_NN {
    tag "Concatenate and sort longcallR-nn VCF files"
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path vcf_files  // channel of VCF files
    val prefix  // Prefix for output VCF file

    output:
    path "${prefix}.sorted.vcf.gz", emit: vcf_file
    path "${prefix}.sorted.vcf.gz.tbi", emit: vcf_index

    script:
    """
    bcftools concat ${vcf_files.join(' ')} | bcftools sort -o ${prefix}.sorted.vcf
    bgzip -c ${prefix}.sorted.vcf > ${prefix}.sorted.vcf.gz
    tabix -p vcf ${prefix}.sorted.vcf.gz
    """
}


process BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR {
    tag "Concatenate and sort longcallR VCF files"
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path vcf_files  // channel of VCF files
    val prefix  // Prefix for output VCF file

    output:
    path "${prefix}.sorted.vcf.gz", emit: vcf_file
    path "${prefix}.sorted.vcf.gz.tbi", emit: vcf_index

    script:
    """
    bcftools concat ${vcf_files.join(' ')} | bcftools sort -o ${prefix}.sorted.vcf
    bgzip -c ${prefix}.sorted.vcf > ${prefix}.sorted.vcf.gz
    tabix -p vcf ${prefix}.sorted.vcf.gz
    """
}

process BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR_EXON {
    tag "Concatenate and sort longcallR exon VCF files"
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path vcf_files  // channel of VCF files
    val prefix  // Prefix for output VCF file

    output:
    path "${prefix}.sorted.vcf.gz", emit: vcf_file
    path "${prefix}.sorted.vcf.gz.tbi", emit: vcf_index

    script:
    """
    bcftools concat ${vcf_files.join(' ')} | bcftools sort -o ${prefix}.sorted.vcf
    bgzip -c ${prefix}.sorted.vcf > ${prefix}.sorted.vcf.gz
    tabix -p vcf ${prefix}.sorted.vcf.gz
    """
}

process BCFTOOLS_CONCAT_SORT_VCF_ONLY_LONGCALLR {
    tag "Concatenate and sort only longcallR VCF files"
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path vcf_files  // channel of VCF files
    val prefix  // Prefix for output VCF file

    output:
    path "${prefix}.sorted.vcf.gz", emit: vcf_file
    path "${prefix}.sorted.vcf.gz.tbi", emit: vcf_index

    script:
    """
    bcftools concat ${vcf_files.join(' ')} | bcftools sort -o ${prefix}.sorted.vcf
    bgzip -c ${prefix}.sorted.vcf > ${prefix}.sorted.vcf.gz
    tabix -p vcf ${prefix}.sorted.vcf.gz
    """
}

