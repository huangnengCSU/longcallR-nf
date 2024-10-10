process BCFTOOLS_CONCAT_SORT_VCF {
    tag "Concatenate and sort VCF files"
    conda 'bioconda::bcftools bioconda::tabix'
    publishDir "${params.outdir}", mode: 'symlink'

    input:
    path vcf_files  // channel of VCF files
    val prefix  // Prefix for output VCF file

    output:
    path "${prefix}.sorted.vcf.gz"
    path "${prefix}.sorted.vcf.gz.tbi"

    script:
    """
    bcftools concat ${vcf_files.join(' ')} | bcftools sort -o ${prefix}.sorted.vcf
    bgzip -c ${prefix}.sorted.vcf > ${prefix}.sorted.vcf.gz
    tabix -p vcf ${prefix}.sorted.vcf.gz
    """
}