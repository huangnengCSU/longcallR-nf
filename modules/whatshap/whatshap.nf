process WHATSHAP {
    tag "Whatshap for ${params.sample_name}"
    publishDir "${params.outdir}/whatshap", mode: 'symlink'

    input:
    path vcf
    path vcf_index
    path bam
    path bam_index
    path ref
    path ref_fai
    path annotation

    output:
    path "whatshap_output/${params.sample_name}/${params.sample_name}.*", emit: whatshap_outputs_ch
}