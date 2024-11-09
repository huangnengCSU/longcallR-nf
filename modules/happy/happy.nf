// process HAPPY {
//     tag "Evaluate haplotype phasing"
//     conda 'python==2.7.* bioconda::hap.py'
//     publishDir "${params.outdir}/happy", mode: 'symlink'

//     input:
//     path query_vcf_file
//     path truth_vcf_file
//     path ref_file
//     path ref_fai
//     path bed_file
//     val prefix  // report prefix

//     output:
//     path "${prefix}.*", emit: report_files

//     script:
//     """
//     hap.py ${truth_vcf_file} ${query_vcf_file} -r ${ref_file} -f ${bed_file} --pass-only --engine=vcfeval --threads=${params.threads} -o ${prefix}
//     """
// }