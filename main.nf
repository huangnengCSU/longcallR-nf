nextflow.enable.dsl=2


// project parameters
params.sample_name = null
params.reads = null
params.ref = null
params.contigs = (1..22).collect { "chr${it}" } // Creates ['chr1', 'chr2', ..., 'chr22']
params.platform = 'ont'
params.datatype = 'cDNA'
params.outdir = null
params.threads = 4
params.memory = '200GB'
params.gpu_device = 0

// LongcallR_dp parameters
params.min_depth = 6
params.min_af = 0.1
params.min_bq = 0


// LongcallR_nn parameters
params.max_depth = 200
params.batch_size = 256
params.no_cuda = false

// LongcallR parameters

// Isoquant parameters
params.annotation = null


// Include the modules
include { MINIMAP2_ALIGN } from './modules/minimap2/main.nf'
include { LONGCALLR_NN_CALL } from './modules/longcallR-nn/main.nf'
include { BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR_NN } from './modules/bcftools/main.nf'
include { BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR } from './modules/bcftools/main.nf'
include { INSTALL_LONGCALLR } from './modules/longcallR/main.nf'
include { LONGCALLR_CALL_PHASE } from './modules/longcallR/main.nf'
include { SAMTOOLS_MERGE_SORT_INDEX } from './modules/samtools/main.nf'
include { ISOQUANT } from './modules/isoquant/main.nf'

// Define the workflow
workflow {

    // Check if the required parameters are provided
    if (!params.reads || !params.ref) {
        error "You must provide both --reads and --ref parameters."
    }

    // Run the Minimap2 alignment process
    MINIMAP2_ALIGN(params.ref, params.reads)

    // Uncomment and add input parameters to run the LONGCALLR_NN_CALL process
    contigs_ch = Channel.from(params.contigs)
    ch_align_bam = MINIMAP2_ALIGN.out.ch_align_bam
    ch_align_bam_bai = MINIMAP2_ALIGN.out.ch_align_bam_bai
    ch_ref = Channel.fromPath(params.ref)
    ch_ref_fai = MINIMAP2_ALIGN.out.ch_ref_fai

    // Run the LONGCALLR_NN_CALL process, collect the output in a channel
    LONGCALLR_NN_CALL(contigs_ch, ch_align_bam, ch_align_bam_bai, ch_ref, ch_ref_fai)
    longcallR_nn_vcfs_ch = LONGCALLR_NN_CALL.out.longcallR_nn_vcfs_ch
    longcallR_nn_vcfs_ch.collect()  // Wait for all chromosomes to finish
    BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR_NN(longcallR_nn_vcfs_ch, "${params.sample_name}_longcallR_nn")
    ch_longcallR_nn_vcf = BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR_NN.out.vcf_file
    ch_longcallR_nn_vcf_index = BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR_NN.out.vcf_index

    // Install longcallR and run the LONGCALLR_CALL_PHASE process
    INSTALL_LONGCALLR()
    longcallr_binary_ch = INSTALL_LONGCALLR.out.longcallr_binary
    
    LONGCALLR_CALL_PHASE(longcallr_binary_ch, ch_align_bam, ch_align_bam_bai, ch_ref, ch_ref_fai, ch_longcallR_nn_vcf, ch_longcallR_nn_vcf_index, contigs_ch)

    longcallR_vcfs_ch = LONGCALLR_CALL_PHASE.out.longcallR_vcfs_ch
    longcallR_phased_bams_ch = LONGCALLR_CALL_PHASE.out.longcallR_phased_bams_ch

    longcallR_vcfs_ch.collect()  // Wait for all chromosomes to finish
    longcallR_phased_bams_ch.collect()  // Wait for all chromosomes to finish

    BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR(longcallR_vcfs_ch, "${params.sample_name}_longcallR")
    ch_longcallR_vcf = BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR.out.vcf_file
    ch_longcallR_vcf_index = BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR.out.vcf_index

    SAMTOOLS_MERGE_SORT_INDEX(longcallR_phased_bams_ch, "${params.sample_name}_longcallR")
    ch_longcallR_bam = SAMTOOLS_MERGE_SORT_INDEX.out.bam_file
    ch_longcallR_bam_index = SAMTOOLS_MERGE_SORT_INDEX.out.bam_index

    ISOQUANT(ch_longcallR_bam, ch_longcallR_bam_index, ch_ref, ch_ref_fai, params.annotation)
    ch_isoquant_outputs = ISOQUANT.out.isoquant_outputs_ch
}
