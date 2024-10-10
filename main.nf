nextflow.enable.dsl=2
params.sample_name = null
params.reads = null
params.ref = null
params.threads = 1
// platform choices: ont, pb
params.platform = 'ont'
// datatype choices: cDNA, dRNA, isoseq, masseq
params.datatype = 'cDNA'
params.outdir = null
params.min_depth = 6
params.min_af = 0.1
params.min_bq = 0
params.contigs = (1..22).collect { "chr${it}" } // Creates ['chr1', 'chr2', ..., 'chr22']

// Include the modules
include { MINIMAP2_ALIGN } from './modules/minimap2/main.nf'
include { LONGCALLR_NN_CALL } from './modules/longcallR-nn/main.nf'
include { CONCAT_SORT_VCF } from './modules/bcftools/main.nf'

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
    CONCAT_SORT_VCF(longcallR_nn_vcfs_ch, "${params.sample_name}_longcallR_nn")
}
