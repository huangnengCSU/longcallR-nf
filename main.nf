nextflow.enable.dsl=2

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
// include { LONGCALLR_NN_CALL } from './modules/longcallR-nn/main.nf'

// Define the workflow
workflow {

    // Check if the required parameters are provided
    if (!params.reads || !params.ref) {
        error "You must provide both --reads and --ref parameters."
    }

    // Run the Minimap2 alignment process
    MINIMAP2_ALIGN(params.ref, params.reads)

    // Uncomment and add input parameters to run the LONGCALLR_NN_CALL process
    // contigs_ch = Channel.value(params.contigs)
    // LONGCALLR_NN_CALL(MINIMAP2_ALIGN.ch_align_bam, params.ref, contigs_ch)
}
