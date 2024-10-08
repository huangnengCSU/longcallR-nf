nextflow.enable.dsl=2

params.reads = null
params.ref = null
params.threads = 1
// platform choices: ont, pb
params.platform = 'ont'
// type choices: cDNA, dRNA, isoseq, masseq
params.type = 'cDNA'
params.outdir = null

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
    MINIMAP2_ALIGN(params.ref, params.reads, params.threads, params.platform, params.type)

    // Uncomment and add input parameters to run the LONGCALLR_NN_CALL process
    // LONGCALLR_NN_CALL()
}
