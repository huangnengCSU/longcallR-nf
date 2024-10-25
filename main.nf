nextflow.enable.dsl=2


// project parameters
params.sample_name = null
params.reads = null
params.bam = null
params.ref = null
params.contigs = null
params.platform = 'ont'
params.datatype = 'cDNA'
params.outdir = null
params.threads = 4
params.threads_per_job = 4 < params.threads ? 4 : params.threads
params.num_jobs = (params.threads/params.threads_per_job) as int
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


// process parameters
def contigs = params.contigs ? 
              params.contigs.split(",") as List : 
              (1..22).collect { "chr${it}" } // Default to ['chr1', ..., 'chr22']
def reads = params.reads ? 
            (params.reads instanceof String ? params.reads.split(",") as List : params.reads) : 
            null


// Include the modules
include { MINIMAP2_ALIGN } from './modules/minimap2/minimap2.nf'
include { INSTALL_LONGCALLR_DP } from './modules/longcallR-dp/longcallR-dp.nf'
include { LONGCALLR_DP } from './modules/longcallR-dp/longcallR-dp.nf'
include { LONGCALLR_NN_CALL } from './modules/longcallR-nn/longcallR-nn.nf'
include { BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR_NN } from './modules/bcftools/bcftools.nf'
include { BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR } from './modules/bcftools/bcftools.nf'
include { BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR_EXON } from './modules/bcftools/bcftools.nf'
include { INSTALL_LONGCALLR } from './modules/longcallR/longcallR.nf'
include { LONGCALLR_CALL_PHASE } from './modules/longcallR/longcallR.nf'
include { LONGCALLR_CALL_PHASE_EXON } from './modules/longcallR/longcallR.nf'
include { SAMTOOLS_MERGE_SORT_INDEX } from './modules/samtools/samtools.nf'
include { SAMTOOLS_MERGE_SORT_INDEX_EXON } from './modules/samtools/samtools.nf'
include { SAMTOOLS_INDEX } from './modules/samtools/samtools.nf'
include { SAMTOOLS_FAIDX } from './modules/samtools/samtools.nf'
include { ISOQUANT } from './modules/isoquant/isoquant.nf'

// Define the workflow
workflow {

    // Check if the required parameters are provided
    if (!params.reads && !params.bam) {
        error "You must provide either --reads or --bam parameter."
    }

    if (params.reads != null && params.bam == null) {
        MINIMAP2_ALIGN(params.ref, reads)
        ch_align_bam = MINIMAP2_ALIGN.out.ch_align_bam
        ch_align_bam_bai = MINIMAP2_ALIGN.out.ch_align_bam_bai
        ch_ref = Channel.fromPath(params.ref)
        ch_ref_fai = MINIMAP2_ALIGN.out.ch_ref_fai
    } else if (params.reads == null && params.bam != null) {
        ch_align_bam = Channel.fromPath(params.bam)
        SAMTOOLS_INDEX(params.bam)
        ch_align_bam_bai = SAMTOOLS_INDEX.out.bam_index
        ch_ref = Channel.fromPath(params.ref)
        SAMTOOLS_FAIDX(params.ref)
        ch_ref_fai = SAMTOOLS_FAIDX.out.fasta_index
    } else {
        error "You must provide either --reads or --bam parameter."
    }

    ch_contigs = Channel.from(contigs)

    // longcallR_dp
    INSTALL_LONGCALLR_DP()
    ch_longcallr_dp_binary = INSTALL_LONGCALLR_DP.out.longcallr_dp_binary
    LONGCALLR_DP(ch_longcallr_dp_binary, ch_contigs, ch_align_bam, ch_align_bam_bai, ch_ref, ch_ref_fai)
    ch_longcallR_dp_features_contigs = LONGCALLR_DP.out.longcallR_dp_features_contigs_ch
    ch_longcallR_dp_features_contigs.collect().set { collected_longcallR_dp_features_contigs }  // Wait for all chromosomes to finish


    // longcallR_nn
    LONGCALLR_NN_CALL(collected_longcallR_dp_features_contigs, ch_ref, ch_ref_fai)
    ch_longcallR_nn_vcfs = LONGCALLR_NN_CALL.out.longcallR_nn_vcfs_ch
    ch_longcallR_nn_vcfs.collect().set { collected_longcallR_nn_vcfs }  // Wait for all chromosomes to finish
    BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR_NN(collected_longcallR_nn_vcfs, "${params.sample_name}_longcallR_nn")
    ch_longcallR_nn_vcf = BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR_NN.out.vcf_file
    ch_longcallR_nn_vcf_index = BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR_NN.out.vcf_index

    // longcallR
    INSTALL_LONGCALLR()
    longcallr_binary_ch = INSTALL_LONGCALLR.out.longcallr_binary    
    LONGCALLR_CALL_PHASE(longcallr_binary_ch, ch_align_bam, ch_align_bam_bai, ch_ref, ch_ref_fai, ch_longcallR_nn_vcf, ch_longcallR_nn_vcf_index, ch_contigs)
    ch_longcallR_vcfs = LONGCALLR_CALL_PHASE.out.longcallR_vcfs_ch
    ch_longcallR_phased_bams = LONGCALLR_CALL_PHASE.out.longcallR_phased_bams_ch
    ch_longcallR_vcfs.collect().set { collected_longcallR_vcfs }  // Wait for all chromosomes to finish
    ch_longcallR_phased_bams.collect().set { collected_longcallR_phased_bams }  // Wait for all chromosomes to finish
    BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR(collected_longcallR_vcfs, "${params.sample_name}_longcallR")
    ch_longcallR_vcf = BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR.out.vcf_file
    ch_longcallR_vcf_index = BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR.out.vcf_index
    SAMTOOLS_MERGE_SORT_INDEX(collected_longcallR_phased_bams, "${params.sample_name}_longcallR")
    ch_longcallR_bam = SAMTOOLS_MERGE_SORT_INDEX.out.bam_file
    ch_longcallR_bam_index = SAMTOOLS_MERGE_SORT_INDEX.out.bam_index

    // longcallR exon
    LONGCALLR_CALL_PHASE_EXON(longcallr_binary_ch, ch_align_bam, ch_align_bam_bai, ch_ref, ch_ref_fai, ch_longcallR_nn_vcf, ch_longcallR_nn_vcf_index, params.annotation, ch_contigs)
    ch_longcallR_exon_vcfs = LONGCALLR_CALL_PHASE_EXON.out.longcallR_exon_vcfs_ch
    ch_longcallR_exon_phased_bams = LONGCALLR_CALL_PHASE_EXON.out.longcallR_exon_phased_bams_ch
    ch_longcallR_exon_vcfs.collect().set { collected_longcallR_exon_vcfs }  // Wait for all chromosomes to finish
    ch_longcallR_exon_phased_bams.collect().set { collected_longcallR_exon_phased_bams }  // Wait for all chromosomes to finish
    BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR_EXON(collected_longcallR_exon_vcfs, "${params.sample_name}_longcallR_exon")
    ch_longcallR_exon_vcf = BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR_EXON.out.vcf_file
    ch_longcallR_exon_vcf_index = BCFTOOLS_CONCAT_SORT_VCF_LONGCALLR_EXON.out.vcf_index
    SAMTOOLS_MERGE_SORT_INDEX_EXON(collected_longcallR_exon_phased_bams, "${params.sample_name}_longcallR_exon")
    ch_longcallR_exon_bam = SAMTOOLS_MERGE_SORT_INDEX_EXON.out.bam_file
    ch_longcallR_exon_bam_index = SAMTOOLS_MERGE_SORT_INDEX_EXON.out.bam_index

    // isoquant
    ISOQUANT(ch_longcallR_bam, ch_longcallR_bam_index, ch_ref, ch_ref_fai, params.annotation)
    ch_isoquant_outputs = ISOQUANT.out.isoquant_outputs_ch
}
