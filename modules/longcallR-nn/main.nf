process LONGCALLR_NN_CALL {
    publishDir "${params.outdir}/longcallR_nn_call", mode: 'symlink'
    conda 'bioconda::bcftools bioconda::tabix'

    inputs:
    path bam
    path ref
    val contigs_list

    outputs:
    path "${basename}.longcallR_nn.sort.vcf", emit: ch_longcallR_nn_vcf
    path "${basename}.longcallR_nn.sort.vcf.gz", emit: ch_longcallR_nn_vcf_gz
    path "${basename}.longcallR_nn.sort.vcf.gz.tbi", emit: ch_longcallR_nn_vcf_gz_tbi

    script:
    """
    basename="${bam.baseName}"
    FEATURE_DIR="${basename}_features"
    mkdir -p ${FEATURE_DIR}
    PREDICTION_DIR="${basename}_predictions"
    mkdir -p ${PREDICTION_DIR}

    threads_per_job=2
    max_parallel_jobs=$(( ${params.threads} / ${threads_per_job} ))

    # Run predictions for each contig in parallel
    parallel -j ${max_parallel_jobs} \
    "longcallR_dp \
    --mode predict \
    --bam-path '${bam}' \
    --ref-path '${ref}' \
    --min-depth '${params.min_depth}' \
    --min-alt-freq '${params.min_af}' \
    --min-baseq '${params.min_bq}' \
    --threads '${threads_per_job}' \
    --contigs {1} \
    --output '${FEATURE_DIR}/{1}'" ::: ${contigs_list[@]}

    // Validate params.platform and params.datatype parameters
    if (params.platform == 'ont') {
        if (params.datatype == 'cDNA') {
            config = 'config/wtc11_cdna.yaml'
            model = 'models/cdna_wtc11_nopass_resnet50_sgd.epoch30.chkpt'
        } else if (params.datatype == 'dRNA') {
            config = 'config/gm12878_drna.yaml'
            model = 'models/drna_gm12878_nopass_resnet50_sgd.epoch30.chkpt'
        } else {
            error "Error: For 'ont' platform, params.datatype must be 'cDNA' or 'dRNA'. Given params.datatype: ${params.datatype}"
        }
    } else if (params.platform == 'pb') {
        if (params.datatype == 'isoseq') {
            config = 'config/hg002_isoseq.yaml'
            model = 'models/hg002_baylor_isoseq_nopass_resnet50_sgd.epoch30.chkpt'
        } else if (params.datatype == 'masseq') {
            config = 'config/hg002_na24385_masseq.yaml'
            model = 'models/hg002_na24385_mix_nopass_resnet50_sgd.epoch30.chkpt'
        } else {
            error "Error: For 'pb' platform, params.datatype must be 'isoseq' or 'masseq'. Given params.datatype: ${params.datatype}"
        }
    } else {
        error "Error: Platform must be 'ont' or 'pb'. Given params.platform: ${params.platform}"
    }

    for ctg in ${contigs_list[@]}; do
        longcallR_nn call \
        -config ${config} \
        -model ${model} \
        -data ${FEATURE_DIR}/${ctg} \
        -ref ${ref} \
        -output ${PREDICTION_DIR}/${basename}_longcallR_nn_${ctg}.vcf \
        -max_depth 200 \
        -batch_size 500
    done

    bcftools concat ${PREDICTION_DIR}/${basename}_longcallR_nn_*.vcf -o ${basename}.longcallR_nn.vcf
    bcftools sort ${basename}.longcallR_nn.vcf -o ${basename}.longcallR_nn.sort.vcf
    bgzip ${basename}.longcallR_nn.sort.vcf > ${basename}.longcallR_nn.sort.vcf.gz
    tabix -p vcf ${basename}.longcallR_nn.sort.vcf.gz
    """
}
