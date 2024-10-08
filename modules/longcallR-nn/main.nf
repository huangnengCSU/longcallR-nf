process LONGCALLR_NN_CALL {
    publishDir "${params.outdir}/longcallR_nn_call", mode: 'symlink'
    apptainer 'apbiobio/longcallr-nn:0.1.0'

    inputs:
    path bam from minimap2_align_bam
    path ref from params.ref
    val read_basename = params.reads.baseName
    val contigs from params.contigs
    val min_depth from params.min_depth
    val min_af from params.min_af
    val min_bq from params.min_bq
    val threads from params.threads
    val platform from params.platform
    val type from params.type

    outputs:
    path("${read_basename}.longcallR_nn.sort.vcf") into longcallR_nn_vcf
    path("${read_basename}.longcallR_nn.sort.vcf.gz") into longcallR_nn_vcf_gz
    path("${read_basename}.longcallR_nn.sort.vcf.gz.tbi") into longcallR_nn_vcf_gz_tbi

    script:
    """
    # Create output directory based on fasta basename
    FEATURE_DIR="${read_basename}_features/contig"
    mkdir -p "${read_basename}_features"
    PREDICTION_DIR="${read_basename}_predictions"
    mkdir -p ${PREDICTION_DIR}



    # Calculate number of parallel jobs that won't exceed the total threads
    threads_per_job=2
    max_parallel_jobs=$(( ${threads} / ${threads_per_job} ))

    # Run predictions for each contig in parallel
    parallel -j ${max_parallel_jobs} \
    "longcallR_nn predict \
    --bam-path '${bam}' \
    --ref-path '${ref}' \
    --min-depth '${min_depth}' \
    --min-alt-freq '${min_af}' \
    --min-baseq '${min_bq}' \
    --threads '${threads_per_job}' \
    --contig {1} \
    --output '${FEATURE_DIR}_{1}'" ::: ${contigs[@]}

    # Run longcallR_nn call for each contig in serial
    # Validate platform and type parameters
    case "${platform}" in
        ont)
            case "${type}" in
                cDNA)
                    config="config/wtc11_cdna.yaml"
                    model="models/cdna_wtc11_nopass_resnet50_sgd.epoch30.chkpt"
                    ;;
                dRNA)
                    config="config/gm12878_drna.yaml"
                    model="models/drna_gm12878_nopass_resnet50_sgd.epoch30.chkpt"
                    ;;
                *)
                    echo "Error: For 'ont' platform, type must be 'cDNA' or 'dRNA'. Given type: ${type}"
                    exit 1
                    ;;
            esac
            ;;
        pb)
            case "${type}" in
                isoseq)
                    config="config/hg002_isoseq.yaml"
                    model="models/hg002_baylor_isoseq_nopass_resnet50_sgd.epoch30.chkpt"
                    ;;
                masseq)
                    config="config/hg002_na24385_masseq.yaml"
                    model="models/hg002_na24385_mix_nopass_resnet50_sgd.epoch30.chkpt"
                    ;;
                *)
                    echo "Error: For 'pb' platform, type must be 'isoseq' or 'masseq'. Given type: ${type}"
                    exit 1
                    ;;
            esac
            ;;
        *)
            echo "Error: Platform must be 'ont' or 'pb'. Given platform: ${platform}"
            exit 1
            ;;
    esac

    for ctg in ${contigs[@]}; do
        longcallR_nn call \
        -config ${config} \
        -model ${model} \
        -data '${FEATURE_DIR}_${ctg}' \
        -ref ${ref} \
        -output '${PREDICTION_DIR}/${read_basename}_longcallR_nn_${ctg}.vcf' \
        -max_depth 200 \
        -batch_size 500
    done

    # Concatenate all VCF files
    bcftools concat ${PREDICTION_DIR}/${read_basename}_longcallR_nn_*.vcf -o ${read_basename}.longcallR_nn.vcf
    bcftools sort ${read_basename}.longcallR_nn.vcf -o ${read_basename}.longcallR_nn.sort.vcf
    bgzip ${read_basename}.longcallR_nn.sort.vcf > ${read_basename}.longcallR_nn.sort.vcf.gz
    tabix -p vcf ${read_basename}.longcallR_nn.sort.vcf.gz
    """
}
