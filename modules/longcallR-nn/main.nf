process LONGCALLR_NN_CALL {
    publishDir "${params.outdir}/longcallR_nn_call", mode: 'symlink'

    input:
        val contig  // Channel of contigs
        path bam   // Channel of BAM file
        path bam_bai    // Channel of BAM index file
        path ref    // Channel of reference genome
        path ref_fai    // Channel of reference genome index file


    output:
    stdout
    path "${params.sample_name}_features"
    path "${params.sample_name}_predictions"
    path "${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_*.vcf", emit: longcallR_nn_vcfs_ch

    script:
    """
    mkdir -p ${params.sample_name}_features
    mkdir -p ${params.sample_name}_predictions

    longcallR_dp \
    --mode predict \
    --bam-path '${bam}' \
    --ref-path '${ref}' \
    --min-depth '${params.min_depth}' \
    --min-alt-freq '${params.min_af}' \
    --min-baseq '${params.min_bq}' \
    --threads '${params.threads}' \
    --contigs '${contig}' \
    --output '${params.sample_name}_features/${contig}'


    source /tools/miniconda3/bin/activate pytorch
    if [ ${params.platform} == "ont" ]; then
        if [ ${params.datatype} == "cDNA" ]; then
            longcallR_nn call \
            -config /tools/longcallR-nn/config/wtc11_cdna.yaml \
            -model /tools/longcallR-nn/models/cdna_wtc11_nopass_resnet50_sgd.epoch30.chkpt \
            -data ${params.sample_name}_features/${contig} \
            -ref ${ref} \
            -output ${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_${contig}.vcf \
            -max_depth 200 \
            -batch_size 500 \
            --no_cuda
        elif [ ${params.datatype} == "dRNA" ]; then
            longcallR_nn call \
            -config /tools/longcallR-nn/config/gm12878_drna.yaml \
            -model /tools/longcallR-nn/models/drna_gm12878_nopass_resnet50_sgd.epoch30.chkpt \
            -data ${params.sample_name}_features/${contig} \
            -ref ${ref} \
            -output ${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_${contig}.vcf \
            -max_depth 200 \
            -batch_size 500 \
            --no_cuda
        else
            echo "Error: For 'ont' platform, datatype must be 'cDNA' or 'dRNA'. Given datatype: ${params.datatype}"
            exit 1
        fi
    elif [ ${params.platform} == "pb" ]; then
        if [ ${params.datatype} == "isoseq" ]; then
            longcallR_nn call \
            -config /tools/longcallR-nn/config/hg002_isoseq.yaml \
            -model /tools/longcallR-nn/models/hg002_baylor_isoseq_nopass_resnet50_sgd.epoch30.chkpt \
            -data ${params.sample_name}_features/${contig} \
            -ref ${ref} \
            -output ${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_${contig}.vcf \
            -max_depth 200 \
            -batch_size 500 \
            --no_cuda
        elif [ ${params.datatype} == "masseq" ]; then
            longcallR_nn call \
            -config /tools/longcallR-nn/config/hg002_na24385_masseq.yaml \
            -model /tools/longcallR-nn/models/hg002_na24385_mix_nopass_resnet50_sgd.epoch30.chkpt \
            -data ${params.sample_name}_features/${contig} \
            -ref ${ref} \
            -output ${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_${contig}.vcf \
            -max_depth 200 \
            -batch_size 500 \
            --no_cuda
        else
            echo "Error: For 'pb' platform, datatype must be 'isoseq' or 'masseq'. Given datatype: ${params.datatype}"
            exit 1
        fi
    else
        echo "Error: Platform must be 'ont' or 'pb'. Given platform: ${params.platform}"
        exit 1
    fi
    """
}
