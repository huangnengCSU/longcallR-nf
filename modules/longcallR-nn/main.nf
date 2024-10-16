process INSTALL_LONGCALLR_DP {
    tag "Install longcallR_dp"

    output:
    path "longcallR-nn/longcallR_dp/target/release/longcallR_dp", emit: longcallr_dp_binary

    script:
    """
    git clone https://github.com/huangnengCSU/longcallR-nn.git
    cd longcallR-nn/longcallR_dp
    cargo build --release
    """
}

process LONGCALLR_NN_CALL {
    tag "${params.sample_name}_longcallR_nn_call"
    conda 'python=3.9 bioconda::longcallr_nn==0.0.1'
    publishDir "${params.outdir}/longcallR_nn_call", mode: 'symlink'

    input:
        path longcallr_dp_binary
        val contig  // Channel of contigs
        path bam   // Channel of BAM file
        path bam_bai    // Channel of BAM index file
        path ref    // Channel of reference genome
        path ref_fai    // Channel of reference genome index file


    output:
    path "models/*", emit: models_ch
    path "${params.sample_name}_features"
    path "${params.sample_name}_predictions"
    path "${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_*.vcf", emit: longcallR_nn_vcfs_ch

    script:
    """
    mkdir -p ${params.sample_name}_features
    mkdir -p ${params.sample_name}_predictions
    mkdir -p models

    ./${longcallr_dp_binary} \
    --mode predict \
    --bam-path '${bam}' \
    --ref-path '${ref}' \
    --min-depth '${params.min_depth}' \
    --min-alt-freq '${params.min_af}' \
    --min-baseq '${params.min_bq}' \
    --threads '${params.threads}' \
    --contigs '${contig}' \
    --output '${params.sample_name}_features/${contig}'

    longcallR_nn download -d models

    if [ ${params.no_cuda} == "true" ]; then
        if [ ${params.platform} == "ont" ]; then
            if [ ${params.datatype} == "cDNA" ]; then
                longcallR_nn call \
                -config models/ont_cdna_config.yaml \
                -model models/ont_cdna_model.chkpt \
                -data ${params.sample_name}_features/${contig} \
                -ref ${ref} \
                -output ${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_${contig}.vcf \
                -max_depth ${params.max_depth} \
                -batch_size ${params.batch_size} \
                --no_cuda
            elif [ ${params.datatype} == "dRNA" ]; then
                longcallR_nn call \
                -config models/ont_drna_config.yaml \
                -model models/ont_drna_model.chkpt \
                -data ${params.sample_name}_features/${contig} \
                -ref ${ref} \
                -output ${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_${contig}.vcf \
                -max_depth ${params.max_depth} \
                -batch_size ${params.batch_size} \
                --no_cuda
            else
                echo "Error: For 'ont' platform, datatype must be 'cDNA' or 'dRNA'. Given datatype: ${params.datatype}"
                exit 1
            fi
        elif [ ${params.platform} == "pb" ]; then
            if [ ${params.datatype} == "isoseq" ]; then
                longcallR_nn call \
                -config models/pb_isoseq_config.yaml \
                -model models/pb_isoseq_model.chkpt \
                -data ${params.sample_name}_features/${contig} \
                -ref ${ref} \
                -output ${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_${contig}.vcf \
                -max_depth ${params.max_depth} \
                -batch_size ${params.batch_size} \
                --no_cuda
            elif [ ${params.datatype} == "masseq" ]; then
                longcallR_nn call \
                -config models/pb_masseq_config.yaml \
                -model models/pb_masseq_model.chkpt \
                -data ${params.sample_name}_features/${contig} \
                -ref ${ref} \
                -output ${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_${contig}.vcf \
                -max_depth ${params.max_depth} \
                -batch_size ${params.batch_size} \
                --no_cuda
            else
                echo "Error: For 'pb' platform, datatype must be 'isoseq' or 'masseq'. Given datatype: ${params.datatype}"
                exit 1
            fi
        else
            echo "Error: Platform must be 'ont' or 'pb'. Given platform: ${params.platform}"
            exit 1
        fi
    else
        export CUDA_VISIBLE_DEVICES=${params.gpu_device}
        if [ ${params.platform} == "ont" ]; then
            if [ ${params.datatype} == "cDNA" ]; then
                longcallR_nn call \
                -config models/ont_cdna_config.yaml \
                -model models/ont_cdna_model.chkpt \
                -data ${params.sample_name}_features/${contig} \
                -ref ${ref} \
                -output ${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_${contig}.vcf \
                -max_depth ${params.max_depth} \
                -batch_size ${params.batch_size}
            elif [ ${params.datatype} == "dRNA" ]; then
                longcallR_nn call \
                -config models/ont_drna_config.yaml \
                -model models/ont_drna_model.chkpt \
                -data ${params.sample_name}_features/${contig} \
                -ref ${ref} \
                -output ${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_${contig}.vcf \
                -max_depth ${params.max_depth} \
                -batch_size ${params.batch_size}
            else
                echo "Error: For 'ont' platform, datatype must be 'cDNA' or 'dRNA'. Given datatype: ${params.datatype}"
                exit 1
            fi
        elif [ ${params.platform} == "pb" ]; then
            if [ ${params.datatype} == "isoseq" ]; then
                longcallR_nn call \
                -config models/pb_isoseq_config.yaml \
                -model models/pb_isoseq_model.chkpt \
                -data ${params.sample_name}_features/${contig} \
                -ref ${ref} \
                -output ${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_${contig}.vcf \
                -max_depth ${params.max_depth} \
                -batch_size ${params.batch_size}
            elif [ ${params.datatype} == "masseq" ]; then
                longcallR_nn call \
                -config models/pb_masseq_config.yaml \
                -model models/pb_masseq_model.chkpt \
                -data ${params.sample_name}_features/${contig} \
                -ref ${ref} \
                -output ${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_${contig}.vcf \
                -max_depth ${params.max_depth} \
                -batch_size ${params.batch_size}
            else
                echo "Error: For 'pb' platform, datatype must be 'isoseq' or 'masseq'. Given datatype: ${params.datatype}"
                exit 1
            fi
        else
            echo "Error: Platform must be 'ont' or 'pb'. Given platform: ${params.platform}"
            exit 1
        fi
    fi
    """
}
