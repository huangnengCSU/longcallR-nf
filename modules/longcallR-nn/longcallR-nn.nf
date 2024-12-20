process LONGCALLR_NN_CALL {
    tag "${params.sample_name}_longcallR_nn_call"
    publishDir "${params.outdir}/longcallR_nn_call", mode: 'symlink'

    input:
        each contig_feature
        path ref    // Channel of reference genome
        path ref_fai    // Channel of reference genome index file


    output:
    path "models/*"
    path "${params.sample_name}_predictions/*"
    path "${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_${contig_feature.baseName}.vcf", emit: longcallR_nn_vcfs_ch

    script:
    """
    mkdir -p ${params.sample_name}_predictions
    mkdir -p models

    longcallR_nn download -d models

    if [ ${params.no_cuda} == "true" ]; then
        if [ ${params.platform} == "ont" ]; then
            if [ ${params.datatype} == "cDNA" ]; then
                longcallR_nn call \
                -config models/ont_cdna_config.yaml \
                -model models/ont_cdna_model.chkpt \
                -data ${contig_feature} \
                -ref ${ref} \
                -output ${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_${contig_feature.baseName}.vcf \
                -max_depth ${params.max_depth} \
                -batch_size ${params.batch_size} \
                --no_cuda
            elif [ ${params.datatype} == "dRNA" ]; then
                longcallR_nn call \
                -config models/ont_drna_config.yaml \
                -model models/ont_drna_model.chkpt \
                -data ${contig_feature} \
                -ref ${ref} \
                -output ${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_${contig_feature.baseName}.vcf \
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
                -data ${contig_feature} \
                -ref ${ref} \
                -output ${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_${contig_feature.baseName}.vcf \
                -max_depth ${params.max_depth} \
                -batch_size ${params.batch_size} \
                --no_cuda
            elif [ ${params.datatype} == "masseq" ]; then
                longcallR_nn call \
                -config models/pb_masseq_config.yaml \
                -model models/pb_masseq_model.chkpt \
                -data ${contig_feature} \
                -ref ${ref} \
                -output ${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_${contig_feature.baseName}.vcf \
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
                -data ${contig_feature} \
                -ref ${ref} \
                -output ${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_${contig_feature.baseName}.vcf \
                -max_depth ${params.max_depth} \
                -batch_size ${params.batch_size}
            elif [ ${params.datatype} == "dRNA" ]; then
                longcallR_nn call \
                -config models/ont_drna_config.yaml \
                -model models/ont_drna_model.chkpt \
                -data ${contig_feature} \
                -ref ${ref} \
                -output ${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_${contig_feature.baseName}.vcf \
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
                -data ${contig_feature} \
                -ref ${ref} \
                -output ${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_${contig_feature.baseName}.vcf \
                -max_depth ${params.max_depth} \
                -batch_size ${params.batch_size}
            elif [ ${params.datatype} == "masseq" ]; then
                longcallR_nn call \
                -config models/pb_masseq_config.yaml \
                -model models/pb_masseq_model.chkpt \
                -data ${contig_feature} \
                -ref ${ref} \
                -output ${params.sample_name}_predictions/${params.sample_name}_longcallR_nn_${contig_feature.baseName}.vcf \
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
