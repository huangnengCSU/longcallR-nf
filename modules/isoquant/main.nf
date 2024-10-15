process ISOQUANT {
    tag "Isoquant for ${params.sample_name}"
    conda "bioconda::isoquant=3.5.0"
    publishDir "${params.outdir}/isoquant", mode: 'symlink'

    input:
    path bam
    path bam_bai
    path ref
    path ref_fai
    path annotation

    output:
    path "isoquant_output/${params.sample_name}/${params.sample_name}.*", emit: isoquant_outputs_ch


    script:
    """
    if [ ${params.platform} == "ont" ]; then
        isoquant.py \
        --reference ${ref} \
        --genedb ${annotation} \
        --bam ${bam} \
        --complete_genedb \
        --data_type nanopore \
        --fl_data \
        --transcript_quantification unique_only \
        --matching_strategy precise \
        --splice_correction_strategy default_ont \
        --model_construction_strategy default_ont \
        --min_mapq 10 \
        --inconsistent_mapq_cutoff 10 \
        --output isoquant_output \
        --prefix ${params.sample_name} \
        -t ${params.threads}
    elif [ ${params.platform} == "pb" ]; then
        isoquant.py \
        --reference ${ref} \
        --genedb ${annotation} \
        --bam ${bam} \
        --complete_genedb \
        --data_type pacbio_ccs \
        --fl_data \
        --transcript_quantification unique_only \
        --matching_strategy precise \
        --splice_correction_strategy default_pacbio \
        --model_construction_strategy fl_pacbio \
        --min_mapq 10 \
        --inconsistent_mapq_cutoff 10 \
        --output isoquant_output \
        --prefix ${params.sample_name} \
        -t ${params.threads}
    else
        echo "Error: Platform must be 'ont' or 'pb'. Given params.platform: ${params.platform}"
        exit 1
    fi
    """

}