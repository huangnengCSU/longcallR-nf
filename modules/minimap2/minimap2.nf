#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process MINIMAP2_ALIGN {
    tag "Minimap2 alignment for ${params.sample_name}"
    publishDir "${params.outdir}/minimap2_align", mode: 'symlink'

    input:
    path ref
    path reads
    val prev_ch
    
    output:
    path "${params.sample_name}.sort.bam", emit: ch_align_bam
    path "${params.sample_name}.sort.bam.bai", emit: ch_align_bam_bai
    path "${ref}.fai", emit: ch_ref_fai

    script:
    """
    # Log the inputs
    echo "Running Minimap2 alignment:"
    echo "Reference: ${ref}"
    echo "Reads: ${reads}"
    echo "Threads: ${params.threads}"
    echo "Platform: ${params.platform}"
    echo "Type: ${params.datatype}"

    # minimap2 -ax splice:hq -uf ${ref} ${reads} -t ${params.threads} --secondary=no > ${params.sample_name}.sam

    if [ ! -f "${ref}.fai" ]; then
        echo "Index for reference ${ref} not found. Generating .fai index..."
        samtools faidx ${ref}
    else
        echo ".fai index for ${ref} already exists."
    fi


    if [ ${params.platform} == "ont" ]; then
        if [ ${params.datatype} == "cDNA" ]; then
            minimap2 -ax splice ${ref} ${reads} -t ${params.threads} --secondary=no | samtools view -bSh -F 2308 - > ${params.sample_name}.bam
        elif [ ${params.datatype} == "dRNA" ]; then
            minimap2 -ax splice -uf -k14 ${ref} ${reads} -t ${params.threads} --secondary=no | samtools view -bSh -F 2308 - > ${params.sample_name}.bam
        else
            echo "Error: For 'ont' params.platform, params.datatype must be 'cDNA' or 'dRNA'. Given params.datatype: ${params.datatype}"
            exit 1
        fi
    elif [ ${params.platform} == "pb" ]; then
        if [ ${params.datatype} == "isoseq" ]; then
            minimap2 -ax splice:hq ${ref} ${reads} -t ${params.threads} --secondary=no | samtools view -bSh -F 2308 - > ${params.sample_name}.bam
        elif [ ${params.datatype} == "masseq" ]; then
            minimap2 -ax splice:hq -uf ${ref} ${reads} -t ${params.threads} --secondary=no | samtools view -bSh -F 2308 - > ${params.sample_name}.bam
        else
            echo "Error: For 'pb' params.platform, params.datatype must be 'isoseq' or 'masseq'. Given params.datatype: ${params.datatype}"
            exit 1
        fi
    else
        echo "Error: Platform must be 'ont' or 'pb'. Given params.platform: ${params.platform}"
        exit 1
    fi

    samtools --version
    samtools sort -@ ${params.threads} -o ${params.sample_name}.sort.bam ${params.sample_name}.bam
    samtools index -@ ${params.threads} ${params.sample_name}.sort.bam
    """
}
