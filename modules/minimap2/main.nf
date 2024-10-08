#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process MINIMAP2_ALIGN {
    tag "Minimap2 alignment for ${reads.baseName}"
    conda 'bioconda::minimap2 bioconda::samtools'
    publishDir "${params.outdir}/minimap2_align", mode: 'symlink'

    input:
    path ref
    path reads
    
    output:
    path "${reads.baseName}.sort.bam", emit: ch_align_bam
    path "${reads.baseName}.sort.bam.bai", emit: ch_align_bam_bai

    script:
    """
    # Log the inputs
    echo "Running Minimap2 alignment:"
    echo "Reference: ${ref}"
    echo "Reads: ${reads}"
    echo "Threads: ${params.threads}"
    echo "Platform: ${params.platform}"
    echo "Type: ${params.datatype}"

    # minimap2 -ax splice:hq -uf ${ref} ${reads} -t ${params.threads} --secondary=no > ${reads.baseName}.sam


    if [ ${params.platform} == "ont" ]; then
        if [ ${params.datatype} == "cDNA" ]; then
            minimap2 -ax splice ${ref} ${reads} -t ${params.threads} --secondary=no | samtools view -bSh -F 2308 - > ${reads.baseName}.bam
        elif [ ${params.datatype} == "dRNA" ]; then
            minimap2 -ax splice -uf -k14 ${ref} ${reads} -t ${params.threads} --secondary=no | samtools view -bSh -F 2308 - > ${reads.baseName}.bam
        else
            echo "Error: For 'ont' params.platform, params.datatype must be 'cDNA' or 'dRNA'. Given params.datatype: ${params.datatype}"
            exit 1
        fi
    elif [ ${params.platform} == "pb" ]; then
        if [ ${params.datatype} == "isoseq" ] || [ ${params.datatype} == "masseq" ]; then
            minimap2 -ax splice:hq -uf ${ref} ${reads} -t ${params.threads} --secondary=no | samtools view -bSh -F 2308 - > ${reads.baseName}.bam
        else
            echo "Error: For 'pb' params.platform, params.datatype must be 'isoseq' or 'masseq'. Given params.datatype: ${params.datatype}"
            exit 1
        fi
    else
        echo "Error: Platform must be 'ont' or 'pb'. Given params.platform: ${params.platform}"
        exit 1
    fi


    samtools sort -@ ${params.threads} -o ${reads.baseName}.sort.bam ${reads.baseName}.bam
    samtools index -@ ${params.threads} ${reads.baseName}.sort.bam
    """
}
