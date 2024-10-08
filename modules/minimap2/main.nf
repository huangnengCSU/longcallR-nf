#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process MINIMAP2_ALIGN {
    tag "Minimap2 alignment for ${reads.baseName}"
    conda 'bioconda::minimap2 bioconda::samtools'
    publishDir "${params.outdir}/minimap2_align", mode: 'symlink'

    input:
    path ref
    path reads
    val threads
    val platform
    val type
    
    output:
    path "${reads.baseName}.sort.bam*"

    script:
    """
    # Log the inputs
    echo "Running Minimap2 alignment:"
    echo "Reference: ${ref}"
    echo "Reads: ${reads}"
    echo "Threads: ${threads}"
    echo "Platform: ${platform}"
    echo "Type: ${type}"

    # minimap2 -ax splice:hq -uf ${ref} ${reads} -t ${threads} --secondary=no > ${reads.baseName}.sam


    if [ ${platform} == "ont" ]; then
        if [ ${type} == "cDNA" ]; then
            minimap2 -ax splice ${ref} ${reads} -t ${threads} --secondary=no | samtools view -bSh -F 2308 - > ${reads.baseName}.bam
        elif [ ${type} == "dRNA" ]; then
            minimap2 -ax splice -uf -k14 ${ref} ${reads} -t ${threads} --secondary=no | samtools view -bSh -F 2308 - > ${reads.baseName}.bam
        else
            echo "Error: For 'ont' platform, type must be 'cDNA' or 'dRNA'. Given type: ${type}"
            exit 1
        fi
    elif [ ${platform} == "pb" ]; then
        if [ ${type} == "isoseq" ] || [ ${type} == "masseq" ]; then
            minimap2 -ax splice:hq -uf ${ref} ${reads} -t ${threads} --secondary=no | samtools view -bSh -F 2308 - > ${reads.baseName}.bam
        else
            echo "Error: For 'pb' platform, type must be 'isoseq' or 'masseq'. Given type: ${type}"
            exit 1
        fi
    else
        echo "Error: Platform must be 'ont' or 'pb'. Given platform: ${platform}"
        exit 1
    fi


    samtools sort -@ ${threads} -o ${reads.baseName}.sort.bam ${reads.baseName}.bam
    samtools index -@ ${threads} ${reads.baseName}.sort.bam
    """
}
