process INSTALL_LONGCALLR {
    tag "Install longcallR"

    output:
    path "longcallR/target/release/longcallR", emit: longcallr_binary

    script:
    """
    # Clone the longcallR repository and install dependencies
    git clone https://github.com/huangnengCSU/longcallR.git
    cd longcallR
    # Build longcallR using Cargo (Rust's package manager)
    cargo build --release
    """
}

process LONGCALLR_CALL_PHASE {
    tag "longcallR SNP calling and phasing"
    conda 'bioconda::samtools bioconda::bcftools bioconda::tabix'
    publishDir "${params.outdir}/longcallR", mode: 'symlink'

    input:
    path longcallr_binary
    path bam_file
    path bam_index
    path ref_file
    path ref_index
    path vcf_file
    path vcf_index
    each contig  // Channel of contigs

    output:
    path "${params.sample_name}_longcallR_${contig}.vcf", emit: longcallR_vcfs_ch
    path "${params.sample_name}_longcallR_${contig}.phased.bam", emit: longcallR_phased_bams_ch

    script:
    """
    if [[ -z "${vcf_file}" ]]; then
        if [ ${params.platform} == "pb" ]; then
            if [ ${params.datatype} == "isoseq" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -o ${params.sample_name}_longcallR_${contig} --platform "hifi" --preset "hifi-isoseq" -x ${contig} -t ${params.threads}
            elif [ ${params.datatype} == "masseq" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -o ${params.sample_name}_longcallR_${contig} --platform "hifi" --preset "hifi-masseq" -x ${contig} -t ${params.threads}
            else
                echo "Error: For 'pb' params.platform, params.datatype must be 'isoseq' or 'masseq'. Given params.datatype: ${params.datatype}"
                exit 1
            fi
        elif [ ${params.platform} == "ont" ]; then
            if [ ${params.datatype} == "cDNA" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -o ${params.sample_name}_longcallR_${contig} --platform "ont" --preset "ont-cdna" -x ${contig} -t ${params.threads}
            elif [ ${params.datatype} == "dRNA" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -o ${params.sample_name}_longcallR_${contig} --platform "ont" --preset "ont-drna" -x ${contig} -t ${params.threads}
            else
                echo "Error: For 'ont' params.platform, params.datatype must be 'cDNA' or 'dRNA'. Given params.datatype: ${params.datatype}"
                exit 1
            fi
        fi
    else
        if [ ${params.platform} == "pb" ]; then
            if [ ${params.datatype} == "isoseq" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -o ${params.sample_name}_longcallR_${contig} --input-vcf ${vcf_file} --platform "hifi" --preset "hifi-isoseq" -x ${contig} -t ${params.threads}
            elif [ ${params.datatype} == "masseq" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -o ${params.sample_name}_longcallR_${contig} --input-vcf ${vcf_file} --platform "hifi" --preset "hifi-masseq" -x ${contig} -t ${params.threads}
            else
                echo "Error: For 'pb' params.platform, params.datatype must be 'isoseq' or 'masseq'. Given params.datatype: ${params.datatype}"
                exit 1
            fi
        elif [ ${params.platform} == "ont" ]; then
            if [ ${params.datatype} == "cDNA" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -o ${params.sample_name}_longcallR_${contig} --input-vcf ${vcf_file} --platform "ont" --preset "ont-cdna" -x ${contig} -t ${params.threads}
            elif [ ${params.datatype} == "dRNA" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -o ${params.sample_name}_longcallR_${contig} --input-vcf ${vcf_file} --platform "ont" --preset "ont-drna" -x ${contig} -t ${params.threads}
            else
                echo "Error: For 'ont' params.platform, params.datatype must be 'cDNA' or 'dRNA'. Given params.datatype: ${params.datatype}"
                exit 1
            fi
        fi
    fi
    """
}