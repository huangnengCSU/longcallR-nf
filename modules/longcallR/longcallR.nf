process INSTALL_LONGCALLR {
    tag "Install longcallR"

    input:
    val prev_ch

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

process ONLY_LONGCALLR_CALL_PHASE {
    tag "only use longcallR to call SNPs and phase"
    conda 'bioconda::samtools bioconda::bcftools bioconda::tabix'
    publishDir "${params.outdir}/only_longcallR", mode: 'symlink'

    input:
    path longcallr_binary
    path bam_file
    path bam_index
    path ref_file
    path ref_index
    each contig  // Channel of contigs
    val prev_ch

    // Control the number of concurrent jobs with `maxForks`
    maxForks params.num_jobs
    // Set how many threads each job can use
    cpus params.threads_per_job

    output:
    path "${params.sample_name}_only_longcallR_${contig}.vcf", emit: only_longcallR_vcfs_ch
    path "${params.sample_name}_only_longcallR_${contig}.phased.bam", emit: only_longcallR_phased_bams_ch

    script:
    """
    if [ ${params.platform} == "pb" ]; then
        if [ ${params.datatype} == "isoseq" ]; then
            ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -o ${params.sample_name}_only_longcallR_${contig} --platform "hifi" --preset "hifi-isoseq" -x ${contig} -t ${params.threads_per_job}
        elif [ ${params.datatype} == "masseq" ]; then
            ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -o ${params.sample_name}_only_longcallR_${contig} --platform "hifi" --preset "hifi-masseq" -x ${contig} -t ${params.threads_per_job}
        else
            echo "Error: For 'pb' params.platform, params.datatype must be 'isoseq' or 'masseq'. Given params.datatype: ${params.datatype}"
            exit 1
        fi
    elif [ ${params.platform} == "ont" ]; then
        if [ ${params.datatype} == "cDNA" ]; then
            ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -o ${params.sample_name}_only_longcallR_${contig} --platform "ont" --preset "ont-cdna" -x ${contig} -t ${params.threads_per_job}
        elif [ ${params.datatype} == "dRNA" ]; then
            ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -o ${params.sample_name}_only_longcallR_${contig} --platform "ont" --preset "ont-drna" -x ${contig} -t ${params.threads_per_job}
        else
            echo "Error: For 'ont' params.platform, params.datatype must be 'cDNA' or 'dRNA'. Given params.datatype: ${params.datatype}"
            exit 1
        fi
    fi
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

    // Control the number of concurrent jobs with `maxForks`
    maxForks params.num_jobs
    // Set how many threads each job can use
    cpus params.threads_per_job

    output:
    path "${params.sample_name}_longcallR_${contig}.vcf", emit: longcallR_vcfs_ch
    path "${params.sample_name}_longcallR_${contig}.phased.bam", emit: longcallR_phased_bams_ch

    script:
    """
    if [[ -z "${vcf_file}" ]]; then
        if [ ${params.platform} == "pb" ]; then
            if [ ${params.datatype} == "isoseq" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -o ${params.sample_name}_longcallR_${contig} --platform "hifi" --preset "hifi-isoseq" -x ${contig} -t ${params.threads_per_job}
            elif [ ${params.datatype} == "masseq" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -o ${params.sample_name}_longcallR_${contig} --platform "hifi" --preset "hifi-masseq" -x ${contig} -t ${params.threads_per_job}
            else
                echo "Error: For 'pb' params.platform, params.datatype must be 'isoseq' or 'masseq'. Given params.datatype: ${params.datatype}"
                exit 1
            fi
        elif [ ${params.platform} == "ont" ]; then
            if [ ${params.datatype} == "cDNA" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -o ${params.sample_name}_longcallR_${contig} --platform "ont" --preset "ont-cdna" -x ${contig} -t ${params.threads_per_job}
            elif [ ${params.datatype} == "dRNA" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -o ${params.sample_name}_longcallR_${contig} --platform "ont" --preset "ont-drna" -x ${contig} -t ${params.threads_per_job}
            else
                echo "Error: For 'ont' params.platform, params.datatype must be 'cDNA' or 'dRNA'. Given params.datatype: ${params.datatype}"
                exit 1
            fi
        fi
    else
        if [ ${params.platform} == "pb" ]; then
            if [ ${params.datatype} == "isoseq" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -o ${params.sample_name}_longcallR_${contig} --input-vcf ${vcf_file} --platform "hifi" --preset "hifi-isoseq" -x ${contig} -t ${params.threads_per_job}
            elif [ ${params.datatype} == "masseq" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -o ${params.sample_name}_longcallR_${contig} --input-vcf ${vcf_file} --platform "hifi" --preset "hifi-masseq" -x ${contig} -t ${params.threads_per_job}
            else
                echo "Error: For 'pb' params.platform, params.datatype must be 'isoseq' or 'masseq'. Given params.datatype: ${params.datatype}"
                exit 1
            fi
        elif [ ${params.platform} == "ont" ]; then
            if [ ${params.datatype} == "cDNA" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -o ${params.sample_name}_longcallR_${contig} --input-vcf ${vcf_file} --platform "ont" --preset "ont-cdna" -x ${contig} -t ${params.threads_per_job}
            elif [ ${params.datatype} == "dRNA" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -o ${params.sample_name}_longcallR_${contig} --input-vcf ${vcf_file} --platform "ont" --preset "ont-drna" -x ${contig} -t ${params.threads_per_job}
            else
                echo "Error: For 'ont' params.platform, params.datatype must be 'cDNA' or 'dRNA'. Given params.datatype: ${params.datatype}"
                exit 1
            fi
        fi
    fi
    """
}

process LONGCALLR_CALL_PHASE_EXON {
    tag "longcallR SNP calling and phasing only for exon region"
    conda 'bioconda::samtools bioconda::bcftools bioconda::tabix'
    publishDir "${params.outdir}/longcallR_exon", mode: 'symlink'

    input:
    path longcallr_binary
    path bam_file
    path bam_index
    path ref_file
    path ref_index
    path vcf_file
    path vcf_index
    path annotation
    each contig  // Channel of contigs
    val prev_ch

    // Control the number of concurrent jobs with `maxForks`
    maxForks params.num_jobs
    // Set how many threads each job can use
    cpus params.threads_per_job

    output:
    path "${params.sample_name}_longcallR_exon_${contig}.vcf", emit: longcallR_exon_vcfs_ch
    path "${params.sample_name}_longcallR_exon_${contig}.phased.bam", emit: longcallR_exon_phased_bams_ch

    script:
    """

    gzip -fdc ${annotation} > ${annotation.baseName}

    if [[ -z "${vcf_file}" ]]; then
        if [ ${params.platform} == "pb" ]; then
            if [ ${params.datatype} == "isoseq" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -a ${annotation.baseName} -o ${params.sample_name}_longcallR_exon_${contig} --platform "hifi" --preset "hifi-isoseq" -x ${contig} -t ${params.threads_per_job}
            elif [ ${params.datatype} == "masseq" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -a ${annotation.baseName} -o ${params.sample_name}_longcallR_exon_${contig} --platform "hifi" --preset "hifi-masseq" -x ${contig} -t ${params.threads_per_job}
            else
                echo "Error: For 'pb' params.platform, params.datatype must be 'isoseq' or 'masseq'. Given params.datatype: ${params.datatype}"
                exit 1
            fi
        elif [ ${params.platform} == "ont" ]; then
            if [ ${params.datatype} == "cDNA" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -a ${annotation.baseName} -o ${params.sample_name}_longcallR_exon_${contig} --platform "ont" --preset "ont-cdna" -x ${contig} -t ${params.threads_per_job}
            elif [ ${params.datatype} == "dRNA" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -a ${annotation.baseName} -o ${params.sample_name}_longcallR_exon_${contig} --platform "ont" --preset "ont-drna" -x ${contig} -t ${params.threads_per_job}
            else
                echo "Error: For 'ont' params.platform, params.datatype must be 'cDNA' or 'dRNA'. Given params.datatype: ${params.datatype}"
                exit 1
            fi
        fi
    else
        if [ ${params.platform} == "pb" ]; then
            if [ ${params.datatype} == "isoseq" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -a ${annotation.baseName} -o ${params.sample_name}_longcallR_exon_${contig} --input-vcf ${vcf_file} --platform "hifi" --preset "hifi-isoseq" -x ${contig} -t ${params.threads_per_job}
            elif [ ${params.datatype} == "masseq" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -a ${annotation.baseName} -o ${params.sample_name}_longcallR_exon_${contig} --input-vcf ${vcf_file} --platform "hifi" --preset "hifi-masseq" -x ${contig} -t ${params.threads_per_job}
            else
                echo "Error: For 'pb' params.platform, params.datatype must be 'isoseq' or 'masseq'. Given params.datatype: ${params.datatype}"
                exit 1
            fi
        elif [ ${params.platform} == "ont" ]; then
            if [ ${params.datatype} == "cDNA" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -a ${annotation.baseName} -o ${params.sample_name}_longcallR_exon_${contig} --input-vcf ${vcf_file} --platform "ont" --preset "ont-cdna" -x ${contig} -t ${params.threads_per_job}
            elif [ ${params.datatype} == "dRNA" ]; then
                ./${longcallr_binary} -b ${bam_file} -f ${ref_file} -a ${annotation.baseName} -o ${params.sample_name}_longcallR_exon_${contig} --input-vcf ${vcf_file} --platform "ont" --preset "ont-drna" -x ${contig} -t ${params.threads_per_job}
            else
                echo "Error: For 'ont' params.platform, params.datatype must be 'cDNA' or 'dRNA'. Given params.datatype: ${params.datatype}"
                exit 1
            fi
        fi
    fi
    """
}