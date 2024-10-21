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

process LONGCALLR_DP {
    tag "${params.sample_name}_longcallR_dp"
    publishDir "${params.outdir}/longcallR_dp", mode: 'symlink'

    input:
        path longcallr_dp_binary
        each contig  // Channel of contigs
        path bam   // Channel of BAM file
        path bam_bai    // Channel of BAM index file
        path ref    // Channel of reference genome
        path ref_fai    // Channel of reference genome index file

    output:
    // path "${params.sample_name}_features", emit: longcallR_dp_features
    path "${params.sample_name}_features/*", emit: longcallR_dp_features_contigs_ch

    // Control the number of concurrent jobs with `maxForks`
    maxForks params.num_jobs
    // Set how many threads each job can use
    cpus params.threads_per_job

    script:
    """
    mkdir -p ${params.sample_name}_features
    
    ./${longcallr_dp_binary} \
    --mode predict \
    --bam-path '${bam}' \
    --ref-path '${ref}' \
    --min-depth '${params.min_depth}' \
    --min-alt-freq '${params.min_af}' \
    --min-baseq '${params.min_bq}' \
    --threads '${params.threads_per_job}' \
    --contigs '${contig}' \
    --output '${params.sample_name}_features/${contig}'
    """
}