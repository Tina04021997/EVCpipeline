nextflow.enable.dsl=2

process FilterMutectCalls {
    scratch true
    label 'process_low'
    publishDir("${params.MUTECT2_dir}", mode: 'copy')
    errorStrategy = { task.attempt <= maxRetries ? 'retry' : 'ignore'}
    maxRetries 3

    input:
    tuple val(patient), path(unfiltered_vcf)
    tuple val(patient), path(contamination_table)
    tuple val(patient), path(segments_table)
    tuple val(patient), path(read_orientation_model_tar)
    tuple val(patient), path(merged_stats)

    output:
    tuple val(patient), path("*vcf"), path("*idx"), path("*stats"), emit: MUTECT2_final_out
    val (patient), emit: Mutect2_out

    script:
    """
    /tscc/projects/ps-lalexandrov/shared/EVC_nextflow/gatk-4.6.0.0/gatk FilterMutectCalls -R ${params.ref} -V ${unfiltered_vcf} --contamination-table ${contamination_table} --ob-priors ${read_orientation_model_tar} -O ${patient}_mutect2_filtered.vcf --stats ${merged_stats} --filtering-stats ${patient}_mutect2_filtered.stats --tumor-segmentation ${segments_table}
    """
}
