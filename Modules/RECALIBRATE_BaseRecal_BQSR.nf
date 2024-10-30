nextflow.enable.dsl=2

process RECALIBRATE_BaseRecal_BQSR {
    scratch true
    label 'RECALIBRATE'
    publishDir("${params.recal_dir}", mode: 'copy')
    errorStrategy = { task.attempt <= maxRetries ? 'retry' : 'ignore'}
    maxRetries 3
        

    input:
    tuple val(patient), val(meta), path(bam), path(bai)
    each chunk

    output:
    tuple val(patient), val(meta.status), val(meta), path("*bam"), emit: MergeBam_input
    path("*bai"), emit: bai

    script:
    if (meta.type == "exome")
        """
        /tscc/projects/ps-lalexandrov/shared/EVC_nextflow/gatk-4.6.0.0/gatk BaseRecalibrator \
        -I ${bam} \
        -R ${params.ref} \
        --known-sites ${params.database_dir}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
        --known-sites ${params.database_dir}/Homo_sapiens_assembly38.known_indels.vcf.gz \
        -L ${params.database_dir}/interval_list_20_exome/${chunk}-scattered.interval_list \
        -O ${meta.patient}_${meta.status}_${chunk}_recal.table

        /tscc/projects/ps-lalexandrov/shared/EVC_nextflow/gatk-4.6.0.0/gatk ApplyBQSR \
        -R ${params.ref} \
        -I ${bam} \
        -L ${params.database_dir}/interval_list_20_exome/${chunk}-scattered.interval_list \
        --bqsr-recal-file ${meta.patient}_${meta.status}_${chunk}_recal.table \
        -O ${meta.patient}_${meta.status}_recalibrated_${chunk}.bam
        """
    else
        """
        /tscc/projects/ps-lalexandrov/shared/EVC_nextflow/gatk-4.6.0.0/gatk BaseRecalibrator \
        -I ${bam} \
        -R ${params.ref} \
        --known-sites ${params.database_dir}/resources_broad_hg38_v0_Homo_sapiens_assembly38.dbsnp138.vcf \
        --known-sites ${params.database_dir}/Homo_sapiens_assembly38.known_indels.vcf.gz \
        -L ${params.database_dir}/interval_list_20/${chunk}-scattered.interval_list \
        -O ${meta.patient}_${meta.status}_${chunk}_recal.table

        /tscc/projects/ps-lalexandrov/shared/EVC_nextflow/gatk-4.6.0.0/gatk ApplyBQSR \
        -R ${params.ref} \
        -I ${bam} \
        -L ${params.database_dir}/interval_list_20/${chunk}-scattered.interval_list \
        --bqsr-recal-file ${meta.patient}_${meta.status}_${chunk}_recal.table \
        -O ${meta.patient}_${meta.status}_recalibrated_${chunk}.bam
        """

}
