nextflow.enable.dsl=2

process FASTQC {
    scratch true
    label 'process_medium'
    publishDir("${params.FASTQC_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    val(meta)

    output:
    path("*.html"), emit: html
    path("*.zip") , emit: zip
    tuple val(meta), path("*.html"), emit: fastqc_out

    script:
    """
    ${params.FASTQC}/fastqc -t 8 -o ./ ${meta.fastq_1} ${meta.fastq_2}
    """
}