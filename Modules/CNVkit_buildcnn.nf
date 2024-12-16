nextflow.enable.dsl=2

process CNVkit_buildcnn {
    conda "${params.cnvkit_env}"
    scratch true
    label 'process_medium'
    publishDir("${params.cnvkit_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(patient), val(meta), path(bam), path(bai)

    output:
    path("reference.cnn"), emit: CNVkit_ref_cnn
    path("*coverage.cnn"), emit: CNVkit_buildcnns

    script:

    if (meta.type == "exome")
        """
        cnvkit.py batch \
        ${bam} \
        -n \
        --targets ${params.database_dir}/GRCh38_exome.bed \
        --fasta ${params.ref} \
        --output-reference reference.cnn \
        -p 8 
        """
    else
        """
        cnvkit.py batch \
        ${bam} \
        -n \
        --method wgs \
        --fasta ${params.ref} \
        --output-reference reference.cnn \
        -p 8
        """


}