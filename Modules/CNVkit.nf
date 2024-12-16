nextflow.enable.dsl=2

process CNVkit {
    conda "${params.cnvkit_env}"
    scratch true
    label 'process_high'
    publishDir("${params.cnvkit_dir}", mode: 'copy')
    errorStrategy 'retry'
    maxRetries 3

    input:
    val(map)
    path(cnn)

    output:
    tuple val(meta), path("*.bed"), emit: CNVkit_bed
    tuple val(meta), path("*.cn*"), emit: CNVkit_cn_files
    tuple val(meta), path("*.pdf"), emit: CNVkit_pdf
    tuple val(meta), path("*.png"), emit: CNVkit_png

    script:

    if (map.type == "exome")
        """
        cnvkit.py batch \
        ${map.tumor} \
        --targets ${params.database_dir}/GRCh38_exome.bed \
        -r ${cnn} \
        -p 16 \
        --scatter --diagram

        cnvkit.py call \
        "${map.patient}_tumor_recal.cns" \
        -o ${map.patient}_calls.cns
        """
    else
        """
        cnvkit.py batch \
        ${map.tumor} \
        -r ${cnn} \
        --method wgs \
        -p 16 \
        --scatter --diagram

        cnvkit.py call \
        "${map.patient}_tumor_recal.cns" \
        -o ${map.patient}_calls.cns
        """
}