process MuSE2 {
    publishDir("${params.muse2_dir}", mode: 'copy')
    scratch true
    label 'MuSE2'
    errorStrategy 'retry'
    maxRetries 3

    input:
    val(map)

    output:
    path("*.vcf"), emit: MuSE2_vcf
    val(map), emit:MuSE2_out


    script:
    if (map.type == "exome")
        """
        module load cpu
        module load curl
        module load bzip2
        module load htslib

        ${params.MuSE2} call -f ${params.ref} -n 8 -O ${map.patient} ${map.tumor} ${map.normal}

        ${params.MuSE2} sump -I ${map.patient}.MuSE.txt -n 8 -E -O ${map.patient}.vcf -D ${params.database_dir}/af-only-gnomad.hg38_no_alt.vcf.gz
        """
    else
        """
        module load cpu
        module load curl
        module load bzip2
        module load htslib

        ${params.MuSE2} call -f ${params.ref} -n 8 -O ${map.patient} ${map.tumor} ${map.normal}

        ${params.MuSE2} sump -I ${map.patient}.MuSE.txt -n 8 -G -O ${map.patient}.vcf -D ${params.database_dir}/af-only-gnomad.hg38_no_alt.vcf.gz
        """
}