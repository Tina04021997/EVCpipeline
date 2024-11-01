nextflow.enable.dsl=2

params.sample = "sample.csv"
params.ref="/tscc/projects/ps-lalexandrov/shared/EVC_nextflow/GRCh38_ref/GRCh38.d1.vd1.fa"
params.bam_dir="$projectDir/RESULTS/BAM"
params.report_dir="$projectDir/RESULTS/REPORT"
params.bed="/tscc/projects/ps-lalexandrov/shared/EVC_nextflow/GRCh38_ref/hg38_chr.bed.gz"
params.database_dir="/tscc/projects/ps-lalexandrov/shared/EVC_nextflow/Databases"
params.mkdup_temp_dir="/your/own/temp/MKDUP/folder/in/restricted/"
params.FASTQC="/tscc/projects/ps-lalexandrov/shared/EVC_nextflow/FastQC"
params.FASTQC_dir="$projectDir/RESULTS/FASTQC"
params.mkdup_dir="$projectDir/RESULTS/MKDUP"
params.recal_dir="$projectDir/RESULTS/RECALIBRATE"
params.SAGE_java="/tscc/projects/ps-lalexandrov/shared/EVC_nextflow/SAGE/sage_v3.3.jar"
params.SAGE_ref_dir="/tscc/projects/ps-lalexandrov/shared/EVC_nextflow/SAGE/v5_34/ref/38"
params.SAGE_dir="$projectDir/RESULTS/SAGE"
params.strelka_dir="$projectDir/RESULTS/STRELKA"
params.MuSE2="/tscc/projects/ps-lalexandrov/shared/EVC_nextflow/MuSE/MuSE"
params.muse2_dir="$projectDir/RESULTS/MuSE2"
params.MUTECT2_dir="$projectDir/RESULTS/Mutect2"
params.mosdepth_dir="$projectDir/RESULTS/mosdepth"
params.conpair_dir="$projectDir/RESULTS/Conpair"
params.conpair="/tscc/projects/ps-lalexandrov/shared/EVC_nextflow/Conpair-0.2"
params.jre="/tscc/projects/ps-lalexandrov/shared/EVC_nextflow/jre1.8.0_401"

params.bwamem2_env = "/tscc/projects/ps-lalexandrov/shared/EVC_nextflow/yml/bwamem2.yml"
params.mkdup_env = "/tscc/projects/ps-lalexandrov/shared/EVC_nextflow/yml/mkdup.yml"
params.conpair_env = "/tscc/projects/ps-lalexandrov/shared/EVC_nextflow/yml/conpair.yml"
params.samtools_env = "/tscc/projects/ps-lalexandrov/shared/EVC_nextflow/yml/samtools.yml"
params.strelka_env = "/tscc/projects/ps-lalexandrov/shared/EVC_nextflow/yml/strelka_env.yml"
params.mosdepth_env = "/tscc/projects/ps-lalexandrov/shared/EVC_nextflow/yml/mosdepth_env.yml"
params.summary_env = "/tscc/projects/ps-lalexandrov/shared/EVC_nextflow/yml/py_summary.yml"

include { FASTQC } from './Modules/FASTQC'
include { BWA_MEM } from './Modules/BWA_MEM'
include { CHECK_BAM_BWA } from './Modules/CHECK_BAM_BWA'
include { COMBINE_REPORTS_BWA } from './Modules/COMBINE_REPORTS_BWA'
include { RECALIBRATE_BaseRecal_BQSR } from './Modules/RECALIBRATE_BaseRecal_BQSR.nf'
include { RECALIBRATE_MergeBam } from './Modules/RECALIBRATE_MergeBam.nf'
include { RECALIBRATE_SortBam } from './Modules/RECALIBRATE_SortBam.nf'
include { CHECK_BAM_RECAL } from './Modules/CHECK_BAM_RECAL'
include { COMBINE_REPORTS_RECAL } from './Modules/COMBINE_REPORTS_RECAL'
include { MKDUP } from './Modules/MKDUP'
include { CHECK_BAM_MKDUP } from './Modules/CHECK_BAM_MKDUP'
include { COMBINE_REPORTS_MKDUP } from './Modules/COMBINE_REPORTS_MKDUP'
include { MOSDEPTH } from './Modules/MOSDEPTH'
include { CONPAIR } from './Modules/CONPAIR'
include { SAGE } from './Modules/SAGE'
include { STRELKA } from './Modules/STRELKA'
include { MuSE2 } from './Modules/MuSE2'

include { MUTECT2_CALLING } from './Modules/Mutect2/MUTECT2_CALLING' 
include { GETpileUP } from './Modules/Mutect2/GETpileUP' 
include { GETpileUP_Merge } from './Modules/Mutect2/GETpileUP_Merge'
include { LearnReadOrientationModel } from './Modules/Mutect2/LearnReadOrientationModel' 
include { MergeMutectStats } from './Modules/Mutect2/MergeMutectStats'
include { MergeVcfs } from './Modules/Mutect2/MergeVcfs' 
include { CalculateContamination } from './Modules/Mutect2/CalculateContamination' 
include { FilterMutectCalls } from './Modules/Mutect2/FilterMutectCalls'

include { SUMMARY } from './Modules/SUMMARY.nf'

include { SAVE_CSV_FASTQC } from './Modules/SAVE_CSV/SAVE_CSV_FASTQC'
include { SAVE_CSV_BWA_MEM } from './Modules/SAVE_CSV/SAVE_CSV_BWA_MEM'
include { SAVE_CSV_MKDUP } from './Modules/SAVE_CSV/SAVE_CSV_MKDUP'
include { SAVE_CSV_RECAL } from './Modules/SAVE_CSV/SAVE_CSV_RECAL'
include { SAVE_CSV_SAGE } from './Modules/SAVE_CSV/SAVE_CSV_SAGE'
include { SAVE_CSV_STRELKA } from './Modules/SAVE_CSV/SAVE_CSV_STRELKA'
include { SAVE_CSV_MuSE2 } from './Modules/SAVE_CSV/SAVE_CSV_MuSE2'
include { SAVE_CSV_Mutect2 } from './Modules/SAVE_CSV/SAVE_CSV_Mutect2'
include { SAVE_CSV_CONPAIR } from './Modules/SAVE_CSV/SAVE_CSV_CONPAIR'
include { SAVE_CSV_MOSDEPTH } from './Modules/SAVE_CSV/SAVE_CSV_MOSDEPTH'

workflow {
    chunk = Channel.of(1..20)

    Channel.fromPath(params.sample)
    | splitCsv( header:true )
    | map { row ->
        meta = row.subMap('patient','status', 'fastq_1', 'fastq_2', 'type')
    } | set { sample_sheet }

    //1
    FASTQC(sample_sheet)
    SAVE_CSV_FASTQC(FASTQC.out.fastqc_out,params.report_dir,params.FASTQC_dir)

    //2
    BWA_MEM(sample_sheet).set { BWA_MEM_out }
    SAVE_CSV_BWA_MEM(BWA_MEM_out.bam, params.report_dir, params.bam_dir)
    CHECK_BAM_BWA(BWA_MEM_out.bam).set{ CHECK_BAM_BWA_out }
    COMBINE_REPORTS_BWA(CHECK_BAM_BWA_out.individual_reports_bwa.collect())

    //3
    MKDUP(BWA_MEM_out).set{ MKDUP_out }
    SAVE_CSV_MKDUP(MKDUP_out.mkdup_bam, params.report_dir, params.mkdup_dir)
    CHECK_BAM_MKDUP(MKDUP_out.mkdup_bam).set{ CHECK_BAM_MKDUP_out }
    COMBINE_REPORTS_MKDUP(CHECK_BAM_MKDUP_out.individual_reports_mkdup.collect())

    //4
    RECALIBRATE_BaseRecal_BQSR(MKDUP_out.pair_mutect, chunk).set{ RECALIBRATE_BaseRecal_BQSR_out }

    RECALIBRATE_BaseRecal_BQSR_out.MergeBam_input.groupTuple(by:[0,1]).set { BaseRecal_BQSR_out_pair }
    BaseRecal_BQSR_out_pair.map{
        [patient:it[0], status:it[1], meta:it[2][0], bam:it[3]]
    }.set{BaseRecal_BQSR_out_MAP}

    RECALIBRATE_MergeBam(BaseRecal_BQSR_out_MAP).set{ RECALIBRATE_MergeBam_out }
    RECALIBRATE_SortBam(RECALIBRATE_MergeBam_out.MergeBam_input).set{ RECALIBRATE_out }

    SAVE_CSV_RECAL(RECALIBRATE_out.recal_bam, params.report_dir, params.recal_dir)
    CHECK_BAM_RECAL(RECALIBRATE_out.pair_recal).set{ CHECK_BAM_RECAL_out }
    COMBINE_REPORTS_RECAL(CHECK_BAM_RECAL_out.individual_reports_recal.collect())

    RECALIBRATE_out.pair_recal.filter{it[1].status == 'normal'}.set{normal}
    RECALIBRATE_out.pair_recal.filter{it[1].status == 'tumor'}.set{tumor}
    normal.cross(tumor){it[0]}.map{
        normal, tumor ->
        [patient:normal[0], type:normal[1].type,normal:normal[2],tumor:tumor[2]]
    }.set{ RECALIBRATE_out_MAP }

    //5,6,7,8,9
    SAGE(RECALIBRATE_out_MAP)
    STRELKA(RECALIBRATE_out_MAP)
    MuSE2(RECALIBRATE_out_MAP)
    CONPAIR(RECALIBRATE_out_MAP)
    MOSDEPTH(RECALIBRATE_out_MAP)

    SAVE_CSV_SAGE(SAGE.out.SAGE_out,params.report_dir, params.SAGE_dir)
    SAVE_CSV_STRELKA(STRELKA.out.STRELKA_out,params.report_dir, params.strelka_dir)
    SAVE_CSV_MuSE2(MuSE2.out.MuSE2_out,params.report_dir, params.muse2_dir)
    SAVE_CSV_CONPAIR(CONPAIR.out.CONPAIR_out,params.report_dir, params.conpair_dir)
    SAVE_CSV_MOSDEPTH(MOSDEPTH.out.MOSDEPTH_out,params.report_dir, params.mosdepth_dir)

    //10
    GETpileUP(RECALIBRATE_out.GETpileUP_input, chunk).set { GETpileUP_out }
    MUTECT2_CALLING(RECALIBRATE_out_MAP, chunk).set { MUTECT2_CALLING_out }
    MUTECT2_CALLING_out.LearnReadOrientationModel_input.groupTuple(by:0).set { LearnReadOrientationModel_input_pair }
    MUTECT2_CALLING_out.MergeMutectStats_input.groupTuple(by:0).set { MergeMutectStats_input_pair }
    MUTECT2_CALLING_out.MergeVcfs_input.groupTuple(by:0).set { MergeVcfs_input_pair }
    LearnReadOrientationModel(LearnReadOrientationModel_input_pair).set { LearnReadOrientationModel_out }
    MergeMutectStats(MergeMutectStats_input_pair).set { MergeMutectStats_out }
    MergeVcfs(MergeVcfs_input_pair).set { MergeVcfs_out }

    GETpileUP_out.GETpileUP_Merge_input.map{patient, status, table, chunk -> [patient,status, [table, chunk]]}.groupTuple(by:[0,1],sort: {it[1]}).set { GETpileUP_out_pair }
    GETpileUP_out_pair.map{
    [patient:it[0], status:it[1], mix:it[2]]
    }.set{ GETpileUP_out_pair_MAP }

    GETpileUP_Merge(GETpileUP_out_pair_MAP).set { GETpileUP_Merge_out }

    GETpileUP_Merge_out.CalculateContamination_input.filter{it[1]== 'normal'}.set{ GETpileUP_normal }
    GETpileUP_Merge_out.CalculateContamination_input.filter{it[1]== 'tumor'}.set{ GETpileUP_tumor }
    GETpileUP_normal.cross(GETpileUP_tumor){it[0]}.map{
        normal, tumor ->
        [patient:normal[0], normal:normal[2],tumor:tumor[2]]
    }.set{ GETpileUP_out_MAP }

    CalculateContamination(GETpileUP_out_MAP).set { CalculateContamination_out }
    FilterMutectCalls(MergeVcfs_out.MUTECT2_vcf, CalculateContamination_out.MUTECT2_contamination_table, CalculateContamination_out.MUTECT2_segments_table, LearnReadOrientationModel_out.MUTECT2_read_orientation, MergeMutectStats_out.MUTECT2_stats)
    SAVE_CSV_Mutect2(FilterMutectCalls.out.Mutect2_out,params.report_dir,params.MUTECT2_dir)

    SAVE_CSV_MOSDEPTH.out.concat(SAVE_CSV_CONPAIR.out, SAVE_CSV_Mutect2.out, SAVE_CSV_MuSE2.out, SAVE_CSV_SAGE.out, SAVE_CSV_STRELKA.out, SAVE_CSV_RECAL.out, SAVE_CSV_FASTQC.out, SAVE_CSV_MKDUP.out, SAVE_CSV_BWA_MEM.out).collect().set{saved_csv}
    SUMMARY(saved_csv).view()
    
}
