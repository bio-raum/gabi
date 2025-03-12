
/*
Import modules
*/
include { MOSDEPTH }            from './../../modules/mosdepth'
include { SAMTOOLS_MERGE }      from './../../modules/samtools/merge'
include { SAMTOOLS_INDEX }      from './../../modules/samtools/index'
include { SAMTOOLS_STATS }      from './../../modules/samtools/stats'

/*
Import Subworkflows
*/
include { ALIGN_SHORT_READS }   from './../align_short_reads'
include { ALIGN_LONG_READS }    from './../align_long_reads'

workflow COVERAGE {

    take:
    assembly
    short_reads
    ont_reads 
    pacbio_reads

    main:

    ch_bam      = Channel.from([])
    ch_versions = Channel.from([])
    ch_summary_by_platform = Channel.from([])

    short_reads.mix(ont_reads).mix(pacbio_reads).map { m,r ->
        [ m.sample_id,m,r]
    }.combine(
        assembly.map { m,a ->
            [ m.sample_id, a]
        }, by: 0
    ).map { i,m,r,a ->
        [ m,r,a]
    }.set { ch_reads_with_assembly }

    /*
    Align short reads using BWA
    */
    ALIGN_SHORT_READS(
        ch_reads_with_assembly.filter{ m,r,a -> 
            m.platform == "ILLUMINA"
        }
    )
    ch_bam = ch_bam.mix(ALIGN_SHORT_READS.out.bam)

    /*
    Align long reads using Minimap2
    */
    ALIGN_LONG_READS(
        ch_reads_with_assembly.filter{ m,r,a -> 
            m.platform == "NANOPORE" || m.platform == "PACBIO"
        }
    )
    ch_bam = ch_bam.mix(ALIGN_LONG_READS.out.bam)

    /*
    We need both platform-specific as well as 
    sample-level coverage reports. So we combine all bams
    for one sample to compute a global coverage in addition
    to the platform-level reports
    */
    bam_mapped = ch_bam.map { meta, bam ->
        def new_meta = [:]
        new_meta.sample_id = meta.sample_id
        def groupKey = meta.sample_id
        tuple( groupKey, new_meta, bam)
    }.groupTuple(by: [0,1]).map { g ,new_meta ,bam -> [ new_meta, bam ] }
            
    bam_mapped.branch {
        single:   it[1].size() == 1
        multiple: it[1].size() > 1
    }.set { bam_to_merge }

    /*
    Merge BAM files across technologies
    */
    SAMTOOLS_MERGE(
        bam_to_merge.multiple
    )

    bam_to_merge.single.mix(SAMTOOLS_MERGE.out.bam).map { m,b -> 
        m.platform = "ALL"
        tuple(m,b)
    }.set { ch_bam_all }

    /*
    Index the BAM files
    */
    SAMTOOLS_INDEX(
        ch_bam.mix(ch_bam_all)
    )

    /*
    Calculate coverage
    */
    MOSDEPTH(
        SAMTOOLS_INDEX.out.bam
    )

    /*
    Compute BAM stats for Illumina reads
    */
    SAMTOOLS_STATS(
        SAMTOOLS_INDEX.out.bam.filter { m,b,i -> m.platform == "ILLUMINA"}
    )

    MOSDEPTH.out.global_txt.branch { m,r ->
        illumina: m.platform == "ILLUMINA"
        ont: m.platform == "NANOPORE"
        pacbio: m.platform == "PACBIO"
    }.set { reports_by_platform }

    emit:
    versions    = ch_versions
    report      = MOSDEPTH.out.global_txt
    report_illumina = reports_by_platform.illumina
    report_ont  = reports_by_platform.ont
    report_pacbio = reports_by_platform.pacbio
    bam         = SAMTOOLS_INDEX.out.bam
    bam_stats   = SAMTOOLS_STATS.out.stats
    summary     = MOSDEPTH.out.summary_txt
    summary_by_platform = ch_summary_by_platform

}