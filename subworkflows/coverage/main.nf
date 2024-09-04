
/*
Import modules
*/
include { MOSDEPTH }            from './../../modules/mosdepth'
include { SAMTOOLS_MERGE }      from './../../modules/samtools/merge'
include { SAMTOOLS_INDEX }      from './../../modules/samtools/index'

/*
Import Subworkflows
*/
include { ALIGN_SHORT_READS }   from './../align_short_reads'
include { ALIGN_LONG_READS }    from './../align_long_reads'

ch_bam      = Channel.from([])
ch_qc       = Channel.from([])
ch_versions = Channel.from([])

workflow COVERAGE {

    take:
    assembly
    short_reads
    ont_reads 
    pacbio_reads

    main:

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

    // Index the BAM files
    SAMTOOLS_INDEX(
        ch_bam
    )

    // Calculate coverage
    MOSDEPTH(
        SAMTOOLS_INDEX.out.bam
    )

    emit:
    versions    = ch_versions
    report      = MOSDEPTH.out.global_txt

}