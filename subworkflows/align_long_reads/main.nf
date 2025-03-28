include { MINIMAP2_ALIGN }    from './../../modules/minimap2/align'

workflow ALIGN_LONG_READS {

    take:
    ch_reads_with_assembly

    main:

    ch_versions = Channel.from([])
    ch_bam      = Channel.from([])

    // Align reads with minimap
    MINIMAP2_ALIGN(
        ch_reads_with_assembly,
        "bam"
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions)
    ch_bam      = ch_bam.mix(MINIMAP2_ALIGN.out.bam)

    emit:
    versions = ch_versions
    bam = ch_bam
}