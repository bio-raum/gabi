/*
include Modules
*/
include { BWAMEM2_INDEX }   from './../../modules/bwamem2/index'
include { BWAMEM2_MEM }     from './../../modules/bwamem2/mem'

workflow ALIGN_SHORT_READS {

    take:
    ch_assembly_with_reads

    main:

    ch_bam = channel.from([])
    ch_versions = channel.from([])

    // Index the assembly
    BWAMEM2_INDEX(
        ch_assembly_with_reads.map { m,r,a ->
            [ m,a]
        }
    )
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)

    // Join reads with index
    ch_assembly_with_reads.map { m,r,a ->
        [ m,r ]
    }.join(
        BWAMEM2_INDEX.out.index
    ).set { ch_reads_with_index }

    // Align reads against index
    BWAMEM2_MEM(
        ch_reads_with_index
    )
    ch_versions = ch_versions.mix(BWAMEM2_MEM.out.versions)
    ch_bam      = ch_bam.mix(BWAMEM2_MEM.out.bam)

    emit:
    versions    = ch_versions
    bam         = ch_bam    

}