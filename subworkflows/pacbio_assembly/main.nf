include { FLYE as FLYE_PACBIO }     from '../../modules/flye'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_PAF } from '../../modules/minimap2/align'
include { KMC }                     from '../../modules/kmc'
include { RACON }                   from '../../modules/racon'

ch_versions = Channel.from([])

workflow PACBIO_ASSEMBLY {

    take:
    reads // [ meta, hifi_reads ]

    main:

    // FLYE long read assembler
    FLYE_PACBIO(
        reads
    )
    ch_versions = ch_versions.mix(FLYE_PACBIO.out.versions)

    ch_flye_with_reads = reads.join(FLYE_PACBIO.out.fasta)

    // Align long reads to initial FLYE assembly
    MINIMAP2_ALIGN_PAF(
        ch_flye_with_reads,
        "paf"
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_PAF.out.versions)

    ch_flye_with_alignment = ch_flye_with_reads.join(MINIMAP2_ALIGN_PAF.out.paf)
    // Use read alignments to polish FLYE assembly
    RACON(
        ch_flye_with_alignment
    )
    ch_versions = ch_versions.mix(RACON.out.versions)

    emit:
    assembly = RACON.out.improved_assembly
    versions = ch_versions

}
