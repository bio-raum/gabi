include { FLYE }                    from '../../modules/flye'
include { MEDAKA_CONSENSUS }        from '../../modules/medaka/consensus'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_PAF } from '../../modules/minimap2/align'
include { RACON }                   from '../../modules/racon'
include { POLYPOLISH_POLISH }       from '../../modules/polypolish/polish'
include { BWAMEM2_INDEX as BWAMEM2_INDEX_POLYPOLISH } from '../../modules/bwamem2/index'
include { BWAMEM2_MEM_POLYPOLISH }  from '../../modules/bwamem2/mem_polypolish'

ch_versions = Channel.from([])

workflow ONT_ASSEMBLY {

    take:
    lreads
    sreads

    main:

    // FLYE long read assembler
    FLYE(
        reads
    )
    ch_versions = ch_versions.mix(FLYE.out.versions)

    // Align long reads to initial FLYE assembly
    MINIMAP2_ALIGN_PAF(
        lreads.join(FLYE.out.fasta),
        "paf"
    )
    ch_versions = ch_versions.mix(MINIMAP2_ALIGN_PAF.out.versions)

    // Use read alignments to polish FLYE assembly
    RACON(
        lreads.join(FLYE.out.fasta).join(MINIMAP2_ALIGN_PAF.out.paf)
    )
    ch_versions = ch_versions.mix(RACON.out.versions)

    // Re-polish initial consensus contigs with Medaka
    MEDAKA_CONSENSUS(
        lreads.join(RACON.out.improved_assembly)
    )
    ch_versions = ch_versions.mix(MEDAKA_CONSENSUS.out.versions)

    // Join polished Medaka assembly with optional short reads
    MEDAKA_CONSENSUS.out.consensus.map { m,c ->
        tuple(m.sample_id,m,c)
    }.join(
        sreads.map {m,r ->
            tuple(m.sample_id,r)
        }, remainder: true
    ).map { key, m, p, r ->
        tuple(m,p,r)
    }.set { medaka_with_short_reads }

    medaka_with_short_reads.branch {
        with: it.last()
        without: !it.last()
    }.set { polished_with_short_reads }

    // Create BWA index
    BWAMEM2_INDEX(
        polished_with_short_reads.map { m,a,r -> 
            tuple(m,a)
        }
    )
    ch_versions = ch_versions.mix(BWAMEM2_INDEX.out.versions)

    // Align short reads and create one SAM file per direction
    BWAMEM2_MEM_POLYPOLISH(
        polished_with_short_reads.map { m,a,r ->
            tuple(m,r)
        }.join(
            BWAMEM2_INDEX.out.index
        )
    )
    ch_versions = ch_versions.mix(BWAMEM2_MEM_POLYPOLISH.out.versions)

    // Run Polypolish if short reads are available
    POLYPOLISH_POLISH(
        polished_with_short_reads.map { m,a,r ->
            tuple(m,a)
        }.join(
            BWAMEM2_MEM_POLYPOLISH.out.sam
        )
    )
    ch_versions = ch_versions.mix(POLYPOLISH_POLISH.out.versions)

    ch_polished_assembly = polished_with_short_reads.without.mix(POLYPOLISH_POLISH.out.fasta)

    emit:
    assembly = ch_polished_assembly
    versions = ch_versions

}