include { FLYE as FLYE_ONT }        from '../../modules/flye'
include { MEDAKA_CONSENSUS }        from '../../modules/medaka/consensus'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_PAF } from '../../modules/minimap2/align'
include { HOMOPOLISH }              from '../../modules/homopolish'
include { DNAAPLER }                from '../../modules/dnaapler'
include { POLYPOLISH_POLISH }       from '../../modules/polypolish/polish'
include { BWAMEM2_INDEX as BWAMEM2_INDEX_POLYPOLISH } from '../../modules/bwamem2/index'
include { BWAMEM2_MEM_POLYPOLISH }  from '../../modules/bwamem2/mem_polypolish'

/* 
This workflow is inspired by https://github.com/rpetit3/dragonflye
Since Dragonflye isn't regularly maintained, GABI re-implements the
basic (slightly simplified) logic into a subworkflow instead
*/
workflow ONT_ASSEMBLY {

    take:
    reads // [ meta, [short_reads], ont_reads ]
    homopolish_db

    main:

    ch_versions = Channel.from([])

    // Get long reads
    reads.map { m,s,o ->
        tuple(m,o)
    }.set { lreads }

    // Get short reads if they exist 
    reads.map { m,s,o ->
        tuple(m,s)
    }.filter { it.last() }
    .set { sreads }

    // FLYE long read assembler
    FLYE_ONT(
        lreads
    )
    ch_versions = ch_versions.mix(FLYE_ONT.out.versions)

    // Re-polish initial consensus contigs with Medaka
    MEDAKA_CONSENSUS(
        lreads.join(FLYE_ONT.out.fasta)
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
    }.branch {
        with: it.last()
        without: !it.last()
    }.set { polished_with_short_reads }

    // Homopolish to remove homopolymer errors when no short reads
    // are available ; skippable if users chooses to  
    if (!params.skip_homopolish) {
            HOMOPOLISH(
            polished_with_short_reads.without.map { m,p,r ->
                tuple(m,p)
            },
            homopolish_db
        )
        ch_versions = ch_versions.mix(HOMOPOLISH.out.versions)
        ch_homopolished = HOMOPOLISH.out.polished
    } else {
        ch_homopolished = MEDAKA_CONSENSUS.out.consensus
    }

    // Create BWA index
    BWAMEM2_INDEX_POLYPOLISH(
        polished_with_short_reads.with.map { m,a,r -> 
            tuple(m,a)
        }
    )
    ch_versions = ch_versions.mix(BWAMEM2_INDEX_POLYPOLISH.out.versions)

    // Align short reads and create one SAM file per mate
    BWAMEM2_MEM_POLYPOLISH(
        polished_with_short_reads.with.map { m,a,r ->
            tuple(m,r)
        }.join(
            BWAMEM2_INDEX_POLYPOLISH.out.index
        )
    )
    ch_versions = ch_versions.mix(BWAMEM2_MEM_POLYPOLISH.out.versions)

    // Run Polypolish if short reads are available
    POLYPOLISH_POLISH(
        polished_with_short_reads.with.map { m,a,r ->
            tuple(m,a)
        }.join(
            BWAMEM2_MEM_POLYPOLISH.out.sam
        )
    )
    ch_versions = ch_versions.mix(POLYPOLISH_POLISH.out.versions)

    // Combine shot-read polished assemblies with homopolish assemblies for which we had no short reads
    ch_polished_assembly = ch_homopolished.mix(POLYPOLISH_POLISH.out.fasta)

    // Consistently orient chromosomes
    DNAAPLER(
        ch_polished_assembly
    )
    ch_versions = ch_versions.mix(DNAAPLER.out.versions)

    emit:
    assembly = DNAAPLER.out.fasta
    versions = ch_versions

}
