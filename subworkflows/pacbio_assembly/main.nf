include { FLYE as FLYE_PACBIO }     from '../../modules/flye'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_PAF } from '../../modules/minimap2/align'
include { KMC }                     from '../../modules/kmc'
include { DNAAPLER }                from '../../modules/dnaapler'
include { POLYPOLISH_POLISH }       from '../../modules/polypolish/polish'
include { BWAMEM2_INDEX as BWAMEM2_INDEX_POLYPOLISH } from '../../modules/bwamem2/index'
include { BWAMEM2_MEM_POLYPOLISH }  from '../../modules/bwamem2/mem_polypolish'

workflow PACBIO_ASSEMBLY {

    take:
    reads // [ meta, hifi_reads ]

    main:

    ch_versions = channel.from([])

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
    FLYE_PACBIO(
        lreads
    )
    ch_versions = ch_versions.mix(FLYE_PACBIO.out.versions)

    FLYE_PACBIO.out.fasta.map { m, a ->
        tuple(m.sample_id, m, a)
    }.join(
        sreads.map {m,r ->
            tuple(m.sample_id,r)
        }, remainder: true
    ).map { key, m, p, r ->
        tuple(m,p,r)
    }.branch {
        with: it.last()
        without: !it.last()
    }.set { assembly_with_short_reads }

    // Create BWA index
    BWAMEM2_INDEX_POLYPOLISH(
        assembly_with_short_reads.with.map { m,a,r -> 
            tuple(m,a)
        }
    )
    ch_versions = ch_versions.mix(BWAMEM2_INDEX_POLYPOLISH.out.versions)

    // Align short reads and create one SAM file per mate
    BWAMEM2_MEM_POLYPOLISH(
        assembly_with_short_reads.with.map { m,a,r ->
            tuple(m,r)
        }.join(
            BWAMEM2_INDEX_POLYPOLISH.out.index
        )
    )
    ch_versions = ch_versions.mix(BWAMEM2_MEM_POLYPOLISH.out.versions)

    // Run Polypolish if short reads are available
    POLYPOLISH_POLISH(
        assembly_with_short_reads.with.map { m,a,r ->
            tuple(m,a)
        }.join(
            BWAMEM2_MEM_POLYPOLISH.out.sam
        )
    )
    ch_versions = ch_versions.mix(POLYPOLISH_POLISH.out.versions)

    ch_polished_assemblies = assembly_with_short_reads.without.map { m,a,s -> [m, a] }.mix(POLYPOLISH_POLISH.out.fasta)

    // Consistently orient chromosomes
    DNAAPLER(
        ch_polished_assemblies
    )
    ch_versions = ch_versions.mix(DNAAPLER.out.versions)

    emit:
    assembly = DNAAPLER.out.fasta
    versions = ch_versions

}
