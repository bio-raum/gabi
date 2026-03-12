include { FLYE as FLYE_PACBIO }             from '../../modules/flye'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_PAF } from '../../modules/minimap2/align'
include { KMC }                             from '../../modules/kmc'
include { DNAAPLER }                        from '../../modules/dnaapler'
include { POLYPOLISH_POLISH }               from '../../modules/polypolish/polish'
include { BWAMEM2_INDEX as BWAMEM2_INDEX_POLYPOLISH } from '../../modules/bwamem2/index'
include { HOMOPOLISH as HOMOPOLISH_PACBIO } from '../../modules/homopolish'
include { BWAMEM2_MEM_POLYPOLISH }          from '../../modules/bwamem2/mem_polypolish'
// include { AUTOCYCLER_WORKFLOW }          from '../../modules/autocycler/full'
include { AUTOCYCLER_WORKFLOW }             from './../autocycler_workflow'

workflow PACBIO_ASSEMBLY {

    take:
    reads // [ meta, illumina, pacbio ] where illumina reads are optional
    homopolish_db

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

    if (params.autocycler) {
        read_type = params.pacbio_hifi ? "pacbio_hifi" : "pacbio_clr"
        // Autocycler bash workflow
        AUTOCYCLER_WORKFLOW(
            lreads,
            read_type
        )
        ch_versions = ch_versions.mix(AUTOCYCLER_WORKFLOW.out.versions)
        ch_long_read_assembly = AUTOCYCLER_WORKFLOW.out.fasta
    } else {
        // FLYE long read assembler
        FLYE_PACBIO(
            lreads
        )
        ch_versions = ch_versions.mix(FLYE_PACBIO.out.versions)
        ch_long_read_assembly = FLYE_PACBIO.out.fasta
    }

    ch_long_read_assembly.map { m, a ->
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

    // Run homopolish only on CLR reads
    if (params.homopolish & !params.pacbio_hifi) {
        HOMOPOLISH_PACBIO(
            assembly_with_short_reads.without.map { m,a,r ->
                tuple(m,a)
            },
            homopolish_db
        )
        ch_versions = ch_versions.mix(HOMOPOLISH_PACBIO.out.versions)
        ch_homopolished = HOMOPOLISH_PACBIO.out.polished
    } else {
        ch_homopolished = assembly_with_short_reads.without.map { m,a,r -> tuple(m,a) }
    }

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

    ch_polished_assemblies = ch_homopolished.mix(POLYPOLISH_POLISH.out.fasta)

    // Consistently orient chromosomes
    DNAAPLER(
        ch_polished_assemblies
    )
    ch_versions = ch_versions.mix(DNAAPLER.out.versions)

    emit:
    assembly = DNAAPLER.out.fasta
    versions = ch_versions

}
