include { FLYE }            from '../../modules/flye'
include { DORADO_POLISH }   from '../../modules/dorado/polish' 

ch_versions = Channel.from([])

workflow ONT_ASSEMBLY {

    take:
    lreads
    sreads

    main:

    FLYE(
        reads
    )
    ch_versions = ch_versions.mix(FLYE.out.versions)

    DORADO_POLISH(
        FLYE.out.fasta
    )
    ch_versions = ch_versions.mix(DORADO_POLISH.out.versions)

    DORADO_POLISH.out.polished.map { m,p ->
        tuple(m.sample_id,m,p)
    }.join(
        reads.map {m,r ->
            tuple(m.sample_id,m,r)
        }, remainder: true
    ).map { key, m, p, r ->
        tuple(m,p,r)
    }.set { dorado_with_short_reads }

    dorado_with_short_reads.branch {
        with: it.last()
        without: !it.last()
    }.set { polished_with_short_reads }

    POLYPOLISH(
        polished_with_short_reads.with
    )

    ch_polished_assembly = polished_with_short_reads.without.mix(POLYPOLISH.out.polished)
    // if short reads available


    emit:
    assembly = ch_polished_assembly

}