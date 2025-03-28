include { KRAKEN2_KRAKEN2 }     from './../../modules/kraken2/kraken2'
include { BRACKEN_BRACKEN }     from './../../modules/bracken/bracken'

workflow TAXONOMY_PROFILING {
    take:
    reads
    kraken2_db

    main:

    ch_versions = Channel.from([])

    // Kraken2 raw estimates
    KRAKEN2_KRAKEN2(
        reads,
        kraken2_db,
        false,
        false
    )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2.out.versions)

    // Bracken deconvoluted estimates
    BRACKEN_BRACKEN(
        KRAKEN2_KRAKEN2.out.report,
        kraken2_db
    )
    ch_versions = ch_versions.mix(BRACKEN_BRACKEN.out.versions)

    emit:
    report = BRACKEN_BRACKEN.out.reports
    report_txt = BRACKEN_BRACKEN.out.txt
    versions = ch_versions
}

