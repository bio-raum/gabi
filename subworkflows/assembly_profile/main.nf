include { KRAKEN2_KRAKEN2 as KRAKEN2_KRAKEN2_ASSEMBLY }     from './../../modules/kraken2/kraken2'
include { BRACKEN_BRACKEN as BRACKEN_BRACKEN_ASSEMBLY }     from './../../modules/bracken/bracken'

workflow ASSEMBLY_PROFILE {
    take:
    reads
    kraken2_db

    main:

    ch_versions = Channel.from([])

    // Kraken2 raw estimates
    KRAKEN2_KRAKEN2_ASSEMBLY(
        reads,
        kraken2_db,
        false,
        false
    )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2_ASSEMBLY.out.versions)

    // Bracken deconvoluted estimates
    BRACKEN_BRACKEN_ASSEMBLY(
        KRAKEN2_KRAKEN2_ASSEMBLY.out.report,
        kraken2_db
    )
    ch_versions = ch_versions.mix(BRACKEN_BRACKEN_ASSEMBLY.out.versions)

    emit:
    report = BRACKEN_BRACKEN_ASSEMBLY.out.reports
    report_txt = BRACKEN_BRACKEN_ASSEMBLY.out.txt
    versions = ch_versions
}

