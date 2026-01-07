include { KRAKEN2_KRAKEN2 as KRAKEN2_KRAKEN2_ASSEMBLY }     from './../../modules/kraken2/kraken2'
include { BRACKEN_BRACKEN as BRACKEN_BRACKEN_ASSEMBLY }     from './../../modules/bracken/bracken'
include { CHECKM2_PREDICT }                                 from './../../modules/checkm2/predict'
include { TAXONKIT_LINEAGE }                                from './../../modules/taxonkit/lineage'
include { TAXONKIT_REFORMAT }                               from './../../modules/taxonkit/reformat'
include { HELPER_FORMAT_TAXONKIT }                          from './../../modules/helper/format_taxonkit'

workflow ASSEMBLY_PROFILE {
    take:
    assembly
    kraken2_db
    checkm_db
    taxdump

    main:

    ch_versions = Channel.from([])
    ch_reports = Channel.from([])

    // Kraken2 raw estimates
    KRAKEN2_KRAKEN2_ASSEMBLY(
        assembly,
        kraken2_db,
        false,
        true
    )
    ch_versions = ch_versions.mix(KRAKEN2_KRAKEN2_ASSEMBLY.out.versions)

    // Run Taxonkit lineage
    TAXONKIT_LINEAGE(
        KRAKEN2_KRAKEN2_ASSEMBLY.out.classified,
        taxdump
    )

    TAXONKIT_REFORMAT(
        TAXONKIT_LINEAGE.out.tsv,
        taxdump
    )

    HELPER_FORMAT_TAXONKIT(
        KRAKEN2_KRAKEN2_ASSEMBLY.out.classified.join(TAXONKIT_REFORMAT.out.tsv)
    )
    ch_reports =  ch_reports.mix(HELPER_FORMAT_TAXONKIT.out.txt)

    // Bracken deconvoluted estimates
    BRACKEN_BRACKEN_ASSEMBLY(
        KRAKEN2_KRAKEN2_ASSEMBLY.out.report,
        kraken2_db
    )
    ch_versions = ch_versions.mix(BRACKEN_BRACKEN_ASSEMBLY.out.versions)
    ch_reports = ch_reports.mix(BRACKEN_BRACKEN_ASSEMBLY.out.reports)

    CHECKM2_PREDICT(
        assembly,
        checkm_db
    )
    ch_versions = ch_versions.mix(CHECKM2_PREDICT.out.versions)
    ch_reports = ch_reports.mix(CHECKM2_PREDICT.out.checkm2_tsv)

    emit:
    report = ch_reports
    report_txt = BRACKEN_BRACKEN_ASSEMBLY.out.txt
    versions = ch_versions
}

