include { GABI_SUMMARY }    from './../../modules/helper/gabi_summary'
include { GABI_QC }         from './../../modules/helper/gabi_qc'
include { GABI_REPORT }     from './../../modules/helper/gabi_report'

ch_versions = Channel.from([])

workflow REPORT {
    take:
    reports
    template
    refs
    yml

    main:

    GABI_SUMMARY(
        reports,
        yml.collect()
    )
    ch_versions = ch_versions.mix(GABI_SUMMARY.out.versions)

    GABI_QC(
        GABI_SUMMARY.out.json,
        refs.collect()
    )
    ch_versions = ch_versions.mix(GABI_QC.out.versions)

    GABI_REPORT(
        GABI_QC.out.json.collect(),
        template,
        yml
    )
    ch_versions = ch_versions.mix(GABI_REPORT.out.versions)

    emit:
    html = GABI_REPORT.out.html
    json = GABI_QC.out.json
    versions = ch_versions
}
