include { MOBSUITE_RECON }  from './../../modules/mobsuite/recon'
include { ABRICATE_RUN as ABRICATE_PLASMIDFINDER }  from './../../modules/abricate/run'

workflow PLASMIDS {
    take:
    ch_assemblies

    main:

    ch_versions = channel.from([])

    MOBSUITE_RECON(
        ch_assemblies
    )
    ch_versions = ch_versions.mix(MOBSUITE_RECON.out.versions)

    ABRICATE_PLASMIDFINDER(
        ch_assemblies.map { m,a -> tuple(m,a, "plasmidfinder")}
    )
    ch_versions = ch_versions.mix(ABRICATE_PLASMIDFINDER.out.versions)

    emit:
    reports = MOBSUITE_RECON.out.mobtyper_results
    plasmids = MOBSUITE_RECON.out.plasmids
    chromosome = MOBSUITE_RECON.out.chromosome
    versions = ch_versions
}
