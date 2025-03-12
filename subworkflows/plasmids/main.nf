include { MOBSUITE_RECON }  from './../../modules/mobsuite/recon'

workflow PLASMIDS {
    take:
    ch_assemblies

    main:

    ch_versions = Channel.from([])

    MOBSUITE_RECON(
        ch_assemblies
    )
    ch_versions = ch_versions.mix(MOBSUITE_RECON.out.versions)

    emit:
    reports = MOBSUITE_RECON.out.mobtyper_results
    plasmids = MOBSUITE_RECON.out.plasmids
    chromosome = MOBSUITE_RECON.out.chromosome
    versions = ch_versions
}
