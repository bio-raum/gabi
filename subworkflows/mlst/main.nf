include { MLST }                            from './../../modules/mlst'

workflow MLST_TYPING {
    take:
    assembly

    main:

    ch_versions = Channel.from([])

    assembly.branch { m, a ->
        annotated: m.taxon != 'unknown'
        unknown: m.taxon == 'unknown'
    }.set { ch_assembly_filtered }

    /* ----------------------------------------
    RUN ALL THE TOOLS
    ------------------------------------------- */

    /*
    Run Thorsten Seemanns MLST tool with the built-in best-match database
    As more than one database may belong to a given label, all databes will be run
    */
    MLST(
        ch_assembly_filtered.annotated
    )
    ch_versions = ch_versions.mix(MLST.out.versions)

    emit:
    versions = ch_versions
    report = MLST.out.json
}
