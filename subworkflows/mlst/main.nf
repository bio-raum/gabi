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

       /*
    We use the previously attempted taxonomic classification to
    choose the appropriate MLST schema(s), if any
    */
    ch_assembly_filtered.annotated.map { m, a ->
        def (genus,species) = m.taxon.toLowerCase().split(' ')
        def dbs = null
        if (params.mlst["${genus}_${species}"]) {
            dbs = params.mlst["${genus}_${species}"]
        } else if (params.mlst[genus]) {
            dbs = params.mlst[genus]
        } else {
            dbs = [ null ]
        }
        tuple(m, a, dbs)
    }.flatMap { m,a,dbs ->
        dbs.collect { [ m,a,it ]}
    }.branch { m, a, db ->
        fail: db == null
        pass: db
    }.set { assembly_with_mlst_db }


    /* ----------------------------------------
    RUN ALL THE TOOLS
    ------------------------------------------- */

    /*
    Run Thorsten Seemanns MLST tool with the built-in best-match database
    As more than one database may belong to a given label, all databes will be run
    */
    MLST(
        assembly_with_mlst_db.pass
    )
    ch_versions = ch_versions.mix(MLST.out.versions)

    emit:
    versions = ch_versions
    report = MLST.out.json
}
