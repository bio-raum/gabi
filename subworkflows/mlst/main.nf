include { CHEWBBACA_ALLELECALL }            from './../../modules/chewbbaca/allelecall'
include { CHEWBBACA_ALLELECALL as CHEWBBACA_ALLELECALL_SINGLE }            from './../../modules/chewbbaca/allelecall'
include { CHEWBBACA_JOINPROFILES }          from './../../modules/chewbbaca/joinprofiles'
include { CHEWBBACA_ALLELECALLEVALUATOR }   from './../../modules/chewbbaca/allelecallevaluator'
include { CHEWBBACA_REMOVEGENES }           from './../../modules/chewbbaca/removegenes'
include { MLST }                            from './../../modules/mlst'

ch_versions = Channel.from([])

workflow MLST_TYPING {
    take:
    assembly

    main:

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

    /*
    We use the previously attempted taxonomic classification
    to choose the appropriate Chewbbaca cgMLST schema, if any
    Assemblies are grouped by taxon to create a multi-sample
    call matrix per species
    */

    ch_assembly_filtered.annotated.map { m, a ->
        def (genus,species) = m.taxon.toLowerCase().split(' ')
        def chewie_db = null
        if (params.chewbbaca[genus]) {
            chewie_db = params.chewbbaca[genus]
            m.db_name = genus
        } else if (params.chewbbaca["${genus}_${species}"]) {
            chewie_db = params.chewbbaca["${genus}_${species}"]
            m.db_name = "${genus}_${species}"
        } else {
            chewie_db = null
        }
        tuple(m, a, chewie_db)
    }.branch { m, a, db ->
        fail: db == null
        pass: db
    }.set { assembly_with_chewie_db }

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

    if (!params.skip_cgmlst) {
        /*
        Inform users about to-be-skipped samples due to a lack of a matching cgMLST database
        */
       
        assembly_with_chewie_db.fail.subscribe { m, s, d ->
            log.warn "${m.sample_id} - could not match a Chewbbaca cgMLST database to ${m.taxon}."
        }

        /*
        Perform cgMLST calling with Chewbbaca
        Part one consists of a joint allele calling approach in which all samples belonging to the same species are jointly call
        In addition, each sample is called invidivually to support downstream analysis of samples from across runs
        */
        CHEWBBACA_ALLELECALL_SINGLE(
            assembly_with_chewie_db.pass
        )
        ch_versions = ch_versions.mix(CHEWBBACA_ALLELECALL_SINGLE.out.versions)

        /*
        Single samples: filter wgMLST down to cgMLST, if we have a list of cgMLST loci
        */
        CHEWBBACA_ALLELECALL_SINGLE.out.profile.map { m, a ->
            def (genus,species) = m.taxon.toLowerCase().split(' ')
            def filter_list = null
            if (params.chewbbaca_filters[genus]) {
                filter_list = params.chewbbaca_filters[genus]
            } else if (params.chewbbaca_filters["${genus}_${species}"]) {
                filter_list = params.chewbbaca_filters["${genus}_${species}"]
            }
            tuple(m, a, filter_list)
        }.branch { m, a, f ->
            fail: f == null
            pass: f
        }.set { chewie_report_with_filter }

        /*
        Remove genes that are not part of the cgMLST scheme, if any
        */
        CHEWBBACA_REMOVEGENES(
            chewie_report_with_filter.pass
        )
    
        ch_single_calls = chewie_report_with_filter.fail.mix(CHEWBBACA_REMOVEGENES.out.profile)

        /* Join assemblies and databases to generate
        [ meta, [ assemblies ], db ] and filter out all
        cases where # assemblies is < 3 (no point to compute relationships)
        */
        assembly_with_chewie_db.pass.map { m, a, d ->
            def meta = [:]
            meta.sample_id = m.db_name
            meta.db_name = m.db_name
            tuple(meta, a, d)
        }.groupTuple(by: [0, 2])
        .set { ch_assemblies_chewie_grouped }

        ch_assemblies_chewie_grouped.filter { m, a, d -> (a.size() >= 3) }
        .set { ch_assemblies_chewie_call }

        CHEWBBACA_ALLELECALL(
            ch_assemblies_chewie_call
        )
        ch_versions = ch_versions.mix(CHEWBBACA_ALLELECALL.out.versions)
        CHEWBBACA_ALLELECALLEVALUATOR(
            CHEWBBACA_ALLELECALL.out.report.map { m, r ->
                def chewie_db = params.chewbbaca[m.db_name]
                tuple(m, r, chewie_db)
            }
        )
        ch_versions = ch_versions.mix(CHEWBBACA_ALLELECALLEVALUATOR.out.versions)
        }

    emit:
    versions = ch_versions
    report = MLST.out.json
}
