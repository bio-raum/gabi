
include { KRAKEN2_DOWNLOAD }                                from './../modules/kraken2/download'
include { CONFINDR_INSTALL  }                               from './../modules/helper/confindr_install'
include { BUSCO_DOWNLOAD as BUSCO_INSTALL }                 from './../modules/busco/download'
include { AMRFINDERPLUS_UPDATE as AMRFINDERPLUS_INSTALL }   from './../modules/amrfinderplus/update'
include { STAGE_FILE as DOWNLOAD_SOURMASH_DB }              from './../modules/helper/stage_file'
include { STAGE_FILE as DOWNLOAD_SOURMASH_NR_DB }           from './../modules/helper/stage_file'
include { GUNZIP as GUNZIP_GENOME }                         from './../modules/gunzip'
include { GUNZIP as GUNZIP_HOMOPOLISH_DB }                  from './../modules/gunzip'
include { BIOBLOOM_MAKER }                                  from './../modules/biobloom/maker'
include { UNTAR as UNTAR_TAXDUMP }                          from './../modules/untar'
include { CHECKM2_DATABASEDOWNLOAD }                        from './../modules/checkm2/databasedownload'


workflow BUILD_REFERENCES {
    main:

    kraken_db_url       = Channel.fromPath(params.references['kraken2'].url)
    confindr_db_url     = Channel.fromPath(params.references['confindr'].url)
    sourmash_db_url     = params.references['sourmashdb'].url
    sourmash_nr_db_url  = params.references['sourmashdb_nr'].url
    taxdb_url           = params.references['taxdb'].url
    homopolish_db       = Channel.fromPath(file(params.references['homopolish_db'].url)).map { f -> [ [target: 'Homopolish'], f] }
    ch_busco_lineage    = Channel.from(['bacteria_odb10'])
    host_genome         = Channel.fromPath(file(params.references['host_genome'].url)).map { f -> [ [target: 'Host'], f] }

    /* 
    Decompress taxdump
    */
    UNTAR_TAXDUMP(
        taxdb_url
    )

    /*
    Download the checkM database
    */
    CHECKM2_DATABASEDOWNLOAD(14897628)
    
    /*
    Decompress homopolish database
    */
    GUNZIP_HOMOPOLISH_DB(
        homopolish_db
    )

    /*
    Download Horse genome from EnsEMBL and build index
    */
    GUNZIP_GENOME(
        host_genome
    )

    BIOBLOOM_MAKER(
        GUNZIP_GENOME.out.gunzip.map { m,f -> f }
    )
    
    /*
    Download SourmashDB
    */
    DOWNLOAD_SOURMASH_DB(
        sourmash_db_url
    )

    DOWNLOAD_SOURMASH_NR_DB(
        sourmash_nr_db_url
    )
    /*
    Download the latest version of the AMRfinderplus DB
    This is not ideal since we cannot select specific versions -  but it works
    since we use a frozen version, and the last release of the DB for that version
    */
    AMRFINDERPLUS_INSTALL()

    /*
    Download the default Busco lineages
    */
    BUSCO_INSTALL(
        ch_busco_lineage
    )

    /*
    Download the Kraken MiniDB
    This should be good enough for our purposes
    */
    KRAKEN2_DOWNLOAD(
        kraken_db_url
    )

    /*
    Download a ConfindR database
    */
    CONFINDR_INSTALL(
        confindr_db_url
    )

    workflow.onComplete = {
        log.info 'Installation complete - deleting staged files. '
        workDir.resolve("stage-${workflow.sessionId}").deleteDir()
    }
}