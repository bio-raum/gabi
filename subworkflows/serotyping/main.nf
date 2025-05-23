include { ECTYPER }     from './../../modules/ectyper'
include { SEQSERO2 }    from './../../modules/seqsero2'
include { LISSERO }     from './../../modules/lissero'
include { SISTR }       from './../../modules/sistr'
include { STECFINDER }  from './../../modules/stecfinder'
include { SCCMEC }      from './../../modules/sccmec'
include { BTYPER3 }     from './../../modules/btyper3'


workflow SEROTYPING {
    take:
    assembly // [ meta, assembly ]

    main:

    ch_versions = Channel.from([])
    ch_reports = Channel.from([])

    assembly.branch { m, a ->
        ecoli: m.taxon ==~ /^Escherichia.*/
        salmonella: m.taxon ==~ /^Salmonella.*/
        listeria: m.taxon ==~ /^Listeria.*/
        staphylococcus: m.taxon ==~ /^Staphylococcus.*/
        bacillus: m.taxon ==~ /^Bacillus.*/
    }.set { assembly_by_taxon }

    /* ~~~~~~~~~~~~~
    Bacillus Btyper3
    ~~~~~~~~~~~~~~*/
    BTYPER3(
        assembly_by_taxon.bacillus
    )
    ch_versions = ch_versions.mix(BTYPER3.out.versions)
    ch_reports = ch_reports.mix(BTYPER3.out.tsv)

    /*
    Run Ectyper - Serotyping of E. coli
    */
    ECTYPER(
        assembly_by_taxon.ecoli
    )
    ch_versions = ch_versions.mix(ECTYPER.out.versions)
    ch_reports = ch_reports.mix(ECTYPER.out.tsv)

    /*
    Run Stecfinder for E. coli
    */
    STECFINDER(
        assembly_by_taxon.ecoli
    )
    ch_versions = ch_versions.mix(STECFINDER.out.versions)
    ch_reports = ch_reports.mix(STECFINDER.out.tsv)
    
    /*
    Run SeqSero2 - Serotyping for Salmonella
    */
    SEQSERO2(
        assembly_by_taxon.salmonella
    )
    ch_versions = ch_versions.mix(SEQSERO2.out.versions)
    ch_reports = ch_reports.mix(SEQSERO2.out.tsv)

    /*
    Run SISTR - Serotyping for Salmonella
    */
    SISTR(
        assembly_by_taxon.salmonella
    )
    ch_versions = ch_versions.mix(SISTR.out.versions)
    ch_reports = ch_reports.mix(SISTR.out.tsv)

    /*
    Run Sccmec for Staphylococcus
    */
    SCCMEC(
        assembly_by_taxon.staphylococcus
    )
    ch_versions = ch_versions.mix(SCCMEC.out.versions)
    ch_reports = ch_reports.mix(SCCMEC.out.tsv)

    /*
    Run LisSero - Serotyping L. monocytogenes
    */
    LISSERO(
        assembly_by_taxon.listeria
    )
    ch_versions = ch_versions.mix(LISSERO.out.versions)
    ch_reports  = ch_reports.mix(LISSERO.out.tsv)

    emit:
    versions = ch_versions
    reports = ch_reports
    }
