/*
Include Modules
*/
include { AMRFINDERPLUS_RUN }               from './../../modules/amrfinderplus/run'
include { AMRFINDERPLUS_UPDATE }            from './../../modules/amrfinderplus/update'
include { HAMRONIZATION_AMRFINDERPLUS }     from './../../modules/hamronization/amrfinderplus'
include { HAMRONIZATION_ABRICATE }          from './../../modules/hamronization/abricate'
include { HAMRONIZATION_SUMMARIZE }         from './../../modules/hamronization/summarize'
include { HAMRONIZATION_SUMMARIZE as HAMRONIZATION_SUMMARIZE_HTML } from './../../modules/hamronization/summarize'
include { ABRICATE_RUN }                    from './../../modules/abricate/run'
include { ABRICATE_RUN as ABRICATE_RUN_ECOLI_VIRULENCE } from './../../modules/abricate/run'

workflow AMR_PROFILING {
    take:
    assembly
    db              // The AMRfinder database to run
    abricate_dbs    // A list of abricate databases to run (should be generic!)

    main:

    ch_versions = channel.from([])
    multiqc_files = channel.from([])
    ch_hamronization_input = channel.from([])

    assembly.branch { m, a ->
        ecoli: m.taxon ==~ /^Escherichia.*/
        salmonella: m.taxon ==~ /^Salmonella.*/
        listeria: m.taxon ==~ /^Listeria.*/
        campylobacter: m.taxon ==~ /^Campylobacter.*/
        bacillus: m.taxon ==~ /^Bacillus.*/
    }.set { assembly_by_taxon }

    /*
    Run AMRFinderPlus and make JSON report
    */

    // if no local DB is defined, we download it on the flye
    if (!params.reference_base) {
        AMRFINDERPLUS_UPDATE()
        ch_amrfinderplus_db = AMRFINDERPLUS_UPDATE.out.db
        ch_versions = ch_versions.mix(AMRFINDERPLUS_UPDATE.out.versions)
    } else {
        ch_amrfinderplus_db = channel.from(db)
    }

    assembly.map { m,a ->
        def organism = create_organism(m.taxon)
        if (organism) {
            m.organism = organism
        }
        tuple(m,a)
    }.set  { assembly_with_organism }

    AMRFINDERPLUS_RUN(
        assembly_with_organism,
        ch_amrfinderplus_db.collect()
    )
    ch_versions = ch_versions.mix(AMRFINDERPLUS_RUN.out.versions)

    HAMRONIZATION_AMRFINDERPLUS(
        AMRFINDERPLUS_RUN.out.report,
        'json',
        AMRFINDERPLUS_RUN.out.tool_version,
        AMRFINDERPLUS_RUN.out.db_version
    )
    ch_hamronization_input = ch_hamronization_input.mix(HAMRONIZATION_AMRFINDERPLUS.out.json)
    ch_versions = HAMRONIZATION_AMRFINDERPLUS.out.versions

    /*
    Run Abricate and make JSON report
    */

    assembly_with_db = assembly.combine(abricate_dbs)
    ABRICATE_RUN(
        assembly_with_db
    )
    ch_versions = ch_versions.mix(ABRICATE_RUN.out.versions)

    /*  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Taxon-specific abricate analyses
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~  */ 
    // E. coli - here we use a specific database!
    ABRICATE_RUN_ECOLI_VIRULENCE(
        assembly_by_taxon.ecoli.map { m,a -> [ m, a, 'ecoli_vf']}
    )
    ch_versions = ch_versions.mix(ABRICATE_RUN_ECOLI_VIRULENCE.out.versions)

    // Join basic Abricate results
    HAMRONIZATION_ABRICATE(
        ABRICATE_RUN.out.report.mix(ABRICATE_RUN_ECOLI_VIRULENCE.out.report).groupTuple(),
        'json',
        '1.0.1',
        '2021-Mar-27'
    )
    ch_versions = ch_versions.mix(HAMRONIZATION_ABRICATE.out.versions)
    ch_hamronization_input = ch_hamronization_input.mix(HAMRONIZATION_ABRICATE.out.json)

    /*
    Summarize reports across tools
    */
    HAMRONIZATION_SUMMARIZE(
        ch_hamronization_input.map { m, j -> j }.collect(),
        params.arg_hamronization_summarizeformat
    )
    ch_versions = ch_versions.mix(HAMRONIZATION_SUMMARIZE.out.versions)

    HAMRONIZATION_SUMMARIZE_HTML(
        ch_hamronization_input.map { m, j -> j }.collect(),
        "interactive"
    )

    emit:
    report = HAMRONIZATION_SUMMARIZE.out.json
    amrfinder_report = AMRFINDERPLUS_RUN.out.report
    abricate_report = ABRICATE_RUN.out.report
    versions = ch_versions
    qc = multiqc_files
}

// Map taxon to AMRfinder compatible organism name, if any
def create_organism(taxon) {

    def taxa_for_join = [ "Acinetobacter baumannii", 
    "Vibrio cholerae", "Klebsiella pneumoniae", "Burkholderia cepacia",
    "Burkholderia pseudomallei", "Citrobacter freundii", "Clostridioides difficile",
    "Enterobacter asburiae", "Enterobacter cloacae", "Enterococcus faecalis", 
    "Enterococcus faecium", "Klebsiella oxytoca", "Neisseria gonorrhoeae", 
    "Neisseria meningitidis", "Pseudomonas aeruginosa", "Serratia marcescens",
    "Staphylococcus aureus", "Staphylococcus pseudintermedius", "Streptococcus agalactiae",
    "Streptococcus pneumoniae", "Streptococcus pyogenes", "Vibrio parahaemolyticus",
    "Vibrio vulnificus"
    ]

    def organism = null

    if (taxon.contains("Escherichia")) {
        organism = "Escherichia"
    } else if (taxon.contains("Shigella")) {
        organism = "Salmonella"
    } else if (taxon.contains("Salmonella")) {
        organism = "Salmonella"
    } else if (taxon.contains("Campylobacter")) {
        organism = "Campylobacter"
    } else if (taxa_for_join.contains(taxon)) {
        organism = taxon.split(" ").join("_")
    }

    return organism

}