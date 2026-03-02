/*
----------------------
Import Modules
----------------------
*/
include { STAGE_FILE as STAGE_SAMPLESHEET } from './../modules/helper/stage_file'
include { INPUT_CHECK }                 from './../modules/input_check'
include { MULTIQC }                     from './../modules/multiqc'
include { MULTIQC as MULTIQC_ILLUMINA } from './../modules/multiqc'
include { MULTIQC as MULTIQC_NANOPORE } from './../modules/multiqc'
include { MULTIQC as MULTIQC_PACBIO }   from './../modules/multiqc'
include { RENAME_CTG as RENAME_EXTERNAL_CTG } from './../modules/rename_ctg'
include { RENAME_CTG as RENAME_PLASMID_CTG } from './../modules/rename_ctg'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from './../modules/custom/dumpsoftwareversions'

/*
-------------------
Import Subworkflows
-------------------
*/
include { GROUP_READS }                 from './../subworkflows/group_reads'
include { QC }                          from './../subworkflows/qc'
include { AMR_PROFILING }               from './../subworkflows/amr_profiling'
include { TAXONOMY_PROFILING }          from './../subworkflows/taxonomy_profiling'
include { ASSEMBLY_PROFILE }            from './../subworkflows/assembly_profile'
include { ASSEMBLY_QC }                 from './../subworkflows/assembly_qc'
include { PLASMIDS }                    from './../subworkflows/plasmids'
include { ANNOTATE }                    from './../subworkflows/annotate'
include { MLST_TYPING }                 from './../subworkflows/mlst'
include { REPORT }                      from './../subworkflows/report'
include { FIND_REFERENCES }             from './../subworkflows/find_references'
include { SEROTYPING }                  from './../subworkflows/serotyping'
include { COVERAGE }                    from './../subworkflows/coverage'
include { VARIANTS }                    from './../subworkflows/variants'
include { ILLUMINA_ASSEMBLY }           from './../subworkflows/illumina_assembly'
include { ONT_ASSEMBLY }                from './../subworkflows/ont_assembly'
include { PACBIO_ASSEMBLY }             from './../subworkflows/pacbio_assembly'
include { ABRICATE_RUN }                from './../modules/abricate/run/main.nf'

workflow GABI {

    main:

    ch_versions     = channel.from([])
    multiqc_files   = channel.from([])
    ch_assemblies   = channel.from([])
    ch_report       = channel.from([])
    ch_multiqc_illumina = channel.from([])  
    ch_multiqc_nanopore = channel.from([])
    ch_multiqc_pacbio   = channel.from([])

    samplesheet = params.input ? channel.fromPath(file(params.input, checkIfExists:true)) : channel.value([])

    refDir = file(params.reference_base + "/gabi/${params.reference_version}")
    if (!refDir.exists()) {
        log.info 'The required reference directory was not found on your system, exiting!'
        System.exit(1)
    }

    ch_multiqc_config = params.multiqc_config   ? channel.fromPath(params.multiqc_config, checkIfExists: true).collect()    : []
    ch_multiqc_logo   = params.multiqc_logo     ? channel.fromPath(params.multiqc_logo, checkIfExists: true).collect()      : []

    ch_report_template = params.template        ? channel.fromPath(params.template, checkIfExists: true).collect()          : []
    ch_report_refs     = params.report_refs     ? channel.fromPath(params.report_refs, checkIfExists: true).collect()          : []

    ch_prokka_proteins = params.prokka_proteins ? channel.fromPath(params.prokka_proteins, checkIfExists: true).collect()   : []
    ch_prokka_prodigal = params.prokka_prodigal ? channel.fromPath(params.prokka_prodigal, checkIfExists:true).collect()    : []

    abricate_dbs    = channel.from(params.abricate_dbs)
    amrfinder_db    = params.reference_base ? file(params.references['amrfinderdb'].db, checkIfExists:true)   : []
    kraken2_db      = params.reference_base ? file(params.references['kraken2'].db, checkIfExists:true)       : []
    homopolish_db   = params.reference_base ? file(params.references['homopolish_db'].db, checkIfExists:true) : []
    checkm_db       = params.reference_base ? file(params.references['checkmdb'].db, checkIfExists: true)     : []
    taxdb           = params.reference_base ? file(params.references['taxdb'].db, checkIfExists:true)         : []

    // Sourmash DB choice - either the full thing or a smaller "nr" one to speed up searches at the cost of some precision
    if (params.fast_ref) {
        sourmashdb      = params.reference_base ? file(params.references['sourmashdb_nr'].db, checkIfExists:true)    : []
    } else {
        sourmashdb      = params.reference_base ? file(params.references['sourmashdb'].db, checkIfExists:true)    : []
    }

    busco_db_path   = params.reference_base ? file(params.references['busco'].db, checkIfExists:true)         : []
    busco_lineage   = params.busco_lineage

    confindr_db     = params.confindr_db ? params.confindr_db : file(params.references['confindr'].db, checkIfExists: true)

    ch_bloom_filter = params.reference_base ? channel.from([ file(params.references["host_genome"].db + ".bf", checkIfExists: true), file(params.references["host_genome"].db + ".txt", checkIfExists: true)]).collect() : []

    STAGE_SAMPLESHEET(samplesheet)

    INPUT_CHECK(samplesheet)

    // Check if the pre-assembled genomes are correctly named
    INPUT_CHECK.out.assemblies.branch { m, a ->
        valid: a.getBaseName() == m.sample_id
        invalid: a.getBaseName() != m.sample_id
    }.set { assembly_by_status }
    
    ch_assemblies = ch_assemblies.mix(assembly_by_status.valid)

    // If we pass existing assemblies instead of raw reads:
    // rename to sample_id if not already the case
    RENAME_EXTERNAL_CTG(
        assembly_by_status.invalid,
        'fasta'
    )
    ch_assemblies = ch_assemblies.mix(RENAME_EXTERNAL_CTG.out)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUB: Run read trimming and contamination check(s)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    QC(
        INPUT_CHECK.out.reads,
        confindr_db,
        ch_bloom_filter
    )
    ch_versions         = ch_versions.mix(QC.out.versions)
    ch_illumina_trimmed = QC.out.illumina
    ch_ont_trimmed      = QC.out.ont
    ch_pacbio_trimmed   = QC.out.pacbio
    multiqc_files       = multiqc_files.mix(QC.out.qc)
    ch_report           = ch_report.mix(QC.out.confindr_reports, QC.out.fastp_json, QC.out.nanoplot_stats)

    ch_multiqc_illumina = ch_multiqc_illumina.mix(QC.out.qc_illumina)
    ch_multiqc_nanopore = ch_multiqc_nanopore.mix(QC.out.qc_nanopore)
    ch_multiqc_pacbio   = ch_multiqc_pacbio.mix(QC.out.qc_pacbio)


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUB: See which samples are Illumina-only, ONT-only, Pacbio-only
    or have a mix of both for hybrid assembly
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    GROUP_READS(
        ch_illumina_trimmed,
        ch_ont_trimmed,
        ch_pacbio_trimmed
    )
    ch_hybrid_reads     = GROUP_READS.out.hybrid_reads
    ch_short_reads_only = GROUP_READS.out.illumina_only
    ch_ont_reads_only   = GROUP_READS.out.ont_only
    ch_pb_reads_only    = GROUP_READS.out.pacbio_only
    ch_pb_hybrid_reads  = GROUP_READS.out.pacbio_hybrid
    ch_dragonflye       = GROUP_READS.out.dragonflye

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUB: Predict taxonomy from read data
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_reads_for_taxonomy = ch_ont_trimmed.mix(ch_illumina_trimmed,ch_pacbio_trimmed)

    TAXONOMY_PROFILING(
        ch_reads_for_taxonomy,
        kraken2_db
    )
    ch_versions         = ch_versions.mix(TAXONOMY_PROFILING.out.versions)
    ch_report           = ch_report.mix(TAXONOMY_PROFILING.out.report)

    ch_multiqc_illumina = ch_multiqc_illumina.mix(TAXONOMY_PROFILING.out.report_txt.filter{m,r -> m.platform == "ILLUMINA"}.map {m,r -> r })
    ch_multiqc_nanopore = ch_multiqc_nanopore.mix(TAXONOMY_PROFILING.out.report_txt.filter{m,r -> m.platform == "NANOPORE"}.map {m,r -> r })
    ch_multiqc_pacbio   = ch_multiqc_pacbio.mix(TAXONOMY_PROFILING.out.report_txt.filter{m,r -> m.platform == "PACBIO"}.map {m,r -> r })

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUB: Assemble reads based on what data is available
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    /*
    Illumina short-reads only
    */
    ILLUMINA_ASSEMBLY(
        ch_short_reads_only
    )
    ch_versions = ch_versions.mix(ILLUMINA_ASSEMBLY.out.versions)
    ch_assemblies = ch_assemblies.mix(ILLUMINA_ASSEMBLY.out.assembly)

    /*
    ONT assembly including polishing with and without short reads
    */    
    ONT_ASSEMBLY(
        ch_dragonflye,
        homopolish_db
    )
    ch_versions = ch_versions.mix(ONT_ASSEMBLY.out.versions)
    ch_assemblies = ch_assemblies.mix(ONT_ASSEMBLY.out.assembly)

    /*
    Option: Pacbio HiFi reads
    Flye
    */
    PACBIO_ASSEMBLY(
        ch_pb_hybrid_reads
    )
    ch_versions     = ch_versions.mix(PACBIO_ASSEMBLY.out.versions)
    ch_assemblies   = ch_assemblies.mix(PACBIO_ASSEMBLY.out.assembly)
    
    // Find empty assemblies and stop them
    ch_assemblies.branch { m,f ->
        fail: f.countFasta() < 1
        pass: f.countFasta() > 0
    }.set { ch_assemblies_size }

    ch_assemblies_size.fail.subscribe { m, f ->
        log.warn "${m.sample_id} - assembly is empty, stopping sample"
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Tag and optionally remove highly fragmented assemblies
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_assemblies_size.pass.branch { m, f ->
        fail: f.countFasta() > params.max_contigs
        pass: f.countFasta() <= params.max_contigs
    }.set { ch_assemblies_status }

    ch_assemblies_status.fail.subscribe { m, f ->
        log.warn "${m.sample_id} - assembly is highly fragmented!"
    }

    if (params.skip_failed) {
        ch_assemblies_filtered = ch_assemblies_status.pass
    } else {
        ch_assemblies_filtered = ch_assemblies_size.pass
    }
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Clean the meta data object to remove stuff we don't need anymore
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ch_assemblies_filtered.map { m, f ->
        def newMeta = [:]
        newMeta.sample_id = m.sample_id
        tuple(newMeta, f)
    }.set { ch_assemblies_clean }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Calculate coverage against assembled genome
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    COVERAGE(
        ch_assemblies_clean,
        ch_illumina_trimmed,
        ch_ont_trimmed,
        ch_pacbio_trimmed
    )
    ch_versions   = ch_versions.mix(COVERAGE.out.versions)

    ch_multiqc_illumina = ch_multiqc_illumina.mix(COVERAGE.out.report_illumina.map{ m,r -> r })
    ch_multiqc_nanopore = ch_multiqc_nanopore.mix(COVERAGE.out.report_ont.map{ m,r -> r })
    ch_multiqc_pacbio   = ch_multiqc_pacbio.mix(COVERAGE.out.report_pacbio.map{ m,r -> r })

    ch_report = ch_report.mix(
        COVERAGE.out.summary,
        COVERAGE.out.report,
        COVERAGE.out.bam_stats
    )

    /*
    Taxonomically profile the assembly to check composition
    */
    ASSEMBLY_PROFILE(
        ch_assemblies_clean.map { m,a ->
            m.single_end = true
            tuple(m, a)
        },
        kraken2_db,
        checkm_db,
        taxdb
    )
    ch_versions = ch_versions.mix(ASSEMBLY_PROFILE.out.versions)
    ch_report = ch_report.mix(ASSEMBLY_PROFILE.out.report)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUB: Identify and analyse plasmids from draft assemblies
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    PLASMIDS(
        ch_assemblies_clean
    )
    ch_versions = ch_versions.mix(PLASMIDS.out.versions)
    ch_report = ch_report.mix(PLASMIDS.out.reports)

    RENAME_PLASMID_CTG(
        PLASMIDS.out.chromosome,
        'chromosomes.fasta'        
    )
    ch_assembly_without_plasmids = RENAME_PLASMID_CTG.out

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUB: Map reads to chromosome assembly to check 
    for polymorphic positions as indication of read or assembly
    errors
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (!params.skip_variants && !params.skip_optional) {
        VARIANTS(
            ch_illumina_trimmed.mix(ch_ont_trimmed).map { m,r ->
                tuple(m.sample_id,m,r)
            }.join(
                ch_assembly_without_plasmids.map { m,a ->
                    tuple(m.sample_id,a)
                }
            ).map { s,m,r,a ->
                tuple(m,r,a)
            }
        )
        ch_versions         = ch_versions.mix(VARIANTS.out.versions)
        ch_multiqc_illumina = ch_multiqc_illumina.mix(VARIANTS.out.stats.filter { m,s -> m.platform == "ILLUMINA"}.map { m,s -> s})
        ch_multiqc_nanopore = ch_multiqc_nanopore.mix(VARIANTS.out.stats.filter { m,s -> m.platform == "NANOPORE"}.map { m,s -> s})
        ch_report           = ch_report.mix(VARIANTS.out.stats)
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUB: Find the appropriate reference genome+annotation for each assembly
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    FIND_REFERENCES(
        ch_assembly_without_plasmids,
        sourmashdb
    )
    ch_versions     = ch_versions.mix(FIND_REFERENCES.out.versions)
    ch_report       = ch_report.mix(FIND_REFERENCES.out.gbk)
    //ch_assemblies_without_plasmids_with_reference_and_gbk = FIND_REFERENCES.out.assembly_with_ref
    //ch_assemblies_without_plasmids_with_taxa = FIND_REFERENCES.out.assembly_with_tax

    // Assembly with plasmids and the detected reference + gbk/gff
    ch_assemblies_clean.map { m, s ->
        tuple(m.sample_id,s)
    }.join(
        FIND_REFERENCES.out.reference.map { m, r, g, k ->
            tuple(m.sample_id,m,r,g,k)
        }
    ).map { d, s, m, r, g, k ->
        tuple(m,s,r,g,k)
    }.set { ch_assemblies_clean_with_reference_and_gbk }

    // as well as a channel with the clean assembly incl Plasmids and taxon information
    ch_assemblies_clean_with_reference_and_gbk.map { m, s, r, g, k ->
        tuple(m,s)
    }.set { ch_assemblies_clean_with_taxa }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUB: Perform serotyping of assemblies
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (!params.skip_serotyping && !params.skip_optional) {
        SEROTYPING(
            ch_assemblies_clean_with_taxa
        )
        ch_versions     = ch_versions.mix(SEROTYPING.out.versions)
        ch_report       = ch_report.mix(SEROTYPING.out.reports)
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUB: Perform MLST typing of assemblies
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (!params.skip_mlst && !params.skip_optional) {
        MLST_TYPING(
            ch_assemblies_clean_with_taxa
        )
        ch_versions = ch_versions.mix(MLST_TYPING.out.versions)
        ch_report = ch_report.mix(MLST_TYPING.out.report)
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUB: Predict gene models
    We use taxonomy-enriched meta hashes to add
    genus/species to the Prokka output(s)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (!params.skip_annotation && !params.skip_optional) {
        ANNOTATE(
            ch_assemblies_clean_with_taxa,
            ch_prokka_proteins,
            ch_prokka_prodigal
        )
        ch_versions = ch_versions.mix(ANNOTATE.out.versions)
        multiqc_files = multiqc_files.mix(ANNOTATE.out.qc.map { m, r -> r })
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUB: Identify antimocrobial resistance genes
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (!params.skip_amr && !params.skip_optional) {
        AMR_PROFILING(
            ch_assemblies_clean_with_taxa,
            amrfinder_db,
            abricate_dbs
        )
        ch_versions = ch_versions.mix(AMR_PROFILING.out.versions)
        ch_report   = ch_report.mix(AMR_PROFILING.out.amrfinder_report, AMR_PROFILING.out.abricate_report)
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUB: Gauge quality of the assembly
    This does not include the plasmids. 
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    ASSEMBLY_QC(
        ch_assemblies_clean_with_reference_and_gbk,
        busco_lineage,
        busco_db_path
    )
    ch_versions     = ch_versions.mix(ASSEMBLY_QC.out.versions)
    ch_assembly_qc  = ASSEMBLY_QC.out.quast
    multiqc_files   = multiqc_files.mix(ASSEMBLY_QC.out.qc.map { m, r -> r })
    ch_report       = ch_report.mix(ch_assembly_qc)
    ch_report       = ch_report.mix(ASSEMBLY_QC.out.busco_json)

    /*
    Gather all version information
    */
    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUB: Make summary report
    This is optonal in case of unforseen
    issues.
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (!params.skip_report) {
        // standardize the meta hash to enable downstream grouping
        ch_report.map { m, r ->
            def meta = [:]
            meta.sample_id = m.sample_id
            tuple(meta, r)
        }.groupTuple().map { meta,r ->
            tuple(meta.sample_id,r)
        }.join(
            FIND_REFERENCES.out.taxon.map { m ->
                tuple(m.sample_id,m)
            }
        ).map { sid,r,meta ->
            tuple (meta,r)
        }.set { ch_reports_grouped }

        REPORT(
            ch_reports_grouped,
            ch_report_template,
            ch_report_refs,
            CUSTOM_DUMPSOFTWAREVERSIONS.out.yml
        )
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Generate QC reports
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    multiqc_files = multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml)

    MULTIQC(
        multiqc_files.collect(),
        ch_multiqc_config,
        ch_multiqc_logo
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Platform-specific MultiQC reports
    since different technologies are difficult to
    display jointly (scale etc)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    MULTIQC_ILLUMINA(
        ch_multiqc_illumina.collect(),
        ch_multiqc_config,
        ch_multiqc_logo
    )
    MULTIQC_NANOPORE(
        ch_multiqc_nanopore.collect(),
        ch_multiqc_config,
        ch_multiqc_logo
    )
    MULTIQC_PACBIO(
        ch_multiqc_pacbio.collect(),
        ch_multiqc_config,
        ch_multiqc_logo
    )
    
    emit:
    qc = MULTIQC.out.report
}
