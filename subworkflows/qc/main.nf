/*
Subworkflows
*/
include { QC_ILLUMINA }     from './../qc_illumina'
include { QC_NANOPORE }     from './../qc_nanopore'
include { QC_PACBIO }       from './../qc_pacbio'

/*
Modules
*/
include { CONFINDR2MQC_SUMMARY } from './../../modules/helper/confindr2mqc_summary'

workflow QC {
    take:
    reads
    confindr_db
    bloomfilter

    main:


    ch_versions         = channel.from([])
    multiqc_files       = channel.from([])
    ch_confindr_reports = channel.from([])
    ch_confindr_json    = channel.from([])
    ch_qc               = channel.from([])

    // Divide reads up into their sequencing technologies
    reads.branch { meta, fastq ->
        illumina: meta.platform == 'ILLUMINA'
        ont: meta.platform == 'NANOPORE'
        pacbio: meta.platform == 'PACBIO'
        torrent: meta.platform == 'TORRENT'
    }.set { ch_reads }

    ch_reads.torrent.subscribe { m, r ->
        log.warn "Torrent data not yet supported, skipping ${m.sample_id}..."
    }

    /*
    Trim and QC Illumina reads
    */
    QC_ILLUMINA(
        ch_reads.illumina,
        confindr_db,
        bloomfilter
    )
    ch_illumina_trimmed = QC_ILLUMINA.out.reads
    ch_confindr_reports = ch_confindr_reports.mix(QC_ILLUMINA.out.confindr_report)
    ch_confindr_json    = ch_confindr_json.mix(QC_ILLUMINA.out.confindr_json)
    ch_versions         = ch_versions.mix(QC_ILLUMINA.out.versions)
    multiqc_files       = multiqc_files.mix(QC_ILLUMINA.out.qc)

    /*
    Trim and QC nanopore reads
    */
    QC_NANOPORE(
        ch_reads.ont,
        confindr_db
    )
    ch_ont_trimmed      = QC_NANOPORE.out.reads
    ch_versions         = ch_versions.mix(QC_NANOPORE.out.versions)
    ch_confindr_reports = ch_confindr_reports.mix(QC_NANOPORE.out.confindr_report)
    ch_confindr_json    = ch_confindr_json.mix(QC_NANOPORE.out.confindr_json)
    /*
    Trim and QC Pacbio HiFi reads
    */
    QC_PACBIO(
        ch_reads.pacbio,
        confindr_db
    )
    ch_pacbio_trimmed   = QC_PACBIO.out.reads
    ch_confindr_reports = ch_confindr_reports.mix(QC_PACBIO.out.confindr_report)
    ch_confindr_json    = ch_confindr_json.mix(QC_PACBIO.out.confindr_json)
    ch_versions         = ch_versions.mix(QC_PACBIO.out.versions)
    multiqc_files       = multiqc_files.mix(QC_PACBIO.out.qc)

    /*
    Summarize all ConfindR reports from the previously
    generated JSON format to find samples that have failed
    in any of their contributing reads (Illumina and Pacbio only)
    */
    CONFINDR2MQC_SUMMARY(
        ch_confindr_json.map { m, j -> j }.collect()
    )
    ch_qc = ch_qc.mix(CONFINDR2MQC_SUMMARY.out.json)

    emit:
    fastp_json      = QC_ILLUMINA.out.fastp_json
    nanoplot_stats  = QC_NANOPORE.out.nanoplot_stats
    confindr_reports = ch_confindr_reports
    qc_illumina     = QC_ILLUMINA.out.qc.mix(QC_ILLUMINA.out.confindr_qc)
    qc_nanopore     = QC_NANOPORE.out.qc.mix(QC_NANOPORE.out.confindr_qc)
    qc_pacbio       = QC_PACBIO.out.qc.mix(QC_PACBIO.out.confindr_qc)
    illumina        = ch_illumina_trimmed
    ont             = ch_ont_trimmed
    pacbio          = ch_pacbio_trimmed
    versions        = ch_versions
    qc              = ch_qc
    }
