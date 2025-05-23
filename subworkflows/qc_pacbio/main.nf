
include { RASUSA }                          from './../../modules/rasusa'
include { CAT_FASTQ  }                      from './../../modules/cat_fastq'
include { FASTQC }                          from './../../modules/fastqc'

/*
subworkflows
*/
include { CONTAMINATION }               from './../contamination'
include { DOWNSAMPLE_READS }            from './../downsample_reads'

workflow QC_PACBIO {
    take:
    reads
    confindr_db

    main:

    ch_versions = Channel.from([])
    multiqc_files = Channel.from([])

    // Merge Nanopore reads per sample
    reads.groupTuple().branch { meta, fastq ->
        single: fastq.size() == 1
            return [ meta, fastq.flatten()]
        multi: fastq.size() > 1
            return [ meta, fastq.flatten()]
    }.set { ch_reads_pb }

    CAT_FASTQ(
        ch_reads_pb.multi
    )

    // The trimmed ONT reads, concatenated by sample
    ch_pb_trimmed = ch_reads_pb.single.mix(CAT_FASTQ.out.reads)

    FASTQC(
        reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    multiqc_files = multiqc_files.mix(FASTQC.out.zip.map { m, z -> z })

    CONTAMINATION(
        ch_pb_trimmed,
        confindr_db
    )
    ch_versions = ch_versions.mix(CONTAMINATION.out.versions)
    ch_reads_decont = CONTAMINATION.out.reads

    if (params.max_coverage) {

        DOWNSAMPLE_READS(
            ch_reads_decont
        )
        ch_versions = ch_versions.mix(DOWNSAMPLE_READS.out.versions)
        ch_processed_reads = DOWNSAMPLE_READS.out.reads
        
    } else {
        ch_processed_reads = ch_reads_decont
    }

    emit:
    confindr_report = CONTAMINATION.out.report
    confindr_json   = CONTAMINATION.out.confindr_json
    confindr_qc = CONTAMINATION.out.qc
    reads = ch_processed_reads
    qc = multiqc_files
    versions = ch_versions
    }
