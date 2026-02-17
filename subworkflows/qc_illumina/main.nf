include { FASTP }                       from './../../modules/fastp'
include { CAT_FASTQ }                   from './../../modules/cat_fastq'
include { FASTQC }                      from './../../modules/fastqc'
include { BIOBLOOM_CATEGORIZER }        from './../../modules/biobloom/categorizer'

/*
subworkflows
*/
include { CONTAMINATION }               from './../contamination'
include { DOWNSAMPLE_READS }            from './../downsample_reads'

workflow QC_ILLUMINA {
    take:
    reads
    confindr_db
    bloomfilter

    main:

    ch_versions = channel.from([])
    multiqc_files = channel.from([])

    // Split trimmed reads by sample to find multi-lane data set
    reads.map {m, fastq ->
        def newMeta = [:]
        newMeta.sample_id = m.sample_id
        newMeta.platform = m.platform
        newMeta.single_end = m.single_end
        tuple(newMeta,fastq)
    }.groupTuple().branch { meta, fastq ->
        single: fastq.size() == 1
            return [ meta, fastq.flatten()]
        multi: fastq.size() > 1
            return [ meta, fastq.flatten()]
    }.set { ch_reads_illumina }

    // Concatenate samples with multiple PE files
    CAT_FASTQ(
        ch_reads_illumina.multi
    )

    ch_reads_merged = ch_reads_illumina.single.mix(CAT_FASTQ.out.reads)
     
    // Short read trimming and QC
    // this will also standardize the read names to sample_id_R1/2_trimmed.fastq.gz
    // This is only acceptable because we concatenate the reads beforehand!
    FASTP(
        ch_reads_merged
    )
    ch_versions = ch_versions.mix(FASTP.out.versions)
    multiqc_files = multiqc_files.mix(FASTP.out.json.map{ m,j -> j})

    FASTQC(
        FASTP.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions)
    multiqc_files = multiqc_files.mix(FASTQC.out.zip.map { m, z -> z })
    
    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Clean reads against a bloom filter to remove any
    potential host contaminations - currently: horse, from
    blood medium used during growth of campylobacter
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (params.remove_host) {
        BIOBLOOM_CATEGORIZER(
            FASTP.out.reads,
            bloomfilter
        )
        ch_illumina_clean = BIOBLOOM_CATEGORIZER.out.reads
        ch_versions = ch_versions.mix(BIOBLOOM_CATEGORIZER.out.versions)
    } else {
        ch_illumina_clean = FASTP.out.reads
    }

    // Run the contamination subworkflow
    CONTAMINATION(
        ch_illumina_clean,
        confindr_db
    )
    ch_versions = ch_versions.mix(CONTAMINATION.out.versions)
    ch_reads_decont = CONTAMINATION.out.reads
    
    // Downsample reads if a genome size is given
    if (params.max_coverage) {

        // perform downsampling of reads
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
    fastp_json      = FASTP.out.json
    confindr_qc     = CONTAMINATION.out.qc
    reads           = ch_processed_reads
    versions        = ch_versions
    qc              = multiqc_files
    }
