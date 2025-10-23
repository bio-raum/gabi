/*
Include Modules
*/

include { PORECHOP_ABI }                from './../../modules/porechop/abi'
include { CAT_FASTQ  }                  from './../../modules/cat_fastq'
include { NANOPLOT }                    from './../../modules/nanoplot'
include { CHOPPER }                     from './../../modules/chopper'
include { SEQKIT_REPLACE }              from './../../modules/seqkit/replace'

// Subworkflows
include { CONTAMINATION }               from './../contamination'
include { DOWNSAMPLE_READS }            from './../downsample_reads'

workflow QC_NANOPORE {
    take:
    reads
    confindr_db

    main:

    ch_versions = Channel.from([])
    multiqc_files = Channel.from([])

    if (!params.skip_porechop) {
        // Nanopore adapter trimming
        PORECHOP_ABI(
            reads
        )
        ch_versions         = ch_versions.mix(PORECHOP_ABI.out.versions)
        multiqc_files       = multiqc_files.mix(PORECHOP_ABI.out.log.map { m, l -> l })
        ch_porechop_reads   = PORECHOP_ABI.out.reads
    } else {
        ch_porechop_reads   = reads
    }

    // Merge Nanopore reads per sample
    ch_porechop_reads.map { m, fastq ->
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
    }.set { ch_reads_ont }

    CAT_FASTQ(
        ch_reads_ont.multi
    )

    // The trimmed ONT reads, concatenated by sample
    ch_ont_trimmed = ch_reads_ont.single.mix(CAT_FASTQ.out.reads)

    // Filter the reads by size and quality
    CHOPPER(
        ch_ont_trimmed
    )
    ch_versions = ch_versions.mix(CHOPPER.out.versions)

    CHOPPER.out.fastq.branch { m,r ->
        pass: r.countFastq() >= params.ont_min_reads
        fail: r.countFastq() < params.ont_min_reads
    }.set { ch_chopped_reads }

    // Stop a sample if the number of ONT reads is under a threshold
    ch_chopped_reads.fail.subscribe { m,r ->
        log.warn "Stopping ONT read set ${m.sample_id} - not enough reads surviving.\nConsider adjusting ont_min_length, ont_min_reads and ont_min_q."
    }

    // Run contamination check
    CONTAMINATION(
        ch_chopped_reads.pass,
        confindr_db
    )
    ch_versions = ch_versions.mix(CONTAMINATION.out.versions)

    // Generate a plot of the trimmed reads
    NANOPLOT(
        ch_chopped_reads.pass
    )
    ch_versions = ch_versions.mix(NANOPLOT.out.versions)
    multiqc_files = multiqc_files.mix(NANOPLOT.out.txt.map { m, r -> r })

    if (params.max_coverage) {

        // Replace tabs in ONT fastq headers, else KMC will not work
        SEQKIT_REPLACE(
            ch_chopped_reads.pass
        )
        ch_versions = ch_versions.mix(SEQKIT_REPLACE.out.versions)

        // Perform downsampling of reads
        DOWNSAMPLE_READS(
            SEQKIT_REPLACE.out.fastx
        )

        ch_versions = ch_versions.mix(DOWNSAMPLE_READS.out.versions)
        ch_processed_reads = DOWNSAMPLE_READS.out.reads

    } else {
        ch_processed_reads = ch_chopped_reads.pass
    }

    emit:
    confindr_report = CONTAMINATION.out.report
    confindr_json   = CONTAMINATION.out.confindr_json
    confindr_qc     = CONTAMINATION.out.qc
    reads           = ch_processed_reads
    qc              = multiqc_files
    nanoplot_stats  = NANOPLOT.out.txt
    versions        = ch_versions
}

