/*
Include Modules
*/

include { PORECHOP_ABI }                from './../../modules/porechop/abi'
include { RASUSA }                      from './../../modules/rasusa'
include { CAT_FASTQ  }                  from './../../modules/cat_fastq'
include { NANOPLOT }                    from './../../modules/nanoplot'
include { CHOPPER }                     from './../../modules/chopper'

// Subworkflows
include { CONTAMINATION }               from './../contamination'


ch_versions = Channel.from([])
multiqc_files = Channel.from([])

workflow QC_NANOPORE {
    take:
    reads
    confindr_db

    main:

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
    ch_porechop_reads.map { m,r ->
        def newMeta = [:]
        newMeta.sample_id = m.sample_id
        newMeta.platform = m.platform
        newMeta.single_end = m.single_end
        tuple(newMeta,r)
    }.groupTuple().branch { meta, reads ->
        single: reads.size() == 1
            return [ meta, reads.flatten()]
        multi: reads.size() > 1
            return [ meta, reads.flatten()]
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

    if (params.genome_size) {
        ch_chopped_reads.pass.countFastq()

        RASUSA(
            ch_chopped_reads.pass
        )
        ch_versions = ch_versions.mix(RASUSA.out.versions)
        ch_processed_reads = RASUSA.out.reads
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
