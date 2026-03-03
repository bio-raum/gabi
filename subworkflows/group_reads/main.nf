/*
This subworkflow bins reads by sample id into singular or mixed emissions
for downstream assembly with the appropriate tool chain
*/
workflow GROUP_READS {
    take:
    reads

    main:

    /*
    ------------------------------
    Setting up the channels
    ------------------------------
    */

    illumina = reads.filter {m,r -> m.platform == "ILLUMINA"}
    ont = reads.filter { m,r -> m.platform == "NANOPORE"}
    pacbio = reads.filter { m,r -> m.platform == "PACBIO"}

    ch_short_reads_for_cross                    = illumina.map { m, r -> [m.sample_id, m, r] }
    ch_ont_reads_for_cross                      = ont.map { m, r -> [m.sample_id, m, r] }
    ch_pb_reads_for_cross                       = pacbio.map { m, r -> [m.sample_id, m, r] }

    ch_short_reads_cross_grouped                = ch_short_reads_for_cross.groupTuple()
    ch_ont_reads_cross_grouped                  = ch_ont_reads_for_cross.groupTuple()
    ch_pb_reads_cross_grouped                   = ch_pb_reads_for_cross.groupTuple()

    /*
    -----------------------------------------------
    Grouping PACBIO reads with optional short reads
    -----------------------------------------------
    */

    // Get the PB only samples
    ch_pb_reads_cross_grouped_joined            = ch_pb_reads_cross_grouped.join(ch_short_reads_cross_grouped, remainder: true) // [ sample_id, [ meta, pb], [meta, ill]]
    ch_pb_reads_cross_grouped_joined_filtered   = ch_pb_reads_cross_grouped_joined.filter { !(it.last()) } // remove those entries that had no joined short reads
    ch_pb_reads_only                            = ch_pb_reads_cross_grouped_joined_filtered.transpose().map { [ it[1], it[2]] }

    // Combine short reads with PB reads based on sample id
    ch_pb_reads_cross_grouped_joined            = ch_short_reads_cross_grouped.join(ch_pb_reads_cross_grouped, remainder: true) // [ sample_id, [ meta, ill], [meta, pb]]

    // Channel were no matching PB reads are available
    ch_pb_reads_cross_grouped_joined_filtered      = ch_pb_reads_cross_grouped_joined.filter { !(it.last()) } // [ sample_id, [meta, ill], null ]

    // Channel with pacbio reads and optional illumina reads for polishing
    // [ meta, [ illumina], pacbio ]

    ch_reads_with_pacbio                        = ch_pb_reads_cross_grouped_joined.filter { it.last() } // [ sample_id, [meta, ill], [meta, pb]]
    ch_reads_with_pacbio_no_short               = ch_reads_with_pacbio.filter { !it[1] }
    ch_reads_with_pacbio.filter { it[1] }.transpose().map { [ it[3], it[2], it[4] ] }.map { m, i, n ->
        def newMeta = [:]
        newMeta.sample_id = m.sample_id
        newMeta.platform = 'ILLUMINA_AND_PACBIO'
        tuple(newMeta, i, n)
    }.set { ch_reads_with_pacbio_and_short }

    /*
    -----------------------------------------------
    Group ONT reads with optional short reads
    -----------------------------------------------
    */

    // Get the ONT only samples
    ch_ont_reads_cross_grouped_joined           = ch_ont_reads_cross_grouped.join(ch_short_reads_cross_grouped, remainder: true)
    ch_ont_reads_cross_grouped_joined_filtered  = ch_ont_reads_cross_grouped_joined.filter { !(it.last()) } // remove those entries that had no joined short reads
    ch_ont_reads_only                           = ch_ont_reads_cross_grouped_joined_filtered.transpose().map { [ it[1], it[2]] }

    // Combine short reads with ONT reads based on sample id
    ch_ont_reads_cross_grouped_joined               = ch_short_reads_cross_grouped.join(ch_ont_reads_cross_grouped, remainder: true)

    // Channel where no matching ONT reads are available
    ch_ont_reads_cross_grouped_joined_filtered      = ch_ont_reads_cross_grouped_joined.filter { !(it.last()) }

    // Channel with nanopore reads and optional illumina reads for polishing
    // [ meta, [ illumina], nanopore ]
    ch_reads_with_nanopore                      = ch_ont_reads_cross_grouped_joined.filter { it.last() }
    ch_reads_with_nanopore_no_short             = ch_reads_with_nanopore.filter { !it[1] }
    ch_reads_with_nanopore.filter { it[1] }.transpose().map { [ it[3], it[2], it[4] ] }.map { m, i, n ->
        def newMeta = [:]
        newMeta.sample_id = m.sample_id
        newMeta.platform = 'ILLUMINA_AND_NANOPORE'
        tuple(newMeta, i, n)
    }.set { ch_reads_with_nanopore_and_short }


    // The paired ONT/Illumina data
    ch_dragonflye                               = ch_reads_with_nanopore_and_short
    // And adding in ONT data without Illumina
    ch_dragonflye                               = ch_dragonflye.mix(ch_reads_with_nanopore_no_short.transpose().map { [ it[2], [], it[3]] })

    // Samples for which we only have short reads
    ch_short_reads_only                         = ch_ont_reads_cross_grouped_joined_filtered.mix(ch_pb_reads_cross_grouped_joined_filtered).transpose().map { [ it[1], it[2]] }.unique()

    /* 
    Now we find short reads that weren't matched to either ONT or PB, or both
    we get this for both ONT and Nanopore (= redundant), so need to reconcile to find truly unmatched short reads
    Structure at this point is: [ meta, [ reads ]]
    We check if any of the long-read/short read pairs contain this short read set, and if so hide it from the short-read only channel
    */

    // repeat the initial strategy of grouping by sample_id
    ch_short_reads_only_with_key                = ch_short_reads_only.map { m, r -> [ m.sample_id, m, r]}
    ch_reads_with_nanopore_and_short_with_key   = ch_reads_with_nanopore_and_short.map { m,s,r -> [ m.sample_id, m,s,r ]}
    ch_reads_with_pacbio_and_short_with_key     = ch_reads_with_pacbio_and_short.map { m,s,r -> [ m.sample_id, m,s,r ]}

    ch_short_reads_with_paired_long_reads       = ch_short_reads_only_with_key.join(ch_reads_with_nanopore_and_short_with_key, remainder: true).join(ch_reads_with_pacbio_and_short_with_key, remainder: true)
    
    // A pure short-read data set will not have joined with either ONT-Ill or PB-Ill
    ch_short_reads_only_final                   = ch_short_reads_with_paired_long_reads.filter { !it[-2] & !it[-1]}.map { [ it[1], it[2]]}

    // Samples with short-reads and matched nanopore reads
    // from [ sample_id, meta1, [illumina_reads ], meta2, [ ont_reads ]]

    emit:
    illumina_only   = ch_short_reads_only_final
    ont_only        = ch_ont_reads_only
    hybrid_reads    = ch_reads_with_nanopore_and_short
    dragonflye      = ch_dragonflye
    pacbio_only     = ch_pb_reads_only
    pacbio_hybrid   = ch_reads_with_pacbio_no_short.transpose().map { [ it[2], [], it[3]] }.mix(ch_reads_with_pacbio_and_short)
    
}
