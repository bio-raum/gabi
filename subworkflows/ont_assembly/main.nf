include { FLYE as FLYE_ONT }                        from '../../modules/flye'
include { MEDAKA_CONSENSUS }                        from '../../modules/medaka/consensus'
include { MINIMAP2_ALIGN as MINIMAP2_ALIGN_PAF }    from '../../modules/minimap2/align'
include { HOMOPOLISH as HOMOPOLISH_ONT }            from '../../modules/homopolish'
include { DNAAPLER }                                from '../../modules/dnaapler'
include { POLYPOLISH_POLISH }                       from '../../modules/polypolish/polish'
include { BWAMEM2_INDEX as BWAMEM2_INDEX_POLYPOLISH } from '../../modules/bwamem2/index'
include { BWAMEM2_MEM_POLYPOLISH }                  from '../../modules/bwamem2/mem_polypolish'
include { AUTOCYCLER_WORKFLOW }                     from './../autocycler_workflow'
include { GENOMESIZE }                              from './../genomesize'

/* 
This workflow is inspired by https://github.com/rpetit3/dragonflye
Since Dragonflye isn't regularly maintained, GABI re-implements the
basic (slightly simplified) logic into a subworkflow instead, with 
some additional steps
*/
workflow ONT_ASSEMBLY {

    take:
    reads // [ meta, [short_reads], ont_reads ]
    homopolish_db

    main:

    ch_versions = channel.from([])

    // Get long reads
    reads.map { m,s,o ->
        tuple(m,o)
    }.set { lreads }

    // Get short reads if they exist 
    reads.map { m,s,o ->
        tuple(m,s)
    }.filter { it.last() }
    .set { sreads }

    // Determine genome size from this read set
    GENOMESIZE(
        lreads
    )
    ch_versions = ch_versions.mix(GENOMESIZE.out.versions)

    if (params.autocycler) {
        AUTOCYCLER_WORKFLOW(
            GENOMESIZE.out.reads_with_genome_size,
            "ont_r10"
        )
        ch_long_read_assembly = AUTOCYCLER_WORKFLOW.out.fasta
        ch_versions = ch_versions.mix(AUTOCYCLER_WORKFLOW.out.versions)
    } else {
        // FLYE long read assembler
        FLYE_ONT(
            GENOMESIZE.out.reads_with_genome_size
        )
        ch_versions = ch_versions.mix(FLYE_ONT.out.versions)
        ch_long_read_assembly = FLYE_ONT.out.fasta
    }

    ch_medaka_polished = channel.from([])
    
    if (params.medaka_model) {

        MEDAKA_CONSENSUS(
            lreads.join(
                ch_long_read_assembly
            )
        )
        ch_versions = ch_versions.mix(MEDAKA_CONSENSUS.out.versions)
        ch_medaka_polished = ch_medaka_polished.mix(MEDAKA_CONSENSUS.out.consensus)

    } else {

        // Check if reads contain information about basecalling model
        lreads.map { m,r ->
            def status = check_ont_model(r)
            tuple(m, r, status )
        }.branch { m, r, s ->
            with_model: s == true
            no_model: s == false
            tuple(m,r)
        }.set { lreads_with_model_status }

        MEDAKA_CONSENSUS(
            lreads_with_model_status.with_model.join(
                ch_long_read_assembly
            )
        )
        ch_versions = ch_versions.mix(MEDAKA_CONSENSUS.out.versions)
        ch_medaka_polished = ch_medaka_polished.mix(MEDAKA_CONSENSUS.out.consensus)

        // Get any assemblies which were not polished
        lreads_with_model_status.no_model.join(
            ch_long_read_assembly
        ).map { m,l,a -> 
            tuple(m, a)
        }.set { ch_skip_polishing }

        ch_medaka_polished = ch_medaka_polished.mix(ch_skip_polishing)

        ch_skip_polishing.subscribe { m, a -> 
            log.warn "Skipping Medaka polishing for ${m.sample_id} - no basecalling model defined in file or from command line!"
        }

    }

    // Join polished Medaka assembly with optional short reads
    ch_medaka_polished.map { m,c ->
        tuple(m.sample_id,m,c)
    }.join(
        sreads.map {m,r ->
            tuple(m.sample_id,r)
        }, remainder: true
    ).map { key, m, p, r ->
        tuple(m,p,r)
    }.branch {
        with: it.last()
        without: !it.last()
    }.set { polished_with_short_reads }

    // Homopolish to remove homopolymer errors when no short reads
    // are available 
    if (params.homopolish) {
            HOMOPOLISH_ONT(
            polished_with_short_reads.without.map { m,p,r ->
                tuple(m,p)
            },
            homopolish_db
        )
        ch_versions = ch_versions.mix(HOMOPOLISH_ONT.out.versions)
        ch_homopolished = HOMOPOLISH_ONT.out.polished
    } else {
        ch_homopolished = polished_with_short_reads.without.map { m,p,r -> tuple(m,p) }
    }

    // Create BWA index
    BWAMEM2_INDEX_POLYPOLISH(
        polished_with_short_reads.with.map { m,a,r -> 
            tuple(m,a)
        }
    )
    ch_versions = ch_versions.mix(BWAMEM2_INDEX_POLYPOLISH.out.versions)

    // Align short reads and create one SAM file per mate
    BWAMEM2_MEM_POLYPOLISH(
        polished_with_short_reads.with.map { m,a,r ->
            tuple(m,r)
        }.join(
            BWAMEM2_INDEX_POLYPOLISH.out.index
        )
    )
    ch_versions = ch_versions.mix(BWAMEM2_MEM_POLYPOLISH.out.versions)

    // Run Polypolish if short reads are available
    POLYPOLISH_POLISH(
        polished_with_short_reads.with.map { m,a,r ->
            tuple(m,a)
        }.join(
            BWAMEM2_MEM_POLYPOLISH.out.sam
        )
    )
    ch_versions = ch_versions.mix(POLYPOLISH_POLISH.out.versions)

    // Combine shot-read polished assemblies with homopolish assemblies for which we had no short reads
    ch_polished_assembly = ch_homopolished.mix(POLYPOLISH_POLISH.out.fasta)

    // Consistently orient chromosomes
    DNAAPLER(
        ch_polished_assembly
    )
    ch_versions = ch_versions.mix(DNAAPLER.out.versions)

    emit:
    assembly = DNAAPLER.out.fasta
    versions = ch_versions

}

/*
Dorado encodes the base calling model into the fastQ headers
Check if it is there, else data has no model information
*/
def check_ont_model(fastq) {

    def has_model = false

    fastq.withInputStream { fq ->
        def isGzip = fastq.name.toString().endsWith('.gz')
        def stream = isGzip ? new java.util.zip.GZIPInputStream(fq) as InputStream : fq as InputStream
        def decoder = new InputStreamReader(stream, 'ASCII')
        def buffered = new BufferedReader(decoder)
        line = buffered.readLine()
        if (line.contains("model")) {
            has_model = true
        }
    }

    return has_model
}