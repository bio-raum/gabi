include { AUTOCYCLER_HELPER }                   from'./../../modules/autocycler/helper'
include { AUTOCYCLER_SUBSAMPLE }                from'./../../modules/autocycler/subsample'
include { AUTOCYCLER_FINISH }                   from'./../../modules/autocycler/finish'
include { GENOMESIZE }                          from './../genomesize'

/* Implements the full autocycler workflow outlined in 
https://github.com/rrwick/Autocycler/blob/main/pipelines/Automated_Autocycler_Bash_script_by_Ryan_Wick/autocycler_full.sh
*/
workflow AUTOCYCLER_WORKFLOW {

    take:
    reads
    read_type

    main:
    ch_versions = channel.from([])

    /*
    Determine the genome size from KMErs
    Note that this will not work perfectly on very poor or too deeply
    sampled data!
    */
    GENOMESIZE(
        reads
    )
    ch_versions = ch_versions.mix(GENOMESIZE.out.versions)

    // Subsample reads for parallel processing
    AUTOCYCLER_SUBSAMPLE(
        GENOMESIZE.out.reads_with_genome_size
    )
    ch_versions = ch_versions.mix(AUTOCYCLER_SUBSAMPLE.out.versions)

    // Subsample reads into 4 random subsets, and pass list index along for file naming
    AUTOCYCLER_SUBSAMPLE.out.reads.flatMap { m, g, chunks ->
        chunks.withIndex().collect {element, index -> [ m, file(element), g, index ]}
    }.set { ch_read_chunks_with_genome_size }

    // Join reads with the proper list of assembly tools
    ch_read_chunks_with_genome_size.map { m, r, g, i ->
        def tools = tool_list(m)
        tuple(m, r, g, i, tools)
    }.flatMap { m, r, g, i, tools -> 
        tools.collect { t -> [ m, r, g, i, t ]}
    }.set { ch_reads_with_size_and_tool }
                
    // Performs the actual assembly with one of the requested tools for the read subsets
    AUTOCYCLER_HELPER(
        ch_reads_with_size_and_tool,
        read_type
    )
    ch_versions = ch_versions.mix(AUTOCYCLER_HELPER.out.versions)

    // Perform consensus finding and produce the final assembly
    AUTOCYCLER_FINISH(
        AUTOCYCLER_HELPER.out.fasta.groupTuple()
    )
    ch_versions = ch_versions.mix(AUTOCYCLER_FINISH.out.versions)

    emit:
    fasta = AUTOCYCLER_FINISH.out.fasta
    versions = ch_versions
    
}

// Not all data types may be assembled with the same tools
def tool_list(meta) {

    def tools = []

    if (params.autocycler_tools) {
        log.info "Using user-provided list of tools for Autocycler"
        tools = params.tools.split(",")
    } else {
        if (meta.platform.contains("NANOPORE")) {
            if (params.onthq) {
                tools = ["flye", "miniasm", "necat", "raven", "plassembler"]
            } else {
                tools = ["flye", "miniasm", "necat", "raven", "plassembler"]
            }
        } else if (meta.platform.contains("PACBIO")) {
            if (params.pacbio_hifi) {
                tools = ["flye", "hifiasm", "plassembler"]
            } else {
                tools = ["flye", "miniasm", "raven", "canu", "plassembler" ]
            }
        } else {
            log.warn "No known sequencing platform attached to reads of sample ${meta.sample_id}"
        }
    }

    return tools
}