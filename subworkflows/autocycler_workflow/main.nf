include { KMC }                                 from '../../modules/kmc'
include { AUTOCYCLER_HELPER }                   from'./../../modules/autocycler/helper'
include { AUTOCYCLER_SUBSAMPLE }                from'./../../modules/autocycler/subsample'
include { AUTOCYCLER_FINISH }                   from'./../../modules/autocycler/finish'

/* Implements the full autocycler workflow outlined in 
https://github.com/rrwick/Autocycler/blob/main/pipelines/Automated_Autocycler_Bash_script_by_Ryan_Wick/autocycler_full.sh
*/
workflow AUTOCYCLER_WORKFLOW {

    take:
    reads
    read_type

    main:
    ch_versions = channel.from([])

    // Determine genome size from Kmers
    // Autocycler genome size check ist way too slow.
    KMC(
        reads
    )
    ch_versions = ch_versions.mix(KMC.out.versions)

    reads.join(
        KMC.out.log
    ).map { m,r,l ->
        def gsize = parse_genome_size(l)
        tuple(m,r,gsize)
    }.set { ch_reads_with_genome_size }

    // Subsample reads for parallel processing
    AUTOCYCLER_SUBSAMPLE(
        ch_reads_with_genome_size
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

def parse_genome_size(aFile) {

    // defaults to 5Mb in case no base count was reported
    def gsize = '5000000'

    aFile.eachLine { line ->
        if (line.contains("unique counted k-mers")) {
            def elements = line.trim().split(/\s+/)
            def raw = elements[-1].toInteger()
            // Capped at 14MB for the largest known bacterial genome - else use 5MB
            if (raw <= 14000000) {
                gsize = raw
            } else {
                log.warn "Genome size estimate exceeds limits for bacterial genomes - capping at 5MB\nMake sure to check reads for contamination."
            }
        }
    }
    return gsize
}

// Not all data types may be assembled with the same tools
def tool_list(meta) {
    def tools = []
    if (meta.platform.contains("NANOPORE")) {
        tools = ["flye", "metamdbg", "miniasm", "necat", "raven"]
    } else if (meta.platform.contains("PACBIO")) {
        if (params.pacbio_hifi) {
            tools = ["flye", "metamdbg", "miniasm", "raven"]
        } else {
            tools = ["flye", "metamdbg", "miniasm", "raven", "canu"]
        }
    } else {
        log.warn "No known sequencing platform attached to reads of sample ${meta.sample_id}"
    }
    return tools
}