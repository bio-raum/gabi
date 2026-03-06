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
    tools = channel.from(["canu", "flye", "metamdbg", "miniasm", "necat", "nextdenovo","raven"])

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

    AUTOCYCLER_SUBSAMPLE.out.reads.flatMap { m, g, chunks ->
        chunks.withIndex().collect {element, index -> [ m, file(element), g, index ]}
    }.set { ch_read_chunks_with_genome_size }

    ch_reads_with_size_and_tool = ch_read_chunks_with_genome_size.combine(tools)
                
    AUTOCYCLER_HELPER(
        ch_reads_with_size_and_tool,
        read_type
    )
    ch_versions = ch_versions.mix(AUTOCYCLER_HELPER.out.versions)

    AUTOCYCLER_HELPER.out.fasta.groupTuple().view()

    AUTOCYCLER_FINISH(
        AUTOCYCLER_HELPER.out.fasta.groupTuple()
    )
    ch_versions = ch_versions.mix(AUTOCYCLER_FINISH.out.versions)

    emit:
    fasta = AUTOCYCLER_FINISH.out.fasta
    versions = ch_versions
    
}

def parse_genome_size(aFile) {

    // defaults to 6Mb in case no base count was reported
    def gsize = '6000000'

    aFile.eachLine { line ->
        if (line.contains("unique counted k-mers")) {
            def elements = line.trim().split(/\s+/)
            def raw = elements[-1].toInteger()
            gsize = raw
        }

    }

    return gsize
}