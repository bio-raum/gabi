include { KMC }     from '../../modules/kmc'
include { RASUSA }  from '../../modules/rasusa'

ch_versions = Channel.from([])

workflow DOWNSAMPLE_READS {

    take:
    reads

    main:

    // Determine genome size from Kmers
    KMC(
        reads
    )
    ch_versions = ch_versions.mix(KMC.out.versions)

    reads.join(
        KMC.out.log
    ).map { m,r,l ->
        gsize = parse_genome_size(l)
        tuple(m,r,gsize)
    }.set { ch_reads_with_genome_size }

    // Downsample reads
    RASUSA(
        ch_reads_with_genome_size
    )
    ch_versions = ch_versions.mix(RASUSA.out.versions)

    emit:
    reads = RASUSA.out.reads
    versions = ch_versions

}

def parse_genome_size(aFile) {

    // defaults to 6Mb in case no base count was reported
    def gsize = '6Mb'

    aFile.eachLine { line ->
        if (line.contains("unique counted k-mers")) {
            def elements = line.trim().split(/\s+/)
            raw = (elements[-1].toInteger()/1000000).round(1)
            gsize = "${raw}Mb"
        }

    }

    return gsize
}
