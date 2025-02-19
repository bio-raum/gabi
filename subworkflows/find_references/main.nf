include { SOURMASH_SKETCH } from './../../modules/sourmash/sketch'
include { SOURMASH_SEARCH } from './../../modules/sourmash/search'
include { DOWNLOAD_GENOME } from './../../modules/helper/download_genome'

ch_versions = Channel.from([])

workflow FIND_REFERENCES {
    take:
    assembly
    sourmashdb

    main:

    // sketch the assembly
    SOURMASH_SKETCH(
        assembly
    )
    ch_versions = ch_versions.mix(SOURMASH_SKETCH.out.versions)

    // search sketch against db
    SOURMASH_SEARCH(
        SOURMASH_SKETCH.out.signature,
        sourmashdb
    )
    ch_versions = ch_versions.mix(SOURMASH_SEARCH.out.versions)

    SOURMASH_SEARCH.out.csv.map { m,c ->
        def newMeta = [:]
        (gbk,taxon) = sourmash_get_acc(c)
        newMeta.sample_id = m.sample_id
        newMeta.gbk = gbk
        newMeta.taxon = taxon
        newMeta.db_name = m.db_name
        tuple(newMeta,c)
    }.set { mash_with_gbk}
            
    mash_with_gbk.map { m, r ->
        m.gbk
    }.unique()
    .set { genome_accessions }

    // Download the best reference genome
    DOWNLOAD_GENOME(
        genome_accessions
    )
    ch_versions = ch_versions.mix(DOWNLOAD_GENOME.out.versions)

    ch_genome_with_gff = DOWNLOAD_GENOME.out.sequence.join(DOWNLOAD_GENOME.out.gff).join(DOWNLOAD_GENOME.out.genbank)

    /*
    We use combine here because several assemblies may
    map to the same reference genome
    */
    mash_with_gbk.map { m, r ->
        tuple(m.gbk, m, r)
    }.combine(
        ch_genome_with_gff, by: 0
    ).map { g, m, r, s, a, k ->
        def meta = [:]
        meta.sample_id = m.sample_id
        meta.taxon = m.taxon
        meta.domain = 'Bacteria'
        meta.db_name = m.db_name
        tuple(meta, s, a, k)
    }.set { meta_with_sequence }

    meta_with_genbank = meta_with_sequence.map{m,s,a,k -> [m,k]}

    /*
    Combine the assembly with the best reference genome and annotation
    */
    assembly.map { m, s ->
        tuple(m.sample_id, s)
    }.join(
        meta_with_sequence.map { m, r, g, k ->
            tuple(m.sample_id, m, r, g, k)
        }
    ).map { i,s, m, r, g, k ->
        tuple(m, s, r, g, k)
    }.set { assembly_with_reference_and_gbk }

    // and we create a channel with taxon-enriched metadata and assembly for other analyses
    assembly_with_reference_and_gbk.map { m,s, r, g, k ->
        tuple(m,s)
    }.set { assembly_with_taxa }

    emit:
    taxon = meta_with_sequence.map {m,s,a,k -> m }
    gbk = meta_with_genbank
    reference = meta_with_sequence
    assembly_with_ref = assembly_with_reference_and_gbk
    assembly_with_tax = assembly_with_taxa
    versions = ch_versions
    
}

def sourmash_get_acc(csv) {
    def gbk = ''
    def taxon = 'unknown'
    lines = file(csv).readLines()
    if (lines.size() > 1 ) {
        def elements = lines[1].trim().split(",")
        gbk_file = elements[3].split(" ")[0]
        taxon = elements[3].split(" ")[1..2].join(" ")
        if (gbk_file.contains('GCF_')) {
            gbk = gbk_file.replace('"','')
        }
    }

    return [ gbk, taxon ]
}