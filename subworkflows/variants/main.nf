include { SNIPPY_RUN }      from './../../modules/snippy/run'
include { MEDAKA_VARIANT }  from './../../modules/medaka/variant'
include { TABIX_TABIX }     from './../../modules/tabix/tabix'
include { TABIX_BGZIP }     from './../../modules/tabix/bgzip'
include { BCFTOOLS_STATS }  from './../../modules/bcftools/stats'

/*
Call variants against the assembled genome
to check how many sites are polymorphic and indicative
of e.g. conatmination or errors
*/
workflow VARIANTS {

    take:
    reads_with_assembly

    main:

    ch_versions = Channel.from([])
    multiqc_files = Channel.from([])
    ch_variants = Channel.from([])

    reads_with_assembly.branch { m,r,a ->
        nanopore: m.platform == "NANOPORE"
        illumina: m.platform == "ILLUMINA"
    }.set { reads_by_platform }

    // Medaka is a stand-alone workflow for ONT reads to perform alignment and variant calling
    MEDAKA_VARIANT(
        reads_by_platform.nanopore
    )
    ch_versions = ch_versions.mix(MEDAKA_VARIANT.out.versions)

    // Compress the Medaka VCF file
    TABIX_BGZIP(
        MEDAKA_VARIANT.out.vcf
    )
    ch_versions = ch_versions.mix(TABIX_BGZIP.out.versions)
    ch_variants = ch_variants.mix(TABIX_BGZIP.out.output)

    // Snippy is a stand-alone workflow for Illumina reads to perform alignment and variant calling
    SNIPPY_RUN(
        reads_by_platform.illumina
    )
    ch_versions = ch_versions.mix(SNIPPY_RUN.out.versions)
    ch_variants = ch_variants.mix(SNIPPY_RUN.out.vcf_gz)

    // Index vcf.gz files
    TABIX_TABIX(
        SNIPPY_RUN.out.vcf_gz.mix(TABIX_BGZIP.out.output)
    )
    ch_versions = ch_versions.mix(TABIX_TABIX.out.versions)

    // Generate stats for VCF files
    BCFTOOLS_STATS(
        ch_variants.join(
            TABIX_TABIX.out.tbi
        )
    )
    ch_versions = ch_versions.mix(BCFTOOLS_STATS.out.versions)
    multiqc_files = multiqc_files.mix(BCFTOOLS_STATS.out.stats.map { m,s -> s })

    emit:
    vcf = ch_variants
    qc  = multiqc_files
    stats = BCFTOOLS_STATS.out.stats
    versions = ch_versions
}