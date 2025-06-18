
include { SHOVILL }                             from './../../modules/shovill'
include { UNICYCLER }                           from './../../modules/unicycler'
include { SHOVILL_FIX_CONTIGS }                 from './../../modules/helper/shovill_fix_contigs'


workflow ILLUMINA_ASSEMBLY {

    take:
    reads

    main:

    ch_versions = Channel.from([])

    if (!params.unicycler) {
         /*
        Shovill
        */
        SHOVILL(
            reads
        )
        ch_versions = ch_versions.mix(SHOVILL.out.versions)

        /*
        Shovill randomizes naming for contigs of equal length which makes results non-reproducible
        This module will apply a reproducible naming schema to the shovill contigs
        */
        SHOVILL_FIX_CONTIGS(
            SHOVILL.out.contigs
        )
        ch_versions = ch_versions.mix(SHOVILL_FIX_CONTIGS.out.versions)
        ch_assembly = SHOVILL_FIX_CONTIGS.out.contigs

    } else {
       
        ch_unicycler_reads = reads.map { m,i -> [ m, i, []]}

        // Use unicycler; module expects optional long reads, which we set to []
        UNICYCLER(
            ch_unicycler_reads
        )
        ch_versions = ch_versions.mix(UNICYCLER.out.versions)
        ch_assembly = UNICYCLER.out.scaffolds
    }

    emit:
    assembly = ch_assembly
    versions = ch_versions
    
}