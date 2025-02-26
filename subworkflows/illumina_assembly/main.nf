
include { SHOVILL }                             from './../../modules/shovill'
include { RENAME_CTG as RENAME_SHOVILL_CTG }    from './../../modules/rename_ctg'

ch_versions = Channel.from([])

workflow ILLUMINA_ASSEMBLY {

    take:
    reads

    main:
    /*
    Shovill
    */
    SHOVILL(
        reads
    )
    ch_versions = ch_versions.mix(SHOVILL.out.versions)

    //Shovill generates generic output names, must rename to sample id
    RENAME_SHOVILL_CTG(
        SHOVILL.out.contigs,
        'fasta'
    )

    emit:
    assembly = RENAME_SHOVILL_CTG.out
    versions = ch_versions
    
}