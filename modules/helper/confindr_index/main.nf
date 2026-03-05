process CONFINDR_INDEX {

    tag "${meta.sample_id}"
    label 'medium_serial'

    conda "${moduleDir}/environment.yml"
    container "mhoeppner/confindr:0.7.4"

    input:
    tuple val(meta), path(fasta)

    output:
    path(fasta)            , emit: fasta
    path("*_kma*")         , emit: kma

    script:
    kma = fasta.getBaseName() + "_kma"

    """
    kma index -i $fasta -o $kma

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        confindr: \$(confindr.py --version 2>&1 | sed -e "s/ConFindr //g")
    END_VERSIONS
    """
}
