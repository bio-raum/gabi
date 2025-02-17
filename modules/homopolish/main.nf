process HOMOPOLISH {
    tag "$meta.sample_id"
    label 'short_parallel'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/homopolish:0.4.1--pyhdfd78af_1' :
        'quay.io/biocontainers/homopolish:0.4.1--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(assembly)
    path(mash)

    output:
    tuple val(meta), path('*polished.fasta')    , emit: polished
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    homopolish polish \\
        $args \\
        -s $mash \\
        -a ${assembly} \\
        -o polished

    cp polished/*.fasta . 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        homopolish: \$( homopolish -v 2>&1 | cut -f4 -d " " )
    END_VERSIONS
    """
}
