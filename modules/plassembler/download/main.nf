process PLASSEMBLER_DOWNLOAD {
    tag "$meta.sample_id"

    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plassembler:1.8.1--pyhdfd78af_0' :
        'quay.io/biocontainers/plassembler:1.8.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("plassembler")    , emit: db
    path 'versions.yml'                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plassembler download -d plassembler

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plassembler: \$( plassembler --version 2>&1 | sed 's/^.*version //' )
    END_VERSIONS
    """
}
