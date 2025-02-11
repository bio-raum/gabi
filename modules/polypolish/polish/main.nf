process POLYPOLISH_POLISH {
    tag "$meta.sample_id"

    label 'medium_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/polypolish:0.6.0--h3ab6199_3' :
        'quay.io/biocontainers/polypolish:0.6.0--h3ab6199_3' }"

    input:
    tuple val(meta), path(assembly), path(sam)

    output:
    tuple val(meta), path('*.polished.fasta')      , emit: fasta
    path 'versions.yml'                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    """
    polypolish polish $args $assembly $sam > ${prefix}.polished.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$( flye --version )
    END_VERSIONS
    """
}
