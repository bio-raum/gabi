process TABIX_BGZIP {
    tag "$meta.sample_id"
    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/htslib:1.20--h5efdd21_2' :
        'quay.io/biocontainers/htslib:1.20--h5efdd21_2' }"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("${output}")    , emit: output
    tuple val(meta), path("${output}.gzi"), emit: gzi, optional: true
    path  "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix   = task.ext.prefix ?: "${meta.sample_id}"
    extension = "vcf"
    output   = "${prefix}.${extension}.gz"
    // Name the index according to $prefix, unless a name has been requested

    """
    bgzip -c $args -@${task.cpus} $input > ${output}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tabix: \$(echo \$(tabix -h 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

}
