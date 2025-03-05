process GABI_QC {
    tag "${meta.sample_id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.27.1--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.27.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(summary)
    path(refs)

    output:
    path('*.qc.json')       , emit: json
    path 'versions.yml'     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.sample_id
    result = prefix + '.qc.json'

    """
    gabi_json_qc.py --input $summary \
    --refs $refs \
    $args \
    --output $result

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version  | sed -e "s/Python //")
    END_VERSIONS
    """
}
