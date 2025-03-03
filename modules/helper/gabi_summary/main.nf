process GABI_SUMMARY {
    tag "${meta.sample_id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dajin2:0.5.5--pyhdfd78af_0' :
        'quay.io/biocontainers/dajin2:0.5.5--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reports, stageAs: '?/*')
    path(yaml)

    output:
    tuple val(meta), path('*.json') , emit: json
    path 'versions.yml'             , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.sample_id
    result = prefix + '.json'

    """
    gabi_json.py --sample ${meta.sample_id} \
    --taxon '${meta.taxon}' \
    --yaml $yaml \\
    $args \
    --output $result

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version  | sed -e "s/Python //")
    END_VERSIONS
    """
}
