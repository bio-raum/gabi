process CHEWBBACA_FILTER_SCHEMA {
    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.23--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.23--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(schema_dir), path(list)

    output:
    tuple val(meta), path(results) , emit: schema

    when:

    script:
    def args = task.ext.args ?: ''
    results = schema_dir.name + "_cgMLST"

    """
    chewbbaca_filter_schema.py --schema $schema_dir --list $list --output $results
    
    """
}
