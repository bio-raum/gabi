process HAMRONIZATION_ABRICATE {
    tag "$meta.sample_id"
    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hamronization:1.1.4--pyhdfd78af_0' :
        'quay.io/biocontainers/hamronization:1.1.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(reports)
    val(format)
    val(software_version)
    val(reference_db_version)

    output:
    tuple val(meta), path('*.json'), optional: true, emit: json
    tuple val(meta), path('*.tsv') , optional: true, emit: tsv
    path 'versions.yml'            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${reports.getBaseName()}"
    """
    hamronize \\
        abricate \\
        ${reports} \\
        $args \\
        --format ${format} \\
        --analysis_software_version ${software_version} \\
        --reference_database_version ${reference_db_version} \\
        --output ${prefix}.${format}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hamronization: \$(echo \$(hamronize --version 2>&1) | cut -f 2 -d ' ' )
    END_VERSIONS
    """
}
