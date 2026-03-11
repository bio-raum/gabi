process CONFINDR_DATABASE_SETUP {

    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    // We need 0.7.4 for actual analysis, but 0.8 allows us to index the db files. 
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/confindr:0.8.2--pyhdfd78af_0':
        'quay.io/biocontainers/confindr:0.8.2--pyhdfd78af_0' }"

    output:
    path('confindr'),   emit: db 
    path 'versions.yml',   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    
    """
    confindr_database_setup \\
        -i \\
        -o confindr \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        confindr: \$(confindr.py --version 2>&1 | sed -e "s/ConFindr //g")
    END_VERSIONS
    """
}
