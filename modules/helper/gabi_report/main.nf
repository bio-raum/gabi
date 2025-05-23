process GABI_REPORT {
    tag "All"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dajin2:0.5.5--pyhdfd78af_0' :
        'quay.io/biocontainers/dajin2:0.5.5--pyhdfd78af_0' }"

    input:
    path(reports)
    path(template)
    path(yml)

    output:
    path('*.html')          , emit: html
    path 'versions.yml'     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: params.run_name
    result = prefix + '.html'

    version = workflow.manifest.version
    call = workflow.commandLine
    wd = workflow.workDir

    """
    gabi_v2.py --template $template \
    --input $yml \
    --version $version \
    --call '$call' \
    --wd $wd \
    $args \
    --output $result

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version  | sed -e "s/Python //")
    END_VERSIONS
    """
}
