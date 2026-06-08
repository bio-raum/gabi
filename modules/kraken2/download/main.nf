process KRAKEN2_DOWNLOAD {
    label 'medium_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    path(url)

    output:
    path("minikraken2*")         , emit: db

    script:
    archive = url.toString().split('/')[-1]

    """
    tar -xvf $archive
    rm *.tgz
    """
}
