process RACON {
    tag "$meta.sample_id"
    label 'short_parallel'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/racon:1.4.20--h9a82719_1' :
        'quay.io/biocontainers/racon:1.4.20--h9a82719_1' }"

    input:
    tuple val(meta), path(reads), path(assembly), path(paf)

    output:
    tuple val(meta), path('*racon.fasta')    , emit: improved_assembly
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    """
    racon $args -t $task.cpus \\
        ${reads} \\
        ${paf} \\
        ${assembly} > \\
        ${prefix}.racon.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        racon: \$( racon --version 2>&1 | sed 's/^.*v//' )
    END_VERSIONS
    """
}
