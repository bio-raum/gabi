process AUTOCYCLER_COMPRESS {
    tag "$meta.sample_id"

    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"
    container "mhoeppner/autocycler:0.6.2"

    input:
    tuple val(meta), val(assemblies, stageAs: 'assemblies/?')

    output:
    tuple val(meta), path('autocycler')                         , emit: compressed
    tuple val(meta), path("autocycler.log")                     , emit: log  
    path 'versions.yml'                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"

    """
    autocycler compress \\
        $args
        -i assemblies \\
        -a ${prefix}_autocycler \\
        -t $task.cpus > autocycler.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | cut -d " " -f2 )
    END_VERSIONS
    """
}
