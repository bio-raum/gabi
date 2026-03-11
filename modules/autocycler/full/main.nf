process AUTOCYCLER_FULL {
    tag "$meta.sample_id"

    label 'extra_long_parallel'

    conda "${moduleDir}/environment.yml"
    container "varunshamanna/autocycler:v0.5.2"

    input:
    tuple val(meta), path(reads)
    val(read_type)

    output:
    tuple val(meta), path('*.assembly.fasta')      , emit: fasta
    path 'versions.yml'                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    """
    autocycler_full.sh \\
        $reads \\
        $task.cpus \\
        $args \\
        $read_type

    cp autocycler_out/consensus_assembly.fasta ${prefix}.assembly.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | cut -d " " -f2 )
    END_VERSIONS
    """
}
