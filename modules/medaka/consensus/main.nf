process MEDAKA_CONSENSUS {
    tag "$meta.sample_id"
    label 'short_parallel'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/medaka:2.0.1--py38h8774169_0' :
        'quay.io/biocontainers/medaka:2.0.1--py38h8774169_0' }"

    input:
    tuple val(meta), path(reads), path(assembly)

    output:
    tuple val(meta), path("*consensus.fasta")   , emit: consensus
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p polish
    medaka_consensus \\
        -t $task.cpus \\
        $args \\
        -i $reads \\
        -d $assembly \\
        -o polish

    cp polish/consensus.fasta ${prefix}.medaka.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        medaka: \$( medaka --version 2>&1 | sed 's/medaka //g' )
    END_VERSIONS
    """
}
