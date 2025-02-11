process DNAAPLER {
    tag "$meta.sample_id"
    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/dnaapler:1.1.0--pyhdfd78af_0':
        'quay.io/biocontainers/dnaapler:1.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(assembly, stageAs: '?/*')

    output:
    tuple val(meta), path('*.fasta') , emit: fasta
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"

    """
    dnaapler $args \\
    --input $assembly \\
    --output reorient \\
    --threads $task.cpus \\
    --prefix $prefix

    cp reorient/*reoriented.fasta ${prefix}.fasta 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dnaapler: \$(dnaapler --version 2>&1 | cut -d ' ' -f 3)
    END_VERSIONS
    """
}
