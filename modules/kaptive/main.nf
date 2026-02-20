process KAPTIVE {
    tag "${meta.sample_id}"

    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kaptive:3.1.0--pyhdfd78af_0' :
        'quay.io/biocontainers/kaptive:3.1.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.txt'), emit: txt
    tuple val(meta), path('*.json'), emit: json
    path('versions.yml'), emit: versions

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.sample_id
    """
    kaptive assembly \\
    $args \\
    $fasta \\
    -j ${prefix}.kaptive.json \\
    -o ${prefix}.kaptive.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kaptive: \$(kaptive -v 2>&1)
    END_VERSIONS

    """
}
