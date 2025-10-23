process SEQKIT_REPLACE {
    tag "${meta.sample_id}"
    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/seqkit:2.10.1--he881be0_0'
        : 'quay.io/biocontainers/seqkit:2.10.1--he881be0_0'}"

    input:
    tuple val(meta), path(fastx, stageAs: '?/*')

    output:
    tuple val(meta), path("*.fastq.gz"), emit: fastx
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
 
    """
    seqkit replace \\
    ${args} \\
    --threads ${task.cpus} \\
    -o ${prefix}.fastq.gz \\
    ${fastx}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqkit: \$( seqkit version | sed 's/seqkit v//' )
    END_VERSIONS
    """

}
