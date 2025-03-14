process RASUSA {
    tag "$meta.sample_id"
    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rasusa:0.3.0--h779adbc_1' :
        'quay.io/biocontainers/rasusa:0.3.0--h779adbc_1' }"

    input:
    tuple val(meta), path(reads), val(gsize)

    output:
    tuple val(meta), path('*.fastq.gz'), emit: reads
    path 'versions.yml'                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    def output   = meta.single_end ? "--output ${prefix}.fastq.gz" : "--output ${prefix}_1.fastq.gz ${prefix}_2.fastq.gz"
    """
    rasusa \\
        $args \\
        --genome-size $gsize \\
        --input $reads \\
        $output
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rasusa: \$(rasusa --version 2>&1 | sed -e "s/rasusa //g")
    END_VERSIONS
    """
}
