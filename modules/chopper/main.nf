process CHOPPER {
    tag "$meta.sample_id"
    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chopper:0.3.0--hd03093a_0':
        'quay.io/biocontainers/chopper:0.3.0--hd03093a_0' }"

    input:
    tuple val(meta), path(fq)

    output:
    tuple val(meta), path('*.chopped.fastq.gz') , emit: fastq
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def args2  = task.ext.args2  ?: ''
    def args3  = task.ext.args3  ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    def result = prefix + ".chopped.fastq.gz"

    if ("$fq" == "${result}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    zcat \\
        $args \\
        $fq | \\
    chopper \\
        --threads $task.cpus \\
        $args2 | \\
    gzip -c \\
        $args3 > $result

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: \$(chopper --version 2>&1 | cut -d ' ' -f 2)
    END_VERSIONS
    """
}
