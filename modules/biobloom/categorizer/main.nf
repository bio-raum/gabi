process BIOBLOOM_CATEGORIZER {
    tag "$meta.sample_id"

    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biobloomtools:2.3.5--hdcf5f25_5' :
        'quay.io/biocontainers/biobloomtools:2.3.5--hdcf5f25_5' }"

    input:
    tuple val(meta), path(reads, stageAs: '?/*')
    tuple path(bf), path(txt)

    output:
    tuple val(meta), path('*fastq.gz')          , emit: reads
    path 'versions.yml'                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    def options = meta.single_end ? "" : "-e"
    
    """
        biobloomcategorizer \\
        $options \\
        $args \\
        -p $prefix \\
        -t $task.cpus \\
        -f $bf \\
        $reads

    if [ -f $prefix'_noMatch_1.fq.gz' ]; then
        ln -s $prefix'_noMatch_1.fq.gz' $prefix'_R1.fastq.gz'
    fi
    if [ -f $prefix'_noMatch_2.fq.gz' ]; then
        ln -s $prefix'_noMatch_2.fq.gz' $prefix'_R2.fastq.gz'
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biobloomcategorizer: \$( biobloomcategorizer --version 2>&1 | head -n1 | cut -f3 -d ' ')
    END_VERSIONS
    """
   
}
