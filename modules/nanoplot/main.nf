process NANOPLOT {
    tag "$meta.sample_id"
    label 'short_parallel'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/nanoplot:1.44.1--pyhdfd78af_0' :
        'quay.io/biocontainers/nanoplot:1.44.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(ontfile)

    output:
    tuple val(meta), path('*.html')                , emit: html
    tuple val(meta), path('*.png') , optional: true, emit: png
    tuple val(meta), path('*.txt')                 , emit: txt
    tuple val(meta), path('*.log')                 , emit: log
    path  'versions.yml'                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def input_file = ("$ontfile".endsWith('.fastq.gz')) ? "--fastq ${ontfile}" :
        ("$ontfile".endsWith('.txt')) ? "--summary ${ontfile}" : ''
    def prefix = meta.sample_id
    
    """
    NanoPlot \\
        $args \\
        -t $task.cpus \\
        $input_file

    mv NanoStats.txt ${prefix}.NanoStats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nanoplot: \$(echo \$(NanoPlot --version 2>&1) | sed 's/^.*NanoPlot //; s/ .*\$//')
    END_VERSIONS
    """
}
