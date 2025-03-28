process FLYE {
    tag "$meta.sample_id"

    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/flye:2.9--py39h6935b12_1' :
        'quay.io/biocontainers/flye:2.9--py39h6935b12_1' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*.assembly.fasta')      , emit: fasta
    tuple val(meta), path('*.assembly_graph.gfa')  , emit: gfa
    tuple val(meta), path('*.assembly_graph.gv')   , emit: gv
    tuple val(meta), path('*.txt')     , emit: txt
    tuple val(meta), path('*.log')     , emit: log
    tuple val(meta), path('*.json')    , emit: json
    path 'versions.yml'                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    """
    flye \\
        $args $reads \\
        --out-dir . \\
        --threads $task.cpus \\

    cp assembly.fasta ${prefix}.assembly.fasta
    cp assembly_graph.gfa ${prefix}.assembly_graph.gfa
    cp assembly_graph.gv ${prefix}.assembly_graph.gv
    mv assembly_info.txt ${prefix}.assembly_info.txt
    mv flye.log ${prefix}.flye.log
    mv params.json ${prefix}.params.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flye: \$( flye --version )
    END_VERSIONS
    """
}
