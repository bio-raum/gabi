process UNICYCLER {
    tag "$meta.sample_id"
    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"
    container "https://depot.galaxyproject.org/singularity/unicycler:0.5.1--py311hc84137b_4"
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/unicycler:0.5.1--py311hc84137b_4' :
    //    'quay.io/biocontainers/unicycler:0.5.1--py311hc84137b_4' }"

    input:
    tuple val(meta), path(shortreads), path(longreads)

    output:
    tuple val(meta), path('*.scaffolds.fa'), emit: scaffolds
    tuple val(meta), path('*.assembly.gfa.gz'), emit: gfa
    tuple val(meta), path('*.log')            , emit: log
    path  "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    def short_reads = shortreads ? ( meta.single_end ? "-s $shortreads" : "-1 ${shortreads[0]} -2 ${shortreads[1]}" ) : ""
    def long_reads  = longreads ? "-l $longreads" : ""
    """
    unicycler \\
        --threads $task.cpus \\
        $args \\
        $short_reads \\
        $long_reads \\
        --out ./

    mv assembly.fasta ${prefix}.scaffolds.fa
    mv assembly.gfa ${prefix}.assembly.gfa
    gzip -n ${prefix}.assembly.gfa
    mv unicycler.log ${prefix}.unicycler.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        unicycler: \$(echo \$(unicycler --version 2>&1) | sed 's/^.*Unicycler v//; s/ .*\$//')
    END_VERSIONS
    """

}
