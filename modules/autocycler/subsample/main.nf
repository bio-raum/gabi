process AUTOCYCLER_SUBSAMPLE {
    tag "$meta.sample_id"

    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"
    container "mhoeppner/autocycler:0.6.2"

    input:
    tuple val(meta), path(reads), val(genome_size)

    output:
    tuple val(meta), val(genome_size), path('subsampled*/*.fastq')    , emit: reads
    tuple val(meta), path("autocycler.log")                     , emit: log  
    path 'versions.yml'                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"

    """
    autocycler subsample \\
        $args \\
        --reads $reads \\
        --genome_size $genome_size \\
        --out_dir subsampled_${prefix} > autocycler.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | cut -d " " -f2 )
    END_VERSIONS
    """
}
