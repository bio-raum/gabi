process AUTOCYCLER_HELPER {
    tag "${meta.sample_id}|${idx}|${tool}"

    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"
    container "mjfos2r/autocycler:0.5.2"

    input:
    tuple val(meta), path(reads), val(genome_size), val(idx), val(tool)
    val(read_type)

    output:
    tuple val(meta), path('*.fasta')            , emit: fasta
    tuple val(meta), path("autocycler.log")     , emit: log  
    path 'versions.yml'                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}_${idx}-${tool}"
    
    """
    autocycler helper \\
        $tool \\
        $args \\
        --reads $reads \\
        --genome_size $genome_size \\
        --threads $task.cpus \\
        -d tmp \\
        --out_prefix $prefix \\
        --read_type $read_type \\
        --min_depth_rel 0.1 > autocycler.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | cut -d " " -f2 )
    END_VERSIONS
    """
}
