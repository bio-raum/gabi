process PLASSEMBLER_RUN {
    tag "$meta.sample_id"

    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/plassembler:1.8.1--pyhdfd78af_0' :
        'quay.io/biocontainers/plassembler:1.8.1--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(sreads), path(lreads),path(flye_folder)
    val(db)

    output:
    tuple val(meta), path('*.assembly.fasta')   , emit: fasta    
    path 'versions.yml'                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    def short_reads = (sreads.size() == 1) ? "-1 $sreads" : "-1 ${sreads[0]} -2 ${sreads[1]}"

    """

    export TERM=xterm-256color
    
    plassembler run \
    -d $db \
    -l $lreads \
    $short_reads \
    -t ${task.cpus} \
    --flye_directory $flye_folder \
    --keep_chromosome \
    $args

    cat plassembler.output/chromosome.fasta plassembler.output/plassembler_plasmids.fasta > ${prefix}.assembly.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        plassembler: \$( plassembler --version 2>&1 | sed 's/^.*version //' )
    END_VERSIONS
    """
}
