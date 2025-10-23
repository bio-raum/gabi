process KMC {

    tag "${meta.sample_id}"

    label 'short_parallel'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kmc:3.2.4--haf24da9_3' :
        'quay.io/biocontainers/kmc:3.2.4--haf24da9_3' }"

    input:
    tuple val(meta), path(reads, stageAs: '?/*')

    output:
    tuple val(meta), path('*kmc.txt')   , emit: log
    path('versions.yml')                , emit: versions

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.sample_id
    def ci = (meta.platform == "NANOPORE") ? 10 : 2

    """
    for i in $reads ; do echo \$i >> files.txt ; done;

    mkdir -p kmc
    $args \\
    kmc -sm \\
    -m${task.memory.toGiga()-1} \\
    -ci${ci} \\
    -k21 \\
    -t${task.cpus} \\
    @files.txt \\
    ${prefix}_kmc.out \\
    kmc > ${prefix}_${meta.platform}_kmc.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kmc: \$(kmc -h 2>&1 | head -n1 | cut -f5 -d ' ')
    END_VERSIONS

    """
}
