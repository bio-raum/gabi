process MINIMAP2_ALIGN {
    tag "$meta.sample_id"
    label 'medium_parallel'

    // Note: the versions here need to match the versions used in the mulled container below and minimap2/index
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' :
        'quay.io/biocontainers/mulled-v2-66534bcbb7031a148b13e2ad42583020b9cd25c4:3161f532a5ea6f1dec9be5667c9efc2afdac6104-0' }"

    input:
    tuple val(meta), path(reads), path(reference)
    val(format)

    output:
    tuple val(meta), path("*.bam")                       , optional: true, emit: bam
    tuple val(meta), path("*.paf")                       , optional: true, emit: paf
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    def bam_index = "${prefix}.bam"
    def algorithm = meta.platform == "NANOPORE" ? "-x map-ont" : "-x map-hifi"
    def rg = "-R \"@RG\\tID:${prefix}_${meta.platform}\\tPL:${meta.platform}\\tSM:${meta.sample_id}\""

    if (format == "paf") {

        """
        minimap2 \\
            $args \\
            -t $task.cpus \\
            $algorithm \\
            $reference \\
            $reads \\
             2>&1 1> ${prefix}.paf

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            minimap2: \$(minimap2 --version 2>&1)
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """

    } else {
        """
        minimap2 \\
            $args \\
            -a \\
            -t $task.cpus \\
            $rg \\
            $algorithm \\
            $reference \\
            $reads \\
            | samtools sort -@ ${task.cpus} -o ${bam_index} -
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            minimap2: \$(minimap2 --version 2>&1)
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }

}
