process BWAMEM2_MEM_POLYPOLISH {
    tag "$meta.sample_id"
    label 'short_parallel'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2d15960ccea84e249a150b7f5d4db3a42fc2d6c3-0' :
        'quay.io/biocontainers/mulled-v2-e5d375990341c5aef3c9aff74f96f66f65375ef6:2d15960ccea84e249a150b7f5d4db3a42fc2d6c3-0' }"

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path("*.sam")  , emit: sam , optional:true
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"
    def extension = "bam"
    def rg = "-R \"@RG\\tID:${prefix}_${meta.platform}\\tPL:${meta.platform}\\tSM:${meta.sample_id}\""
    def R1 = reads[0]
    if (meta.single_end) {

        """
        INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

        bwa-mem2 \\
            mem \\
            $args \\
            -t $task.cpus \\
            $rg \\
            \$INDEX \\
            $R1 > ${prefix}_R1.sam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """

    } else {
        def R2 = reads[1]
        """
        INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`

        bwa-mem2 \\
            mem \\
            $args \\
            -t $task.cpus \\
            $rg \\
            \$INDEX \\
            $R1 > ${prefix}_R1.sam
        
        bwa-mem2 \\
            mem \\
            $args \\
            -t $task.cpus \\
            $rg \\
            \$INDEX \\
            $R2 > ${prefix}_R2.sam

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bwamem2: \$(echo \$(bwa-mem2 version 2>&1) | sed 's/.* //')
            samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        END_VERSIONS
        """
    }

}
