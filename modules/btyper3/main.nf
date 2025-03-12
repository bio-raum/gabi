process BTYPER3 {
    tag "${meta.sample_id}"

    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/btyper3:3.4.0--pyhdfd78af_0' :
        'quay.io/biocontainers/btyper3:3.4.0--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*.tsv'), emit: tsv
    path('versions.yml'), emit: versions

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.sample_id
    """
    btyper3 -i $fasta \
    $args \
    -o btyper

    cp btyper/btyper3_final_results/*final_results.txt ${prefix}.btyper3.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        btyper3: \$(btyper3 --version 2>&1 | cut -f2 -d " ")
    END_VERSIONS

    """
}
