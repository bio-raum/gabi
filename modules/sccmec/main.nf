process SCCMEC {
    tag "${meta.sample_id}"

    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sccmec:1.2.0--hdfd78af_0' :
        'quay.io/biocontainers/sccmec:1.2.0--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path('*sccmec.tsv'), emit: tsv
    path('versions.yml'), emit: versions

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.sample_id
    """
    sccmec --input $fasta \
    $args \
    --prefix $prefix

    cp ${prefix}.tsv ${prefix}.sccmec.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sccmec: \$(sccmec --version 2>&1 | grep sccmec_targets | cut -f4 -d " ")
    END_VERSIONS

    """
}
