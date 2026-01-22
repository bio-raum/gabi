process TAXONKIT_REFORMAT {
    tag "$meta.sample_id"
    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/taxonkit:0.18.0--h9ee0642_0':
        'quay.io/biocontainers/taxonkit:0.18.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(report)
    path taxdb

    output:
    tuple val(meta), path("*.reformat.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"

    """
    taxonkit \\
        reformat \\
        $args \\
        --format "k__{k}|p__{p}|c__{c}|o__{o}|f__{f}|g__{g}|s__{s}" --miss-rank-repl "unclassified" \\
        --data-dir $taxdb \\
        --threads $task.cpus \\
        --out-file ${prefix}.reformat.tsv \\
        $report
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        taxonkit: \$( taxonkit version | sed 's/.* v//' )
    END_VERSIONS
    """

}
