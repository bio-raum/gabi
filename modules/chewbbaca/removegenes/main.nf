process CHEWBBACA_REMOVEGENES {

    tag "${meta.sample_id}"

    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/chewbbaca:3.3.4--pyhdfd78af_0' :
        'quay.io/biocontainers/chewbbaca:3.3.4--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(report), val(filters)

    output:
    tuple val(meta), path(results)  , emit: report
    tuple val(meta), path(results)  , emit: profile
    path('versions.yml')            , emit: versions

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.sample_id
    results = "${prefix}_results_alleles_cgMLST.tsv"

    """
    chewBBACA.py RemoveGenes \\
    -i $report \\
    -g $filters \\
    -o $results \\

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chewBBACA: \$(chewBBACA.py --version 2>&1 | sed -e "s/chewBBACA version: //g")
    END_VERSIONS

    """
}
