process SHOVILL_FIX_CONTIGS {
    tag "${meta.sample_id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython.convert:1.3.3--pyh5e36f6f_0' :
        'quay.io/biocontainers/biopython.convert:1.3.3--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(contigs, stageAs: '?/*')

    output:
    tuple val(meta), path('*.fasta') , emit: contigs
    path 'versions.yml'             , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.sample_id

    """
    shovill_fix_contigs.py --input $contigs \\
    --output ${prefix}.fasta $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version  | sed -e "s/Python //")
    END_VERSIONS
    """
}
