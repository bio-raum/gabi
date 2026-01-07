process HELPER_FORMAT_TAXONKIT {
    tag "${meta.sample_id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04' :
        'ubuntu:20.04' }"

    input:
    tuple val(meta), path(kraken), path(taxkit)

    output:
    tuple val(meta), path('*.taxonkit.txt')  , emit: txt
    path 'versions.yml'     , emit: versions

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.sample_id
    result = prefix + '.taxonkit.txt'

    """
    cut -f3 $taxkit > tmp
    awk 'BEGIN{{FS="\t"}}; FNR==NR{{ a[FNR""] = \$2 FS \$3 FS \$4; next }}{{ print a[FNR""] FS \$0 }}' $kraken tmp > $result

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(echo \$(awk --version 2>&1) | sed 's/^.*(AWK) //; s/, API.*\$//')
    END_VERSIONS
    """
}
