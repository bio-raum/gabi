process FASTP {
    tag "${meta.sample_id}"

    label 'short_parallel'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.4--hadf994f_2' :
        'quay.io/biocontainers/fastp:0.23.4--hadf994f_2' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path('*trimmed.fastq.gz'), emit: reads
    tuple val(meta), path("*.json"), emit: json
    path('versions.yml'), emit: versions

    script:

    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.sample_id

    def options = meta.single_end ? "--trim_tail1 1" : "--trim_tail1 1 --trim_tail2 1"

    r1 = reads[0]

    suffix = '_trimmed.fastq.gz'

    json = prefix + '.fastp.json'
    html = prefix + '.fastp.html'

    r1_trim = prefix + "_R1" + suffix

    if (meta.single_end) {
        """
        fastp --in1 ${r1} \
        --out1 $r1_trim \
        $options \
        -w ${task.cpus} \
        -j $json \
        -h $html $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS
        """
    } else {
        r2 = reads[1]
        r2_trim = prefix + "_R2" + suffix

        """
        fastp --in1 ${r1} --in2 ${r2} \
        --out1 $r1_trim \
        --out2 $r2_trim \
        --detect_adapter_for_pe \
        $options \
        -w ${task.cpus} \
        -j $json \
        -h $html \
        $args

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            fastp: \$(fastp --version 2>&1 | sed -e "s/fastp //g")
        END_VERSIONS

        """
    }
}
