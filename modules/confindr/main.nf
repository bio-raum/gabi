process CONFINDR {
    tag "$meta.sample_id"
    label 'medium_parallel'

    conda "${moduleDir}/environment.yml"
    // This is a work-around since the official biocontainer container for confindr produces incongruent results to the conda package!
    // Something broken there - until fixed, this solution remains in place.
    container "mhoeppner/confindr:0.7.4"

    input:
    tuple val(meta), path(reads, stageAs: 'input_dir/*')
    path db

    output:
    tuple val(meta), path('confindr_results/*contamination.csv'),   emit: csv, optional: true
    tuple val(meta), path('confindr_results/*confindr_log.txt'),    emit: log
    tuple val(meta), path('confindr_results/*confindr_report.csv'), emit: report
    tuple val(meta), path('confindr_results/*_rmlst.csv'),          emit: rmlst, optional: true
    path 'versions.yml',                                            emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}_${meta.platform}"
    def db_options = db ? "-d ${db}" : ''
    def options = meta.platform == "NANOPORE" ? "-dt Nanopore -bf 0.1 -q 20 -b 5" : "-dt Illumina -bf 0.05 -q 20 -b 2 -fid _R1 -rid _R2"
    
    """
    confindr.py \\
        -Xmx ${task.memory.toGiga()}G \\
        --threads $task.cpus \\
        -i input_dir \\
        -o confindr_results \\
        $args $options $db_options

    mv confindr_results/confindr_log.txt confindr_results/${prefix}_confindr_log.txt
    mv confindr_results/confindr_report.csv confindr_results/${prefix}_confindr_report.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        confindr: \$(confindr.py --version 2>&1 | sed -e "s/ConFindr //g")
    END_VERSIONS
    """
}
