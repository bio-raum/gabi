
/*
This module combines the final Autocycler steps into one
This is for efficency purposes as all steps are quite short
*/
process AUTOCYCLER_FINISH {
    tag "$meta.sample_id"

    label 'short_parallel'

    conda "${moduleDir}/environment.yml"
    container "mjfos2r/autocycler:0.5.2"

    input:
    tuple val(meta), path(assemblies, stageAs: 'assemblies/*')

    output:
    tuple val(meta), path('*.assembly.fasta')   , emit: fasta
    tuple val(meta), path("autocycler_out")     , emit: results
    tuple val(meta), path("autocycler.log")     , emit: log  
    tuple val(meta), path("*.summary.tsv")      , optional: true, emit: tsv
    path 'versions.yml'                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample_id}"

    """
    shopt -s nullglob
    # Give circular contigs from Plassembler extra clustering weight
    for f in assemblies/*plassembler*.fasta; do
        sed -i 's/circular=True/circular=True Autocycler_cluster_weight=3/' "\$f"
    done

    # Give contigs from Canu and Flye extra consensus weight
    for f in assemblies/*canu*.fasta assemblies/*flye*.fasta; do
        sed -i 's/^>.*\$/& Autocycler_consensus_weight=2/' "\$f"
    done
    shopt -u nullglob

    autocycler compress \\
        $args \\
        -i assemblies \\
        -a autocycler_out \\
        -t $task.cpus > autocycler.log

    autocycler cluster -a autocycler_out >> autocycler.log

    for c in autocycler_out/clustering/qc_pass/cluster_*; do
        autocycler trim -c "\$c" 2>> autocycler.stderr
        autocycler resolve -c "\$c" 2>> autocycler.stderr
    done 

    autocycler combine \\
    -a autocycler_out \\
    -i autocycler_out/clustering/qc_pass/cluster_*/5_final.gfa

    cp autocycler_out/consensus_assembly.fasta ${prefix}.assembly.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        autocycler: \$(autocycler --version | cut -d " " -f2 )
    END_VERSIONS
    """
}
