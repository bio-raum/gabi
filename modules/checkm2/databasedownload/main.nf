def downloadZenodoApiEntry(zenodo_id) {
    // Download metadata from Zenodo API, setting "Accept: application/json" header
    def api_url  = "https://zenodo.org/api/records/${zenodo_id}"
    def conn     = new URL(api_url).openConnection()
    conn.setRequestProperty('Accept', 'application/json')
    conn.setRequestProperty('User-Agent', "Nextflow ${nextflow.version ?: ''}".trim())

    def api_text = conn.getInputStream().getText('UTF-8')
    def parser   = new groovy.json.JsonSlurper()

    return parser.parseText(api_text)
}

process CHECKM2_DATABASEDOWNLOAD {
    label 'short_serial'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ngsfetch:0.1.1--pyh7e72e81_0':
        'quay.io/biocontainers/ngsfetch:0.1.1--pyh7e72e81_0' }"

    input:
    val(db_zenodo_id)

    output:
    tuple val(meta), path("checkm2_db_v${db_version}.dmnd"), emit: database
    path("versions.yml")                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    zenodo_id  = db_zenodo_id ?: 14897628  // Default to version 3 if no ID provided
    api_data   = downloadZenodoApiEntry(zenodo_id)
    db_version = api_data.metadata.version
    checksum   = api_data.files[0].checksum.replaceFirst(/^md5:/, "md5=")
    meta       = [id: 'checkm2_db', version: db_version]
    """
    # Automatic download is broken when using singularity/apptainer (https://github.com/chklovski/CheckM2/issues/73)
    # So it's necessary to download the database manually
    aria2c \
        ${args} \
        --checksum ${checksum} \
        https://zenodo.org/records/${zenodo_id}/files/checkm2_database.tar.gz

    tar -xzf checkm2_database.tar.gz
    db_path=\$(find -name *.dmnd)
    mv \$db_path checkm2_db_v${db_version}.dmnd

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aria2: \$(echo \$(aria2c --version 2>&1) | grep 'aria2 version' | cut -f3 -d ' ')
    END_VERSIONS
    """
}
