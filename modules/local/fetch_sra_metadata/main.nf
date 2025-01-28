process FETCH_SRA_METADATA {
    tag "Fetch SRA metadata"
    label "process_single"

    conda "${moduleDir}/environment.yml"

    input:
    val bioproject_id
    val email
    path sra_metadata_file
    path trimming_config

    output:
    path "*_updated.csv", optional: true
    path "*_to_process.csv", emit: new_samples

    script:
    def trim_config = trimming_config.name != "NO_FILE" ? "--trimming_config ${trimming_config}" : ""
    """
    fetch_sra_metadata.py \\
        --bioproject_ids "${bioproject_id}" \\
        --email ${email} \\
        --metadata ${sra_metadata_file} \\
        ${trim_config}
    """

    stub:
    """
    touch ${sra_metadata_file}
    """
}
