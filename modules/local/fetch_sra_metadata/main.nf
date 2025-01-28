process FETCH_SRA_METADATA {
    tag "Fetch SRA metadata"
    label "process_single"

    conda "${moduleDir}/environment.yml"

    input:
    val bioproject_id
    val email
    path sra_metadata_file

    output:
    path "*_updated.csv", optional: true
    path "*_to_process.csv", emit: new_samples

    script:
    """
    fetch_sra_metadata.py \\
        --bioproject_ids "${bioproject_id}" \\
        --email ${email} \\
        --metadata ${sra_metadata_file}
    """

    stub:
    """
    touch ${sra_metadata_file}
    """
}
