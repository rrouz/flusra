process FETCH_SRA_METADATA {
    tag "Fetch SRA metadata"
    label "process_high"

    conda "${moduleDir}/environment.yml"

    input:
    val bioproject_id
    val email
    path sra_metadata_file

    output:
    path "*_new.txt", emit: new_sra_metadata_file, optional: true
    path "*_updated.csv", emit: updated_sra_metadata_file, optional: true

    script:
    """
    fetch_sra_metadata.py \\
        --bioproject_id ${bioproject_id} \\
        --email ${email} \\
        --metadata ${sra_metadata_file}
    """

    stub:
    """
    touch ${sra_metadata_file}
    """
}
