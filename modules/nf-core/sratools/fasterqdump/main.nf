process SRATOOLS_FASTERQDUMP {
    tag "$sra"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    val sra

    output:
    tuple val(sra), path('*.fastq'), emit: reads
    path "versions.yml"            , emit: versions

    script:
    """
    fasterq-dump \\
        --threads $task.cpus \\
        ${sra}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sratools: \$(fasterq-dump --version 2>&1 | grep -Eo '[0-9.]+')
    END_VERSIONS
    """
}
