process SAMTOOLS_INDEX {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.bai")  , emit: bai
    path  "versions.yml"            , emit: versions

    script:
    """
    samtools \\
        index \\
        -@ ${task.cpus-1} \\
        $input

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${input}.bai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
