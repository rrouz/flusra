process MINIMAP2_ALIGN {
    tag "$sra"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sra), path(reads)
    val reference

    output:
    tuple val(sra), path("*.bam") , emit: bam
    path  "versions.yml"           , emit: versions

    script:
    """
    minimap2 \\
        -ax sr \\
        -t $task.cpus \\
        $reference \\
        $reads \\
        | samtools view \\
            -bS \\
            | samtools sort \\
                -o ${sra}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
