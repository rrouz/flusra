process FREYJA_VARIANTS {
    tag "$sra"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sra), path(bamFile)
    path reference

    output:
    tuple val(sra), path("${sra}_variants.tsv"), path("${sra}_depth.tsv") , emit: variants
    path  "versions.yml", emit: versions

    script:
    """
    samtools index $bamFile

    freyja variants \\
        $bamFile \\
        --variants ${sra}_variants.tsv \\
        --depths ${sra}_depth.tsv \\
        --ref $reference

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja: \$(echo \$(c --version 2>&1) | sed 's/^.*version //' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${sra}_variant.tsv
    touch ${sra}_depth.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
