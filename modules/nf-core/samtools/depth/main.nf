process SAMTOOLS_DEPTH {
    tag "$meta - $referenceGene"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    each referenceGene
    tuple val(meta), path(bamFile)

    output:
    path "*_depth.tsv"        , emit: tsv
    path "versions.yml" , emit: versions

    script:
    def gene = referenceGene.split("\\|")[0]
    """
    samtools \\
        depth \\
        --threads ${task.cpus-1} \\
        -aa \\
        $bamFile \\
        > ${meta}_${gene}_depth.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta}_${gene}_depth.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        ivar: \$(echo \$(ivar version 2>&1) | sed 's/^.*version //; s/Please.*\$//')
    END_VERSIONS
    """
}
