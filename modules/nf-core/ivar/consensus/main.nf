process IVAR_CONSENSUS {
    tag "$meta - $referenceGene"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    input:
    each referenceGene
    tuple val(meta), path(bamFile)

    output:
    tuple path("*.fa"),   emit: consensus
    path  "versions.yml", emit: versions

    script:
    def gene = referenceGene.split("\\|")[0]
    """
    samtools index $bamFile

    samtools mpileup -r \"$referenceGene\" -A -d 0 -aa -Q 0 $bamFile | ivar consensus -p ${meta}_${gene}_cns -t 0.5 -m 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        ivar: \$(echo \$(ivar version 2>&1) | sed 's/^.*version //; s/Please.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta}_${gene}_cns.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        ivar: \$(echo \$(ivar version 2>&1) | sed 's/^.*version //; s/Please.*\$//')
    END_VERSIONS
    """
}
