process IVAR_VARIANTS {
    tag "$meta - $referenceGene"
    label 'process_high'

    conda "${moduleDir}/environment.yml"

    input:
    each referenceGene
    tuple val(meta), path(bamFile)
    path reference

    output:
    path "*_variants.tsv",         emit: variants
    path  "versions.yml", emit: versions

    script:
    def gene = referenceGene.split("\\|")[0]
    """
    samtools index $bamFile

    samtools mpileup --reference $reference -r \"$referenceGene\" -A -d 0 -aa -Q 0 $bamFile | ivar variants -p ${meta}_${gene}_variants -t 0.01 -m 1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        ivar: \$(echo \$(ivar version 2>&1) | sed 's/^.*version //; s/Please.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta}_${gene}_variants.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        ivar: \$(echo \$(ivar version 2>&1) | sed 's/^.*version //; s/Please.*\$//')
    END_VERSIONS
    """
}
