process IVAR_CONSENSUS {
    tag "$sra - $referenceGene"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    each referenceGene
    tuple val(sra), path(bamFile)
    path reference
    val consensus_threshold
    val consensus_min_depth

    output:
    path  "*_cns.fa",         emit: consensus
    path  "versions.yml", emit: versions

    script:
    def gene = referenceGene.split("\\|")[0]
    """
    samtools index $bamFile

    samtools mpileup -r \"$referenceGene\" -A -d 0 -aa -Q 0 $bamFile | ivar consensus -p ${sra}_${gene}_cns -t $consensus_threshold -m $consensus_min_depth

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        ivar: \$(echo \$(ivar version 2>&1) | sed 's/^.*version //; s/Please.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${sra}_${gene}_cns.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        ivar: \$(echo \$(ivar version 2>&1) | sed 's/^.*version //; s/Please.*\$//')
    END_VERSIONS
    """
}
