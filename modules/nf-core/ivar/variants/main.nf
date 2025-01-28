process IVAR_VARIANTS {
    tag "$meta.id - $referenceGene"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"

    input:
    each referenceGene
    tuple val(meta), path(bamFile)
    path reference
    val gff_files
    val variant_threshold
    val variant_min_depth

    output:
    path "*_variants.tsv",         emit: variants
    path  "versions.yml", emit: versions

    script:
    def gene = referenceGene.split("\\|")[0]
    def gff_file_arg = ""
    if (!gff_files.isEmpty()) {
        if (gff_files.containsKey(gene)) {
            if (new File(gff_files[gene]).exists()) {
                gff_file_arg = "-g ${gff_files[gene]}"
            } else {
                throw new Exception("GFF file for gene $gene does not exist at ${gff_files[gene]}")
            }
        }
    }
    """
    samtools index $bamFile

    samtools mpileup --reference $reference -r \"$referenceGene\" -A -d 0 -aa -Q 0 $bamFile | ivar variants -p ${meta.id}_${gene}_variants -t $variant_threshold -m $variant_min_depth -r $reference $gff_file_arg

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        ivar: \$(echo \$(ivar version 2>&1) | sed 's/^.*version //; s/Please.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_${gene}_variants.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
        ivar: \$(echo \$(ivar version 2>&1) | sed 's/^.*version //; s/Please.*\$//')
    END_VERSIONS
    """
}
