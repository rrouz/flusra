process FREYJA_DEMIX {
    tag "$meta.id"
    label 'process_high', 'process_low_memory'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(variant), path(depth)
    path barcode
    val demix_autoadapt
    val demix_depthcutoff

    output:
    path "*.demixed",     emit: demixed
    path  "versions.yml", emit: versions

    script:
    def autoadapt = demix_autoadapt ? "--autoadapt" : ""
    def depth_cutoff = demix_depthcutoff ? "--depthcutoff $demix_depthcutoff" : ""
    """
    freyja demix \\
        $variant \\
        $depth \\
        $autoadapt \\
        $depth_cutoff \\
        --barcodes $barcode \\
        --output ${meta.id}.demixed


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja: \$(echo \$(freyja --version 2>&1) | sed 's/^.*version //' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}.demixed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja: \$(echo \$(freyja --version 2>&1) | sed 's/^.*version //' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
