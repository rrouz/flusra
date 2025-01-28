process FREYJA_DEMIX {
    tag "$meta.id"
    label 'process_high', 'process_low_memory'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(meta), path(variant), path(depth)
    path barcode

    output:
    path "*.demixed",     emit: demixed
    path  "versions.yml", emit: versions

    script:
    """
    freyja demix \\
        $variant \\
        $depth \\
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
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
