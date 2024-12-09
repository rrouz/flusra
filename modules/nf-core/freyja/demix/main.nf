process FREYJA_DEMIX {
    tag "$sra"
    label 'process_low'

    conda "${moduleDir}/environment.yml"

    input:
    tuple val(sra), path(variant), path(depth)
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
        --output ${sra}.demixed


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        freyja: \$(echo \$(freyja --version 2>&1) | sed 's/^.*version //' )
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${sra}.demixed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
