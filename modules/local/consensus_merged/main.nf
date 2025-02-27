process CONSENSUS_MERGED {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), val(genes), path(consensus_files)

    output:
    tuple val(meta), path("${meta.id}.fa"), emit: merged_fasta
    path "versions.yml",                     emit: versions

    script:
    """
    # Define segment order array
    segments=("PB2" "PB1" "PA" "HA" "NP" "NA" "MP" "NS")

    # Process each segment in order
    for segment in "\${segments[@]}"; do
        for consensus in ${consensus_files}; do
            if [[ \$consensus == *"_\${segment}_cns.fa" ]]; then
                cat \$consensus >> "${meta.id}.fa"
                echo >> "${meta.id}.fa"
            fi
        done
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(echo \$(bash --version 2>&1) | sed 's/^.*version //; s/ Release.*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.id}.fa"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bash: \$(echo \$(bash --version 2>&1) | sed 's/^.*version //; s/ Release.*\$//')
    END_VERSIONS
    """
}