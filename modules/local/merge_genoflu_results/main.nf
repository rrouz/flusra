process MERGE_GENOFLU_RESULTS {
    tag "merge_results"
    label 'process_medium'
    
    input:
    path(tsv_files)

    output:
    path "genoflu_results.tsv", emit: merged_results
    path "versions.yml", emit: versions

    script:
    """
    export OUTPUT_DIR="${params.genoflu_outdir ?: params.outdir + '/genoflu'}"
    merge_genoflu.py \\
        --input '${tsv_files}'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genoflu: \$(genoflu.py --version 2>&1 || echo "NA")
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch genoflu_results.tsv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genoflu: \$(genoflu.py --version 2>&1 || echo "NA")
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}