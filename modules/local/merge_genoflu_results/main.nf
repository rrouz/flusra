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
    export OUTPUT_DIR="${params.outdir}/genoflu"
    merge_genoflu.py \\

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
