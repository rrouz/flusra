process MERGE_GENOFLU_RESULTS {
    tag "merge_results"
    label 'process_medium'
    
    conda "${moduleDir}/environment.yml"
    
    input:
    path tsv_files
    path existing_results

    output:
    path "${params.genoflu_results}", emit: merged_results
    path "versions.yml", emit: versions

    script:
    """    
    merge_genoflu.py \\
        --output_file "${params.genoflu_results}"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genoflu: \$(genoflu.py --version 2>&1 || echo "NA")
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch "${params.genoflu_results}"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genoflu: \$(genoflu.py --version 2>&1 || echo "NA")
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}