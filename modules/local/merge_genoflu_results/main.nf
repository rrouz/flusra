process MERGE_GENOFLU_RESULTS {
    tag "merge_results"
    label 'process_medium'
    
    publishDir path: "${params.outdir}/genoflu", mode: 'copy', pattern: 'genoflu_results.tsv'

    input:
    path(tsv_files) 

    output:
    path "genoflu_results.tsv", emit: merged_results
    path "versions.yml",        emit: versions

    script:
    """
    if [ -z "${tsv_files}" ]; then
        echo "Error: No TSV files provided to merge"
        exit 1
    fi

    head -n 1 \$(ls -1 *_stats.tsv | head -1) > genoflu_results.tsv
    
    for f in *_stats.tsv; do
        tail -n +2 \$f >> genoflu_results.tsv
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genoflu: \$(genoflu.py --version 2>&1)
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch genoflu_results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genoflu: \$(genoflu.py --version 2>&1)
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}