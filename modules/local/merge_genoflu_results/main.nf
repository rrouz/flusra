process MERGE_GENOFLU_RESULTS {
    tag "merge_results"
    label 'process_medium'
    

    input:
    path(tsv_files) 

    output:
    path "genoflu_results.tsv", emit: merged_results
    path "versions.yml",        emit: versions

    script:
    """
    OUTPUT_DIR="${params.outdir}/genoflu"
    EXISTING_FILE="\${OUTPUT_DIR}/genoflu_results.tsv"
    
    if [ ! \$(ls -1 *_stats.tsv 2>/dev/null | wc -l) -gt 0 ]; then
        if [ -f "\${EXISTING_FILE}" ]; then
            cp "\${EXISTING_FILE}" ./genoflu_results.tsv
        else
            echo "sample_id\tgenomic_segment\tsubtype\tclade\tlineage\tpercent_identity" > genoflu_results.tsv
        fi
    elif [ -f "\${EXISTING_FILE}" ]; then
        cp "\${EXISTING_FILE}" ./existing_results.tsv
        
        head -n 1 \$(ls -1 *_stats.tsv | head -1) > header.tmp
        
        head -n 1 \$(ls -1 *_stats.tsv | head -1) > new_results.tmp
        for stats_file in *_stats.tsv; do
            if [ -s "\$stats_file" ]; then
                tail -n +2 "\$stats_file" >> new_results.tmp
            fi
        done
        
        tail -n +2 ./existing_results.tsv | cut -f1 | sort > existing_ids.tmp
        tail -n +2 new_results.tmp | cut -f1 | sort > new_ids.tmp
        
        comm -23 new_ids.tmp existing_ids.tmp > unique_new_ids.tmp
        
        if [ -s unique_new_ids.tmp ]; then
            cat header.tmp > unique_new_results.tmp
            
            while read sample_entry; do
                grep -P "^\${sample_entry}\\t" new_results.tmp >> unique_new_results.tmp
            done < unique_new_ids.tmp
            
            cat existing_results.tsv > genoflu_results.tsv
            tail -n +2 unique_new_results.tmp >> genoflu_results.tsv
        else
            cp existing_results.tsv genoflu_results.tsv
        fi
        
        rm -f *.tmp
    else
        head -n 1 \$(ls -1 *_stats.tsv | head -1) > genoflu_results.tsv
        
        for stats_file in *_stats.tsv; do
            if [ -s "\$stats_file" ]; then
                tail -n +2 "\$stats_file" >> genoflu_results.tsv
            fi
        done
    fi
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genoflu: \$(genoflu.py --version 2>&1 || echo "NA")
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    """
    touch genoflu_results.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        genoflu: \$(genoflu.py --version 2>&1 || echo "NA")
        python: \$(python3 --version | sed 's/Python //g')
    END_VERSIONS
    """
}
