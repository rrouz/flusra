#!/usr/bin/env python3

import pandas as pd
import os
import glob

def merge_genoflu_results():
    # Using glob to find all individual genoflu run files
    stats_files = [f for f in glob.glob("*_stats.tsv") if os.path.isfile(f) and os.path.getsize(f) > 0]
    
    results_file = "genoflu_results.tsv"
    output_dir = os.environ.get("OUTPUT_DIR", "./")
    existing_file_path = os.path.join(output_dir, results_file)
    
    # If genoflu results file exists
    if os.path.isfile(existing_file_path):
        existing_results = pd.read_csv(existing_file_path, sep='\t')
        new_results_list = []
        for stats_file in stats_files:
            df = pd.read_csv(stats_file, sep='\t')
            new_results_list.append(df)

        if new_results_list:
            new_results = pd.concat(new_results_list)

            # Use the sample ID column (first column) as the identifier for duplicates
            sample_id_col = existing_results.columns[0]
            unique_new_results = new_results[~new_results[sample_id_col].isin(existing_results[sample_id_col])]
            combined_results = pd.concat([existing_results, unique_new_results])
            combined_results.to_csv(results_file, sep='\t', index=False)
    
    else:
        # If no existing genoflu results
        new_results_list = []
        for stats_file in stats_files:
            df = pd.read_csv(stats_file, sep='\t')
            new_results_list.append(df)
        
        if new_results_list:
            combined_results = pd.concat(new_results_list)
            combined_results.to_csv(results_file, sep='\t', index=False)

if __name__ == "__main__":
    merge_genoflu_results()