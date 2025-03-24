#!/usr/bin/env python3

import pandas as pd
import os
import glob
import argparse

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='Merge GenofLu results files.')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Name of the output merged TSV file')
    return parser.parse_args()

def merge_genoflu_results(output_file):
    # Using glob to find all individual genoflu run files
    stats_files = [f for f in glob.glob("*_stats.tsv") if os.path.isfile(f) and os.path.getsize(f) > 0]
    new_results_list = []
    
    for stats_file in stats_files:
        df = pd.read_csv(stats_file, sep='\t')
        new_results_list.append(df)
    
    if new_results_list:
        new_results = pd.concat(new_results_list)
        
        # If genoflu results file exists
        if os.path.isfile(output_file):
            existing_results = pd.read_csv(output_file, sep='\t')
            
            # Use the sample ID column (first column) as the identifier for duplicates
            sample_id_col = existing_results.columns[0]
            unique_new_results = new_results[~new_results[sample_id_col].isin(existing_results[sample_id_col])]
            combined_results = pd.concat([existing_results, unique_new_results])
        else:
            # If no existing file, use just the new results
            combined_results = new_results
        
        combined_results.to_csv(output_file, sep='\t', index=False)

def main():
    args = parse_args()
    merge_genoflu_results(args.output_file)

if __name__ == "__main__":
    main()
