#!/usr/bin/env python3

import pandas as pd
import os
import glob

def merge_genoflu_results():
    # Using glob to find all individual genoflu runs
    stats_files = [f for f in glob.glob("*_stats.tsv") if os.path.isfile(f) and os.path.getsize(f) > 0]
    
    results_file = "genoflu_results.tsv"
    
    output_dir = os.environ.get("OUTPUT_DIR", "./")
    existing_file_path = os.path.join(output_dir, results_file)
    
    # If genoflu results file exists
    if os.path.isfile(existing_file_path):
        with open(existing_file_path, "r") as f:
            existing_header = f.readline().strip()
            existing_lines = {}
            existing_ids = set()
            
            for line in f:
                line = line.rstrip('\n')
                if line:
                    parts = line.split('\t')
                    if parts:
                        sample_id = parts[0]
                        existing_lines[sample_id] = line
                        existing_ids.add(sample_id)
        
        new_results = {}
        new_ids = set()
        
        # Get header from the first genoflu file
        with open(stats_files[0], "r") as f:
            new_header = f.readline().strip()
        
        # Read all genoflu files
        for stats_file in stats_files:
            with open(stats_file, "r") as f:
                f.readline()  # Skip header
                for line in f:
                    line = line.strip()
                    if line:
                        parts = line.split('\t')
                        if parts and parts[0]:
                            sample_id = parts[0]
                            new_results[sample_id] = line
                            new_ids.add(sample_id)
        
        # Find new genoflu sample ids
        unique_new_ids = new_ids - existing_ids
        
        # Write combined results
        with open(results_file, "w") as f:
            f.write(existing_header + "\n")
            
            # Write existing lines
            for sample_id in existing_ids:
                f.write(existing_lines[sample_id] + "\n")
            
            # Append unique new lines
            for sample_id in unique_new_ids:
                f.write(new_results[sample_id] + "\n")
    
    else:
        # Create new genoflu results
        with open(results_file, "w") as out_f:
            with open(stats_files[0], "r") as f:
                header = f.readline().strip()
                out_f.write(header + "\n")
            
            for stats_file in stats_files:
                with open(stats_file, "r") as f:
                    f.readline()
                    for line in f:
                        line = line.strip()
                        if line:
                            out_f.write(line + "\n")

if __name__ == "__main__":
    merge_genoflu_results()