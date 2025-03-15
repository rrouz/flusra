#!/usr/bin/env python3

import os
import glob
import sys

def merge_genoflu_results():
    stats_files = [f for f in glob.glob("*_stats.tsv") if os.path.isfile(f) and os.path.getsize(f) > 0]
    
    results_file = "genoflu_results.tsv"
    
    output_dir = os.environ.get("OUTPUT_DIR", "./")
    existing_file_path = os.path.join(output_dir, "genoflu_results.tsv")
    
    if not stats_files:
        if os.path.isfile(existing_file_path):
            with open(existing_file_path, "r") as src, open(results_file, "w") as dst:
                dst.write(src.read())
        else:
            with open(results_file, "w") as f:
                f.write("sample_id\tgenomic_segment\tsubtype\tclade\tlineage\tpercent_identity\n")
        return
    
    if os.path.isfile(existing_file_path):
        existing_results = {}
        existing_ids = set()
        existing_header = ""
        
        with open(existing_file_path, "r") as f:
            existing_header = f.readline().strip()
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split("\t")
                    if parts:
                        sample_id = parts[0]
                        existing_results[sample_id] = line
                        existing_ids.add(sample_id)
        
        new_results = {}
        new_ids = set()
        new_header = ""
        
        with open(stats_files[0], "r") as f:
            new_header = f.readline().strip()
        
        for stats_file in stats_files:
            with open(stats_file, "r") as f:
                f.readline()
                for line in f:
                    line = line.strip()
                    if line:
                        parts = line.split("\t")
                        if parts:
                            sample_id = parts[0]
                            new_results[sample_id] = line
                            new_ids.add(sample_id)
        
        unique_new_ids = new_ids - existing_ids
        
        with open(results_file, "w") as f:
            f.write(existing_header or new_header)
            f.write("\n")
            
            for sample_id in existing_ids:
                f.write(existing_results[sample_id])
                f.write("\n")
            
            for sample_id in unique_new_ids:
                f.write(new_results[sample_id])
                f.write("\n")
    
    else:
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