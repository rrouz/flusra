#!/usr/bin/env python3

import argparse
import time
import pandas as pd
import xml.etree.ElementTree as ET
from tqdm import tqdm
from Bio import Entrez
import urllib.error
import http.client

def get_new_srps(search_term, email):
    """Fetch SRA metadata for a BioProject."""
    Entrez.email = email

    # Search for SRA submissions
    handle = Entrez.esearch(db="sra", idtype='acc', retmax=4000, sort='recently_added', term=search_term)
    record = Entrez.read(handle)
    handle.close()

    # Fetch metadata
    retries = 2
    for attempt in range(retries):
        try:
            handle = Entrez.efetch(db="sra", id=record['IdList'], rettype="gb", retmode='text')
            returned_meta = handle.read().decode('UTF-8')
            handle.close()
            break
        except (urllib.error.HTTPError, http.client.IncompleteRead) as e:
            print(f"Error fetching data: {e}, retrying...")
            time.sleep(10)
            if attempt == retries - 1:
                raise

    # Write fetched metadata to a file for debugging purposes
    with open("NCBI_metadata.xml", "w") as f:
        f.write(returned_meta)

    # Parse XML metadata
    try:
        root = ET.fromstring(returned_meta)
    except ET.ParseError:
        print('Error parsing XML, retrying...')
        time.sleep(10)
        root = ET.fromstring(returned_meta)

    allDictVals = {}
    for root0 in tqdm(root.findall('.//EXPERIMENT_PACKAGE')):
        try:
            isolate = root0.find(".//SAMPLE").attrib['alias']
            collection_date = root0.find(".//SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='collection_date']/VALUE").text
            us_state = ''
            farm_id = ''
            run = root0.find(".//RUN").attrib['accession']
            bioproject = root0.find(".//STUDY_REF/IDENTIFIERS/EXTERNAL_ID").text
            biosample = root0.find(".//SAMPLE/IDENTIFIERS/EXTERNAL_ID").text
            center_name = root0.find(".//SUBMISSION").attrib.get('center_name', '')
            host = root0.find(".//SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE[TAG='host']/VALUE").text
            allDictVals[isolate] = {
                'Date': collection_date,
                'US State': us_state,
                'FarmID': farm_id,
                'Run': run,
                'BioProject': bioproject,
                'BioSample': biosample,
                'Center Name': center_name,
                'Host': host
            }
        except AttributeError as e:
            print(f"Error processing entry: {e}, skipping...")

    metadata = pd.DataFrame.from_dict(allDictVals, orient='index')
    metadata.reset_index(inplace=True)
    metadata.rename(columns={'index': 'isolate'}, inplace=True)
    return metadata

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Fetch SRA metadata for a BioProject.')
    parser.add_argument('-b', '--bioproject_id', type=str, required=True, help='BioProject ID to monitor')
    parser.add_argument('-e', '--email', type=str, required=True, help='Email address for Entrez')
    parser.add_argument('-m', '--metadata', type=str, required=True, help='Path to old metadata file')
    return parser.parse_args()

def main():
    args = parse_args()
    search_term = f"{args.bioproject_id}[BioProject]"
    new_metadata = get_new_srps(search_term, args.email)
    
    prev_metadata = pd.read_csv(args.metadata)
    
    # Identify new SRAs not present in the previous metadata
    new_sras = new_metadata[~new_metadata['Run'].isin(prev_metadata['Run'])]
    
    if new_sras.empty:
        print("No new SRA runs found.")
    else:
        print(f"Found {len(new_sras)} new SRA runs.")

        # Combine the old and new metadata
        combined_metadata = pd.concat([prev_metadata, new_sras])
        
        # Save the updated metadata
        combined_metadata.to_csv(args.metadata.replace('.csv', '_updated.csv'), index=False)
        new_sras.to_csv(args.metadata.replace('.csv', '_new.csv'), index=False)

if __name__ == "__main__":
    main()
