#!/usr/bin/env python3

import argparse
import time
import pandas as pd
import xml.etree.ElementTree as ET
from tqdm import tqdm
from Bio import Entrez
import urllib.error
import http.client
import yaml

def get_new_srps(search_term, email):
    """Fetch SRA metadata for a BioProject."""
    Entrez.email = email

    # Search for SRA submissions
    print(f"Searching for SRA submissions with the term: {search_term}")
    handle = Entrez.esearch(db="sra", idtype='acc', retmax=10000, sort='recently_added', term=search_term)
    record = Entrez.read(handle)
    handle.close()

    # Fetch metadata
    retries = 2
    for attempt in range(retries):
        try:
            handle = Entrez.efetch(db="sra", id=record['IdList'], rettype="full", retmode='text', retmax=10000)
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
    
    fields = {
        'Run' : {
            'find' : './/RUN',
            'attrib' : 'accession'
        },
        'Assay Type' : {
            'find' : './/LIBRARY_STRATEGY',
            'text' : True
        },
        'AvgSpotLen' : {
            'find' : './/Statistics/Read[@index="0"]',
            'attrib' : 'average'
        },
        'Bases' : {
            'find' : './/RUN_SET',
            'attrib' : 'bases'
        },
        'BioProject' : {
            'find' : './/EXTERNAL_ID[@namespace="BioProject"]',
            'text' : True
        },
        'BioSample' : {
            'find' : './/EXTERNAL_ID[@namespace="BioSample"]',
            'text' : True
        },
        'BioSample Accession' : {
            'find' : './/SAMPLE_DESCRIPTOR',
            'attrib' : 'accession'
        },
        'BioSampleModel' : {
            'find' : './/SAMPLE_ATTRIBUTE[TAG="BioSampleModel"]/VALUE',
            'text' : True
        },
        'Bytes' : {
            'find' : './/RUN',
            'attrib' : 'size'
        },
        'Center Name' : {
            'find' : './/SUBMISSION',
            'attrib' : 'center_name'
        },
        'Collection_Date' : {
            'find' : './/SAMPLE_ATTRIBUTE[TAG="collection_date"]/VALUE',
            'text' : True
        },
        'Consent' : {
            'find' : './/SRAFile',
            'attrib' : 'cluster'

        },
        'DATASTORE filetype' : {
            'find' : './/CloudFile',
            'attrib' : 'filetype'
        },
        'DATASTORE provider' : {
            'find' : './/CloudFile',
            'attrib' : 'provider'
        },
        'DATASTORE region' : {
            'find' : './/CloudFile',
            'attrib' : 'location'
        },
        'Experiment' : {
            'find' : './/EXPERIMENT',
            'attrib' : 'accession'
        },
        'geo_loc_name_country' : {
            'find' : './/SAMPLE_ATTRIBUTE[TAG="geo_loc_name"]/VALUE',
            'text' : True
        },
        'geo_loc_name_country_continent' : {
            'find' : './/SAMPLE_ATTRIBUTE[TAG="geo_loc_name_country_continent"]/VALUE',
            'default' : 'North America'
        },
        'geo_loc_name' : {
            'find' : './/SAMPLE_ATTRIBUTE[TAG="geo_loc_name"]/VALUE',
            'text' : True
        },
        'Host' : {
            'find' : './/SAMPLE_ATTRIBUTE[TAG="host"]/VALUE',
            'text' : True
        },
        'Instrument' : {
            'find' : './/INSTRUMENT_MODEL',
            'text' : True
        },
        'isolate' : {
            'find' : './/SAMPLE_ATTRIBUTE[TAG="isolate"]/VALUE',
            'text' : True
        },
        'Library Name' : {
            'find' : './/LIBRARY_NAME',
            'text' : True
        },
        'LibraryLayout' : {
            'find' : './/PAIRED',
            'text' : False
        },
        'LibrarySelection' : {
            'find' : './/LIBRARY_SELECTION',
            'text' : True
        },
        'LibrarySource' : {
            'find' : './/LIBRARY_SOURCE',
            'text' : True
        },
        'Organism' : {
            'find' : './/SAMPLE_NAME/SCIENTIFIC_NAME',
            'text' : True
        },
        'Platform' : {
            'find' : './/PLATFORM',
            'text' : False
        },
        'ReleaseDate' : {
            'find' : './/RUN',
            'attrib' : 'published'
        },
        'create_date' : {
            'find' : './/SRAFile',
            'attrib' : 'date'
        },
        'version' : {
            'find' : './/SRAFile',
            'attrib' : 'version'
        },
        'Sample Name' : {
            'find' : './/SAMPLE',
            'attrib' : 'alias'
        },
        'SRA Study' : {
            'find' : './/STUDY',
            'attrib' : 'accession'
        },
        'serotype' : {
            'find' : './/SAMPLE_ATTRIBUTE[TAG="serotype"]/VALUE',
            'text' : True
        },
        'isolation_source' : {
            'find' : './/SAMPLE_ATTRIBUTE[TAG="isolation_source"]/VALUE',
            'text' : True
        }
    }

    allDictVals = {}

    for root0 in tqdm(root.findall('.//EXPERIMENT_PACKAGE')):
        # Initialize a dictionary to store the extracted data
        try:
            data = {field: '' for field in fields}
            for field in fields:
                if isinstance(fields[field], dict):
                    try:
                        if 'DATASTORE' in field:
                            value = []
                            for child in root0.findall(fields[field]['find']):
                                val = child.attrib[fields[field]['attrib']]
                                if val:
                                    value.append(val)
                            data[field] = ','.join(list(set(value)))
                            continue
                        if 'attrib' in fields[field]:
                            data[field] = root0.find(fields[field]['find']).attrib[fields[field]['attrib']]
                        elif 'text' in fields[field] and fields[field]['text']:
                            data[field] = root0.find(fields[field]['find']).text
                        elif 'text' in fields[field] and not fields[field]['text']:
                            data[field] = root0.find(fields[field]['find']).tag
                    except AttributeError:
                        data[field] = ''
                else:
                    data[field] = ''
        except AttributeError:
            print(f"Error parsing data for {root0.__dict__}")

        allDictVals[data['Run']] = {k: v for k, v in data.items() if k != 'Run'}        

    metadata = pd.DataFrame.from_dict(allDictVals, orient='index', columns=[k for k in fields.keys() if k != 'Run'])
    metadata.reset_index(inplace=True)
    metadata.rename(columns={'index': 'Run'}, inplace=True)
    return metadata

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Fetch SRA metadata for a BioProject.')
    parser.add_argument('-b', '--bioproject_ids', type=str, required=True, help='BioProject ID to monitor. Multiple IDs should be separated by commas.')
    parser.add_argument('-e', '--email', type=str, required=True, help='Email address for Entrez')
    parser.add_argument('-m', '--metadata', type=str, required=True, help='Path to old metadata file')
    parser.add_argument('-t', '--trimming_config', type=str, required=False, help='Path to trimming yaml file')
    return parser.parse_args()

def main():
    args = parse_args()
     # Process BioProject IDs
    bioproject_ids = args.bioproject_ids.split(',')
    bioproject_ids = [id.strip() for id in bioproject_ids if id.strip()]
    if len(bioproject_ids) == 1:
        search_term = f"{bioproject_ids[0]}[BioProject]"
    else:
        search_term = f"{' OR '.join([f'{id}[BioProject]' for id in bioproject_ids])}"

    new_metadata = get_new_srps(search_term, args.email)
    
    prev_metadata = pd.read_csv(args.metadata)
    
    # Identify new SRA runs not in the previous metadata
    new_sras = new_metadata[~new_metadata['Run'].isin(prev_metadata['Run'])]
    
    if new_sras.empty:
        print("No new SRA runs found.")
    else:
        print(f"Found {len(new_sras)} new SRA run(s).")

        # Combine the old and new metadata
        combined_metadata = pd.concat([prev_metadata, new_sras])
        
        # Save updated metadata
        combined_metadata.to_csv(args.metadata.replace('.csv', '_updated.csv'), index=False)
        # Extract runs where 'isolation_source' contains "milk"
        new_sras['is_milk'] = new_sras['isolation_source'].str.contains('milk', case=False, na=False)

        new_sras['process_flag'] = True

        save_columns = ['Run', 'process_flag', 'is_milk']

        if args.trimming_config:
            with open(args.trimming_config, 'r') as f:
                trimming_config = yaml.safe_load(f)
            # if 'Assay Type' in new_sras contains 'amp' then check if project in trimming_config and set global trimming to that value
            new_sras['global_trimming'] = new_sras[['BioProject', 'Assay Type']].apply(
                lambda row: trimming_config.get(row['BioProject'], trimming_config['default']).get('global_trimming', None) if 'amp' in str(row['Assay Type']).lower() else None,
                axis=1
            )
            save_columns.append('global_trimming')
            
        new_sras[save_columns].to_csv(args.metadata.replace('.csv', '_to_process.tsv'), index=False, sep='\t')


if __name__ == "__main__":
    main()
