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

FIELDS = {
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
            'find' : './/LIBRARY_LAYOUT',
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


def read_xml_data(xml_data):
# Parse XML metadata
    try:
        root = ET.fromstring(xml_data)
    except ET.ParseError:
        print('Error parsing XML, retrying...')
        time.sleep(10)
        root = ET.fromstring(xml_data)

    allDictVals = {}

    for root0 in tqdm(root.findall('.//EXPERIMENT_PACKAGE')):
        # Initialize a dictionary to store the extracted data
        try:
            data = {field: '' for field in FIELDS}
            for field in FIELDS:
                if isinstance(FIELDS[field], dict):
                    try:
                        if 'DATASTORE' in field:
                            value = []
                            for child in root0.findall(FIELDS[field]['find']):
                                val = child.attrib[FIELDS[field]['attrib']]
                                if val:
                                    value.append(val)
                            data[field] = ','.join(list(set(value)))
                            continue
                        if 'attrib' in FIELDS[field]:
                            data[field] = root0.find(FIELDS[field]['find']).attrib[FIELDS[field]['attrib']]
                        elif 'text' in FIELDS[field] and FIELDS[field]['text']:
                            data[field] = root0.find(FIELDS[field]['find']).text
                        elif 'text' in FIELDS[field] and not FIELDS[field]['text']:
                            data[field] = root0.find(FIELDS[field]['find'])[0].tag
                    except AttributeError:
                        data[field] = ''
                else:
                    data[field] = ''
        except AttributeError:
            print(f"Error parsing data for {root0.__dict__}")

        allDictVals[data['Run']] = {k: v for k, v in data.items() if k != 'Run'}        

    metadata = pd.DataFrame.from_dict(allDictVals, orient='index', columns=[k for k in FIELDS.keys() if k != 'Run'])
    metadata.reset_index(inplace=True)
    metadata.rename(columns={'index': 'Run'}, inplace=True)
    return metadata

def get_new_srps(search_term, email):
    """Fetch SRA metadata for a BioProject."""
    Entrez.email = email
    retstart = 0
    retmax = 1000
    all_ids = []
    print(f"Searching for SRA submissions with the term: {search_term}")
    while True:
        print(f"Fetching records {retstart} to {retstart + retmax}")
        retries = 3
        for attempt in range(retries):
            try:
                with Entrez.esearch(db="sra", idtype='acc', retstart=retstart, retmax=retmax, sort='recently_added', term=search_term) as handle:
                    record = Entrez.read(handle)
                break
            except (urllib.error.HTTPError, http.client.IncompleteRead) as e:
                print(f"Error fetching data for records {retstart} to {retstart + retmax}: {e}, retrying...")
                time.sleep(10)
                if attempt == retries - 1:
                    raise
        if not record['IdList']:
            break
        all_ids.extend(record['IdList'])
        retstart += retmax

    print(f"Found {len(all_ids)} SRA submissions.")

    # Fetch metadata in batches
    metadata_batches = []
    for start in range(0, len(all_ids), retmax):
        batch_ids = all_ids[start:start+retmax]
        print(f"Fetching metadata for batch starting at {start} of size {len(batch_ids)}")
        retries = 2
        for attempt in range(retries):
            try:
                with Entrez.efetch(db="sra", id=batch_ids, rettype="full", retmode='text') as handle:
                    batch_meta = handle.read().decode('UTF-8')
                break
            except (urllib.error.HTTPError, http.client.IncompleteRead) as e:
                print(f"Error fetching data for batch starting at {start}: {e}, retrying...")
                time.sleep(10)
                if attempt == retries - 1:
                    raise
        if batch_meta.strip():
            # Optionally, keep a local copy of the XML batch
            with open(f'batch_meta_{start}.xml', 'w') as f:
                f.write(batch_meta)
            metadata_batches.append(read_xml_data(batch_meta))
        else:
            break

    metadata = pd.concat(metadata_batches, ignore_index=True) if metadata_batches else pd.DataFrame()
    return metadata

def check_retracted_runs(old_metadata, new_metadata):
    """Check for retracted SRA runs."""
    retracted_runs = old_metadata[~old_metadata['Run'].isin(new_metadata['Run'])]['Run']
    return retracted_runs

def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Fetch SRA metadata for a BioProject.')
    parser.add_argument('-b', '--bioproject_ids', type=str, required=True, help='BioProject ID to monitor. Multiple IDs should be separated by commas.')
    parser.add_argument('-e', '--email', type=str, required=True, help='Email address for Entrez')
    parser.add_argument('-m', '--metadata', type=str, required=True, help='Path to old metadata file')
    parser.add_argument('-t', '--trimming_config', type=str, required=False, help='Path to trimming yaml file')
    parser.add_argument('-r', '--check_retracted', action='store_true', help='Check for retracted SRA runs')
    return parser.parse_args()

def main():
    args = parse_args()

    # Process BioProject IDs
    bioproject_ids = [bid.strip() for bid in args.bioproject_ids.split(',') if bid.strip()]
    search_term = " OR ".join(f"{bid}[BioProject]" for bid in bioproject_ids)

    new_metadata = get_new_srps(search_term, args.email)
    prev_metadata = pd.read_csv(args.metadata)

    new_sras = new_metadata.loc[~new_metadata['Run'].isin(prev_metadata['Run'])]
    combined_metadata = prev_metadata.copy()

    if new_sras.empty:
        print("No new SRA runs found.")
    else:
        print(f"Found {len(new_sras)} new SRA run(s).")
        combined_metadata = pd.concat([prev_metadata, new_sras], ignore_index=True)
        new_sras = new_sras.assign(
            is_milk=new_sras['isolation_source'].str.contains('milk', case=False, na=False),
            process_flag=True
        )
        save_columns = ['Run', 'process_flag', 'is_milk']

        if args.trimming_config:
            with open(args.trimming_config, 'r') as f:
                trimming_config = yaml.safe_load(f)
            new_sras['global_trimming'] = new_sras.apply(
                lambda row: trimming_config.get(row['BioProject'], trimming_config['default']).get('global_trimming')
                            if 'amp' in str(row['Assay Type']).lower() else None,
                axis=1
            )
            save_columns.append('global_trimming')

        new_sras[save_columns].to_csv(
            args.metadata.replace('.csv', '_to_process.tsv'),
            index=False,
            sep='\t'
        )

    if args.check_retracted:
        # Check for retracted runs and update metadata accordingly
        retracted_runs = check_retracted_runs(prev_metadata, new_metadata)
        combined_metadata['is_retracted'] = combined_metadata['Run'].isin(retracted_runs)

        # Ensure the retraction detection date column exists
        if 'retraction_detection_date_utc' not in combined_metadata.columns:
            combined_metadata['retraction_detection_date_utc'] = pd.NaT

        # Update detection date for newly retracted runs
        now_str = pd.Timestamp.now(tz='UTC').strftime('%Y-%m-%d %H:%M:%S')
        mask = combined_metadata['is_retracted'] & combined_metadata['retraction_detection_date_utc'].isna()
        combined_metadata.loc[mask, 'retraction_detection_date_utc'] = now_str

    if combined_metadata.equals(prev_metadata):
        print("No new updates found.")
    else:
        updated_path = args.metadata.replace('.csv', '_updated.csv')
        combined_metadata.to_csv(updated_path, index=False)
        print(f"Updated metadata saved to {updated_path}")


if __name__ == "__main__":
    main()
