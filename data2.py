import sqlite3
import requests

# Define the API endpoints for each source
api_endpoints = {
    'The Cancer Genome Interpreter': 'https://www.cancergenomeinterpreter.org/webservices/v4.0/predictions',
    'The Clinical Proteomic Tumor Analysis Consortium (CPTAC)': 'https://cptac-data-portal.georgetown.edu/api/study/InterpretationResource/run/c1e7e1cb-b5a3-11e8-bf9d-0a7c1eab007a',
    'The Human Phenotype Ontology (HPO)': 'https://hpo.jax.org/api/hpo/search/',
    'The International Cancer Genome Consortium (ICGC)': 'https://dcc.icgc.org/api/v1/specimen/search?filters=%7B%22projectCode%22:%7B%22is%22:%5B%22LIRI-JP%22%5D%7D,%22donor.primarySite%22:%7B%22is%22:%5B%22intestine%22%5D%7D%7D&size=1&sort=donor.ageAtDiagnosis:desc',
    'The National Library of Medicine (NLM) Sequence Read Archive (SRA)': 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=sra&term=intervertebral+disc+degeneration&retmax=10',
}

# Create a new SQLite database and connect to it
conn = sqlite3.connect('clinical_data.db')

# Create a new table to store the API records
c = conn.cursor()
c.execute('''CREATE TABLE api_records
             (source text, record text)''')

# Loop through each API endpoint and extract the records
for source, endpoint in api_endpoints.items():
    print(f"Extracting records from {source}...")
    try:
        # Make the API request and retrieve the records
        response = requests.get(endpoint)
        records = response.json()
        
        # If the source is the NLM Sequence Read Archive (SRA), extract the accession IDs and append them to the endpoint
        if 'SRA' in source:
            accessions = records['esearchresult']['idlist']
            endpoint = f'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=sra&id={",".join(accessions)}'
            response = requests.get(endpoint)
            records = response.text
        
        # Insert the records into the SQLite database
        c.execute("INSERT INTO api_records (source, record) VALUES (?, ?)", (source, str(records)))
        print(f"{len(records)} records extracted from {source}.")
    except Exception as e:
        print(f"Error retrieving data from {source}: {e}")

# Commit the changes to the database and close the connection
conn.commit()
conn.close()

print("Extraction completed successfully!")
