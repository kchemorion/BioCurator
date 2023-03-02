import requests
import sqlite3

# Dictionary containing the API endpoints for each source
api_endpoints = {
    'The NCBI ClinVar': 'https://www.ncbi.nlm.nih.gov/clinvar/',
    'The Cancer Imaging Archive (TCIA)': 'https://services.cancerimagingarchive.net/services/v4/',
    'The Cancer Genome Interpreter': 'https://www.cancergenomeinterpreter.org/api/v1/',
    'The Clinical Proteomic Tumor Analysis Consortium (CPTAC)': 'https://api.gdc.cancer.gov/',
    'The Database of Single Nucleotide Polymorphisms (dbSNP)': 'https://api.ncbi.nlm.nih.gov/variation/v0/',
    'The Human Genome Diversity Project (HGDP)': 'https://api.nhs.uk/',
    'The Human Microbiome Project (HMP)': 'https://api.ncbi.nlm.nih.gov/microbiome/v1/',
    'The Human Phenotype Ontology (HPO)': 'https://hpo.jax.org/api/hpo/search/',
    'The International Cancer Genome Consortium (ICGC)': 'https://dcc.icgc.org/api/v1/',
    'The International HapMap Project': 'https://www.ebi.ac.uk/gwas/api/search/',
    'The International Mouse Phenotyping Consortium (IMPC)': 'https://www.mousephenotype.org/api/',
    'The International Parkinson Disease Genomics Consortium (IPDGC)': 'https://www.parkinsonsroadmap.org/api/v1/',
    'The Mouse Genome Informatics (MGI) database': 'http://api.mousephenotype.org/webservices/mpws/',
    'The Neuroimaging Informatics Tools and Resources Clearinghouse (NITRC)': 'https://www.nitrc.org/rest/api/',
    'The Protein Data Bank (PDB)': 'https://data.rcsb.org/graphql',
    'The Schizophrenia Research Forum (SRF)': 'https://www.schizophreniaforum.org/api/',
    'The Simons Foundation Autism Research Initiative (SFARI)': 'https://api.sfari.org/',
    'The Structural Genomics Consortium (SGC)': 'https://www.thesgc.org/api/',
    'The United Kingdom Biobank (UKB)': 'https://biobank.ndph.ox.ac.uk/api/',
    'The Wellcome Sanger Institute': 'https://api.sanger.ac.uk/',
}

# Connect to SQLite database
conn = sqlite3.connect('data.db')
c = conn.cursor()

# Create table to store API records if it doesn't exist
c.execute('''CREATE TABLE IF NOT EXISTS api_records
             (source text, record text)''')

# Extract data from each source and insert it into the SQLite database
for source, endpoint in api_endpoints.items():
    print(f"Extracting records from {source}...")
    try:
        response = requests.get(endpoint)
        records = response.json()
        for record in records:
            # Insert the record into the SQLite database
            c.execute("INSERT INTO api_records (source, record) VALUES (?, ?)", (source, str(record)))
        print(f"{len(records)} records extracted from {source}.")
    except Exception as e:
        print(f"Error retrieving data from {source}: {e}")

# Commit the changes to the database and close the connection
conn.commit()
conn.close()

print("Extraction completed successfully!")
