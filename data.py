import requests
import pandas as pd
import sqlite3

# Define API endpoints for each source
sources = {
    'Osteoarthritis Initiative (OAI)': 'https://data-archive.nimh.nih.gov/oai/api/v1/',
    'Framingham Heart Study (FHS)': 'https://api.gbhi.org/fhs/',
    'National Health and Nutrition Examination Survey (NHANES)': 'https://wwwn.cdc.gov/Nchs/Nhanes/',
    'National Institute of Neurological Disorders and Stroke (NINDS) Common Data Elements (CDE) for Spinal Cord Injury': 'https://www.commondataelements.ninds.nih.gov/api/v1/',
    'Spine Patient Outcomes Research Trial (SPORT)': 'https://www.sporttrial.org/api/v1/',
    'National Institute of Arthritis and Musculoskeletal and Skin Diseases (NIAMS) Multidisciplinary Clinical Research Centers (MCRC) for Low Back Pain': 'https://www.niams.nih.gov/api/v1/',
    'Mayo Clinic Lumbosacral Plexus Injury and Herniated Lumbar Disk Study': 'https://www.mayo.edu/api/v1/',
    'Rotterdam Study': 'https://www.erasmus-epidemiology.nl/api/v1/',
    'Cochrane Back Review Group': 'https://www.cochrane.org/api/v1/',
    'North American Spine Society (NASS) Clinical Guidelines for Diagnosis and Treatment of Degenerative Lumbar Spinal Stenosis': 'https://www.spine.org/api/v1/',
    'Canadian Spine Outcomes and Research Network (CSORN)': 'https://www.csorn.org/api/v1/',
    'Health and Retirement Study (HRS)': 'https://hrs.isr.umich.edu/api/v1/',
    'Medical Expenditure Panel Survey (MEPS)': 'https://meps.ahrq.gov/mepsweb/api/v1/',
    'American Association of Neurological Surgeons (AANS) Spine Section Clinical Guidelines': 'https://www.aans.org/api/v1/',
    'International Society for the Study of the Lumbar Spine (ISSLS)': 'https://www.issls.org/api/v1/',
    'Back Pain Outcomes Using Longitudinal Data (BOLD) Project': 'https://www.boldproject.org/api/v1/',
    'International Collaboration on Repair Discoveries (ICORD) Spinal Cord Injury Registry': 'https://www.icord.org/api/v1/',
    'Interdisciplinary Collaborative Research Grant in Biomechanics and Motion Analysis': 'https://www.ubcicbr.ca/api/v1/',
    'The Genotype-Tissue Expression (GTEx) Project': 'https://gtexportal.org/rest/v1/',
    'The Cancer Genome Atlas (TCGA)': 'https://api.gdc.cancer.gov/',
    'Gene Expression Omnibus (GEO)': 'https://www.ncbi.nlm.nih.gov/geo/',
    'Human Protein Atlas': 'https://www.proteinatlas.org/api/',
    'Protein Data Bank (PDB)': 'https://www.rcsb.org/pdb/rest/',
    'UniProt': 'https://www.uniprot.org/',
    'The Human Genome Project': 'https://api.ncbi.nlm.nih.gov/datasets/v1alpha/',
    'ExAC Browser': 'http://exac.hms.harvard.edu/rest/',
    'Ensembl': 'https://rest.ensembl.org/',
    'PharmGKB': 'https://api.pharmgkb.org/v1/',
    'The Genotype-Tissue Expression (GTEx) Project': 'https://gtexportal.org/rest/v1/',
    'The Cancer Genome Atlas (TCGA)': 'https://api.gdc.cancer.gov/',
    'Gene Expression Omnibus (GEO)': 'https://www.ncbi.nlm.nih.gov/geo/',
    'Human Protein Atlas': 'https://www.proteinatlas.org/api/',
    'Protein Data Bank (PDB)': 'https://www.rcsb.org/pdb/rest/',
    'UniProt': 'https://www.uniprot.org/',
    'The Human Genome Project': 'https://api.ncbi.nlm.nih.gov/datasets/v1alpha/',
    'ExAC Browser': 'http://exac.hms.harvard.edu/rest/',
    'Ensembl': 'https://rest.ensembl.org/',
    'PharmGKB': 'https://api.pharmgkb.org/v1/',
    'PubChem': 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/',
    'ClinicalTrials.gov': 'https://clinicaltrials.gov/api/',
    'DrugBank': 'https://go.drugbank.com/',
    'Database of Genotypes and Phenotypes (dbGaP)': 'https://www.ncbi.nlm.nih.gov/gap/phegeni/api/',
    'Human Gene Mutation Database (HGMD)': 'https://portal.biobase-international.com/hgmd/pro/start.php/api',
    'National Center for Biotechnology Information (NCBI) Gene': 'https://www.ncbi.nlm.nih.gov/gene/',
    'The Human Protein Atlas - Tissue Atlas': 'https://www.proteinatlas.org/api/',
    'The Human Protein Atlas - Pathology Atlas': 'https://www.proteinatlas.org/api/',
    'The Human Protein Atlas - Cell Atlas': 'https://www.proteinatlas.org/api/'
}

# Extract 10 records of clinical data from each source
for source_name, api_endpoint in sources.items():
    try:
        # Make API request to retrieve clinical or biological data
        response = requests.get(api_endpoint)
        data = response.json()

        # Extract 10 records of data and convert to DataFrame
        records = data[:10]
        df = pd.DataFrame.from_records(records)

        # Store data in SQLite database
        with sqlite3.connect('data.db') as conn:
            df.to_sql(source_name, conn, if_exists='replace', index=False)
            print(f"Data from {source_name} stored in database")
    except Exception as e:
        print(f"Error retrieving data from {source_name}: {str(e)}")
