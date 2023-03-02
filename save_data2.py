import sqlite3
from collections import defaultdict
from Bio import Entrez
from Bio import SeqIO
from io import StringIO

# Set API key
Entrez.email = "francis.chemorion@insilicotrials.com"
Entrez.api_key = "2c16d55f288b716a61602429308eeb1bd208"

# Define search term
search_term = "intervertebral disc degeneration"

# Define list of sources to search
sources = [
    "pubmed",
    "pmc",
    "mesh",
    "books",
    "omim",
    "journals",
    "popset",
    "probe",
    "protein",
    "nuccore",
    "nucgss",
    "nucest",
    "structure",
    "taxonomy",
    "biocollections",
    "biosample",
    "gtr",
    "snp",
    "gap",
    "dbvar",
]

# Define a dictionary to keep track of row count per source
row_count = defaultdict(int)

#create database if it does not exist:
folder_path = "."
if not os.path.exists(folder_path):
    os.makedirs(folder_path)

db_path = os.path.join(folder_path, "database.db")
conn = sqlite3.connect(db_path)

# Define a function to create tables
def create_tables(cursor, source):
    if source == "pubmed":
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS pubmed (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                pmid VARCHAR(255) UNIQUE,
                title TEXT,
                abstract TEXT
            )
            """
        )
    elif source == "pmc":
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS pmc (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                pmcid VARCHAR(255) UNIQUE,
                title TEXT,
                abstract TEXT
            )
            """
        )
    elif source == "mesh":
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS mesh (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                mesh_id VARCHAR(255) UNIQUE,
                term TEXT,
                definition TEXT,
                parent_terms TEXT
            )
            """
        )
    elif source == "books":
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS books (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                book_id VARCHAR(255) UNIQUE,
                title TEXT,
                abstract TEXT
            )
            """
        )
    elif source == "omim":
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS omim (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                omim_id VARCHAR(255) UNIQUE,
                title TEXT,
                abstract TEXT
            )
            """
        )
    elif source == "journals":
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS journals (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                journal_id VARCHAR(255) UNIQUE,
                title TEXT,
                abstract TEXT
            )
            """
        )
    elif source == "popset":
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS popset (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                popset_id VARCHAR(255) UNIQUE,
                title TEXT,
                abstract TEXT
            )
            """
        )
    elif source == "probe":
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS probe (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                probe_id VARCHAR(255) UNIQUE,
                title TEXT,
                abstract TEXT
            )
            """
        )
    elif source == "protein":
        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS protein (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                protein_id VARCHAR(255)



# Retrieve data from DBVar database
dbvar_handle = Entrez.esearch(db="dbvar", term=search_term, retmax=max_records, retstart=0, idtype="acc")
dbvar_record = Entrez.read(dbvar_handle)
dbvar_handle.close()
dbvar_id_list = dbvar_record["IdList"]
dbvar_total_count = int(dbvar_record["Count"])
print(f"Total {dbvar_total_count} DBVar records found.")
print("Retrieving data from DBVar database...")

# Retrieve data from DBVar database and save to SQLite database
batch_size = 500
for i in range(0, dbvar_total_count, batch_size):
    dbvar_batch_handle = Entrez.efetch(db="dbvar", id=dbvar_id_list[i : i + batch_size], rettype="fasta", retmode="text")
    dbvar_batch_records = dbvar_batch_handle.read()
    dbvar_batch_handle.close()

    for record in SeqIO.parse(StringIO(dbvar_batch_records), "fasta"):
        dbvar_id = record.id
        title = record.description
        abstract = str(record.seq)
        cursor.execute("INSERT OR IGNORE INTO dbvar (dbvar_id, title, abstract) VALUES (?, ?, ?)", (dbvar_id, title, abstract))
    conn.commit()
    print(f"{len(dbvar_id_list[:i+batch_size])} records retrieved from DBVar database.")
print(f"Finished retrieving data from DBVar database.")

print(f"All tables created successfully.")
