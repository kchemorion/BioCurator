import sqlite3
from collections import defaultdict
from Bio import Entrez
from Bio import SeqIO
from io import StringIO

# your code using SeqIO here

import sys, os
import nltk
sys.modules['nltk'] = nltk
from nltk.tokenize import sent_tokenize, word_tokenize
from nltk.corpus import stopwords

# Set API key
Entrez.email = "francis.chemorion@insilicotrials.com"
Entrez.api_key = "2c16d55f288b716a61602429308eeb1bd208 "

# Define search term
search_term = "intervertebral disc degeneration"

#create database if it does not exist:folder_path = "path/to/folder"
folder_path = "."

if not os.path.exists(folder_path):
    os.makedirs(folder_path)

db_path = os.path.join(folder_path, "database.db")
conn = sqlite3.connect(db_path)

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

# Define a function to create tables
def create_tables(cursor):
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
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS protein (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            protein_id VARCHAR(255) UNIQUE,
            title TEXT,
            abstract TEXT
        )
        """
    )
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS nuccore (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            nuccore_id VARCHAR(255) UNIQUE,
            title TEXT,
            abstract TEXT
        )
        """
    )
    cursor.execute(
        """
        CREATE TABLE IF NOT EXISTS nucgss (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            nucgss_id VARCHAR(255) UNIQUE,
            title TEXT,
            abstract TEXT
        )
        """
    )
    conn.commit()
    print("Table 'nucgss' created successfully.")

    # Retrieve data from Nucleotide GSS database
    nucgss_handle = Entrez.esearch(db="nucgss", term=search_term, retmax=max_records, retstart=0, idtype="acc")
    nucgss_record = Entrez.read(nucgss_handle)
    nucgss_handle.close()
    nucgss_id_list = nucgss_record["IdList"]
    nucgss_total_count = int(nucgss_record["Count"])
    print(f"Total {nucgss_total_count} Nucleotide GSS records found.")
    print("Retrieving data from Nucleotide GSS database...")

    # Retrieve data from Nucleotide GSS database and save to SQLite database
    batch_size = 500
    for i in range(0, nucgss_total_count, batch_size):
        nucgss_batch_handle = Entrez.efetch(db="nucgss", id=nucgss_id_list[i : i + batch_size], rettype="fasta", retmode="text")
        nucgss_batch_records = nucgss_batch_handle.read()
        nucgss_batch_handle.close()

        for record in SeqIO.parse(StringIO(nucgss_batch_records), "fasta"):
            nucgss_id = record.id
            title = record.description
            abstract = str(record.seq)
            cursor.execute("INSERT OR IGNORE INTO nucgss (nucgss_id, title, abstract) VALUES (?, ?, ?)", (nucgss_id, title, abstract))
        conn.commit()
        print(f"{len(nucgss_id_list[:i+batch_size])} records retrieved from Nucleotide GSS database.")
