import os
import requests
import xml.etree.ElementTree as ET
import time
from pathlib import Path
import sqlite3

BASE_DIR = Path(__file__).resolve().parent

def download_articles(query, num_articles):
    # Define the chunk size for fetching article information
    chunk_size = 100

    # Define the PubMed API base URL for searching articles
    search_base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

    # Build the API request URL for searching articles
    search_params = {
        "db": "pubmed",
        "term": query,
        "retmax": num_articles,
        "usehistory": "y",
    }

    # Send the API request for searching articles and retrieve the response
    search_response = requests.get(search_base_url, params=search_params)

    # Check if the response contains an error message
    if 'Error' in search_response.text:
        print("Error: {}".format(search_response.text))
        return

    # Parse the response XML to extract the article IDs
    root = ET.fromstring(search_response.text)
    ids = [elem.text for elem in root.iter("Id")]

    # Define the PubMed API base URL for fetching article information
    fetch_base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    # Initialize the progress counter
    progress = 0

    # Connect to the SQLite database
    db_path = os.path.join(BASE_DIR, 'pubmed_articles.db')
    conn = sqlite3.connect(db_path)
    c = conn.cursor()

    # Create the table if it doesn't exist
    c.execute('''CREATE TABLE IF NOT EXISTS pubmed_articles
                 (id INTEGER PRIMARY KEY,
                  title TEXT,
                  abstract TEXT,
                  authors TEXT,
                  pub_date TEXT,
                  journal TEXT)''')

    # Loop over the article IDs in chunks
    start = 0
    while start < len(ids):
        end = start + chunk_size
        chunk_ids = ids[start:end]

        # Build the API request URL for fetching article information
        fetch_params = {
            "db": "pubmed",
            "id": ",".join(chunk_ids),
            "retmode": "xml",
            "rettype": "abstract,text",  # specify the type of information to retrieve
            "journal": "journal",
            "pub_date": "pub_date",
            "authors":"authors",
        }

        # Send the API request for fetching article information and retrieve the response
        fetch_response = requests.get(fetch_base_url, params=fetch_params)

        # Check if the response contains an error message
        if 'Error' in fetch_response.text:
            print("Error: {}".format(fetch_response.text))
            return

        # Parse the response XML to extract the article information
        root = ET.fromstring(fetch_response.text)
        for article in root.iter('PubmedArticle'):
            pmid = article.find('.//PMID').text
            title = article.find('.//ArticleTitle').text
            abstract = article.find('.//AbstractText').text
            authors_list = [elem.text for elem in article.iter('Author') if elem.text is  not None]
            authors = ", ".join(authors_list) if authors_list else None
            pub_date = article.find('.//PubDate').text
            journal = article.find('.//Title').text

        # Insert the article information into the database
        c.execute('''INSERT INTO pubmed_articles
                     (id, title, abstract, authors, pub_date, journal)
                     VALUES (?, ?, ?, ?, ?, ?)''',
                  (pmid, title, abstract, authors, pub_date, journal))

    # Commit the changes to the database
    conn.commit()

    print("Downloaded {} - {}".format(start, end))

    # Wait for a few seconds before sending the next request
    time.sleep(3)
    start = end

    # Increment the progress counter
    progress += chunk_size

    # Close the database connection
    conn.close()

    print("Download complete. Total number of articles downloaded: {}".format(len(ids)))

if __name__ == "__main__":
    query = "intervertebral disc"
    num_articles = 100
    download_articles(query, num_articles)