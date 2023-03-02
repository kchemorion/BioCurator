import os
import requests
import xml.etree.ElementTree as ET
import time
from pathlib import Path

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

    # Parse the response XML to extract the article IDs
    root = ET.fromstring(search_response.text)
    ids = [elem.text for elem in root.iter("Id")]

    # Define the PubMed API base URL for fetching article information
    fetch_base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

    # Initialize the progress counter
    progress = 0

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

        # Save the article information to a file
        file_path = os.path.join(BASE_DIR, 'pubmed_articles.xml')
        with open(file_path, 'ab') as f:
            f.write(fetch_response.content)
            print("Downloaded {} - {}".format(start, end))

        # Wait for a few seconds before sending the next request
        time.sleep(3)
        start = end

        # Increment the progress counter
        progress += chunk_size

    print("Download complete! {} articles downloaded.".format(progress))

if __name__ == "__main__":
    query = "intervertebral disc"
    num_articles = 100
    download_articles(query, num_articles)
