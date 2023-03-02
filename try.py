import requests
import xml.etree.ElementTree as ET
import time

# Define the search query and the number of articles to retrieve
query = "intervertebral disc degeneration"
num_articles = 1000

# Define the chunk size for fetching article summaries
chunk_size = 100

# Define the PubMed API base URL
base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

# Build the API request URL
params = {
    "db": "pubmed",
    "term": query,
    "retmax": num_articles,
    "usehistory": "y",
}

# Send the API request and retrieve the response
response = requests.get(base_url, params=params)

# Parse the response XML to extract the article IDs
root = ET.fromstring(response.text)
ids = [elem.text for elem in root.iter("Id")]

# Define the PubMed API base URL for fetching article summaries
fetch_base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Loop over the article IDs in chunks
start = 0
while start < len(ids):
    end = start + chunk_size
    chunk_ids = ids[start:end]

    # Build the API request URL for fetching article summaries
    fetch_params = {
        "db": "pubmed",
        "id": ",".join(chunk_ids),
        "retmode": "xml",
    }

    # Send the API request for fetching article summaries and retrieve the response
    fetch_response = requests.get(fetch_base_url, params=fetch_params)

    # Save the article summaries to a file
    with open("pubmed_articles.xml", "a") as f:
        f.write(fetch_response.text)
        print("Downloaded {} - {}".format(start, end))

    # Wait for a few seconds before sending the next request
    time.sleep(3)
    start = end