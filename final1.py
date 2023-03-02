import os
import requests
import xml.etree.ElementTree as ET
import time
import sqlite3
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

    # Check if the XML file already exists and delete it if it does
    file_path = os.path.join(BASE_DIR, 'pubmed_articles.xml')
    if os.path.exists(file_path):
        os.remove(file_path)

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
            "authors": "authors",
            "pubmedcentral": "pmc",
            "version": "version"
        }

        # Send the API request for fetching article information and retrieve the response
        fetch_response = requests.get(fetch_base_url, params=fetch_params)

        # Save the article information to a file
        with open(file_path, 'ab') as f:
            f.write(fetch_response.content)
            print("Downloaded {} - {}".format(start, end))

        # Wait for a few seconds before sending the next request
        time.sleep(3)
        start = end

        # Increment the progress counter
        progress += chunk_size

    print("Download complete! {} articles downloaded.".format(progress))
    
    # Create a new SQLite database and table to store article information
    conn = sqlite3.connect(os.path.join(BASE_DIR, 'articles.db'))
    c = conn.cursor()
    c.execute('''CREATE TABLE IF NOT EXISTS pubmed_articles
                 (id INTEGER PRIMARY KEY,
                  title TEXT,
                  abstract TEXT,
                  full_text TEXT,
                  authors TEXT,
                  pub_date TEXT,
                  journal TEXT)''')

    # Parse the downloaded article information and insert it into the database
    tree = ET.parse(file_path)
    root = tree.getroot()

    for article in root.iter("PubmedArticle"):
        # Extract the article information
        pmid = int(article.find('.//PMID').text)
        title = article.find('.//ArticleTitle').text
        abstract_elem = article.find('.//AbstractText')
        full_text_elem = article.find('.//FullTextArticle')
        abstract = abstract_elem.text if abstract_elem is not None else None
            # Extract the full article text
    article_text = ''
    for elem in article.iter():
        if elem.text is not None and elem.tag == 'AbstractText':
            article_text += elem.text.strip() + ' '
        elif elem.text is not None and elem.tag == 'p':
            article_text += elem.text.strip() + ' '
        elif elem.text is not None and elem.tag == 'sec':
            article_text += elem.text.strip() + ' '

    authors = ", ".join([elem.findtext('LastName') + ' ' + elem.findtext('Initials') if elem.findtext('LastName') is not None else '' for elem in article.iter('Author')])
    pub_date_elem = article.find('.//PubDate')
    if pub_date_elem is not None:
        pub_date_year = pub_date_elem.findtext('Year')
        pub_date_month = pub_date_elem.findtext('Month')
        pub_date_day = pub_date_elem.findtext('Day')
        pub_date = f"{pub_date_month} {pub_date_day}, {pub_date_year}" if pub_date_year is not None and pub_date_month is not None and pub_date_day is not None else None
    else:
        pub_date = None
    journal_elem = article.find('.//Journal')
    if journal_elem is not None:
        journal = journal_elem.findtext('Title')
    else:
        journal = None

    # Insert the article information into the database
    c.execute('''INSERT INTO pubmed_articles (id, title, abstract, article_text, authors, pub_date, journal)
                VALUES (?, ?, ?, ?, ?, ?, ?)''',
            (pmid, title, abstract, article_text, authors, pub_date, journal))
    conn.commit()
    conn.close()
print("Article information saved to database.")


if __name__ == "__main__":
    query = "intervertebral disc"
    num_articles = 100
    download_articles(query, num_articles)