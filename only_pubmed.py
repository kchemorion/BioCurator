import requests
import xml.etree.ElementTree as ET
import sqlite3
from collections import defaultdict
from nltk.tokenize import sent_tokenize, word_tokenize
from nltk.corpus import stopwords

# Define the search query and the number of articles to retrieve
query = "intervertebral disc degeneration"
num_articles = 100

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

# Define the Entrez email and API key
email = "francis.chemorion@insilicotrials.com"
api_key = "698506c26b829c306fde3933018c5d57dd08"

# Initialize the Entrez API
from Bio import Entrez
Entrez.email = email
Entrez.api_key = api_key

# Retrieve article data from PubMed and save to a SQLite database
conn = sqlite3.connect("pubmed_articles.db")
c = conn.cursor()
c.execute('''CREATE TABLE IF NOT EXISTS articles 
             (id INTEGER PRIMARY KEY, title TEXT, journal TEXT, abstract TEXT, full_text TEXT, authors TEXT, pub_types TEXT, tokens TEXT)''')

count = 0
for article_id in ids:
    # Build the API request URL for fetching article data
    fetch_params = {
        "db": "pubmed",
        "id": article_id,
        "retmode": "xml",
    }

    # Send the API request for fetching article data and retrieve the response
    handle = Entrez.efetch(db="pubmed", id=article_id, retmode="xml")

    # Parse the response XML to extract the article data
    record = Entrez.read(handle)
    if "PubmedArticle" in record:
        article = record["PubmedArticle"][0]["MedlineCitation"]["Article"]

    # Extract the article data
    for article in root.findall(".//PubmedArticle"):
        article_data = {}

        # Extract the article metadata
        article_id = article.findtext(".//PMID")
        article_data["title"] = article.findtext(".//ArticleTitle")
        article_data["journal"] = article.findtext(".//Journal/Title")
        article_data["pub_date"] = article.findtext(".//PubDate/MedlineDate")
        article_data["abstract"] = article.findtext(".//AbstractText")

        # Extract the article authors
        authors = article.findall(".//AuthorList/Author")
        author_names = []
        for author in authors:
            last_name = author.findtext(".//LastName")
            fore_name = author.findtext(".//ForeName")
            author_names.append(last_name + " " + fore_name)
        article_data["authors"] = ", ".join(author_names)

        # Extract the article keywords
        keywords = article.findall(".//KeywordList/Keyword")
        keyword_list = [keyword.text for keyword in keywords]
        article_data["keywords"] = ", ".join(keyword_list)

        # Extract the article text
        article_text = ""
        for text_elem in article.findall(".//AbstractText") + article.findall(".//ArticleTitle") + article.findall(".//AuthorList/Author/ForeName") + article.findall(".//AuthorList/Author/LastName") + article.findall(".//KeywordList/Keyword"):
            if text_elem is not None:
                article_text += text_elem.text.strip() + " "
        article_data["article_text"] = article_text.strip()

        # Save the article data to the database
        cursor.execute('''INSERT INTO articles (article_id, title, journal, pub_date, abstract, authors, keywords, article_text) 
                          VALUES (?, ?, ?, ?, ?, ?, ?, ?)''', (article_id, article_data["title"], article_data["journal"], article_data["pub_date"], article_data["abstract"], article_data["authors"], article_data["keywords"], article_data["article_text"]))
        count += 1

    # Commit the changes to the database
    conn.commit()

# Print the number of records inserted into the database
print(f"{count} records inserted into the database.")

# Print the tables in the database and the number of records in each table
cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
tables = cursor.fetchall()
for table in tables:
    table_name = table[0]
    cursor.execute("SELECT COUNT(*) FROM {};".format(table_name))
    num_rows = cursor.fetchone()[0]
    print("{}: {} rows".format(table_name, num_rows))

# Close the database connection
conn.close()
