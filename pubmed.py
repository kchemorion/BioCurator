import requests
from collections import defaultdict
import xml.etree.ElementTree as ET
from Bio import Entrez
from spacy.matcher import Matcher
from spacy.lang.en import English
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity
import networkx as nx
import numpy as np
import sqlite3
from gensim.models.doc2vec import TaggedDocument, Doc2Vec
from gensim.summarization import TextRank, summarize
from gensim.summarization.summarizer import TextRank



# Set up the connection to the database
conn = sqlite3.connect('ivdd_database.db')
cur = conn.cursor()

# Set up the PubMed API access
from configparser import ConfigParser
from os.path import expanduser

# Read the configuration file
config_file = expanduser("~/.ncbi.cfg")
config = ConfigParser()
config.read(config_file)

email = 'francis.chemorion@insilicotrials.com'
api_key = config.get('NCBI', 'api_key')

# Provide the email and API key
Entrez.email = email
Entrez.api_key = api_key


# Set up the SpaCy NLP pipeline
nlp = English()
nlp.add_pipe('sentencizer')


# Search for articles about intervertebral disc degeneration on PubMed
query = 'intervertebral disc degeneration'
handle = Entrez.esearch(db='pubmed', term=query, retmax=5000)
record = Entrez.read(handle)
id_list = record['IdList']

# Iterate over the article IDs and download each article
for pmid in id_list:
    # Download the MEDLINE record for the article
    handle = Entrez.efetch(db='pubmed', id=pmid, rettype='medline', retmode='text')
    record = handle.read()

    # Preprocess the text
    # Replace newlines and tabs with spaces
    record = record.replace('\n', ' ').replace('\t', ' ')
    # Remove any extra spaces
    record = ' '.join(record.split())

    # Parse the MEDLINE record to extract the relevant metadata
    medline_dict = Entrez.read(Entrez.elink(dbfrom='pubmed', id=pmid, linkname='pubmed_pubmed_refs'))

    try:
        pub_date = medline_dict[0]['PubDate']
    except KeyError:
        pub_date = 'N/A'

    if 'ArticleIds' in medline_dict and 'doi' in medline_dict['ArticleIds']:
        doi = medline_dict['ArticleIds']['doi']
    else:
        doi = None

    if 'MedlineJournalInfo' in medline_dict[0]:
        journal = medline_dict[0]['MedlineJournalInfo']['MedlineTA']
    else:
        journal = None

    if 'AuthorList' in medline_dict:
        author_list = ', '.join([author.get('LastName', '') + ' ' + author.get('Initials', '') for author in medline_dict['AuthorList']])
    else:
        author_list = ''

    if 'KeywordList' in medline_dict[0]:
        keywords = '; '.join(medline_dict[0]['KeywordList'])
    else:
        keywords = ''


    # Download the full text for the article
    fulltext_url = f'https://www.ncbi.nlm.nih.gov/pmc/utils/oa/oa.fcgi?id={pmid}'
    response = requests.get(fulltext_url)
    root = ET.fromstring(response.content)
    article_text = ''
    for child in root.iter():
        if child.tag == '{http://www.ncbi.nlm.nih.gov}sec':
            if child.attrib.get('sec-type', '') == 'methods' or child.attrib.get('sec-type', '') == 'results':
                for text in child.itertext():
                    article_text += text.strip() + ' '
    if not article_text:
        if 'Abstract' in medline_dict[0]:
            article_text = medline_dict[0]['Abstract']
        else:
            article_text = None


    # Perform text preprocessing
    if article_text:
        article_text = article_text.replace('\n', ' ').replace('\t', ' ')
    else:
        article_text = ""

    article_text = ' '.join(article_text.split())

    # Perform Named Entity Recognition
    doc = nlp(article_text)
    entities = []
    for ent in doc.ents:
        entities.append((ent.text, ent.label_, ent.start_char, ent.end_char))

  # Perform Ontology-based curation using standardized ontologies
hpo_matches = defaultdict(int)
do_matches = defaultdict(int)
for entity in entities:
    if entity[1] == 'DISEASE':
        if 'degeneration' in entity[0].lower():
            cur.execute("INSERT INTO diseases (pmid, disease_name) VALUES (?, ?)",
                        (pmid, entity[0]))
            conn.commit()
        # Check if the disease entity matches a term in the Disease Ontology
        # (https://disease-ontology.org/)
        do_query = f"{entity[0]} AND 'Homo sapiens'[Organism]"
        handle = Entrez.esearch(db='doid', term=do_query)
        record = Entrez.read(handle)
        if record['IdList']:
            do_id = record['IdList'][0]
            do_matches[do_id] += 1
            cur.execute("INSERT INTO disease_ontology_matches (pmid, do_id) VALUES (?, ?)",
                        (pmid, do_id))
            conn.commit()
    elif entity[1] == 'GENE_OR_GENE_PRODUCT':
        # Check if the gene entity matches a term in the Human Phenotype Ontology
        # (https://hpo.jax.org/)
        hpo_query = f"{entity[0]}[Gene/Protein] AND 'Homo sapiens'[Organism]"
        handle = Entrez.esearch(db='gene', term=hpo_query)
        record = Entrez.read(handle)
        if record['IdList']:
            gene_id = record['IdList'][0]
            hpo_matches[gene_id] += 1
            cur.execute("INSERT INTO hpo_matches (pmid, gene_id) VALUES (?, ?)",
                        (pmid, gene_id))
            conn.commit()


# Perform Relation Extraction
edges = []
matcher = Matcher(nlp.vocab)
pattern = [{'POS': 'NOUN'}, {'IS_PUNCT': True, 'OP': '?'}, {'POS': 'ADJ'}, {'LOWER': 'disc'}, {'LOWER': 'degeneration'}]
matcher.add('disc_degeneration', None, pattern)
matches = matcher(doc)
for match_id, start, end in matches:
    entity_text = doc[start:end].text.lower()
    edges.append(('disease', entity_text))

# Perform Text Classification
# Use a simple keyword matching approach to classify articles into categories
categories = defaultdict(int)
for entity in entities:
    if entity[1] == 'DISEASE':
        if 'degeneration' in entity[0].lower():
            categories['degeneration'] += 1
        if 'herniation' in entity[0].lower():
            categories['herniation'] += 1
        if 'prolapse' in entity[0].lower():
            categories['prolapse'] += 1
    elif entity[1] == 'GENE_OR_GENE_PRODUCT':
        if 'col2a1' in entity[0].lower():
            categories['col2a1'] += 1
        if 'acan' in entity[0].lower():
            categories['acan'] += 1
        if 'col9a1' in entity[0].lower():
            categories['col9a1'] += 1
        if 'col11a1' in entity[0].lower():
            categories['col11a1'] += 1
        if 'col1a1' in entity[0].lower():
            categories['col1a1'] += 1
    elif entity[1] == 'CHEMICAL':
        if 'hyaluronic acid' in entity[0].lower():
            categories['hyaluronic acid'] += 1
        if 'collagen' in entity[0].lower():
            categories['collagen'] += 1
        if 'proteoglycan' in entity[0].lower():
            categories['proteoglycan'] += 1
        if 'steroids' in entity[0].lower():
            categories['steroids'] += 1
        if 'nonsteroidal anti-inflammatory drugs' in entity[0].lower():
            categories['NSAIDs'] += 1
    elif entity[1] == 'ANATOMICAL_ENTITY':
        if 'disc' in entity[0].lower():
            categories['disc'] += 1
        if 'vertebra' in entity[0].lower():
            categories['vertebra'] += 1
    elif entity[1] == 'PATHWAY':
        if 'mapk' in entity[0].lower() or 'mitogen-activated protein kinase' in entity[0].lower():
            categories['MAPK pathway'] += 1
        if 'tnf' in entity[0].lower() or 'tumor necrosis factor' in entity[0].lower():
            categories['TNF pathway'] += 1
# Insert the categories into the database
for category, count in categories.items():
    cur.execute("INSERT INTO categories (pmid, category, count) VALUES (?, ?, ?)",
                (pmid, category, count))

# Perform Text Similarity Matching
# Use Doc2Vec to identify similar articles based on their content
# Note: This may take a while to run
corpus = []
for pmid in id_list:
    handle = Entrez.efetch(db='pubmed', id=pmid, rettype='medline', retmode='text')
    record = handle.read()
    record = record.replace('\n', ' ').replace('\t', ' ')
    record = ' '.join(record.split())
    corpus.append(record)

# Train a Doc2Vec model on the corpus
documents = [TaggedDocument(doc, [i]) for i, doc in enumerate(corpus)]
model = Doc2Vec(documents, vector_size=100, window=5, min_count=5, workers=4, epochs=20)

# Find similar articles to each article in the corpus
similarities = defaultdict(list)
for i, pmid in enumerate(id_list):
    for j, sim in model.docvecs.most_similar(i):
        similar_pmid = id_list[j]
        similarities[pmid].append((similar_pmid, sim))

# Insert the similarities into the database
for pmid, similar_papers in similarities.items():
    for similar_pmid, similarity_score in similar_papers:
        cur.execute("INSERT INTO similarities (pmid1, pmid2, similarity_score) VALUES (?, ?, ?)",
                    (pmid, similar_pmid, similarity_score))

# Perform Text Summarization
# Use TextRank to generate a summary of each article
summaries = []
for pmid in id_list:
    handle = Entrez.efetch(db='pubmed', id=pmid, rettype='medline', retmode='text')
    record = handle.read()
    record = record.replace('\n', ' ').replace('\t', ' ')
    record = ' '.join(record.split())
    doc = nlp(record)
    sentences = [sent.string.strip() for sent in doc.sents]
    tr = TextRank(sentences)
    summary = tr.summarize()
    summaries.append(summary)

# Insert the summaries into the database
for pmid, summary in zip(id_list, summaries):
    cur.execute("INSERT INTO summaries (pmid, summary) VALUES (?, ?)",
                (pmid, summary))

# Commit changes to the database
conn.commit()

# Perform Event Extraction
events = []
for sentence_index, sentence in enumerate(doc.sents):
    for token_index, token in enumerate(sentence):
        if token.dep_ == 'nsubj' and token.head.pos_ == 'VERB':
            event = {'event_type': token.head.text, 'event_text': token.head.text.lower(),
                     'sentence_index': sentence_index, 'start_index': token.idx, 'end_index': token.idx + len(token.text)}
            events.append(event)

# Perform Text Similarity Matching
similarity_scores = []
for index, row in df.iterrows():
    similarity_score = compute_similarity(doc, nlp(row['text']))
    similarity_scores.append(similarity_score)
df['similarity_score'] = similarity_scores

# Perform Text Summarization
summary = summarize(doc.text)

# Insert the metadata into the database
cur.execute("INSERT INTO articles (pmid, title, abstract, journal, pub_date, doi, author_list, keywords) VALUES (?, ?, ?, ?, ?, ?, ?, ?)",
            (pmid, title, abstract, journal, pub_date, doi, author_list, keywords))
conn.commit()

# Insert the entities into the database
for entity in entities:
    cur.execute("INSERT INTO entities (pmid, entity_text, entity_type, start_index, end_index) VALUES (?, ?, ?, ?, ?)",
                (pmid, entity[0], entity[1], entity[2], entity[3]))
    conn.commit()

# Insert the Disease Ontology matches into the database
for do_id, count in do_matches.items():
    cur.execute("INSERT INTO disease_ontology_counts (pmid, do_id, count) VALUES (?, ?, ?)",
                (pmid, do_id, count))
    conn.commit()

# Insert the HPO matches into the database
for gene_id, count in hpo_matches.items():
    cur.execute("INSERT INTO hpo_counts (pmid, gene_id, count) VALUES (?, ?, ?)",
                (pmid, gene_id, count))
    conn.commit()

# Insert the edges into the database
for edge in edges:
    cur.execute("INSERT INTO edges (pmid, source, target) VALUES (?, ?, ?)",
                (pmid, 'disc', edge[1]))
    conn.commit()

# Insert the categories into the database
for category, count in categories.items():
    cur.execute("INSERT INTO categories (pmid, category, count) VALUES (?, ?, ?)",
                (pmid, category, count))
    conn.commit()

# Insert the events into the database
for event in events:
    cur.execute("INSERT INTO events (pmid, event_type, event_text, sentence_index, start_index, end_index) VALUES (?, ?, ?, ?, ?, ?)",
                (pmid, event['event_type'], event['event_text'], event['sentence_index'], event['start_index'], event['end_index']))
    conn.commit()

# Insert the similarity score into the database
cur.execute("INSERT INTO similarity_scores (pmid, score) VALUES (?, ?)",
            (pmid, df['similarity_score'].mean()))
conn.commit()

# Insert the summary into the database
cur.execute("INSERT INTO summaries (pmid, summary) VALUES (?, ?)",
            (pmid, summary))
conn.commit()


# Close the connection to the database
conn.close()
