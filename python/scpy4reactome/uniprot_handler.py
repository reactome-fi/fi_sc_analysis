"""This script is used to pull out evidence from UniProt functional annotation for proteins and genes to support
the inferrered ligands and TFs in the regulatory networks.
"""

import time
import re
from langchain_core.documents import Document
from langchain_community.vectorstores.faiss import FAISS
from langchain_community.embeddings.huggingface import HuggingFaceEmbeddings
import pandas as pd

UNIPROT_DIR = '/Volumes/ssd/datasets_fi_rf/UniProt/release_2024_3'
HUMAN_UNIPROT_FILE = '{}/{}'.format(UNIPROT_DIR, 'uniprot_sprot_human.dat')
MOUSE_UNIPROT_FILE = '{}/{}'.format(UNIPROT_DIR, 'uniprot_sprot_rodents.dat')
# Used to embedding
EMBEDDING_MODEL_NAME = 'pritamdeka/S-PubMedBert-MS-MARCO'

def check_embedding_similarity(gene: str|list,
                               query_text: str,
                               gene_2_functions: dict,
                               k=10,
                               fetch_k=50) -> pd.DataFrame:
    if isinstance(gene, str):
        gene = [gene]
    all_gene_docs = []
    for current_gene in gene:
        if current_gene not in gene_2_functions.keys():
            continue
        gene_docs = get_function_docs(current_gene, gene_2_functions)
        all_gene_docs.extend(gene_docs)
    # Build the vector database first
    embeddings = HuggingFaceEmbeddings(model_name=EMBEDDING_MODEL_NAME)
    vector_db = FAISS.from_documents(all_gene_docs, embedding=embeddings)
    # Let find the documents that are most likely what we want
    if fetch_k < k:
        fetch_k = k + 10
    matched_doc_with_score = vector_db.similarity_search_with_score(query_text, k=k, fetch_k=fetch_k)
    # Convert it as a dataframe
    df = pd.DataFrame(columns=['gene', 'score(lower better)', 'annotation'])
    row = 0
    for doc_with_score in matched_doc_with_score:
        df.loc[row] = [doc_with_score[0].metadata['gene'],
               doc_with_score[1], 
               doc_with_score[0].page_content]
        row = row + 1
    return df


def load_gene_2_function(species: str = ['human', 'mouse']) -> dict:
    file = None
    if species == 'human':
        file = HUMAN_UNIPROT_FILE
    elif species == 'mouse':
        file = MOUSE_UNIPROT_FILE
    if file is None:
        raise ValueError('{} is not supported: '.format(species))
    gene_2_function = {}
    current_gene = None
    current_functions = None
    current_species = None
    is_in_function = False
    with open(file, 'r') as f:
        for line in f:
            line = line.strip()
            # print(line)
            if line.startswith('//'): # reset
                if current_species == species and current_gene and current_functions:
                    gene_2_function[current_gene] = current_functions
                current_gene = None
                current_functions = None
                # Just in case
                is_in_function = False
            if line.startswith('GN   Name='): # Get the gene name
                index1 = line.index('=')
                line = line[index1+1: ]
                current_gene = re.split(';| ', line)[0]
                # print('Current gene: {}'.format(current_gene))
            if line.startswith('OS'):
                line = line.lower()
                if line.endswith('(mouse).'):
                    current_species = 'mouse'
                elif line.endswith('(human).'):
                    current_species = 'human'
            if line.startswith('CC   -!- FUNCTION'):
                is_in_function = True
                index = line.index(':')
                line = line[index + 1:].strip()
                current_functions = line
            elif line.startswith('CC') and is_in_function:
                if line.startswith('CC   -!-'): # Switch to another CC section
                    is_in_function = False 
                else:
                    line = line[2:].strip()
                    current_functions = current_functions + ' ' + line
            # if len(gene_2_function) == 10:
            #     break

    return gene_2_function


def get_function_docs(gene: str,
                      gene_2_functions: dict) -> list:
    if gene not in gene_2_functions.keys():
        raise ValueError('{} has not function annotation.'.format(gene))
    functions = gene_2_functions[gene]
    tokens = _splitFunctionText(functions)
    documents = []
    for token in tokens:
        metadata = {
            'gene': gene
        }
        doc = Document(page_content=token,
                       metadata=metadata)
        documents.append(doc)
    return documents


def _splitFunctionText(text: str) -> list:
    result = re.split(r'(\(PubMed:\d+\).|\(By similarity\).)', text)
    tokens = [result[i] + result[i+1] for i in range(0, len(result), 2) if i+1 < len(result)]
    tokens.append(result[len(result) - 1])
    # Perform some merging for ECO
    rtn = []
    next_used = False
    for i in range(0, len(tokens) - 1):
        current = tokens[i]
        next = tokens[i + 1]
        next_used = False
        matched = re.match(r'\{ECO:.+?\}.', next.strip())
        if matched:
            current = current + next
            i = i + 1 # Escape 1
            next_used = True
        rtn.append(current)
        if i == len(result):
            break
    if not next_used:
        rtn.append(tokens[len(tokens) - 1])
    return rtn


if __name__ == '__main__':
    time1 = time.time()
    gene2function = load_gene_2_function('human')
    print(len(gene2function))
    genes = ['FN1', 'FGF9', 'MYC', 'DLL1']
    query_string = 'intestinal development'
    # query_string = 'cellular development'

    docs_with_scores = check_embedding_similarity(genes, query_string, gene2function)
    for doc_with_score in docs_with_scores:
        print(doc_with_score)
    time2 = time.time()
    print('Total time: {} seconds'.format(time2 - time1))

