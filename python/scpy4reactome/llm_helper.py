
import networkx as nx
from langchain_core.prompts import ChatPromptTemplate
from langchain_core.runnables import RunnablePassthrough
from pathlib import Path


def get_evaluate_relatedness_prompt_template():
    prompt_template = """
You are an expert in biology and biomedicine. Evaluate the relatedness of the following text to the biological concept "{concept}." 
Use a score between 0 and 10, with 10 indicating the highest degree of relatedness and 0 indicating no relatedness at all. Make sure
the score is returned as the first line with the format like this: Score: 5\n\n.

{text}
"""
    prompt = ChatPromptTemplate.from_template(prompt_template)
    return prompt


def evaluate_relatedness(concept: str,
                         text: str,
                         model: any) -> any:
    llm_parameters = {
        'concept': concept,
        'text': text
    }

    prompt = get_evaluate_relatedness_prompt_template()
    return invoke_llm(model, llm_parameters, prompt)


def invoke_llm(model, llm_parameters, prompt):
    answer = {
        'answer': prompt | model 
    }

    dummy = RunnablePassthrough()
    llm_chain = dummy | answer
    result = llm_chain.invoke(llm_parameters)
    return result


def get_short_descrption_of_term_prompt_template():
    prompt_template = """
You are an expert in biology and biomedicine. Write a short description about {term} with no more than {total_words} words. Focus on
molecular mechanisms.
"""
    prompt = ChatPromptTemplate.from_template(prompt_template)
    return prompt


def write_short_description_of_term(term: str,
                                    model: any,
                                    total_words: int=100) -> any:
    llm_parameters = {
        'term': term,
        'total_words': total_words
    }

    desc_prompt = get_short_descrption_of_term_prompt_template()

    return invoke_llm(model, llm_parameters, desc_prompt)


def get_llm_summary_prompt_template():
    network_summary_prompt_template = """
You are an expert in the field of {study_context}. Relying on the provided lists for reference, use your domain knowledge and write 
some paragraphs in a scholar way to highlight the most important pathways and their functional relationships with transcription 
factors and ligands with no more than {total_words} words. Be specific and don't speculate!

transcription factors: {transcriptaion_factors}

pathways: {pathways}

ligands: {ligands}

functional relationships: {functional_relationships}

"""
    network_summary_prompt = ChatPromptTemplate.from_template(network_summary_prompt_template)
    return network_summary_prompt


def summarize_network_via_llm(network: nx.DiGraph,
                              study_context: str,
                              model: any,
                              total_words: int = 500) -> any:
    tfs, pathways, ligands, fis = extract_information(network)
    # print('TFs: {}\nPathways: {}\nLigands: {}\nFIs: {}'.format(tfs, pathways, ligands, fis))
    llm_parameters = {
        'transcriptaion_factors': ';'.join(tfs),
        'pathways': ';'.join(pathways),
        'ligands': ';'.join(ligands),
        'functional_relationships': ';'.join(fis),
        'total_words': total_words
    }

    # Add a context about that this study is for
    llm_parameters['study_context'] = study_context
    network_summary_prompt = get_llm_summary_prompt_template()

    return invoke_llm(model, llm_parameters, network_summary_prompt)


def extract_information(network: nx.DiGraph | str) -> tuple[list, list, list, list]:
    """Extract transcription factors, pathways, ligands, and their relationships from the network.

    Args:
        network (nx.DiGraph): _description_

    Returns:
        tuple[list, list, list, list]: four lists in the order of tfs, pathways, ligands and relationships.
    """
    # Need to load the network first if it is a file
    if isinstance(network, str):
        path = Path(network)
        if not path.exists():
            raise ValueError('No graphml file exists at: {}'.format(network))
        network = nx.read_graphml(network)
    tfs = []
    pathways = []
    ligands = []
    fis = [] # functional relationships among nodes
    for node, data in network.nodes(data=True):
        type = data['type']
        if type == 'TF':
            tfs.append(node)
        elif type == 'Pathway':
            pathways.append(node)
        elif type == 'Ligand':
            ligands.append(node)
        else:
            raise ValueError('Unknown type: {}'.format(type))
    
    # Extract edge information
    # Mapping action for different edge types
    edge_annotation_actions = {
        'ligand_pathway_activation': lambda _ : "activates",
        'tf_pathway_inhibition': lambda _ : 'inhibits',
        'tf_tf_activation': lambda _ : 'activates',
        'tf_pathway_activation': lambda _ : 'activates',
        'pathway_tf_annotation': lambda _ : 'regulates',
        'pathway_pathway_hierarchy': lambda _ : 'is contained by',
        'ligand_tf_activation': lambda _ : 'activates',
    }
    for src, target, data in network.edges(data=True):
        annotation = data['annotation']
        edge_meaning = edge_annotation_actions[annotation]
        if edge_meaning is None:
            raise ValueError('No meaning defined for edge: {}'.format(annotation))
        fi = '"{}" {} "{}"'.format(src, edge_meaning(annotation), target) # The second element calls the function
        fis.append(fi)

    return tfs, pathways, ligands, fis