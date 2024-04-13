# This script wraps the most important functions from CellPhoneDB

import scanpy as sc
from pathlib import Path
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method
import ktplotspy as kpy
import pandas as pd
import networkx as nx
import numpy as np
from collections import defaultdict
import seaborn as sns


def run_cellphonedb(counts_h5ad_file: str,  # Since CellPhoneDB needs the original file, we cannot pass a AnnData.
                    cell_type_key: str,
                    sample: str,
                    cellphonedb_dir: str,
                    cpdb_file_path: str = '../resources/cellphonedb_07_18_2023.zip',
                    convert_from_human_to_mouse: bool = False,
                    cell_type_prefix: str = None,
                    iterations: int = 1000):
    """Here the default 1,000 samplings are used. This is different from Chris' analysis. He used a much
    larger sampling number, 1,000,000. Not sure if it is really necessary for this type of analysis.

    Args:
        sample (str): the sample value is used for subdirectory cellphonedb_dir.
        hybrid_clusters (list): _description_
        cellphonedb_dir (str): the analysis output directory
        cpdb_file_path (str): the cellphonedb resource file. This should be a zip file downloaded from CellPhoneDB
        base_dir (str): _description_
        cell_type_prefix: to address a bug in cellphonedb when cluster is used for cell_type_key
        convert_from_human_to_mouse: The passed data is for mouse. However, CellPhoneDB can be done for mouse only. 
        Use this flag to do a conversion. 

    Returns:
        _type_: _description_
    """
    # Make sure the directory is correct
    cellpohone_out_dir = Path(cellphonedb_dir, sample)
    if not cellpohone_out_dir.exists():
        cellpohone_out_dir.mkdir()
    cellphonedb_meta_file = Path(
        cellpohone_out_dir, 'meta_{}.csv'.format(sample))

    # Need to create a csv file for the following function
    adata = sc.read_h5ad(counts_h5ad_file)
    if convert_from_human_to_mouse:
        adata = convert_adata_from_mouse_to_human(adata)
        # Save the file
        adata_file_name = Path(cellpohone_out_dir, 'converted_human_data.h5ad')
        adata.write_h5ad(adata_file_name)
        counts_h5ad_file = adata_file_name

    # Do a transformation if cell_type_prefix is provided to address a bug in CellPhoneDB
    if cell_type_prefix:
        adata.obs[cell_type_key] = adata.obs[cell_type_key].map(
            lambda v: '{}{}'.format(cell_type_prefix, v))
    adata.obs.to_csv(cellphonedb_meta_file,
                     columns=[cell_type_key],
                     index_label='cell')

    # Run the random sampling based method. For the time being, using 1,000 permuations
    cpdb_results = cpdb_statistical_analysis_method.call(
        # mandatory: CellPhoneDB database zip file.
        cpdb_file_path=cpdb_file_path,
        # mandatory: tsv file defining barcodes to cell label.
        meta_file_path=cellphonedb_meta_file,
        # mandatory: normalized count matrix.
        counts_file_path=counts_h5ad_file,
        # defines the gene annotation in counts matrix.
        counts_data='hgnc_symbol',
        # optional (default: None): defines cells per microenvironment.
        microenvs_file_path=None,
        # denotes the number of shufflings performed in the analysis.
        iterations=iterations,
        # defines the min % of cells expressing a gene for this to be employed in the analysis.
        threshold=0.1,
        # number of threads to use in the analysis.
        threads=4,
        # debug randome seed. To disable >=0.
        debug_seed=42,
        # Sets the rounding for the mean values in significan_means.
        result_precision=3,
        # P-value threshold to employ for significance.
        pvalue=0.05,
        # To enable subsampling the data (geometri sketching).
        subsampling=False,
        # (mandatory) enable subsampling log1p for non log-transformed data inputs.
        subsampling_log=False,
        # Number of componets to subsample via geometric skectching (dafault: 100).
        subsampling_num_pc=100,
        # Number of cells to subsample (integer) (default: 1/3 of the dataset).
        subsampling_num_cells=1000,
        # Sets the string to employ to separate cells in the results dataframes "cellA|CellB".
        separator='|',
        # Saves all intermediate tables employed during the analysis in pkl format.
        debug=False,
        # Path to save results.
        output_path=cellpohone_out_dir,
        # Replaces the timestamp in the output files by a user defined string in the  (default: None).
        output_suffix=None,
        # Give us a score for the interactions
        score_interactions=True
    )

    # # Output some results
    # kpy.plot_cpdb_heatmap(
    #     adata = adata,
    #     pvals = pvalues,
    #     celltype_key = cellphonedb_celltype_key,
    #     figsize = (5,5),
    #     title = "{}: Number of significant interactions".format(sample),
    #     symmetrical = False,
    #     return_tables=False,
    # )
    # Results should be a dict
    return cpdb_results, adata


def plot_cpdb(adata: sc.AnnData,
              cell_type1: str,
              selected_celltypes: str,
              means: pd.DataFrame,
              pvalues: pd.DataFrame,
              cellphonedb_celltype_key: str,
              title: str,
              genes: str | list = None,
              need_table: bool = False,
              alpha: float = 0.05,
              figsize=(10, 5)):
    ggplot = kpy.plot_cpdb(
        adata=adata,
        cell_type1=cell_type1,
        cell_type2=selected_celltypes,
        means=means,
        pvals=pvalues,
        celltype_key=cellphonedb_celltype_key,
        genes=genes,
        figsize=figsize,
        title=title,
        max_size=6,
        alpha=alpha,
        highlight_size=0.75,
        standard_scale=True,
        special_character_regex_pattern='/',
        # gene_family="costimulatory"
    )
    print(ggplot)
    if need_table:
        # Print out the table
        cpdb_table = kpy.plot_cpdb(
            adata=adata,
            cell_type1=cell_type1,
            cell_type2=selected_celltypes,
            means=means,
            pvals=pvalues,
            celltype_key=cellphonedb_celltype_key,
            special_character_regex_pattern='/',
            standard_scale=True,
            alpha=alpha,
            return_table=True
        )
        # print_cpdb_table(cpdb_table)
        return cpdb_table


def build_cellphonedb_network(pvalues: pd.DataFrame,
                              means: pd.DataFrame,
                              interaction_scores: pd.DataFrame = None,
                              pvalue_cutoff: float = 0.01,
                              mean_cutoff: float = 0,
                              score_cutoff: float = 0,
                              result_start_index: int = 11,
                              save_graphml: str = None):
    """
    Build a networkx object from CellPhoneDB analysis results described by means and pvalues, two DataFrame objects.
    """
    # Make sure pvalues and means have the same shape
    if pvalues.shape != means.shape:
        raise ValueError('pvalues and means don\'t have the same shape: '.format(
            pvalues.shape, means.shape))
    network = nx.DiGraph()
    # Start to build the network by checking each pair in the two data frames
    for i in range(pvalues.shape[0]):
        pvalue_row = pvalues.iloc[i]
        mean_row = means.iloc[i]
        score_row = None
        if interaction_scores is not None:
            score_row = interaction_scores.iloc[i]
        # Extract receptors and ligands
        partner_a = pvalue_row['partner_a'].split(':')[1]
        partner_b = pvalue_row['partner_b'].split(':')[1]
        gene_a = pvalue_row['gene_a']
        if (isinstance(gene_a, float) and not np.isnan(gene_a)) or (isinstance(gene_a, str) and len(gene_a) > 0):  # Is is a float if it is empty
            partner_a = gene_a  # Means this is a gene and we should use gene name
        gene_b = pvalue_row['gene_b']
        if (isinstance(gene_b, float) and not np.isnan(gene_b)) or (isinstance(gene_b, str) and len(gene_b) > 0):
            partner_b = gene_b
        # The following will be used as edge attributes
        receptor_a = 'receptor' if pvalue_row['receptor_a'] else 'ligand'
        receptor_b = 'receptor' if pvalue_row['receptor_b'] else 'ligand'
        direction = receptor_a + '-' + receptor_b
        secreted = pvalue_row['secreted']
        is_integrin = pvalue_row['is_integrin']
        # Start to parsing individual value
        for j in range(result_start_index, pvalues.shape[1]):
            pvalue = pvalue_row[j]
            # print('{}: {}'.format(j, pvalue))
            if pvalue > pvalue_cutoff:
                continue
            mean = mean_row[j]
            if mean <= mean_cutoff:
                continue
            score = None
            if score_row is not None:
                score = score_row[j]
                if score <= score_cutoff:
                    continue
            header = pvalue_row.index[j]
            # Need to parse it into a pair of cell groups
            cell_groups = header.split('|')
            # Naming nodes as cell_group:partner
            src_node = '{}:{}'.format(cell_groups[0], partner_a)
            target_node = '{}:{}'.format(cell_groups[1], partner_b)
            # print('{}, {}, {}, {}, {}'.format(pvalue, mean, direction, secreted, is_integrin))
            if score is not None:
                network.add_edge(src_node, target_node,
                                pvalue=pvalue, 
                                mean=mean, 
                                score=score,
                                direction=direction, 
                                secreted=str(secreted), 
                                is_integrin=str(is_integrin))  # bool cannot be saved in graphml. Have to convert it into a string.
            else:
                network.add_edge(src_node, target_node,
                                pvalue=pvalue, 
                                mean=mean, 
                                score=score,
                                direction=direction, 
                                secreted=str(secreted), 
                                is_integrin=str(is_integrin))
            network.nodes[src_node]['cluster'] = cell_groups[0].split(
                '_')[1] if '_' in cell_groups[0] else cell_groups[0]
            network.nodes[target_node]['cluster'] = cell_groups[1].split(
                '_')[1] if '_' in cell_groups[1] else cell_groups[1]
            network.nodes[src_node]['entity'] = partner_a
            network.nodes[target_node]['entity'] = partner_b
    if save_graphml:
        print('Writing the network to {}.'.format(save_graphml))
        nx.write_graphml(network, save_graphml)
    return network


def convert_adata_from_mouse_to_human(adata: sc.AnnData) -> sc.AnnData:
    adata_copy = adata.copy()
    mousegene2human = load_mouse_genes_to_human_map()
    # Do a map for gene names from mouse to human

    def map_fun(mouse_gene):
        if mouse_gene in mousegene2human.keys():
            return mousegene2human[mouse_gene][0]
        return mouse_gene
    adata_copy.var_names = [map_fun(mouse_gene)
                            for mouse_gene in adata_copy.var_names]
    return adata_copy


def load_mouse_genes_to_human_map(map_file_name: str = '../resources/HOM_MouseHumanSequence.rpt.txt'):
    """Load the map from mouse genes to human genes. A mouse gene may be mapped to more than one
    human genes. In cases like this, human genes are ordered alphabetatically.
    """
    # The mapping file downloaded from http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt.
    map_file = Path(map_file_name)
    if not map_file.exists():
        raise FileNotFoundError("Cannot file file: {}".format(map_file_name))
    map_df = pd.read_csv(map_file, index_col=None, sep='\t')
    # The following code is ported from Java: https://github.com/reactome-fi/corews/blob/
    # d4ad515e53353f0638bcfc7151d9ea72984ef3b9/src/main/java/org/reactome/r3/fi/HumanMouseGeneMapper.java#L36
    id2entries = defaultdict(list)
    for _, row in map_df.iterrows():
        id = row['HomoloGene ID']
        taxon_id = row['NCBI Taxon ID']
        gene = row['Symbol']
        id2entries[id].append('{}\t{}'.format(taxon_id, gene))
    # Map genes from mouse to human
    mousegene2human = dict()
    for entries in id2entries.values():
        mouse_genes = [item.split('\t')[1]
                       for item in entries if item.startswith('10090')]
        human_genes = [item.split('\t')[1]
                       for item in entries if item.startswith('9606')]
        # Just in case
        if len(mouse_genes) == 0 or len(human_genes) == 0:
            continue
        # sort human_genes
        human_genes.sort()
        # Split the mouse genes if any
        for mouse_gene in mouse_genes:
            mousegene2human[mouse_gene] = human_genes
    return mousegene2human


def plot_ligand_src_clustermap(network_with_ligand: nx.DiGraph,
                               cpdb_network: nx.DiGraph,
                               focused_cluster: str | int,
                               figsize = [16, 6]):
    # To avoid the naming confusion between human genes and mouse genes
    # we used all upper case for names
    ligands_in_network = []
    for node, data in network_with_ligand.nodes(data=True):
        if 'type' in data and data['type'] == 'Ligand':
            ligands_in_network.append(node.upper())
    ligands_in_network.sort()
    print('Total ligands in network: {}'.format(len(ligands_in_network)))

    # Build a matrix for these ligands to show their original clusters
    cluster_2_ligands = {}
    for src, target, data in cpdb_network.edges(data=True):
        if cpdb_network.nodes[target]['cluster'] is not str(focused_cluster):
            continue
        src_ligand = src.split(':')[1]
        if not src_ligand.upper() in ligands_in_network:
            continue
        src_cluster = cpdb_network.nodes[src]['cluster']
        if not src_cluster in cluster_2_ligands.keys():
            src_list = np.repeat(0, len(ligands_in_network))
            cluster_2_ligands[src_cluster] = src_list 
        index = ligands_in_network.index(src_ligand.upper())
        cluster_2_ligands[src_cluster][index] += data['score']

    ligand_src_df = pd.DataFrame(columns=ligands_in_network, index=cluster_2_ligands.keys(), dtype=float)
    for cluster in cluster_2_ligands.keys():
        ligand_src_df.loc[cluster] = cluster_2_ligands[cluster]
    
    # Note the transformation
    clustermap = sns.clustermap(data=ligand_src_df.T, figsize=figsize)
    _ = clustermap.ax_heatmap.set_xlabel('Cluster')
    _ = clustermap.ax_heatmap.set_ylabel('Ligand')

    return clustermap