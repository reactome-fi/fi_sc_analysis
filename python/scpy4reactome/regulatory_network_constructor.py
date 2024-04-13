# This script contains functions that are used to contrustruct regulatory networks composed of TFs, pathways, and
# ligands.

import scanpy as sc
import pandas as pd
import pathway_analyzer as pa
import statsmodels.api as sm
from statsmodels.stats.multitest import fdrcorrection
from tqdm import tqdm
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Lasso
from scipy.stats import hypergeom
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import networkx as nx
from scipy.stats import spearmanr
import os
import xml.etree.ElementTree as et


def analyze_tf_pathway_correlations_via_spearman(adata: sc.AnnData,
                                                 pathways: list,
                                                 tfs: list,
                                                 pathway_key: str = pa.AUCELL_KEY,
                                                 tf_key: str = pa.TF_AUCELL_KEY) -> pd.DataFrame:
    result_df = pd.DataFrame(columns=['pathway', 'tf', 'cor', 'cor_p_value'])
    index = 0
    for pathway in tqdm(pathways):
        pathway_score = adata.obsm[pathway_key][pathway]
        for tf in tfs:
            tf_score = adata.obsm[tf_key][tf]
            cor, pvalue = spearmanr(pathway_score, tf_score)
            result_df.loc[index] = [pathway, tf, cor, pvalue]
            index += 1
    # also adding fdr 
    result_df['fdr'] = fdrcorrection(result_df['cor_p_value'])[1]
    return result_df


def analyze_correlation_via_glm(adata: sc.AnnData,
                                pathway: str,
                                tfs: list,
                                pathway_key: pa.AUCELL_KEY,
                                tf_key: pa.TF_AUCELL_KEY):
    """
    Perform a correlation analysis from tfs to a pathway
    :param adata:
    :param pathway:
    :param tfs:
    :return:
    """
    pathway_score = adata.obsm[pathway_key][pathway]
    tfs_scores = adata.obsm[tf_key][tfs]
    x_scores = sm.add_constant(tfs_scores, prepend=False)
    glm_model = sm.GLM(pathway_score, x_scores)
    glm_result = glm_model.fit()
    return glm_result


def analyze_correlations_via_glm(adata: sc.AnnData,
                                 pathways: list,
                                 tfs: list,
                                 pathway_key: str = pa.AUCELL_KEY,
                                 tf_key: str = pa.TF_AUCELL_KEY) -> pd.DataFrame:
    """
    Perform GLM correlation analysis for a set of pathways
    :param adata:
    :param pathways:
    :param tfs:
    :return:
    """
    result_df = pd.DataFrame()
    for pathway in tqdm(pathways):
        # print('Working on pathway: {}'.format(pathway))
        glm_result = analyze_correlation_via_glm(
            adata, pathway, tfs, pathway_key, tf_key)
        result_df[pathway + '_param'] = glm_result.params
        result_df[pathway + '_pvalue'] = glm_result.pvalues
    return result_df


def analyze_correlation_via_lasso(adata: sc.AnnData,
                                  pathway: str,
                                  tfs: list,
                                  pathway_key: pa.AUCELL_KEY,
                                  tf_key: pa.TF_AUCELL_KEY,
                                  alpha: float = 0.01):
    """
    Perform a correlation analysis from tfs to a pathway
    :param adata:
    :param pathway:
    :param tfs:
    :return:
    """
    pathway_score = adata.obsm[pathway_key][pathway]
    tfs_scores = adata.obsm[tf_key][tfs]
    # Standardize features (important for regularized regression)
    scaler = StandardScaler()
    tfs_scores_scaled = scaler.fit_transform(tfs_scores)

    # Create and fit Lasso regression model
    lasso = Lasso(alpha=alpha)  # Adjust alpha for regularization strength
    lasso.fit(tfs_scores_scaled, pathway_score)

    # Since we'd like to get p-values, therefore, use this version for lasso
    # There is no p-values provided here!
    # tfs_scores_scaled = sm.add_constant(tfs_scores_scaled)
    # # Create and fit Lasso regression model using statsmodels
    # # set L1_wt = 1 to make it lasso!
    # lasso_model = sm.OLS(pathway_score, tfs_scores_scaled).fit_regularized(method='elastic_net', alpha=alpha, L1_wt=1)

    return lasso


def _calculate_overlap_p_value(genes1: list,
                               genes2: list,
                               total_genes: int):
    """
    Perform a hypergeomic test to calculate p-values for the overlap.
    :param genes1:
    :param genes2:
    :param total_genes:
    :return:
    """
    if not isinstance(genes1, set):
        genes1 = set(genes1)
    if not isinstance(genes2, set):
        genes2 = set(genes2)
    overlapped = genes1 & genes2
    rv = hypergeom(total_genes, len(genes1), len(genes2))
    return 1 - rv.cdf(len(overlapped) - 1)  # Need to counter itself.


def plot_glm_results(pathway_tf_glm: pd.DataFrame):
    # Plot both pvalues and parameter values from the above dataframe
    # Perform for param
    param_cols = [c for c in pathway_tf_glm.columns if c.endswith('_param')]
    param_df = pathway_tf_glm[param_cols]
    # print(param_df.shape)
    param_df_melt = pd.melt(
        param_df, value_vars=param_df.columns, value_name='param')
    # print(param_df_melt.shape)
    # Want to use absolute values
    param_df_melt['param'] = param_df_melt['param'].map(lambda x: abs(x))
    param_df_melt.sort_values(by='param', inplace=True, ascending=False)
    param_df_melt.index = np.arange(1, len(param_df_melt) + 1)
    # print(param_df_melt.head())

    # For p-values
    pvalue_cols = [c for c in pathway_tf_glm.columns if c.endswith('_pvalue')]
    pvalue_df = pathway_tf_glm[pvalue_cols]
    pvalue_df_melt = pd.melt(
        pvalue_df, value_vars=pvalue_df.columns, value_name='pvalue')
    pvalue_df_melt.sort_values(by='pvalue', inplace=True)
    pvalue_df_melt['-log10(pvalue)'] = pvalue_df_melt['pvalue'].map(
        lambda p: -np.log10(p))
    pvalue_df_melt.index = np.arange(1, len(pvalue_df_melt) + 1)
    # print(pvalue_df_melt.head())

    # Plot
    fig, axs = plt.subplots(nrows=1, ncols=2)
    fig.set_figwidth(30)
    fig.set_figheight(6)
    param_plot = sns.scatterplot(
        x=param_df_melt.index, y=param_df_melt.param, ax=axs[0], s=100, alpha=1.0)
    param_plot.set(xlabel=None)
    pvalue_plot = sns.scatterplot(
        x=pvalue_df_melt.index, y=pvalue_df_melt['-log10(pvalue)'], ax=axs[1])
    pvalue_plot.set(xlabel=None)
    return param_plot, pvalue_plot


def plot_annova(annova_df: pd.DataFrame):
    # The df should be like this:
    #     	F	PR(>F)	fdr_bh
    # Tnr	1790.695490	0.0	0.0
    # Cd34	1930.880181	0.0	0.0
    # Plot both F scores and fdr_bg
    annova_df['-log10(fdr)'] = annova_df['fdr_bh'].map(lambda fdr: -np.log(fdr))
    return _plot_df(annova_df, 'F', '-log10(fdr)')


def _plot_df(df: pd.DataFrame, y1: str, y2: str):
    fig, axs = plt.subplots(nrows=1, ncols=2)
    fig.set_figwidth(30)
    fig.set_figheight(6)
    y1_plot = sns.scatterplot(
        x=df.index, y=df[y1], ax=axs[0], s=100, alpha=1.0)
    y1_plot.set(xlabel=None)
    y2_plot = sns.scatterplot(x=df.index, y=df[y2], ax=axs[1])
    y2_plot.set(xlabel=None)
    return y1_plot, y2_plot


def plot_diff_activity(da_df: pd.DataFrame,
                       y1: str = 'stat',
                       y2: str = '-log10(fdr)',
                       sort_by_fdr: bool = False):
    # The df should be like this:
    #           feature       stat   p-value  median_diff  mean_diff       fdr
    # 266      Cnmd  8287986.0  0.000000     0.163679   0.180981  0.000000
    # 272     Cntn5  8278190.5  0.000000     0.179228   0.194804  0.000000
    da_df['-log10(fdr)'] = da_df['fdr'].map(lambda fdr: -np.log(fdr))
    if sort_by_fdr:
        da_df.sort_values(by='-log10(fdr)', inplace=True, ascending=False)
        da_df.reset_index(drop=True, inplace=True)
    return _plot_df(da_df, y1, y2)


def build_network_for_tfs_pathways(pathways_tf_cors: pd.DataFrame,
                                   cor_cutoff: float = 0.25,
                                   fdr_cutoff: float = 1.0E-6) -> nx.DiGraph:
    """
    Build  network from TFs to pathways. 
    Note: As of April, 2024, we switch to use pariwise correlation to add edges from TFs to pathways
    on the consideration that there are strong co-linear relationships between TF activities since they
    may share quite a lot of target genes.
    :param pathways_tf_cors:
    :return:
    """
    # pathways_tf_cors should be like this:
    # 	pathway	tf	cor	cor_p_value	fdr
    # 0	Metabolism of RNA	Myc	0.814391	0.000000e+00	0.000000e+00
    # 1	Metabolism of RNA	E2f4	0.629553	0.000000e+00	0.000000e+00
    # We will check row by row
    network = nx.DiGraph()
    for _, row in pathways_tf_cors.iterrows():
        cor = row['cor']
        fdr = row['fdr']
        if np.abs(cor) <= cor_cutoff or fdr >= fdr_cutoff:
            continue
        pathway = row['pathway']
        tf = row['tf']
        if not pathway in network:
            network.add_node(pathway, type='Pathway', color='#CCCCFF')
        if not tf in network:
            network.add_node(tf, type='TF', color='#99FFFF')
        color = '#008000' if cor > 0 else '#FF0000'
        annotation = 'tf_pathway_activation' if cor > 0 else 'tf_pathway_inhibition'
        network.add_edge(tf, pathway, color=color,
                         value=np.abs(cor), annotation=annotation)

    return network


def _collect_genes_in_dict(key2genes: dict) -> set:
    all_genes = []
    for genes in key2genes.values():
        all_genes += genes
    return set(all_genes)


def simplify_network_for_tfs_pathways(network: nx.DiGraph,
                                      pathway_gmt_file: str = '../resources/ReactomePathways_Rel_79_122921.gmt',
                                      pathway_hierarchy_file: str = '../resources/HumanReactomePathways_Hierarchy_Rel_79.xml',
                                      tf_file: str = '../resources/dorothea_hs.tsv',
                                      dorothea_evidence: list = [
                                          'A', 'B', 'C', 'D', 'E'],
                                      need_parent_pathway: bool = True,
                                      for_pathway_only: bool = False,
                                      use_direct_interaction: bool = True,
                                      p_value_for_direction_interaction: float = 0.05,
                                      add_tf_links: bool = True,
                                      check_with_tf_cor: bool = True,
                                      adata: sc.AnnData = None,
                                      tf_cor_cutoff: float = 0.25,
                                      add_pathway_to_tf_links: bool = True,
                                      delete_pathway_only_component: bool = True,
                                      delete_leaf_pathway: bool = True) -> tuple[nx.DiGraph, pd.DataFrame]:
    """
    Build a network from TFs to pathways by analyzing the covergage and simplifying connections.
    :param dorothea_evidence: the evidence level for TF/target dict.
    :param use_direct_interaction: true to include directed targetted pathways of TFs based on overlap analysis
    :param p_value_for_direction_interaction: when use_use_direct_interaction is true, use this p_value to selec directed
    interactio
    :param add_pathway_to_tf_links: add links from a pathway to a tf if the pathway contains tf if true.
    :param check_with_tf_cor: true to check correlations for linking tfs
    :param adata: make sure adata exists if check_with_tf_cor is true
    :param delete_pathway_only_component: delete pathway only weakly components if true
    :param delete_leaf_pathway: true to delete pathways that are not impacte by any TF by linked to other pathways.
    """
    pathway2genes = pa._load_reactome_gmt(pathway_gmt_file)
    all_pathway_genes = _collect_genes_in_dict(pathway2genes)
    # For the overlap analysis, we pick up as many targets as possible
    tf2genes = pa.load_dorothea_data(tf_file, dorothea_evidence)
    all_tf_genes = _collect_genes_in_dict(tf2genes)
    overlapped_genes = all_pathway_genes & all_tf_genes
    total_genes = len(overlapped_genes)

    # Make sure these condictions are right
    if check_with_tf_cor:
        if adata is None or pa.TF_AUCELL_KEY not in adata.obsm.keys():
            raise ValueError('Make sure adata is passed and {} is in '
                             'adata.obsm since check_with_tf_cor is true!'.format(pa.TF_AUCELL_KEY))
    # As of August 2, 2022, collect all pathways first regardless if they are connected or not
    pathways_in_network = [x for x, y in network.nodes(
        data=True) if y['type'] == 'Pathway']
    overlap_p_values = pd.DataFrame(columns=['TF', 'Pathway', 'p-Value'])
    row = 0
    # Use the Reactome pathway hierachy for parent/child relationships
    n_network = build_parent_network(
        pathways_in_network, pathway_hierarchy_file)
    # Assign type and color
    nx.set_node_attributes(n_network, 'Pathway', 'type')
    nx.set_node_attributes(n_network, '#CCCCFF', 'color')
    nx.set_edge_attributes(
        n_network, 'pathway_pathway_hierarchy', 'annotation')
    # Assign color
    if not for_pathway_only:
        # Add TFs
        tfs_in_network = [x for x, y in network.nodes(
            data=True) if y['type'] == 'TF']
        for tf in tfs_in_network:
            tf_edges = network.out_edges(tf)
            if len(tf_edges) == 0:
                continue
            n_network.add_node(tf, type='TF', color='#99FFFF')
            tf_edge_pathways = [v for u, v in tf_edges]
            # Get the original edges and make sure only one is added in the new network if multiple pathways are hit in
            # the same branch
            # For each pathway, get its ancestors and remove them from the list
            to_be_removed = set()
            for edge_pathway in tf_edge_pathways:
                ancestors = []
                _collect_ancesctor_pathways(n_network, edge_pathway, ancestors)
                to_be_removed.update(ancestors)
            for edge_pathway in tf_edge_pathways:
                if edge_pathway in to_be_removed:
                    # print('Not added: {}'.format(edge_pathway))
                    continue
                if edge_pathway not in n_network:
                    print('Warning: "{}" is not in the Reactome pathway tree.'.format(
                        edge_pathway))
                    continue
                overlap_p_value = _calculate_overlap_p_value(tf2genes[tf],
                                                             pathway2genes[edge_pathway],
                                                             total_genes)
                overlap_p_values.loc[row] = [tf, edge_pathway, overlap_p_value]
                row += 1
                color_new = network[tf][edge_pathway]['color']
                annotation = network[tf][edge_pathway]['annotation']
                if overlap_p_value < p_value_for_direction_interaction:
                    n_network.add_edge(tf,
                                       edge_pathway,
                                       color=color_new,
                                       annotation=annotation,
                                       value=network[tf][edge_pathway]['value'])  # Value used for weight
                elif not use_direct_interaction:
                    if color_new == '#008000':
                        color_new = '#CCFFCC'
                    else:
                        color_new = '#FFCCCC'  # 90% lighter
                    n_network.add_edge(tf,
                                       edge_pathway,
                                       color=color_new,
                                       annotation=annotation + "_indirect",
                                       value=network[tf][edge_pathway]['value'])  # Value used for weight
    if add_tf_links:
        # Add TFs
        tfs_in_network = [x for x, y in n_network.nodes(
            data=True) if y['type'] == 'TF']
        build_tfs_network(tfs_in_network, tf2genes,
                          check_with_tf_cor, tf_cor_cutoff, adata, n_network)
    if add_pathway_to_tf_links:
        tfs_in_network = [x for x, y in n_network.nodes(
            data=True) if y['type'] == 'TF']
        pathways_in_network = [x for x, y in n_network.nodes(
            data=True) if y['type'] == 'Pathway']
        add_pathway_to_tf_network(
            pathways_in_network, tfs_in_network, pathway2genes, n_network)
    if delete_pathway_only_component:
        _remove_pathway_only_component(n_network)
    if delete_leaf_pathway:
        _remove_leaf_pathways(n_network)
    # Remove nodes that are not connected
    to_be_removed = []
    for node, atts in n_network.nodes(data=True):
        if n_network.degree(node) == 0:
            to_be_removed.append(node)
        elif not need_parent_pathway:
            # Delete pathway nodes that don't link to any TFs
            if atts['type'] == 'Pathway':
                needed = False
                for neighbor in n_network.predecessors(node):
                    if n_network.nodes[neighbor]['type'] == 'TF':
                        needed = True
                        break
                if not needed:
                    # Check other ways
                    for neighbor in n_network.successors(node):
                        if n_network.nodes[neighbor]['type'] == 'TF':
                            needed = True
                            break
                if not needed:
                    to_be_removed.append(node)
    n_network.remove_nodes_from(to_be_removed)
    return n_network, overlap_p_values


def build_tfs_network(tfs: list,
                      tf2targets: dict,
                      check_with_tf_cor: bool = True,
                      tf_cor_cutoff: float = 0.25,
                      adata: sc.AnnData = None,
                      network: nx.DiGraph = None) -> nx.DiGraph:
    """
    Build TF/target network for a list of tfs.
    :param tfs:
    :param tf2targets:
    :return:
    """
    if network is None:
        network = nx.DiGraph()
    # To scale the value so that we have the same value for tf -> pathway
    max_value = None
    all_values = [np.abs(v) for k, v in nx.get_edge_attributes(
        network, 'value').items()]
    if len(all_values) > 0:
        max_value = max(all_values)
    for tf1 in tfs:
        for tf2 in tfs:
            # We may see feedback loops to tf itself or others. Therefore both check
            if tf2 in tf2targets[tf1]:
                if check_with_tf_cor:
                    if adata is not None and pa.TF_AUCELL_KEY in adata.obsm.keys():
                        # Calculate correlation
                        tf1_aucell = adata.obsm[pa.TF_AUCELL_KEY][tf1]
                        tf2_aucell = adata.obsm[pa.TF_AUCELL_KEY][tf2]
                        # Ignore p-value for the time being here
                        cor_value = spearmanr(tf1_aucell, tf2_aucell)[0]
                        if max_value is not None:
                            cor_value *= max_value  # The maximum should be the same
                        if cor_value > tf_cor_cutoff:
                            network.add_edge(tf1, tf2, color='#8FBC8F',
                                             annotation='tf_tf_activation', value=cor_value)  # dark gree for positiveTF edges
                        elif cor_value < -tf_cor_cutoff:
                            network.add_edge(tf1, tf2, color='#A52A2A',
                                             annotation='tf_tf_inhibition', value=-cor_value)  # brown for negative TF edges
                else:
                    # brown for TF edges
                    network.add_edge(
                        tf1, tf2, annotation='tf_tf_interaction', color='#A52A2A')
    return network


def build_parent_network(pathways: list,
                         pathway_hierarchy_file: str) -> nx.DiGraph:
    """
    Basically create a simplified pathway hierarchy for the list of passed pathways.
    :param pathways:
    :return:
    """
    reactome_network = load_pathway_hierarchy(pathway_hierarchy_file)
    # To be returned
    pathway_network = nx.DiGraph(reactome_network)
    selected_nodes = set()
    # Check each pathway in reactome_network and see if it should be included
    for node in reactome_network:
        if include_node(node, reactome_network, pathways):
            selected_nodes.add(node)
    pathway_network.remove_nodes_from(
        set(pathway_network.nodes) - selected_nodes)
    # Remove nodes that are not in the pathways list and link nodes in the pathways list together
    has_changed = True
    to_be_removed = []
    while has_changed:
        has_changed = False
        to_be_removed.clear()
        # To avoid error
        nodes = list(pathway_network.nodes)
        for node in nodes:
            # Check if this node should be removed
            if node not in pathways:
                to_be_removed.append(node)
                # Re-link: our network links from child pathways to parent pathways
                for parent in pathway_network.successors(node):
                    for child in pathway_network.predecessors(node):
                        pathway_network.add_edge(child, parent)
        if len(to_be_removed) > 0:
            pathway_network.remove_nodes_from(to_be_removed)
            has_changed = True
    return pathway_network


def add_pathway_to_tf_network(pathways: list,
                              tfs: list,
                              pathway2genes: dict,
                              network: nx.DiGraph) -> nx.DiGraph:
    """
    Build a network or add new edges bewteen pathways and tfs. If a pathway has a TF annotated, an edge
    will be added between this pathway to TF (from pathway to TF). This is a highly simplified version.
    Manual check should be performed in the future to ensure these links are correct.
    :param pathways:
    :param tfs:
    :param pathway2genes:
    :param netowrk:
    :return:
    """
    # This is a little bit more complicated. We want to push this link as low as possible.
    # Start with TFs
    for tf in tfs:
        selected_pathways = []
        for pathway in pathways:
            pathway_genes = pathway2genes[pathway]
            if tf in pathway_genes:
                selected_pathways.append(pathway)
        # To a prunning to remove other pathways that can be covered by child pathways
        if len(selected_pathways) > 0:
            for pathway in selected_pathways:
                if not _is_in_child_pathways(pathway, selected_pathways, network):
                    # Point from pathway to tf. Dark Salmon
                    network.add_edge(
                        pathway, tf, annotation='pathway_tf_annotation', color='#E9967A')
    return network


def _is_in_child_pathways(pathway, selected_pathways, network) -> bool:
    in_edges = network.in_edges(pathway)
    for u, v in in_edges:
        if network.nodes[u]['type'] == 'TF':
            continue
        if u in selected_pathways:
            return True
        is_included = _is_in_child_pathways(u, selected_pathways, network)
        if is_included:
            return True
    return False


def add_median_to_nodes(network: nx.DiGraph,
                        df: pd.DataFrame):
    """Add median values to the nodes as a ttribute
    Args:
        network (_type_): the targetted network
        df (_type_): a dataframe that should have columns called 'median_diff' and 'feature'.
    """
    node2att = dict()
    for _, row in df.iterrows():
        att = {'value': row['median_diff']}
        node2att[row['feature']] = att
    nx.set_node_attributes(network, node2att)


def _collect_ancesctor_pathways(network: nx.DiGraph,
                                pathway: str,
                                ancestors: list):
    if pathway is None:
        return
    out_edges = network.out_edges(pathway)
    for u, v in out_edges:
        ancestors.append(v)
        _collect_ancesctor_pathways(network, v, ancestors)


def _remove_leaf_pathways(network: nx.DiGraph) -> nx.DiGraph:
    is_changed = True
    while is_changed:
        is_changed = False
        pathways = [n for n, a in network.nodes(
            data=True) if a['type'] == 'Pathway']
        to_be_removed = []
        for pathway in pathways:
            if network.out_degree(pathway) == 1 and network.in_degree(pathway) == 0:
                for u, v in network.out_edges(pathway):
                    if network.nodes[v]['type'] == 'Pathway':
                        to_be_removed.append(u)
        if len(to_be_removed) > 0:
            is_changed = True
            network.remove_nodes_from(to_be_removed)
    return network


def _remove_pathway_only_component(network: nx.DiGraph) -> nx.DiGraph:
    """
    Delete components that don't have any TFs
    :param network:
    :return:
    """
    to_be_removed = []
    has_tf = False
    for component in nx.weakly_connected_components(network):
        has_tf = False
        for node in component:
            if network.nodes[node]['type'] == 'TF':
                has_tf = True
                break
        if not has_tf:
            # component is a set. Make sure to use extend to unecapsulte elements
            to_be_removed.extend(component)
    network.remove_nodes_from(to_be_removed)
    return network


def include_node(node, reactome_network, pathways) -> bool:
    if node in pathways:
        return True
    # If both ancestors and descendents of node are in pathways, return true since
    # this node may be used for connection. Otherwise return false.
    # Check ancesctors first
    if not _include_node(node, reactome_network, pathways, 1):
        return False
    # Check descendents
    if not _include_node(node, reactome_network, pathways, 0):
        return False
    return True


def _include_node(node, reactome_network, pathways, index) -> bool:
    # This should be a width-first search
    if index == 1:  # ancestors check
        out_edges = reactome_network.out_edges(node)
        for out_edge in out_edges:
            if out_edge[index] in pathways:
                return True
        # Width first
        for out_edge in out_edges:
            if _include_node(out_edge[1], reactome_network, pathways, index):
                return True  # Any ancesctor should be fine
    elif index == 0:  # descedent
        in_edges = reactome_network.in_edges(node)
        for in_edge in in_edges:
            if in_edge[index] in pathways:
                return True
        for in_edge in in_edges:
            if _include_node(in_edge[0], reactome_network, pathways, index):
                return True  # Any descend should be good
    return False


def load_pathway_hierarchy(reactome_pathway_hierarchy_xml_file: str) -> nx.DiGraph:
    """
   Load the Reactome mouse pathway hierarchy tree into a network DiaGrah.
   :param file_name:
   :return:
   """
    # Make sure the file is there
    if not os.path.exists(reactome_pathway_hierarchy_xml_file):
        raise ValueError("{} doesn't exist.".format(
            reactome_pathway_hierarchy_xml_file))
    # The file is an XML file. Do parsing here.from
    graph = nx.DiGraph()
    tree = et.parse(reactome_pathway_hierarchy_xml_file)
    root = tree.getroot()
    for child in root:
        _load_pathway_hierarchy(child, graph)
    return graph


def add_ligands_to_network_via_pathways(ligands: list,
                                        network: nx.DiGraph,
                                        adata: sc.AnnData,
                                        pathway_gmt_file: str,
                                        ligand2targets: dict,
                                        pathway_key: str = pa.AUCELL_KEY,
                                        ligand_key: str = pa.LIGAND_AUCELL_KEY,
                                        overlap_fdr_cutoff: float = 0.01,
                                        cor_cutoff: float = 0.20,
                                        cor_pvalue_cutoff: float = 0.001,
                                        use_tf_linked_pathways: bool = True,
                                        graphml: str = None) -> tuple[nx.DiGraph, pd.DataFrame]:
    """Add ligands into a regulatory network composed of TFs and pathways.

    Args:
        ligands (list): a list of ligands to be added
        network (nx.DiGraph): the target network
        adata (sc.AnnData): the data source having all information to calculate correlation between pathways and ligands.
        pathway_gmt_file (str): Reactome GMT file to be loaded for overlap analysis
        ligand2targets (dict): Extract ligand's target for overlap analysis
        overlap_pvalue_cutoff (float, optional): the cutoff p-value to add edge from ligand to pathway. Defaults to 0.01.
        cor_cutoff (float, optional): the correlation cutoff from ligand to pathway from the GLM analysis. Defaults to 0.20.
        cor_pvalue_cutoff (float, optional): this cutoff usually should not be used.
        graphml (str, optional): the target file to export the network. Defaults to None.

    Returns:
        tuple[nx.DiGraph, pd.DataFrame]: the network with ligands added and the dataframe with the information used to add ligands to the network.
        Note: The returned network is the same object of the network passed to this function.
    """
    # Get the genes
    pathway2genes = pa._load_reactome_gmt(pathway_gmt_file)

    # Need to figure out the target pathways in the network
    target_pathways = set()
    # Used to check if the edge linked to there is needed
    if use_tf_linked_pathways:
        for source, target in network.edges():
            source_type = network.nodes[source].get('type')
            target_type = network.nodes[target].get('type')
            if source_type == "TF" and target_type == "Pathway":
                target_pathways.add(target)
            elif source_type == "Pathway" and target_type == "TF":  # Either way should be fine
                target_pathways.add(source)
    else:  # Use all pathways
        for node, data in network.nodes(data=True):
            node_type = data.get('type')
            if node_type == "pathway":
                target_pathways.add(node)

    # Since it is expected to see strong colinear relationships among ligands activities, caused by
    # shared overlapped genes, it is better to use simple pairwise correlation analysis here.
    # More sophsiticated method may be explored but not using GLM without more detailed analysis
    # to avoid co-linear relationships!
    ligands_pathways_cor_df = analyze_ligand_pathway_correlation_via_spearman(adata=adata,
                                                                              target_pathways=target_pathways,
                                                                              pathway_key=pathway_key,
                                                                              pathway2genes=pathway2genes,
                                                                              ligands=ligands,
                                                                              ligand_key=ligand_key,
                                                                              ligand2targets=ligand2targets)
    for ligand in ligands:
        network.add_node(ligand, type='Ligand', color='#21D5F1')
    # The results should be something like this:
    #     	ligand	pathway	cor	cor_p_value	overlap_p_value
    # 0	Nectin2	Signaling by Nuclear Receptors	0.489708	0.000000e+00	4.930094e-04
    # 1	Dsc2	Signaling by Nuclear Receptors	0.462185	0.000000e+00	1.073196e-04
    # We will check row by row
    for _, row in ligands_pathways_cor_df.iterrows():
        cor = row['cor']
        cor_pvalue = row['cor_p_value']
        overlap_fdr = row['overlap_fdr']
        if np.abs(cor) <= cor_cutoff or cor_pvalue >= cor_pvalue_cutoff or overlap_fdr >= overlap_fdr_cutoff:
            continue  # Escape this row
            # Add this relationship
        ligand = row['ligand']
        pathway = row['pathway']
        color = '#86FF00' if cor > 0 else '#FF9600'
        annotation = 'ligand_pathway_activation' if cor > 0 else 'ligand_pathway_inhibition'
        network.add_edge(ligand, pathway, color=color,
                         value=np.abs(cor), annotation=annotation)  # value used for weight in pyviz

    # If a pathway has an edge to a ligand, any edges having the same type to the ligand will be removed
    # since semanatically we may be able to explain the observed relationship between this ligand to all
    # pathways in the branch. The following code is used to do this clean-up using the directed edges
    for ligand in ligands:
        # Remove ligands that don't have any edges
        if network.degree(ligand) == 0:
            network.remove_node(ligand)
            continue
        for pathway in target_pathways:
            # See if ther is an edge between this ligand and this pathway
            if not network.has_edge(ligand, pathway):
                continue
            # activation or inhibition
            edge_type = network.get_edge_data(
                ligand, pathway).get('annotation').split('_')[2]
            # In the current network, child pathways are precedessors
            lower_pathways = network.predecessors(pathway)
            need_remove = False
            for lower_pathway in lower_pathways:
                if not network.has_edge(ligand, lower_pathway):
                    continue
                lower_edge_type = network.get_edge_data(
                    ligand, lower_pathway).get('annotation').split('_')[2]
                if edge_type == lower_edge_type:
                    need_remove = True
                    break
            if need_remove:
                network.remove_edge(ligand, pathway)

    if graphml is not None:
        write_network_to_graphml(network, graphml)

    return network, ligands_pathways_cor_df


def analyze_ligand_pathway_correlation_via_spearman(adata: sc.AnnData,
                                                    target_pathways: list,
                                                    pathway_key: str,
                                                    pathway2genes: dict,
                                                    ligands: list,
                                                    ligand_key: str,
                                                    ligand2targets: dict) -> pd.DataFrame:
    """Perform pairwise correlation analysis. If we see strong co-linear relationships between predictors, we should
    use this pairwise correlation analysis. Otherwise, the relationships among features (e.g. TFs and ligands activities)
    maybe difficult to explain.

    Args:
        adata (sc.AnnData): _description_
        target_pathways (list): _description_
        pathway_key (str): _description_
        pathway2genes (dict): _description_
        ligands (list): _description_
        ligand_key (str): _description_
        ligand2targets (dict): _description_

    Returns:
        pd.DataFrame: _description_
    """
    total_pathway_genes = _collect_genes_in_dict(pathway2genes)
    total_ligand_genes = _collect_genes_in_dict(ligand2targets)
    total_shared_genes = total_pathway_genes & total_ligand_genes

    ligands_pathway_pd = pd.DataFrame(
        columns=['ligand', 'pathway', 'cor', 'cor_p_value', 'overlap_p_value'])

    # Perform correlation analysis
    index = 0
    for pathway in tqdm(target_pathways):
        pathway_score = adata.obsm[pathway_key][pathway]
        pathway_genes_1 = pathway2genes[pathway]
        pathway_genes_1 = set(pathway_genes_1) & total_shared_genes
        for ligand in ligands:
            ligand_score = adata.obsm[ligand_key][ligand]
            cor, pvalue = spearmanr(pathway_score, ligand_score)
            ligand_genes_1 = ligand2targets[ligand]
            ligand_genes_1 = set(ligand_genes_1) & total_shared_genes
            overlap_p_value = _calculate_overlap_p_value(
                pathway_genes_1, ligand_genes_1, len(total_shared_genes))
            ligands_pathway_pd.loc[index] = [
                ligand, pathway, cor, pvalue, overlap_p_value]
            index += 1
    ligands_pathway_pd['overlap_fdr'] = fdrcorrection(ligands_pathway_pd['overlap_p_value'])[1]

    return ligands_pathway_pd


def _load_pathway_hierarchy(elm, graph):
    if elm.tag == 'Reaction':
        return  # Don't include reaction
    if len(list(elm)) == 0:
        return  # Get to the bottom. Nothing to do
    current_name = elm.attrib['displayName']
    # Some of names having trailing space
    current_name = current_name.strip()
    # Escape disease pathways
    if current_name == 'Disease':
        return
    for child in elm:
        if child.tag == 'Reaction':
            continue
        child_name = child.attrib['displayName'].strip()
        graph.add_edge(child_name, current_name)
        _load_pathway_hierarchy(child, graph)


def write_network_to_graphml(network: nx.DiGraph,
                             filename: str):
    nx.write_graphml(network, filename)
