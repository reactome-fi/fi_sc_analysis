"""
This JSONRPC server is used to provide server wrapper aroung some major functions provided by scanpy and its related Python packages.
These functions are provided as json-based REST calls for the ReactomeFIViz app.

"""
# Required by scanpy
# import numpy as np
# import pandas as pd
from builtins import ValueError

import scanpy as sc
from scanpy.plotting import _utils
from anndata import AnnData
from typing import Optional
import numpy as np
from numpy import ndarray
# for Perason correlation
from scipy import stats
# for nnlm
from scipy import optimize
# for pagerank and other stuff
import networkx as nx
import pandas as pd

# As suggested in the tutorial by scanpy doc: Preprocessing and clustering 3k PBMCs
sc.settings.verbosity = 3
sc.logging.print_versions()
sc.set_figure_params(dpi = 80, facecolor = 'white')

# Global level flags to control the use of the serve
# server = SimpleJSONRPCServer("localhost")
# isWaiting = True

# This dict is used for cache loaded AnnaData
cache = dict()
KEY_DATA, KEY_PROCESSED, KEY_MERGED = 'adata', 'processed', 'merged'

def get_loaded_data() :
    if KEY_DATA not in cache.keys() :
        return None
    return cache[KEY_DATA]

def cache_data(adata: AnnData) :
    cache[KEY_DATA] = adata

def get_processed_data() :
    if KEY_PROCESSED not in cache.keys() :
        return None
    return cache[KEY_PROCESSED]

def cache_processed_data(adata: AnnData) :
    cache[KEY_PROCESSED] = adata

def cache_merged_data(adata: AnnData) :
    cache[KEY_MERGED] = adata

def get_merged_data() :
    if KEY_MERGED not in cache.keys() :
        return None
    return cache[KEY_MERGED]

def open_10_genomics_data(dir) :
    """
        Open a 10x genomics data into an annadata object.
    """
    adata = sc.read_10x_mtx(dir, var_names='gene_symbols', cache=True)
    adata.var_names_make_unique()
    return adata

def filter(adata, min_genes = 200, min_cells = 3, copy = False) :
    """
        Perform a simple filtering by remove cells having less than min_genes and genes having less than
        min_cells.
    """
    if (copy) :
        adata = adata.copy()
    sc.pp.filter_cells(adata, min_genes = min_genes)
    sc.pp.filter_genes(adata, min_cells = min_cells)
    # Regardless return the object
    return adata

def preprocess(adata, copy=False, need_scale=True, regression_keys = None):
    """
        Perform a series of steps to pre-process the data for clustering analysis.
        :param copy: true to make a copy of adata
        :param need_scale:  true for including these steps: slicing to highly_variable, regress_out based on total_counts
        and conduct gene-wise z-score
        :return the client to this method should get the returned adata since the referred object will be changed.
    """
    if (copy) :
        adata=adata.copy()
    filter(adata, copy=False) # Regardless we should not copy it
    # Check mitochondrial genes to make sure they don't interfere our analysis
    # To support both human and mouse genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
    # calculate a bunch of stats for as qc metrics
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    # Mark genes with highly_variable flag
    # Default parameters are used in this function call.
    sc.pp.highly_variable_genes(adata)
    if not need_scale :
        return adata
    # Keep the normalized data for future use
    adata.raw = adata
    # Keep the highly variable genes only.
    # Note: A new adata object is created!
    adata = adata[:, adata.var.highly_variable]
    regress(adata, keys=regression_keys)
    return adata

def regress(adata, keys = None) :
    """
    Conduct an optional regress
    :param adata:
    :param keys:
    :return:
    """
    # Default using total_counts though pct_counts_mt is also used i the tutorial
    # TODO: Need to consider to add regress_out for cell cycle score!
    # The following regress_out is based on tutorial: https://scanpy-tutorials.readthedocs.io/en/latest/pbmc3k.html, where
    # total_counts is included even though cell normalization has been performned.
    # Don't procide any regression-out if the user doesn't want. No-regressing produces the best results with
    # mouse intestinal stem cell data.
    if keys is None :
        # keys = ['total_counts']
        # keys = ['pct_counts_mt']
        keys = ['total_counts', 'pct_counts_mt']
    if keys is not None :
        sc.pp.regress_out(adata, keys=keys)
    # gene-wise z-score
    sc.pp.scale(adata, max_value=10)

def cluster(adata, plot=False) :
    """
    This method is a collection of steps: perform PCA, compute neighborhood graph, embed the neighborhood graph
    via UMAP, and clustering the neighborbood graph using the leiden algorithm. The passed adata will hve a bunch
    of new variables added in the unc field.
    :param adata:
    :param plot true to plot the final umap.
    :return: None
    """
    sc.tl.pca(adata) # tl: tools
    # The default parameters are used.
    # TODO: Need to expose the parameters used in the function definition.
    sc.pp.neighbors(adata) # pp: preprocess
    # clustering using the Leiden algorithm. The clustering should be run first before the following paga
    # related statements
    sc.tl.leiden(adata)
    # Since we want to use paga for trajectory inference, it will be nicer to do this right now. The performance
    # penality is not that big based on the mouse ISC data set. These procedures are based on the scanpy tutorial.
    # All default paramters are used here.
    random_state = 1234 # Make sure it can be repeatable
    sc.tl.paga(adata)
    # For some reason, each run may produce a little bit different layout
    sc.pl.paga(adata, plot=False, random_state = random_state) # No need to plot. pl: plot
    # Do umap
    sc.tl.umap(adata, init_pos='paga', random_state=random_state)
    adata.uns['paga']['pos'] = reset_paga_pos(adata)
    # Finding markers genes. Use t-test here for the time being.
    # TODO: This needs to be improved in the future. Two more iterations:
    # 1). Expose the method in the function definition so that the user can choose what method should be used.
    # 2). Use a more sophisticated method: e.g. diffx and something based on ZNBM.
    # Move this to its own place
    # sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
    # tl.umap call should produce the coordinates. If we just need coordinates, there is no need to call the following
    if plot:
        sc.pl.umap(adata, color=['leiden'])
    # There is no need to return adata

def rank_genes_groups(adata,
                      groupby = 'leiden',
                      groups = 'all',
                      reference = 'rest',
                      method = 't-test_overestim_var') :
    """
    Conduct a differential expression analysis for a group of genes.
    :param adata:
    :param groupby:
    :param groups:
    :param reference:
    :param method:
    :return:
    """
    # We want to get all genes
    total_genes = None
    if adata.raw is not None :
        total_genes = adata.raw.n_vars
    else :
        total_genes = adata.n_vars
    # Make sure groups is a list
    if groups is not 'all' and isinstance(groups, str) :
        groups = [groups]
    sc.tl.rank_genes_groups(adata,
                            groupby = groupby,
                            use_raw = True,
                            groups = groups,
                            reference = reference,
                            rankby_abs = True,
                            method = method,
                            n_genes = total_genes) # Get all genes


def reset_paga_pos(adata) -> ndarray:
    """
    The following code is based on scatterplots.py to calculate paga positions directly from
    the umap. The original paga's pos is not matched with umap because of different layout algorithm
    used.
    :param adata:
    :return: the median coordinates of all cells in individual clusters.
    """
    key = 'leiden'
    umap_key = "X_umap"
    clusters = adata.obs[key].cat.categories
    cluster_pos = np.zeros((len(clusters), 2))
    for icluster, cluster in enumerate(clusters) :
        selected = adata.obs[key] == cluster
        selected_umap = adata.obsm[umap_key][selected, :]
        x_pos, y_pos = np.median(selected_umap, axis=0)
        cluster_pos[icluster] = [x_pos, y_pos]
    return cluster_pos


def dpt(adata: AnnData,
        root_cell : str) :
    """
    Perform dpt analysis with the provided root_cell as the root.
    :param adata: pre-processed data
    :param root_cell:
    :return:
    """
    if root_cell not in adata.obs_names :
        raise ValueError(root_cell + " is not in the obs_names.")
    # Get the number index of the root_cell
    root_index = np.flatnonzero(adata.obs_names == root_cell)[0]
    adata.uns['iroot'] = root_index
    sc.tl.dpt(adata)
    if 'dpt_pseudotime' not in adata.obs.keys() :
        raise ValueError("Cannot find dpt_pseudotime in the adata's keys.")
    return adata.obs['dpt_pseudotime']

def project(new_dir_name: str,
            adata_ref: AnnData,
            batch_categories: Optional[str] = None,
)->AnnData:
    """
    This function is used to project adata to adata_ref using ingest function in scanpy. The procedures here are based
    on the tutorial: integration-data-using-inject. Both adata ad adata_ref should have been preprocessed by running the
    preprocess function but with need_scale false.
    :param new_dir_name: the data directory to be projected
    :param adata_ref: the reference data
    :param batch_categories: the name for the reference should be the first element
    :return: return a merged AnnData object with two originally data copied and merged.
    """
    # Check if leiden has been conducted.
    if 'leiden' not in adata_ref.obs :
        cluster(adata_ref)
    adata = open_10_genomics_data(new_dir_name)
    # Make sure we do the same thing as in the original data. But we don't want to keep the original data
    adata = preprocess(adata, copy=False, need_scale=False)
    # Make sure both of them have the same set of variables
    shared_var_names = adata_ref.var_names.intersection(adata.var_names)
    sc.logging.info("shared_var_names: {}".format(len(shared_var_names)))
    # slicing the data to make copies
    adata = adata[:, shared_var_names]
    # Call regress here so that we have almost the same number of genes selected by the adata_ref (aka highly invariable genes)
    regress(adata)
    adata_ref = adata_ref[:, shared_var_names]
    # inject based on the leiden
    sc.tl.ingest(adata, adata_ref, obs = 'leiden')
    # merge two objects
    if not batch_categories :
        batch_categories = ['ref', 'new']
    adata_merged = adata_ref.concatenate(adata, batch_categories = batch_categories)
    return adata_merged

def cytotrace(adata: AnnData) -> list :
    """
    This is a port of the original R implementation of Cytotrace by taking advantage of the existing connectivities
    and preprocessed data.
    :param adata:
    :return: a list of scaled cytotrace scores.
    """
    # Make sure we have all data pre-calcualte
    sc.logging.info("Running cytotrace...")
    gene_count_key = 'n_genes_by_counts'
    if gene_count_key not in adata.obs.keys() :
        raise ValueError(gene_count_key + " is not in the obs names. Run preprocess first.")
    connectivities_key = 'connectivities'
    if connectivities_key not in adata.obsp.keys() :
        raise ValueError(connectivities_key + " is not in the obsp. Run preprocess first.")
    # Calculate the Pearson correlation
    sc.logging.info("Calculating Pearson correlations between gene counts and expressions...")
    gene_count_exp_pc = list()
    # Use raw for correlations
    for gene in adata.var_names :
        pc = stats.pearsonr(adata.raw.obs_vector(gene), adata.obs[gene_count_key])[0]
        gene_count_exp_pc.append(pc)
    # sort to pick up the first 200 genes having the largest positive correlation
    gene_count_exp_pc_sorted_indice = np.argsort(gene_count_exp_pc)
    # Need other way around
    gene_count_exp_pc_sorted_indice = np.flip(gene_count_exp_pc_sorted_indice)
    # Pick top 200 genes
    top_200_genes = adata.var_names[gene_count_exp_pc_sorted_indice[0:200]]
    # Slice the adata.raw since we used the raw data
    top_200_adata = adata.raw[:, top_200_genes]
    # Get the mean
    sc.logging.info("Calculating cell means for chosen genes...")
    cell_means = list()
    for cell in adata.obs_names :
        # cell_mean = np.mean(top_200_adata.var_vector(cell))
        # To calculate geom mean
        cell_vector = top_200_adata.var_vector(cell)
        # Use this to avoid 0
        cell_vector = np.exp(cell_vector)
        cell_mean = stats.gmean(cell_vector)
        cell_mean = np.log(cell_mean)
        cell_means.append(cell_mean)
    # Conduct nnls
    # Need to convert to array as required by nnls
    sc.logging.info("Running nnls...")
    # print(cell_means) # Just for check
    nnls_result = optimize.nnls(adata.obsp[connectivities_key].toarray(), cell_means)
    # nnls_result is tupe. Just need the first element, which is an array
    regressed_scores = adata.obsp[connectivities_key] * nnls_result[0]
    # regressed_scores = adata.obs[gene_count_key]
    page_ranks = run_pagerank(adata, connectivities_key, regressed_scores)
    # page_ranks is a dict. Scale the value from 0 to 1
    page_ranks_min = min(page_ranks.values())
    page_ranks_max = max(page_ranks.values())
    scale = page_ranks_max - page_ranks_min
    cytotrace = [None] * len(page_ranks)
    for key in page_ranks.keys() :
        value = page_ranks[key]
        cytotrace[key] = (value - page_ranks_min) / scale
    # Save to the adata as observation
    adata.obs['cytotrace'] = cytotrace
    sc.logging.info("Done cytotrace.")
    return cytotrace

def infer_cell_root(adata : AnnData,
                    candidate_clusters = None) :
    """
    This function is used to infer a possible cell root for dpt. The algorithm works like this:
    1). Conduct a page rank analysis use n_genes_by_counts as the personalization scores and the graph base
    on connectivities.
    2). If no candidate_clusters have been provided, do the following:
        a). Sort cell clusters based on median page ranks.
        b). Conduct a t-test starting from the top-most clusters based on page ranks and merge cell clusters
    that don't have any significance based on FDR cutoff 0.05
        c). Selected merged clusters as candidate_clusters
    3). Sort cells in the candidate_clsuters and choose the cell having the largest pagerank score to return.
    :param adata: pre-processed and clustered data
    :param candidate_clusters: one or more than one clusters
    :return:
    """
    sc.logging.info("Running infer_cell_root...")
    gene_count_key = 'n_genes_by_counts'
    if gene_count_key not in adata.obs.keys() :
        raise ValueError(gene_count_key + " is not in the obs names. Run preprocess first.")
    connectivities_key = 'connectivities'
    if connectivities_key not in adata.obsp.keys() :
        raise ValueError(connectivities_key + " is not in the obsp. Run preprocess first.")
    leiden_key = 'leiden'
    if leiden_key not in adata.obs.keys() :
        raise ValueError(leiden_key + " is not in the obs names. Run cluster first.")
    run_pagerank(adata, connectivities_key, adata.obs[gene_count_key])
    # Get cells
    if candidate_clusters is None :
        candidate_clusters = choose_clusters_for_cell_root(adata)
    which_cells = adata.obs_vector('leiden').isin(candidate_clusters)
    selected_cells = adata.obs.index[which_cells]
    selected_page_rank = adata.obs_vector('page_rank')[which_cells]
    top_index = np.flip(np.argsort(selected_page_rank))[0]
    selected_cell_index = selected_cells[top_index]
    return (selected_cell_index, selected_page_rank[top_index])

def choose_clusters_for_cell_root(adata: AnnData) :
    sc.logging.info("Choosing clusters for cell root...")
    # all clusters
    clusters = adata.obs_vector('leiden').unique()
    clusters_pagerank_medians = pd.Series()
    page_ranks = adata.obs_vector('page_rank')
    leidens = adata.obs_vector('leiden')
    for cluster in clusters :
        which_cells = leidens.isin([cluster])
        selected_page_rank = page_ranks[which_cells]
        clusters_pagerank_medians[cluster] = np.median(selected_page_rank)
    # sort it
    clusters_pagerank_medians = clusters_pagerank_medians.sort_values(ascending=False)
    # print(clusters_pagerank_medians)
    # Check the first two clusters and see we should merge them based on t-test
    first_clusters, second_clusters = list(), list()
    first_clusters.append(clusters_pagerank_medians.index[0])
    current_next = 1
    second_clusters.append(clusters_pagerank_medians.index[current_next])
    while True :
        which_cells = leidens.isin(first_clusters)
        first_clusters_pageranks = page_ranks[which_cells]
        which_cells = leidens.isin(second_clusters)
        second_cluster_pageranks = page_ranks[which_cells]
        # Mann–Whitney U test is used to avoid the complexity requirement by t-test
        t_test_result = stats.mannwhitneyu(first_clusters_pageranks, second_cluster_pageranks)
        if t_test_result.pvalue < 0.05 : # The default p-value cutoff
            break
        # Merge second_clusters to first_clsutets
        first_clusters.extend(second_clusters)
        current_next += current_next
        # No more. We are done here by collecting all cells
        if current_next == len(clusters_pagerank_medians) :
            break
        second_clusters.clear()
        second_clusters.append(clusters_pagerank_medians.index[current_next])
    sc.logging.info("Found clusters: " + str(first_clusters))
    return first_clusters


def run_pagerank(adata: AnnData,
                 connectivities_key = 'connectivities',
                 scores = None) :
    page_rank_key = 'page_rank'
    if page_rank_key in adata.obs.keys() :
        page_rank = adata.obs_vector('page_rank')
        return dict(zip(np.arange(0, len(page_rank)), page_rank))
    # We will use personalized pagerank for network smooth, which is the same in the original R implementation
    cell_graph = nx.Graph(adata.obsp[connectivities_key])
    if scores is None :
        scores = adata.obs_vector('n_genes_by_counts')
    cell_scores = dict(zip(np.arange(0, len(scores)), scores))
    # Cannot work with pargerank_scipy. So use this method.
    sc.logging.info("Conducting pagerank...")
    page_ranks = nx.pagerank_numpy(cell_graph, alpha=0.85, personalization=cell_scores)
    # values are sorted based on nodes' indices from 0 to total_cells. So we can use it directly.
    adata.obs['page_rank'] = page_ranks.values()
    return page_ranks

# The following is test code and should be removed
dir_17_5 = "/Users/wug/Documents/missy_single_cell/seq_data_v2/17_5_gfp/filtered_feature_bc_matrix"
dir_12_5 = "/Users/wug/Documents/missy_single_cell/seq_data_v2/12_5_gfp/filtered_feature_bc_matrix"
# adata_17_5 = open_10_genomics_data(dir_17_5)
# adata_17_5 = preprocess(adata_17_5)
# cluster(adata_17_5)
# dpt(adata_17_5, 'TTGACCCGTTAGCGGA-1')
# sc.pl.paga_path(adata_17_5, adata_17_5.obs['leiden'].sort_values(), ['n_genes_by_counts'])
# cytotrace(adata_17_5)
# adata_merged = project(dir_12_5, adata_17_5)
# sc.pl.umap(adata_merged, color = ('leiden', 'n_genes_by_counts', 'cytotrace', 'batch'))
# sc.pl.paga(adata_17_5, pos = adata_17_5.uns['paga']['pos'], add_pos = False, color = ('leiden', 'n_genes_by_counts', 'cytotrace'))
# sc.pl.umap(adata_17_5, color = ('leiden', 'n_genes_by_counts'))
# adata_12_5 = open_10_genomics_data(dir_12_5)
# adata_12_5 = preprocess(adata_12_5)
# cluster(adata_12_5)
# sc.pl.umap(adata_12_5, color = ('leiden', 'n_genes_by_counts'))
# adata_merged = project(dir_12_5, adata_17_5)
# sc.pl.umap(adata_merged, color = ('leiden', 'n_genes_by_counts', 'batch'))


