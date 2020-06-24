"""
This JSONRPC server is used to provide server wrapper aroung some major functions provided by scanpy and its related Python packages.
These functions are provided as json-based REST calls for the ReactomeFIViz app.

"""
# Required by scanpy
# import numpy as np
# import pandas as pd
import scanpy as sc
from scanpy.plotting import _utils
from anndata import AnnData
from typing import Optional
import numpy as np
from numpy import ndarray

# As suggested in the tutorial by scanpy doc: Preprocessing and clustering 3k PBMCs
sc.settings.verbosity = 3
sc.logging.print_versions()
sc.set_figure_params(dpi = 80, facecolor = 'white')

# Global level flags to control the use of the serve
# server = SimpleJSONRPCServer("localhost")
# isWaiting = True

# This dict is used for cache loaded AnnaData
cache = dict()
KEY_DATA, KEY_PROCESSED = 'adata', 'processed'

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

def preprocess(adata, copy=False, need_scale=True):
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
    # TODO: We may need to send out the above qc results to ReactomeFIViz for users to investigate 
    # and then pick up sensible parameters for test.
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    # Mark genes with highly_variable flag
    # Default parameters are used in this function call.
    sc.pp.highly_variable_genes(adata)
    if not need_scale :
        return adata
    # Keep the normalized data for future use
    adata.raw = adata
    # Keep the highly variable genes only
    adata = adata[:, adata.var.highly_variable]
    # Default using total_counts though pct_counts_mt is also used i the tutorial
    # TODO: Need to consider to add regress_out for cell cycle score!
    sc.pp.regress_out(adata, keys=['total_counts'])
    # gene-wise z-score
    sc.pp.scale(adata, max_value=10)
    return adata

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
    sc.tl.paga(adata)
    # For some reason, each run may produce a little bit different layout
    sc.pl.paga(adata, plot=False) # No need to plot. pl: plot
    # Do umap
    sc.tl.umap(adata, init_pos='paga')
    adata.uns['paga']['pos'] = reset_paga_pos(adata)
    # Finding markers genes. Use t-test here for the time being.
    # TODO: This needs to be improved in the future. Two more iterations:
    # 1). Expose the method in the function definition so that the user can choose what method should be used.
    # 2). Use a more sophisticated method: e.g. diffx and something based on ZNBM.
    sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
    # tl.umap call should produce the coordinates. If we just need coordinates, there is no need to call the following
    if plot:
        sc.pl.umap(adata, color=['leiden'])
    # There is no need to return adata

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

def project(adata: AnnData,
            adata_ref: AnnData,
            batch_categories: Optional[str] = None,
            copy: bool = False
)->AnnData:
    """
    This function is used to project adata to adata_ref using ingest function in scanpy. The procedures here are based
    on the tutorial: integration-data-using-inject. Both adata ad adata_ref should have been preprocessed by running the
    preprocess function but with need_scale false.
    :param adata: the data to be projected
    :param adata_ref: the reference data
    :param batch_categories: the name for the reference should be the first element
    :param copy: true to copy the passed AnnData objects. Otherwise, the original objects will be sliced.
    :return: return a merged AnnData object with two originally data copied and merged.
    """
    # Check if leiden has been conducted.
    if copy :
        adata = adata.copy()
        adata_ref = adata_ref.copy()
    if 'leiden' not in adata_ref.obs :
        cluster(adata_ref)
    # Make sure both of them have the same set of variables
    shared_var_names = adata_ref.var_names.intersection(adata.var_names)
    # slicing the data
    adata = adata[:, shared_var_names]
    adata_ref = adata_ref[:, shared_var_names]
    # inject based on the leiden
    sc.tl.ingest(adata, adata_ref, obs = 'leiden')
    # merge two objects
    if not batch_categories :
        batch_categories = ['ref', 'new']
    adata_merged = adata_ref.concatenate(adata, batch_categories = batch_categories)
    return adata_merged

# The following is test code and should be removed
dir_17_5 = "/Users/wug/Documents/missy_single_cell/seq_data_v2/17_5_gfp/filtered_feature_bc_matrix"
dir_12_5 = "/Users/wug/Documents/missy_single_cell/seq_data_v2/17_5_gfp/filtered_feature_bc_matrix"
#adata_17_5 = open_10_genomics_data(dir_17_5)
#adata_12_5 = open_10_genomics_data(dir_12_5)
# In python we need to define a function first before we can call it.
#adata = open_10_genomics_data(dir)
#print(adata)
# sc.pl.paga_compare()