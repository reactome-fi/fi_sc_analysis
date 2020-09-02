"""
Provide some functions for mutation information calculation based on Scribe.

"""
from anndata import AnnData
import numpy
from numpy import ndarray
from scipy import stats
# Turn of MI-based relationship calculation from Scribe because of the performance reason.
# from Scribe import information_estimators as iet
import math  # check float isNaN

DEBUG = False
# Minumum values for calculation gene relationships
MIN_VALUE_NUMBER = 12


def calculate_gene_relations(gene_pairs: str,
                             adata: AnnData,
                             cell_time_key: str,
                             layer=None,
                             delay_window=7,
                             mode='spearman') -> dict:
    """
    Calculate relationships for a passed list of gene pairs.
    :param gene_pairs: list of gene pairs with two genes in the pair tab-delimited
    :param adata:
    :param cell_time_key:
    :param layer:
    :param delay_window: the cell time window used to check relationships between two vectors
    :param mode: currently it supports the following types: spearman, pearson, kendal,
    mi (mutual information), rdi (restricted directed information). The default is spearman
    since it runs much faster than others.
    :return:
    """
    # Make sure the basic required information is there
    if DEBUG:
        print(adata)
    if cell_time_key not in adata.obs_keys():
        raise ValueError(cell_time_key + " is not in the data's obs_keys.")
    cell_time_sorted_index = adata.obs_vector(cell_time_key).argsort()
    if delay_window is None:
        delay_window = 7
    if mode is None:
        mode = 'spearman'
    # Get a list of gene pairs
    rtn = {}
    delays = [i for i in range(1, delay_window + 1)]
    for gene_pair in gene_pairs:
        if DEBUG:
            print(gene_pair)
        genes = gene_pair.split("\t")
        # expected two genes
        if len(genes) < 2:
            raise ValueError("gene_pair has less than 2 genes: " + gene_pair)
        # Have to make sure both of these genes are in the var_vectors
        if genes[0] not in adata.var_names or genes[1] not in adata.var_names:
            continue
        value1 = adata.obs_vector(genes[0], layer=layer)[cell_time_sorted_index]
        value2 = adata.obs_vector(genes[1], layer=layer)[cell_time_sorted_index]
        # Make sure we have enough NaN values
        if len(value1[~numpy.isnan(value1)]) < MIN_VALUE_NUMBER or len(value2[~numpy.isnan(value2)]) < MIN_VALUE_NUMBER:
            continue
        rel = None
        if mode in ('spearman', 'kendall', 'pearson'):
            rel = _calculate_cor(value1, value2, delays, mode)
        elif mode == 'mi':
            rel = _calculate_mi(value1, value2, delays)
        elif mode == 'rdi':
            rel = _calculate_rdi(value1, value2, delays)
        if rel is not None:
            rtn[gene_pair] = rel
    return rtn


def _calculate_cor(value1: ndarray,
                   value2: ndarray,
                   delays=(1, 2, 3, 4, 5, 6, 7),
                   mode='spearman') -> [float, float]:
    """
    :param value1:
    :param value2:
    :param delays:
    :param mode:
    :return: two values the first is the correlation and the second pvale.
    """
    cor = None
    for delay in delays:
        x = value1[:-delay]
        y = value2[delay:]
        temp = None
        if mode == 'spearman':
            temp = stats.spearmanr(x, y)
        elif mode == 'kendall':
            temp = stats.kendalltau(x, y)
        elif mode == 'pearson':
            temp = stats.pearsonr(x, y)
        if temp is None:
            continue
        if DEBUG:
            print("{}: {}.".format(delay, temp))
        if math.isnan(temp[0]) or math.isnan(temp[1]):
            continue  # Don't return NaN
        if cor is None or abs(temp[0]) > abs(cor[0]):
            cor = temp
    return cor


def _calculate_mi(value1: ndarray,
                  value2: ndarray,
                  delays=(1, 2, 3)) -> float:
    # Required by scribe to use list
    x_ori = [[i] for i in value1]
    y_ori = [[i] for i in value2]
    mi = 0.0  # This should be the minimum
    for delay in delays:
        # The following code is based on Scribe.causal_network
        x = x_ori[:-delay]
        y = y_ori[delay:]
        # temp = iet.mi(x, y)
        # print(temp)
        # if temp > mi:
        #    mi = temp
    return mi


def _calculate_rdi(value1: ndarray,
                   value2: ndarray,
                   delays=(1, 2, 3)) -> float:
    # Required by scribe to use list
    x_ori = [[i] for i in value1]
    y_ori = [[i] for i in value2]
    rdi = 0.0  # This should be the minimum
    for delay in delays:
        # The following code is based on Scribe.causal_network
        x = x_ori[:-delay]
        y_minux_1 = y_ori[delay - 1:-1]
        y = y_ori[delay:]
        # temp = iet.cmi(x, y, y_minux_1)
        # if temp > rdi :
        #    rdi = temp
    return rdi


# For test
gene_pairs = ['Cps1\tBicc1', 'Prom1\tMuc4', 'Cic\tArhgef2']
# import scanpy as sc
# adata = sc.read("/Users/wug/temp/17_5_gfp_velocity_dynamic.h5ad")
