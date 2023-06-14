"""
This package is used to do pathway enrichment analysis by adoptiong ccommonly used algorithms.
"""
import cProfile
import os
import time
from collections import OrderedDict
from logging import info
from typing import Optional, Union, List

import numpy as np
import pandas as pd
import scanpy as sc
import statsmodels.api as sm
import statsmodels.stats.multitest
from anndata import AnnData
from gseapy.gsea import SingleSampleGSEA
from joblib import Parallel, delayed
from pyscenic import aucell
from pyscenic.genesig import GeneSignature
from seaborn import clustermap
from statsmodels.formula.api import ols

# Pathway analysis keys
SSGSEA_KEY = "X_ssgsea"
AUCELL_KEY = "X_aucell"
TF_SSGSEA_KEY = "X_tf_ssgsea"
TF_AUCELL_KEY = "X_tf_aucell"

def _load_data(adata: Union[AnnData, str],
               need_transpose: bool = True):
    """
    Just a simple utility method to load the expression data if needed.
    :param adata:
    :return:
    """
    # Load the data if necessary
    if isinstance(adata, str) and os.path.isfile(adata):
        adata = sc.read_h5ad(adata)
    if isinstance(adata, AnnData):
        # Convert the adata into a DataFrame
        adata_df = adata.to_df()  # Use the expression matrix and make sure genes for rows and cells for cols
        if need_transpose:
            adata_df = adata_df.T
    elif isinstance(adata, pd.DataFrame):
        adata_df = adata
    else:
        raise Exception("Cannot handle parameter adata.")
    return adata, adata_df


def _load_reactome_gmt(reactome_gmt: Union[dict, str]):
    """
    Another help method to load the Reactome GMT file
    :param reactome_gmt:
    :return:
    """
    # Load the GMT file if the passed is a file name
    if isinstance(reactome_gmt, dict):
        genesets_dict = reactome_gmt
    elif os.path.isfile(reactome_gmt):
        genesets_dict = _read_reactome_gmt(reactome_gmt)
    else:
        raise Exception("Error in parsing reactome_gmt: Need a dict or file name.")
    return genesets_dict


# The following function is modified from VEGA: https://github.com/LucasESBS/vega/blob/c74b8582148126930c4a3405c9cde4e0f2769dda/vega/utils.py#L124
# The function copied here to avoid the dependency to VEGA.
def _read_reactome_gmt(fname):
    """
    Read GMT file into an OrderedDict.

    Args:
        fname (str): Path to gmt file
        sep (str): Separator used to read gmt file.
        min_g (int): Minimum of gene members in gene module.
        max_g (int): Maximum of gene members in gene module.
    Returns:
        OrderedDict: Dictionary of gene_module:genes.
    """
    dict_pathway = OrderedDict()
    with open(fname) as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            val = line.split("\t")
            dict_pathway[val[0]] = val[2:]
    return dict_pathway


class PathwayAnalyzer(object):
    """
    Base class for holding some common functions.
    """

    def __init__(self):
        self.adata_key = "pathway"  # hold the key to be used in anndata
        self.pathway_list_key = "PATHWAY_LIST_KEY"
        self.adata = None

    def scale_pathways(self,
                       pathway_scores: pd.DataFrame) -> pd.DataFrame:
        """
        Perform scale for individual values bewteen 0 and 1.
        :param pathway_scores:
        :return:
        """
        rtn = (pathway_scores - pathway_scores.min()) / (pathway_scores.max() - pathway_scores.min())
        return rtn

    def perform_pathway_anova(self,
                              cluster_key: str = "leiden") -> pd.DataFrame:
        """
        Perform one-way anova for pathway activities across cell clusters.
        The implementation of this function is based on https://www.pythonfordatascience.org/anova-python/.
        """
        self._check_adata()
        pathways_activities = self._get_pathway_activities()
        pathways_list = list(pathways_activities.columns)
        if cluster_key not in self.adata.obs.keys():
            raise Exception("{} is not in obs. Do a clustering first.".format(cluster_key))
        cell_clusters = self.adata.obs[cluster_key]
        # Do annova for each pathway
        results_dict = {}
        for i in range(len(pathways_list)):
            pathway = pathways_list[i]
            # print("Working for {}".format(pathway))
            pathway_scores = pathways_activities.iloc[:, i]
            df = pd.DataFrame({"cluster": cell_clusters,
                               "pathway": pathway_scores})
            model = ols("pathway ~ C(cluster)", data=df).fit()
            aov_table = sm.stats.anova_lm(model, typ=2)
            # Escape in case nothing there
            if np.isnan(aov_table["F"][0]):
                continue
            results_dict[pathway] = [aov_table["F"][0], aov_table["PR(>F)"][0]]
            # if i == 50:
            #     break
        result_df = pd.DataFrame(results_dict).T
        result_df.columns = ["F", "PR(>F)"]
        # Sort based on p-value
        result_df.sort_values(by=['PR(>F)', "F"], inplace=True)
        # Do a multiple test correction
        adjust_method = "fdr_bh"
        adjusted_p_values = statsmodels.stats.multitest.multipletests(result_df['PR(>F)'], method=adjust_method)
        result_df[adjust_method] = adjusted_p_values[1]
        return result_df

    def _check_adata(self):
        # Get all necessary data
        if self.adata is None:
            raise Exception("No adata for analysis!")
        if not isinstance(self.adata, sc.AnnData):
            raise Exception("adata in the object is not an AnnData")

    def _get_pathway_activities(self) -> pd.DataFrame:
        pathways_activities = self.adata.obsm[self.adata_key]
        if isinstance(pathways_activities, pd.DataFrame):
            return pathways_activities
        raise Exception("{} for pathways_activities is not supported".format(type(pathways_activities)))

    def get_pathway_score(self,
                          pathway: str
                          ) -> pd.Series:
        """
        Get the trained pathway score
        :param adata:
        :param pathway:
        :return:
        """
        self._check_adata()
        return self._get_pathway_activities().loc[:, pathway]

    def color_pathway_umap(self,
                           pathways: list or str,
                           color_cluster: bool = True,
                           obs_name_prefix: str = ''
                           ) -> None:
        """
        Color cells in the umap plot based on the inferred pathway scores from trained VEGA.
        :param adata:
        :param pathway:
        :param color_cluster: check if clusters should be plotted too.
        :param obs_name_prefix: For TF, we need to add "tf_" before the obs name to avoid conflict with
        the original gene names.
        :return:
        """
        self._check_adata()
        if isinstance(pathways, str):
            pathways = [pathways]
        pathway_activities = self._get_pathway_activities()
        pathway_list = list(pathway_activities.columns)
        color_names = []
        for pathway in pathways:
            pathway_index = pathway_list.index(pathway)  # Get the first pathway index. There should be one only.
            pathway_score = pathway_activities.iloc[:, pathway_index]
            # Need to register it into the observation for plot
            key_name = obs_name_prefix + pathway
            self.adata.obs[key_name] = pathway_score
            color_names.append(key_name)
        if color_cluster:
            color_names.append("leiden")
        sc.pl.umap(self.adata, color=color_names)


class ReactomeSSGSEA(SingleSampleGSEA):
    """
    A customized ssGSEA application with the hope of improved performance by using pandas index-based
    isIn method. For 100 samples with around 4,000 genes, the original running time is about 300 seconds.
    With this customized code, the running time has decreased to 10 seconds.
    The code was modified from the original gseapy code hosted at https://github.com/zqfang/GSEApy. The original
    author, Dr. Zhuoqing Fang, has the copyright to the original code. The license of the original code is
    BSD 3-Clause "New" or "Revised" License: https://github.com/zqfang/GSEApy/blob/master/LICENSE.
    """
    def __init__(self, data, gene_sets, sample_norm_method='rank',
                 min_size=15, max_size=2000, permutation_num=0, weighted_score_type=0.25,
                 scale=True, ascending=False, processes=1, figsize=(7,6), format='pdf',
                 graph_num=20, no_plot=True, seed=None, verbose=False):
        super(ReactomeSSGSEA, self).__init__(data, gene_sets, outdir="GSEA_SingleSample", sample_norm_method='rank',
                 min_size=15, max_size=2000, permutation_num=0, weighted_score_type=0.25,
                 scale=True, ascending=False, processes=1, figsize=(7,6), format='pdf',
                 graph_num=20, no_plot=True, seed=None, verbose=False)

    def load_gmt(self,
                 gene_list: pd.Index,
                 gmt):
        """
        Re-implement to use isin for fast performance
        :param gene_list:
        :param gmt:
        :return:
        """
        """load gene set dict"""
        if isinstance(gmt, dict):
            genesets_dict = gmt
        else:
            raise Exception("Error parsing gmt parameter for gene sets")

        subsets = list(genesets_dict.keys())
        self.n_genesets = len(subsets)
        for subset in subsets:
            subset_list = genesets_dict.get(subset)
            if isinstance(subset_list, set):
                subset_list = list(subset_list)
                genesets_dict[subset] = subset_list
            # tag_indicator = np.in1d(gene_list, subset_list, assume_unique=True)
            tag_indicator = gene_list.isin(subset_list)
            tag_len = tag_indicator.sum()
            if self.min_size <= tag_len <= self.max_size: continue
            del genesets_dict[subset]

        filsets_num = len(subsets) - len(genesets_dict)
        self._logger.info("%04d gene_sets have been filtered out when max_size=%s and min_size=%s" % (
        filsets_num, self.max_size, self.min_size))

        if filsets_num == len(subsets):
            self._logger.error("No gene sets passed through filtering condition!!!, try new parameters again!\n" + \
                               "Note: check gene name, gmt file format, or filtering size.")
            raise Exception("No gene sets passed through filtering condition")

        self._gmtdct = genesets_dict
        return genesets_dict

    def enrichment_score_tensor(self,
                                data: pd.Series,
                                gene_sets,
                                weighted_score_type,
                                scale=False):
        """ A local implementation to support ssgsea only using Pandas for quick performance based on
            pySCENIC AUCell code. The code is simplied from the original gseapy code.
            Next generation algorithm of GSEA and ssGSEA.

            :param data: ranked gene expression values indexed by gene names
            :param dict gene_sets:  gmt file dict.
            :param float weighted_score_type:     weighting by the correlation.
                                    options: 0(classic), 1, 1.5, 2. default:1 for GSEA and 0.25 for ssGSEA.
            :param bool scale:      If True, normalize the scores by number of genes_mat.
            :return: a tuple contains::

                     | ES: Enrichment score (real number between -1 and +1), for ssGSEA, set scale eq to True.
                     | ESNULL: Enrichment score calculated from random permutation.
                     | Hits_Indices: Indices of genes if genes are included in gene_set.
                     | RES: The running enrichment score for all locations in the gene list.

        """
        # gene_mat -> 1d: prerank, ssSSEA or 2d: GSEA
        keys = sorted(gene_sets.keys())
        cor_mat = data.values
        gene_mat = data.index.values
        if weighted_score_type == 0:
            # don't bother doing calcuation, just set to 1
            cor_mat = np.ones(cor_mat.shape)
        elif weighted_score_type > 0:
            pass
        else:
            raise ValueError("Using negative values of weighted_score_type, not allowed.")

        cor_mat = np.abs(cor_mat)
        # Used for ssGSEA only
        if cor_mat.ndim != 1:
            raise ValueError("The dimension of cor_mat must be equal to 1.")

        # ssGSEA or Prerank
        # genestes->M, genes->N
        N, M = len(gene_mat), len(keys)
        # generate gene hits matrix
        # It is very slow to use np.in1d. Use the index.isin instead for quick performance.
        # tag_indicator = np.vstack([np.in1d(gene_mat, gene_sets[key], assume_unique=True) for key in keys])
        # Only this line is changed. All others are the same as the original code but focusing on ssgsea only.
        tag_indicator = np.vstack([data.index.isin(gene_sets[key]) for key in keys])
        tag_indicator = tag_indicator.astype(int)
        # index of hits
        hit_ind = [np.flatnonzero(tag).tolist() for tag in tag_indicator]
        # generate permutated hits matrix
        nperm = 0  # For ssgsea
        perm_tag_tensor = np.repeat(tag_indicator, nperm + 1).reshape((M, N, nperm + 1))
        # missing hits
        no_tag_tensor = 1 - perm_tag_tensor
        # calculate numerator, denominator of each gene hits
        rank_alpha = (perm_tag_tensor * cor_mat[np.newaxis, :, np.newaxis]) ** weighted_score_type

        axis = 1
        P_GW_denominator = np.sum(rank_alpha, axis=axis, keepdims=True)
        P_NG_denominator = np.sum(no_tag_tensor, axis=axis, keepdims=True)
        REStensor = np.cumsum(rank_alpha / P_GW_denominator - no_tag_tensor / P_NG_denominator, axis=axis)
        # ssGSEA: scale es by gene numbers ?
        # https://gist.github.com/gaoce/39e0907146c752c127728ad74e123b33
        if scale: REStensor = REStensor / len(gene_mat)
        # ssGSEA
        esmatrix = REStensor.sum(axis=axis)
        es, esnull, RES = esmatrix[:, -1], esmatrix[:, :-1], REStensor[:, :, -1]
        return es, esnull, hit_ind, RES

    def run(self):
        """run entry"""
        self._logger.info("Parsing data files for ssGSEA...........................")
        # load data
        data = self.load_data()
        # normalized samples, and rank
        normdat = self.norm_samples(data)
        # filtering out gene sets and build gene sets dictionary
        gmt = self.load_gmt(gene_list=normdat.index, gmt=self.gene_sets)
        self._logger.info("%04d gene_sets used for further statistical testing....."% len(gmt))
        # set cpu numbers
        self._set_cores()
        # start analysis
        self._logger.info("Start to run ssGSEA...Might take a while................")
        # ssGSEA without permutation
        self.runSamples(df=normdat, gmt=gmt)

    def runSamples(self, df, gmt):
        """
        Provide a quicker implementation based on pySCENIC AUCell using pandas. The code
        was simplied from the original implementation with the replacement of a new
        wrapped method.
        :param df:
        :param gmt:
        :return:
        """
        # df.index.values are gene_names
        # Save each sample results to odict
        self.resultsOnSamples = OrderedDict()
        # run ssgsea for gct expression matrix
        # multi-threading
        subsets = sorted(gmt.keys())
        tempdat = []
        names = []
        # name should be cell id and ser is named gene expression ranks
        for name, ser in df.iteritems():
            # prepare input
            dat = ser.sort_values(ascending=self.ascending)
            tempdat.append(dat)
            names.append(name)

        tempes = Parallel(n_jobs=self._processes)(
            delayed(self.enrichment_score_tensor)(
                dat,
                gmt,
                self.weighted_score_type,
                self.scale)
            for dat in tempdat)

        # save results
        for i, temp in enumerate(tempes):
            name = names[i]  # cell id
            self._logger.debug("Calculate Enrichment Score for Sample: %s " % name)
            es, esnull, hit_ind, RES = temp
            # save results
            self.resultsOnSamples[name] = pd.Series(data=es, index=subsets, name=name)


class GSEAPyWrapper(PathwayAnalyzer):
    """
    A wrapper of GSEAPY for doing ssgsea pathway analysis: https://gseapy.readthedocs.io/en/latest/index.html
    """
    def __init__(self):
        super().__init__()
        self.adata_key = SSGSEA_KEY

    def ssgsea(self,
               adata: Union[str, AnnData],
               reactome_gmt: Union[str, dict],
               data_key: str,
               weigted_score_type: float = 0.25,
               n_processes: int = 1
               ) -> Optional[pd.DataFrame]:
        """
        Perform ssgsea analysis
        :param adata: an AnnData object or a hd5 file
        :param reactome_gmt:
        :param data_key
        :param weigted_score_type
        :param n_processes
        :return:
        """
        if data_key is not None:
            self.adata_key = data_key
        adata, adata_df = _load_data(adata)
        genesets_dict = _load_reactome_gmt(reactome_gmt)

        ss = ReactomeSSGSEA(adata_df,
                            genesets_dict,
                            weighted_score_type=weigted_score_type,
                            scale=False,
                            processes=n_processes)
        ss.run()
        # ss = ssgsea_runner.ssgsea(adata_df,
        #                genesets_dict,
        #                weighted_score_type=weigted_score_type,
        #                scale=False,
        #                processes=n_processes,
        #                outdir=None)  # Don't want to keep the output directory
        # Keep the raw enrichment scores
        if isinstance(adata, AnnData):
            adata.obsm[self.adata_key] = pd.DataFrame(ss.resultsOnSamples).T  # Need to revert it backs
            self.adata = adata
        else:
            return pd.DataFrame(ss.resultsOnSamples)


class AUCellWrapper(PathwayAnalyzer):
    """
    The wrapper for doing AUCell based pathway analysis based on pyScenic: https://pyscenic.readthedocs.io/en/latest/
    """
    def __init__(self):
        super(AUCellWrapper, self).__init__()
        self.adata_key = AUCELL_KEY

    def aucell(self,
               adata: Union[AnnData, str],
               reactome_gmt: Union[OrderedDict, str],
               data_key: str,
               filter_with_max_score: float = 0.0001, # Default filter to avoid generating too many zeros.
               need_plot: bool = False,
               need_scale: bool = True
               ) -> Optional[pd.DataFrame]:
        """
        Perform aucell-based pathway analysis. The code here is based on Example 2: Gene signatures from a GMT
        file from https://github.com/aertslab/pySCENIC/blob/master/notebooks/pySCENIC%20-%20AUCell%20example.ipynb.
        :param adata:
        :param reactome_gmt:
        :param filter_with_max_score: pathway cols with this max values will be filtered out
        :param need_plot: true to generate a hierarchical clustering map
        :param need_scale: true for scale for individual pathways
        :return:
        """
        if data_key is not None:
            self.adata_key = data_key
        adata, adata_df = _load_data(adata, need_transpose=False)
        genesets_dict = _load_reactome_gmt(reactome_gmt)
        # Need to convert genesets_dict to list of GeneSignatures
        gene_signatures = [GeneSignature(name=name, gene2weight=genes) for name, genes in genesets_dict.items()]
        # Get the percentiles for genes that should be used
        percentiles = aucell.derive_auc_threshold(adata_df)
        # Choose top 5 percentiles genes. This may need to be parameterized.
        aucs_matrix = aucell.aucell(adata_df, gene_signatures, percentiles[0.05])
        if filter_with_max_score is not None: # Do a filtering
            info("The size of aucx_matrix before filtering: {}".format(aucs_matrix.shape))
            aucs_matrix = aucs_matrix.loc[:,
                          [col for col in aucs_matrix.columns if aucs_matrix.loc[:, col].max() > filter_with_max_score]]
            info("The size of aucx_matrix after filtering: {}".format(aucs_matrix.shape))
        if need_scale:
            aucs_matrix = super().scale_pathways(aucs_matrix)
        if need_plot:
            clustermap(aucs_matrix, figsize=(14, 14)) # Give a larger figure size
        # Store the information
        adata.obsm[self.adata_key] = aucs_matrix
        self.adata = adata
        return aucs_matrix


def reactome_ssgsea(adata: Union[AnnData, str],
                    reactome_gmt: Union[dict, str],
                    data_key: str = SSGSEA_KEY) -> PathwayAnalyzer:
    """
    Perform ssgsea analysis based on Reactome pathways
    :param adata:
    :param reactome_gmt:
    :return:
    """
    gsea_wrapper = GSEAPyWrapper()
    gsea_wrapper.ssgsea(adata, reactome_gmt, data_key)
    return gsea_wrapper


def reactome_aucell(adata: Union[AnnData, str],
                    reactome_gmt: Union[dict, str],
                    data_key: str = AUCELL_KEY,
                    filter_with_max_score: float = 1.0E-4, # Default with 0.0001
                    need_plot: bool = False) -> PathwayAnalyzer:
    """
    Perform aucell-based pathway analysis for Reactome pathways.
    """
    aucell_wrapper = AUCellWrapper()
    aucell_wrapper.aucell(adata,
                          reactome_gmt,
                          data_key,
                          filter_with_max_score=filter_with_max_score,
                          need_plot=need_plot)
    return aucell_wrapper


def fetch_analyzed_pathway_keys(adata : AnnData) -> List[str]:
    """
    Get the keys for analyzed pathway activities.
    :param adata:
    :return:
    """
    rtn_keys = list()
    if AUCELL_KEY in adata.obsm.keys():
        rtn_keys.append(AUCELL_KEY)
    if SSGSEA_KEY in adata.obsm.keys():
        rtn_keys.append(SSGSEA_KEY)
    return rtn_keys


def fetch_cluster_pathway_activities(adata: AnnData,
                                     cluster: str,
                                     pathway_activity_key: str = [SSGSEA_KEY, AUCELL_KEY],
                                     cluster_key: str = "leiden") -> pd.Series:
    """
    Fetch the median of pathway activities for a cluster
    :param adata:
    :param cluster_key: the cluster algorithm name
    :param cluster: the cluster to be selected
    :param pathway_activity_key:
    :return:
    """
    # Make sure pathway has been analyzes
    if pathway_activity_key not in adata.obsm.keys():
        raise ValueError("Cannot find key in adata.obsm for {}.".format(pathway_activity_key))
    # Make sure cluster_key is in there
    if cluster_key not in adata.obs.keys():
        raise ValueError("Cannot find key in adata.obs for {}".format(cluster))
    # Get the cells
    cell_cluster = adata.obs[cluster_key] == str(cluster) # Just in case cluster is an int
    cluster_pathways = adata.obsm[pathway_activity_key].loc[cell_cluster, :]
    cluster_pathway_median = cluster_pathways.median()
    return cluster_pathway_median


def load_dorothea_data(file_name: str,
                       confidence_levels: str = ['A','B']):
    """
    Load the Dorothea TF/target interactions data into a dict so that it can be used for
    pathway activities analysis.
    :param file_name:
    :param confidence_levels:
    :return:
    """
    file = open(file_name)
    tf_genes = dict()
    if isinstance(confidence_levels, str):
        confidence_levels = [confidence_levels]
    # escape the first line
    line = file.readline()
    while True:
        line = file.readline()
        if not line:
            break
        # tf\tconfidence\ttarget\tmor
        tokens = line.split("\t")
        if not tokens[1] in confidence_levels:
            continue
        if not tokens[0] in tf_genes.keys():
            tf_genes[tokens[0]] = []
        tf_genes[tokens[0]].append(tokens[2]) # Don't care about the action mode
    file.close()
    return tf_genes


def test_reactome_ssgsea():
    adata = sc.read_h5ad("/Users/wug/Documents/wgm/work/FIPlugIns/test_data/ScRNASeq/SavedResults/mouse/17_5_gfp.h5ad")
    # Just need the first 10 cells for this test
    adata = adata[:100, :]
    reactome_gmt = "../data/gmt/MouseReactomePathways_Rel_75_122220.gmt"
    time0 = time.time()
    pr = cProfile.Profile()
    pr.enable()
    rtn = reactome_ssgsea(adata, reactome_gmt)
    pr.disable()
    pr.dump_stats("ssgsea_new_062321_1.prof")
    time1 = time.time()
    print("Time used for ssgea: {} seconds".format(time1 - time0))
    return rtn


def test_reactome_vega():
    adata = sc.read_h5ad("/Users/wug/Documents/wgm/work/FIPlugIns/test_data/ScRNASeq/SavedResults/mouse/17_5_gfp.h5ad")
    reactome_gmt = "data/gmt/MouseReactomePathways_Rel_75_122220.gmt"
    return reactome_vega(adata, reactome_gmt, n_epochs=3)


def test_reactome_aucell():
    # Need to use all expression data without filtering
    file_name = "/Users/wug/Documents/missy_single_cell/seq_data_v2/17_5_gfp/filtered_feature_bc_matrix"
    adata = sc.read_10x_mtx(file_name, cache=True)
    # adata = sc.read_h5ad("/Users/wug/Documents/wgm/work/FIPlugIns/test_data/ScRNASeq/SavedResults/mouse/17_5_gfp.h5ad")
    reactome_gmt = "data/gmt/MouseReactomePathways_Rel_75_122220.gmt"
    return reactome_aucell(adata, reactome_gmt, 0.01, False)

# test_reactome_ssgsea()
# test_reactome_aucell()