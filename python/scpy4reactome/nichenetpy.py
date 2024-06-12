# This script wraps some essential functions in NicheNet, an R package, so that we may be able to perform
# NicheNet analysis in the Python-based workflows. The calling of R functions is based on Python package
# rpy2: https://rpy2.github.io.

from rpy2.robjects.packages import importr
from rpy2.robjects.lib.dplyr import DataFrame
from rpy2.robjects import rl, default_converter, conversion
from rpy2.robjects.pandas2ri import converter
from rpy2.robjects.vectors import StrVector

import scanpy as sc
import pandas as pd
import seaborn as sns
import networkx as nx

r_base = importr('base')
nichenet = importr('nichenetr')

# Make sure all nichenet RDS files are in this folder for both human and mouse.
NICHENET_DATA_DIR = '../resources/nichenet/'

class NicheNetWrapper:

    def __init__(self,
                 species: str,
                 nichenet_data_dir: str = NICHENET_DATA_DIR):
        self.species = species
        self.nichenet_data_dir = nichenet_data_dir
        self._init()

    def _init(self):
        """Initialize the network information
        """
        if self.species is None:
            raise ValueError('Species has not defined!')
        if self.species == 'mouse':
            self.lr_network = r_base.readRDS(
                self.nichenet_data_dir + "lr_network_mouse_21122021.rds")
            self.gr_network = r_base.readRDS(
                self.nichenet_data_dir + "gr_network_mouse_21122021.rds")
            self.sig_network = r_base.readRDS(
                self.nichenet_data_dir + "signaling_network_mouse_21122021.rds")
            self.ligand_target_matrix = r_base.readRDS(
                self.nichenet_data_dir + "ligand_target_matrix_nsga2r_final_mouse.rds")
            self.ligand_tf_matrix = r_base.readRDS(
                self.nichenet_data_dir + 'ligand_tf_matrix_nsga2r_final_mouse.rds')
            self.weighted_networks = r_base.readRDS(
                self.nichenet_data_dir + "weighted_networks_nsga2r_final_mouse.rds")
        elif self.species == 'human':
            self.lr_network = r_base.readRDS(
                self.nichenet_data_dir + "lr_network_human_21122021.rds")
            self.gr_network = r_base.readRDS(
                self.nichenet_data_dir + "gr_network_human_21122021.rds")
            self.sig_network = r_base.readRDS(
                self.nichenet_data_dir + "signaling_network_human_21122021.rds")
            self.ligand_target_matrix = r_base.readRDS(
                self.nichenet_data_dir + "ligand_target_matrix_nsga2r_final_human.rds")
            self.ligand_tf_matrix = r_base.readRDS(
                self.nichenet_data_dir + 'ligand_tf_matrix_nsga2r_final_human.rds')
            self.weighted_networks = r_base.readRDS(
                self.nichenet_data_dir + "weighted_networks_nsga2r_final_human.rds")
        else:
            raise ValueError('Species need to be either mouse or human!')

        # There are some duplcations in the lr_network matrix. The following is used to remove duplication
        self.lr_network = DataFrame(
            self.lr_network).distinct(rl('from'), rl('to'))
        # Add a new network
        self.weighted_networks_lr_sig = DataFrame(self.weighted_networks.rx2(
            'lr_sig')).inner_join(self.lr_network, by=StrVector(['from', 'to']))

    def rpy_2_pd_df(self, rpy_df: DataFrame) -> pd.DataFrame:
        with (default_converter + converter).context():
            pd_from_r_df = conversion.get_conversion().rpy2py(rpy_df)
        return pd_from_r_df
    

    def convert_matrix_pd(self, m):
        col2value = {}
        index = 1 # Start with index = 1 since it is in R in rx
        for col in m.colnames:
            col2value[col] = m.rx(True, index)
            index += 1
        df = DataFrame(col2value)
        pd_df = self.rpy_2_pd_df(df)
        pd_df.index = m.rownames
        return pd_df
    

    def get_expressed_genes(self,
                            adata: sc.AnnData,
                            obs_value: str | list = None,
                            obs_key: str = 'leiden',
                            pct: float = 0.10) -> list:
        """Get the genes expressed in at least pct of cells.

        Args:
            adata (sc.AnnData): _description_
            obs_value (str | list, optional): _description_. Defaults to None.
            obs_key (str, optional): _description_. Defaults to 'leiden'.
            pct (float, optional): _description_. Defaults to 0.10.
        """
        adata_sub = None
        if not obs_value:
            adata_sub = adata.copy()
        else:
            if isinstance(obs_value, str):
                obs_value = [obs_value]
            adata_sub = adata[adata.obs[obs_key].isin(obs_value)]
        sc.pp.filter_genes(adata_sub, min_cells=int(pct * adata_sub.n_obs))
        return adata_sub.var_names.to_list()

    def get_background_expressed_genes(self, expressed_genes_receiver: list) -> list:
        return list(set(expressed_genes_receiver).intersection(self.ligand_target_matrix.rownames))

    def list_ligands(self) -> list:
        return list(self.lr_network.rx2('from'))

    def list_receptors(self) -> list:
        return list(self.lr_network.rx2('to'))

    def get_potential_ligands(self,
                              expressed_genes_sender: list,
                              expressed_genes_receiver: list) -> set:
        potential_ligands = set()
        for row in self.lr_network.iter_row():
            # print('{} {}'.format(row[0], row[1]))
            _from = row.rx2('from')[0]
            _to = row.rx2('to')[0]
            if _from in expressed_genes_sender and _to in expressed_genes_receiver:
                potential_ligands.add(_from)
        return potential_ligands

    def predict_ligand_activities(self,
                                  geneset: list,
                                  background_expressed_genes: list,
                                  potential_ligands: list) -> pd.DataFrame:
        """Predict ligand activities. This is the main function of nichenet!

        Args:
            geneset (list): This usually should be the list of DEGs
            background_expressed_genes (list): the list of background expressed genes (all) in the receiver cells.
            potential_ligands (list): the ligands to be checked.

        Raises:
            ValueError: In case the object has not been properly initialized.

        Returns:
            pd.DataFrame: a Pandas DataFrame object is returned.
        """
        if self.ligand_target_matrix is None:
            raise ValueError(
                'ligand_target_matrix is none. Cannot run this function before initializing this matrix.')
        # The dataframe in R
        ligand_activities = nichenet.predict_ligand_activities(geneset=StrVector(geneset),
                                                               background_expressed_genes=StrVector(
                                                                   background_expressed_genes),
                                                               ligand_target_matrix=self.ligand_target_matrix,
                                                               potential_ligands=StrVector(potential_ligands))
        # Convert it into Python
        ligand_activities_df = self.rpy_2_pd_df(ligand_activities)
        # Sort
        ligand_activities_df.sort_values(
            by='aupr_corrected', inplace=True, ascending=False)
        ligand_activities_df.reset_index(inplace=True, drop=True)
        return ligand_activities_df

    def get_weighted_ligand_target_links(self,
                                         ligands: list,
                                         geneset: list,
                                         top_n_targets: int = 250) -> pd.DataFrame:
        active_ligand_target_links_df = None
        for ligand in ligands:
            active_ligand_target_links = nichenet.get_weighted_ligand_target_links(ligand,
                                                                                   geneset=StrVector(geneset),
                                                                                   ligand_target_matrix=self.ligand_target_matrix,
                                                                                   n=top_n_targets)
            if active_ligand_target_links_df is None:
                active_ligand_target_links_df = self.rpy_2_pd_df(active_ligand_target_links)
            else:
                tmp_df = self.rpy_2_pd_df(active_ligand_target_links)
                active_ligand_target_links_df = pd.concat([active_ligand_target_links_df, tmp_df], axis=0)
        return active_ligand_target_links_df

    def plot_ligand_target_clustermap(self,
                                      active_ligand_target_links_df: pd.DataFrame,
                                      weight_cutoff: float = 0.05,
                                      figsize=(16, 6),
                                      xticklabels=1,
                                      xlabelsize=8):
        active_ligand_target_links_df = active_ligand_target_links_df[active_ligand_target_links_df['weight'] > weight_cutoff]
        active_ligand_target_links_df_wide = active_ligand_target_links_df.pivot(index='ligand', columns='target', values='weight')
        # print(active_ligand_target_links_df_wide.shape)
        active_ligand_target_links_df_wide.fillna(0, inplace=True)
        clustermap = sns.clustermap(active_ligand_target_links_df_wide, figsize=figsize, xticklabels=xticklabels)
        if xlabelsize:
            clustermap.ax_heatmap.tick_params(axis='x', labelsize=8)
        return clustermap
        
    
    def get_lr_network_top_df(self,
                              best_upstream_ligands: list,
                              expressed_receptors: list):
        lr_network_df = self.rpy_2_pd_df(self.lr_network)
        lr_network_top = lr_network_df[lr_network_df['from'].isin(best_upstream_ligands) & 
                                       lr_network_df['to'].isin(expressed_receptors)].drop_duplicates(subset=['from', 'to'])
        best_upstream_receptors = lr_network_top['to'].unique().tolist()

        lr_network_top_df_large = self.rpy_2_pd_df(self.weighted_networks_lr_sig)
        lr_network_top_df_large = lr_network_top_df_large[lr_network_top_df_large['from'].isin(best_upstream_ligands) & 
                                                          lr_network_top_df_large['to'].isin(best_upstream_receptors)]
        return lr_network_top_df_large
    
    def plot_ligand_receptor_clustermap(self,
                                        lr_network_top_df_large: pd.DataFrame,
                                        weight_cutoff: float,
                                        figsize = (16, 6),
                                        xticklabels = 1,
                                        xlabelsize = 8):
        lr_network_top_df_large = lr_network_top_df_large[lr_network_top_df_large['weight'] > weight_cutoff]
        lr_network_top_df_wide = lr_network_top_df_large.pivot(index='from', columns='to', values='weight')
        lr_network_top_df_wide.fillna(0, inplace=True)
        cluster = sns.clustermap(lr_network_top_df_wide, figsize=figsize, xticklabels=xticklabels)
        cluster.ax_heatmap.tick_params(axis='x', labelsize=xlabelsize)
        _ = cluster.ax_heatmap.set_xlabel('Receptor')
        _ = cluster.ax_heatmap.set_ylabel('Ligand')
        return cluster
    
    def get_ligand_signaling_path(self,
                                  ligands_all: list,
                                  targets_all: list,
                                  weight_cutoff: float = None,
                                  need_data_sources: bool = False) -> pd.DataFrame:
        """Get the ligands to target pathways, including ligands, receptors, and targets, covering the whole signaling paths.

        Args:
            ligands_all (list): the ligands to be checked
            targets_all (list): the targets to be checked
            weight_cutoff (float): the weight cutoff to select paths
            need_data_sources: flag if the original data sources (e.g. harmonizome) is needed
        Returns:
            pd.DataFrame: two or three (when need_data_sources is true) pd.DataFrames to be returned, the first is sig_network and the second gr_network.
            The third is the data source pd.DataFrame if requested.
        """
        active_signaling_network = nichenet.get_ligand_signaling_path(ligand_tf_matrix=self.ligand_tf_matrix, 
                                                                      ligands_all=StrVector(ligands_all), 
                                                                      targets_all=StrVector(targets_all), 
                                                                      weighted_networks=self.weighted_networks)
        sig_network_pd = self.rpy_2_pd_df(active_signaling_network.rx2['sig'])
        gr_network_pd = self.rpy_2_pd_df(active_signaling_network.rx2['gr'])
        if weight_cutoff:
            sig_network_pd = sig_network_pd[sig_network_pd['weight'] > weight_cutoff]
            gr_network_pd = gr_network_pd[gr_network_pd['weight'] > weight_cutoff]
        if not need_data_sources:
            return sig_network_pd, gr_network_pd
        else:
            data_sources = nichenet.infer_supporting_datasources(signaling_graph_list = active_signaling_network,
                                                            lr_network=self.lr_network,
                                                            sig_network = self.sig_network,
                                                            gr_network = self.gr_network)
            data_sources_pd = self.rpy_2_pd_df(data_sources)
            return sig_network_pd, gr_network_pd, data_sources_pd
    
    def build_network_from_paths(self,
                                 sig_network_df: pd.DataFrame,
                                 gr_network_df: pd.DataFrame,
                                 network_data_sources_df: pd.DataFrame = None,
                                 graphml_file: str = None) -> nx.DiGraph:
        """Build a networkx DiGraph object from the two paths.

        Args:
            sig_network_df (pd.DataFrame): signaling network paths in pd.DataFrame.
            gr_network_df (pd.DataFrame): gene regulatory network paths in pd.DataFrame.
            graphml_file (str, optional): the file name for exporting the built network object. Defaults to None.

        Returns:
            nx.DiGraph: the constructed network DiGraph object.
        """
        edge_attributes = ['weight']
        if network_data_sources_df is not None:
            sig_network_df = pd.merge(sig_network_df, network_data_sources_df, how='inner', on=['from', 'to'])
            gr_network_df = pd.merge(gr_network_df, network_data_sources_df, how='inner', on=['from', 'to'])
            edge_attributes.extend(network_data_sources_df.columns[2:])
        sig_network = nx.from_pandas_edgelist(sig_network_df, 'from', 'to', edge_attributes, create_using=nx.DiGraph())
        gr_network = nx.from_pandas_edgelist(gr_network_df, 'from', 'to', edge_attributes, create_using=nx.DiGraph())
        network = nx.compose(sig_network, gr_network)
        
        if graphml_file:
            nx.write_graphml(network, graphml_file)

        return network

    def get_ligand_2_targets(self,
                             background_genes: list = None,    
                             n_top :int = 250) -> dict:
        """Get a dict from ligands to targets.

        Args:
            background_genes (list, optional): genes that should be considered,
            n_top (int, optional): ranked based on the ligand_target_matrix. Defaults to 250.

        Returns:
            dict: _description_
        """
        ligand2targets = {}

        # Convert it to pd.DataFrame for easy handling in Python
        ligand_target_df = self.convert_matrix_pd(self.ligand_target_matrix)
        if background_genes is not None:
            # Do a filter
            ligand_target_df = ligand_target_df.loc[ligand_target_df.index.isin(background_genes)]
        for ligand in ligand_target_df.columns:
            targets = ligand_target_df[ligand].sort_values(ascending=False)[:n_top]
            ligand2targets[ligand] = targets.index.to_list()

        return ligand2targets
    

    def get_ligand_2_tfs(self,
                         weight_cutoff: float=0.001) -> dict:
        ligand2tfs = {}
        ligand_tf_df = self.convert_matrix_pd(self.ligand_tf_matrix)
        for ligand in ligand_tf_df.columns:
            tfs = ligand_tf_df[ligand_tf_df[ligand] >= weight_cutoff].index.to_list()
            if len(tfs) == 0:
                continue
            ligand2tfs[ligand] = tfs
        return ligand2tfs

    
    def list_all_targets(self) -> list:
        ligand2targets = self.get_ligand_2_targets()
        all_targets = [target for targets in ligand2targets.values() for target in targets]
        all_targets = list(set(all_targets))
        return all_targets