import pandas as pd
from anndata import AnnData
from jsonrpclib.SimpleJSONRPCServer import SimpleJSONRPCServer
import logging as logger
import sys
import inspect
import scvelo as scv
import scanpy as sc
from typing import Optional
from typing import Union

from . import scanpy_wrapper as analyzer
from . import gene_rel_eval as rel
from . import pathway_analyzer as pa


def echo(text):
    return "You sent: " + text


def analyze_tfs(dorothea_file_name: str,
                method: str = "ssgsea",
                need_rtn: str = "false") -> Optional[str]:
    """
    Perform TF activities analysis using a Dorothea file for interactions.
    :param dorothea_file_name:
    :param method:
    :param need_rtn:
    :return:
    """
    logger.info("Perform transcriptional factor analysis using {} via method {}".format(dorothea_file_name, method))
    tf_genes = pa.load_dorothea_data(dorothea_file_name)
    data_key = None
    if method == "ssgsea":
        data_key = pa.TF_SSGSEA_KEY
    elif method == "aucell":
        data_key = pa.TF_AUCELL_KEY
    return _analyze_pathways(tf_genes, method, data_key, need_rtn)


def analyze_pathways(gmt_file_name: str,
                     method: str = "ssgsea",
                     need_rtn: str = "false") -> Optional[str]:
    """
    conduct pathway analysis.
    :param gmt_file_name: the reactome gmt file.
    :param method: one of ssgsea and aucell, which are supported currently
    :param need_rtn: true for returning the analysis results as a json text
    otherwiser, return the key for saved analysis results in the object. Need to use str here
    since json api cannot automatically convert it into a Boolean even though it is typed as it.
    :return: json text of the calculated pathway activities
    """
    logger.info("Perform pathway analysis using {} via method {}".format(gmt_file_name, method))
    return _analyze_pathways(gmt_file_name, method, None, need_rtn)


def _analyze_pathways(gmt_file_name: Union[str, dict],
                      method: str,
                      data_key: str,
                      need_rtn: str):
    if method != "ssgsea" and method != "aucell":
        return "error: Method {} is not supported. Only ssgsea and aucell are supported now.".format(method)
    adata = analyzer.get_processed_data()
    if adata is None:
        return "error: no processed data is available."
    # Is raw is available, use raw so that we have all genes
    # To copy the results back to the main adata
    adata.used = adata
    if adata.raw is not None:
        # Need to find the actual adata
        if isinstance(adata.raw, AnnData):
            adata.used = adata.raw
        else:
            adata.used = adata.raw.to_adata()  # Have to call this and cannot use _adata, which is the same as the adat
    if method == "ssgsea":
        if data_key is None:
            results = pa.reactome_ssgsea(adata.used, gmt_file_name)
        else:
            results = pa.reactome_ssgsea(adata.used, gmt_file_name, data_key)
    elif method == "aucell":
        if data_key is None:
            results = pa.reactome_aucell(adata.used, gmt_file_name)
        else:
            results = pa.reactome_aucell(adata.used, gmt_file_name, data_key)
    if results is None:
        return "error: Cannot perform a pathway analysis."
    # The pathway scores are in a panda DataFrames, which can be converted into json.
    if id(adata) != id(adata.used):  # Keep the analysis results in the original adata is it is not the same
        adata.obsm[results.adata_key] = results.adata.obsm[results.adata_key]
    if need_rtn == "true":
        return results.adata.obsm[results.adata_key].to_dict()  # Use dictionay for easy handling at the Java end
    else:
        return results.adata_key


def anova_pathway(adata_key: str) -> str:
    """
    Perform pathway anova analysis for clusters.
    :param adata_key: the key used to save the pathway analysis results.
    :return:
    """
    # Make sure adata_key is there
    adata = analyzer.get_processed_data()
    if adata_key not in adata.obsm.keys():
        return "error: {} not in the obsm keys. Run pathway analysis first with the required method.".format(adata_key)
    # Use a PathwayAnalyer object to do this
    annova_analyzer = pa.PathwayAnalyzer()
    annova_analyzer.adata = adata
    annova_analyzer.adata_key = adata_key
    results = annova_analyzer.perform_pathway_anova()
    # results is a DataFrame. Though DataFrame has to_json() function, however, the Java end prefers to use
    # dict. Therefore, we have this converting.
    return results.to_dict(orient="index") # Keyed by pathways


def pathway_activities(pathway: str,
                       adata_key: str) -> str:
    """
    Fetch the pathway activities.
    :param pathway:
    :param adata_key:
    :return:
    """
    # Make sure adata_key is there
    adata = analyzer.get_processed_data()
    if adata_key not in adata.obsm.keys():
        return "error: {} not in the obsm keys. Run pathway analysis first with the required method.".format(adata_key)
    # This should return a DataFrame
    results = adata.obsm[adata_key]
    # Convert from Series to a dict for Java
    return results[pathway].to_dict()


def scv_open(file_name):
    logger.info('scv_open(file_name = {})...'.format(file_name))
    adata = analyzer.scv_open(file_name)
    analyzer.cache_data(adata)
    return str(adata)


def scv_preprocess():
    logger.info('{}...'.format(inspect.currentframe().f_code.co_name))
    adata = analyzer.get_loaded_data()
    analyzer.scv_preprocess(adata)
    # This is the same as the loaded data for scv
    analyzer.cache_processed_data(adata)
    return str(adata)


def scv_velocity(mode):
    logger.info('scv_velocity(mode = {})...'.format(mode))
    adata = analyzer.get_processed_data()
    analyzer.scv_velocity(adata, mode=mode)
    return str(adata)


def scv_velocity_plot(gene: str):
    logger.info('scv_velocity_plot(gene = {})...'.format(gene))
    adata = analyzer.get_processed_data()
    if gene not in adata.var_names:
        return "error: " + gene + " cannot be found."
    file_name = gene + '_velocity.pdf'
    scv.pl.velocity(adata, gene, color='leiden', show=False, save=file_name)
    return "scvelo_" + file_name


def scv_rank_velocity_genes():
    logger.info('{}...'.format(inspect.currentframe().f_code.co_name))
    adata = analyzer.get_processed_data()
    scv.tl.rank_velocity_genes(adata, groupby='leiden', n_genes=analyzer.n_rank_genes)
    return adata.uns['rank_velocity_genes']['names'].tolist()


def scv_rank_dynamic_genes():
    logger.info('{}...'.format(inspect.currentframe().f_code.co_name))
    adata = analyzer.get_processed_data()
    # Have to make sure dynamic mode is used for RNA velocity analysis
    if adata.uns['velocity_params']['mode'] != 'dynamical':
        return "Error: The dynamical mode for RNA velocity analysis must be used to rank dynamic genes."
    scv.tl.rank_dynamical_genes(adata, groupby='leiden', n_genes=analyzer.n_rank_genes)
    return adata.uns['rank_dynamical_genes']['names'].tolist()


def scv_embedding(color_key=None):
    logger.info('{}...'.format(inspect.currentframe().f_code.co_name))
    adata = analyzer.get_processed_data()
    if color_key is None:
        color_key = 'leiden'
    elif color_key not in adata.var_names:
        return "error: " + color_key + " cannot be found."
    file_name = color_key + '_umap_embedding.pdf'
    # Just dump the plot to a file and let the client do whatever it needs
    # cannot generate a non-blocking interactive plot here.
    # TODO: Study how to use an async call for the following statement
    scv.pl.velocity_embedding(adata, basis="umap", color=color_key, show=False, save=file_name)
    # For some unknown reason, the actual file name having scvelo prefixed
    return "scvelo_" + file_name


def scv_embedding_grid(color_key=None):
    logger.info('{}...'.format(inspect.currentframe().f_code.co_name))
    adata = analyzer.get_processed_data()
    if color_key is None:
        color_key = 'leiden'
    elif color_key not in adata.var_names:
        return "error: " + color_key + " cannot be found."
    file_name = color_key + '_umap_embedding_grid.pdf'
    # Just dump the plot to a file and let the client do whatever it needs
    # cannot generate a non-blocking interactive plot here.
    # TODO: Study how to use an async call for the following statement
    scv.pl.velocity_embedding_grid(adata, basis="umap", color=color_key, show=False, save=file_name)
    # For some unknown reason, the actual file name having scvelo prefixed
    return "scvelo_" + file_name


def scv_embedding_stream(color_key=None):
    logger.info('{}...'.format(inspect.currentframe().f_code.co_name))
    adata = analyzer.get_processed_data()
    if color_key is None:
        color_key = 'leiden'
    elif color_key not in adata.var_names:
        return "error: " + color_key + " cannot be found."
    file_name = color_key + '_umap_embedding_stream.png'
    # Just dump the plot to a file and let the client do whatever it needs
    # cannot generate a non-blocking interactive plot here.
    # TODO: Study how to use an async call for the following statement
    scv.pl.velocity_embedding_stream(adata, basis="umap", color=color_key, show=False, save=file_name)
    # For some unknown reason, the actual file name having scvelo prefixed
    return "scvelo_" + file_name


def open_data(dir_name,
              method: str):
    logger.info('open_data(dir_name = {}, method = {})...'.format(dir_name, method))
    adata = analyzer.open_10_genomics_data(dir_name, method)
    analyzer.cache_data(adata)
    # Just return a string for the client
    return str(adata)


def open_analyzed_data(file_name: str) -> str:
    logger.info('open_analyzed_data(file_name = {})...'.format(file_name))
    """
    Open a processed adata writted by function write_data below.
    :param file_name:
    :return:
    """
    adata = sc.read(file_name)
    analyzer.cache_processed_data(adata)
    return str(adata)


def write_data(file_name: str) -> str:
    logger.info("write_data(file_name = {})...".format(file_name))
    """
    Write the loaded data into a file in the h5ad format
    :param file_name:
    :return:
    """
    adata = analyzer.get_processed_data()
    if adata is None:
        return "error: no data loaded for writing."
    adata.write(file_name, compression='gzip')
    return str(adata)  # for debug purpose


def project(dir_name,
            method: str,
            scv: bool = False):
    """
    Project a new data onto the existing one.
    :param dir_name: the new project directory or file name
    :param method: the file format. This is the same as in open_data.
    :param scv: A Bool to check if the RNA veclocity is used.
    :return:
    """
    logger.info('project(dir_name = {}, method = {}, scv = {})...'.format(dir_name, method, scv))
    adata = analyzer.get_processed_data()
    if adata is None:
        return "error: no pre-processed reference data is available."
    if scv is not None and scv == 'true':
        scv = True
    else:
        scv = False
    merged_data = analyzer.project(dir_name, adata, method, scv)
    analyzer.cache_merged_data(merged_data)
    # Return the location of UMAP coordinates for new_data.
    merged_new_data = merged_data[merged_data.obs['batch'] == 'new']
    zipped = zip(merged_new_data.obs.index.to_list(), merged_new_data.obsm['X_umap'].tolist(),
                 merged_new_data.obs['leiden'].tolist())
    rtn = dict()
    for cell, umap, leiden in zipped:
        rtn[cell] = (umap[0], umap[1], leiden)
    return rtn


def preprocess_data(regress_out_keys=None,
                    imputation: str = 'magic'):
    """
    Run preprocess steps.
    :param regress_out_keys:
    :param imputation: if it is not null, use 'magic' for preprocess. Currently no other is support
    :return:
    """
    # Convert the gress_out_keys into a list
    logger.info('preprocess_data(regress_out_keys = {}, imputation = {})...'.format(regress_out_keys, imputation))
    if regress_out_keys is not None:
        if len(regress_out_keys) == 0:
            regress_out_keys = None
        else:
            regress_out_keys = str.split(regress_out_keys, ",")
    if imputation is not None:
        if len(imputation) == 0:
            imputation = None
        elif imputation != 'magic':
            return "error: The supported imputation method is 'magic' only!"
    adata = analyzer.get_loaded_data()
    if adata is None:
        return "error: no data loaded. Call open_data first."
    processed = analyzer.preprocess(adata,
                                    copy=True,
                                    need_scale=True,
                                    regressout_keys=regress_out_keys,
                                    imputation=imputation)
    analyzer.cache_processed_data(processed)
    return str(processed)


def cluster_data():
    logger.info('{}...'.format(inspect.currentframe().f_code.co_name))
    # We should use pre-processed data for clustering analysis
    adata = analyzer.get_processed_data()
    if adata is None:
        return "error: no preprocessed data. Call open_data and preproces_data first."
    analyzer.cluster(adata, plot=False)  # Plot should be turn off
    # Expect to see more variables after clustering
    return str(adata)


def get_umap():
    logger.info('{}...'.format(inspect.currentframe().f_code.co_name))
    adata = analyzer.get_processed_data()
    key = 'X_umap'
    if key not in adata.obsm.keys():
        return "error: no clustering data. Call open_data, preprocess_data, and cluster_data first."
    return adata.obsm[key].tolist()


def get_connectivites():
    logger.info('{}...'.format(inspect.currentframe().f_code.co_name))
    adata = analyzer.get_processed_data()
    key = 'connectivities'
    if key not in adata.obsp.keys():
        return "error: no connectivities data. Call open_data, preprocess_data, and cluster_data first."
    # Use the network structure for output the connectivities
    import networkx as nx
    network = nx.Graph(adata.obsp[key])
    rtn = list()
    for edge in network.edges:
        rtn.append((str(edge[0]), str(edge[1]), str(network[edge[0]][edge[1]]['weight'])))
    return rtn


def cytotrace():
    logger.info('{}...'.format(inspect.currentframe().f_code.co_name))
    adata = analyzer.get_processed_data()
    key = "cytotrace"
    if key not in adata.obs.keys():
        analyzer.cytotrace(adata)
    return adata.obs[key].tolist()


def get_cluster():
    logger.info('{}...'.format(inspect.currentframe().f_code.co_name))
    return get_obs('leiden')


def get_obs(obs_name: str):
    """
    Get the values in the obs data frame for individual cells.
    :param obs_name:
    :return:
    """
    logger.info('get_obs(obs_name = {})...'.format(obs_name))
    adata = analyzer.get_processed_data()
    if adata is None:
        return "error: no preprocessed data. Call open_data and preproces_data first."
    if obs_name not in adata.obs.keys():
        return "error: " + obs_name + " is not in the preprocessed data."
    return adata.obs[obs_name].tolist()


def get_obs_names():
    """"
    Get a list of obs_names.
    """
    logger.info('{}...'.format(inspect.currentframe().f_code.co_name))
    adata = analyzer.get_processed_data()
    if adata is None:
        return "error: no preprocessed data. Call open_data and preproces_data first."
    rtn = adata.obs.keys().to_list()
    # # cluster should be handled elsewhere
    # if 'leiden' in rtn :
    #     rtn.remove('leiden')
    return rtn


def get_cell_ids():
    logger.info('{}...'.format(inspect.currentframe().f_code.co_name))
    # Have to use the processed data. Otherwise, cell ids may be too many
    adata = analyzer.get_processed_data()
    if adata is None:
        return "error: no preprocessed data. Call open_data and preproces_data first."
    return adata.obs.index.to_list()


def rank_genes_groups(groups='all',
                      reference='rest',
                      groupby='leiden') -> dict:
    logger.info('rank_genes_groups(groups = {}, reference = {}, groupby = {})...'.format(groups, reference, groupby))
    adata = analyzer.get_processed_data()
    if adata is None:
        return "error: no preprocessed data. Call open_data and preproces_data first."
    analyzer.rank_genes_groups(adata,
                               groups=groups,
                               reference=reference,
                               groupby=groupby)
    key = 'rank_genes_groups'
    if key not in adata.uns.keys():
        return "error: rank_genes_groups() cannot finish."
    # Generate a disc for return
    rtn = dict()
    for key1 in adata.uns[key]:
        if key1 == 'params':
            continue  # Don't need to expose this
        values = adata.uns[key][key1]
        values_converted = list()
        for value in values:
            values_converted.append(value.tolist())
        rtn[key1] = values_converted
    return rtn


def get_paga():
    logger.info('{}...'.format(inspect.currentframe().f_code.co_name))
    adata = analyzer.get_processed_data()
    key = 'paga'
    if key not in adata.uns.keys():
        return "error: no clustering data. Call open_data, preprocess_data, and cluster_data first."
    # Need some process for json converting
    rtn = dict()
    # A list of list of double
    rtn['pos'] = adata.uns[key]['pos'].tolist()
    # Since this is a graph for clusters and the adjacency matrix is not that sparse,
    # using this should be fine. This should be a list of list of double for a n x n
    # matrix (n is the number of clusters)
    edge_key = 'transitions_confidence'  # Directed cluster adjacency matrix from velocity analysis
    if edge_key not in adata.uns['paga'].keys():
        edge_key = 'connectivities'  # undirected cluster adjacency matrix: symmetric
    rtn['connectivities'] = adata.uns[key][edge_key].toarray().tolist()
    return rtn


def dpt(root_cell: str):
    logger.info('dpt(root_cell = {})...'.format(root_cell))
    adata = analyzer.get_processed_data()
    if adata is None:
        return "error: no preprocessed data. Call open_data, preprocess first."
    return analyzer.dpt(adata, root_cell).to_list()


def get_gene_exp(gene: str):
    """
    Get the gene expression value for the passed gene. If there is a raw, use the raw value. Otherwise, use the
    processed data.
    :param gene:
    :return:
    """
    logger.info('get_gene_exp(gene = {})...'.format(gene))
    adata = analyzer.get_processed_data()
    if adata is None:
        adata = analyzer.get_loaded_data()
    if adata is None:
        return "error: no data is loaded. Call open_data first."
    # Check if the query gene is in the var list
    # We will prefer to use the raw if it is there
    var_names = None
    if adata.raw is not None:
        var_names = adata.raw.var_names
    else:
        var_names = adata.var_names
    if gene not in var_names:
        return "error: " + gene + " doesn't have any expression data."
    rtn = None
    if adata.raw is not None:
        rtn = adata.raw.obs_vector(gene)
    else:
        rtn = adata.obs_vector(gene)
    if rtn is None:
        return "error: cannot find expression values for " + gene + "."
    return rtn.tolist()


def infer_cell_root(*args):
    logger.info('infer_cell_root({})...'.format(args))
    adata = analyzer.get_processed_data()
    if adata is None:
        return "error: no preprocessed data. Call open_data, preprocess first."
    # Generate a list
    target_clusters = None
    if len(args) > 0:
        target_clusters = args
    return analyzer.infer_cell_root(adata, target_clusters)


def get_cell_time_keys() -> str:
    logger.info('{}...'.format(inspect.currentframe().f_code.co_name))
    adata = analyzer.get_processed_data()
    # Get the keys for gene relationships calcuations
    rtn = ['latent_time', 'velocity_pseudotime', 'cytotrace', 'dpt_pseudotime']
    return [i for i in rtn if i in adata.obs_keys()]


def calculate_gene_relations(gene_pairs: str,
                             groups: str,
                             cell_time_key: str,
                             layer: str = None,
                             delay_window=None,
                             mode='spearman') -> dict:
    """
    Calculate gene-gene relations for a list of passed gene pairs.
    :param gene_pairs: list of gene pairs with two genes tab-delimited. The two genes are directed.
    :param groups: calculations should be conducted based on cells in the specific groups.
    :param cell_time_key: one of latent_time, velocity_pseudotime, cytotrace, and dpt_pseudotime.
    The first two values are generated for RNA-velocity data analysis.
    :param layer:
    :param delay_window:
    :param mode:
    :return:
    """
    logger.info('calculate_gene_regulations(groups = {}, cell_time_key = {}, '
                'layer = {}, delay_window = {}, '
                'mode = {})...'.format(groups, cell_time_key, layer, delay_window, mode))
    adata = analyzer.get_processed_data()
    if cell_time_key not in adata.obs_keys():
        return "error: " + cell_time_key + " is not in the observation keys."
    # A specific key for caller
    if layer is not None and layer == 'null':
        layer = None
    if layer is not None and layer not in adata.layers:
        return "error: " + layer + " is not in the dataset."
    adata_slice = None
    if groups is None or groups == 'all':
        adata_slice = adata
    else:
        groups = [i for i in groups.split(',')]
        selected = adata.obs_vector('leiden').isin(groups)
        adata_slice = adata[selected, :]
    if delay_window is not None:
        delay_window = int(delay_window)  # Force the string to an int
    # Want to focus on the passed clusters
    gene_pairs = [i for i in gene_pairs.split('\n')]
    return rel.calculate_gene_relations(gene_pairs,
                                        adata_slice,
                                        cell_time_key,
                                        layer,
                                        delay_window,
                                        mode)


def main():
    # server = SimpleJSONRPCServer(('localhost', 8085))
    server.register_function(open_data)
    server.register_function(write_data)
    server.register_function(open_analyzed_data)
    server.register_function(preprocess_data)
    server.register_function(cluster_data)
    server.register_function(get_umap)
    server.register_function(get_cluster)
    server.register_function(get_cell_ids)
    server.register_function(get_connectivites)
    server.register_function(get_paga)
    server.register_function(get_gene_exp)
    server.register_function(get_obs)
    server.register_function(get_obs_names)
    server.register_function(rank_genes_groups)
    server.register_function(cytotrace)
    server.register_function(project)
    server.register_function(dpt)
    server.register_function(infer_cell_root)
    server.register_function(scv_open)
    server.register_function(scv_preprocess)
    server.register_function(scv_velocity)
    server.register_function(scv_embedding)
    server.register_function(scv_embedding_grid)
    server.register_function(scv_embedding_stream)
    server.register_function(scv_velocity_plot)
    server.register_function(scv_rank_velocity_genes)
    server.register_function(scv_rank_dynamic_genes)
    server.register_function(calculate_gene_relations)
    server.register_function(get_cell_time_keys)
    server.register_function(analyze_pathways)
    server.register_function(analyze_tfs)
    server.register_function(anova_pathway)
    server.register_function(pathway_activities)
    server.register_function(echo)
    server.register_function(stop)
    logger.info("Start server...")
    # server.serve_forever()
    start()


def start():
    while (isWaiting):
        server.handle_request()
        # logger.debug("isWaiting: {}".format(isWaiting))
    server.server_close()  # Don't use shutdown(). It will block the call.
    logger.info("Server stopped.")


"""
To call the following method in Java, we need a complicated JSON object like the following:
    @Test
    public void testPythonServer() throws IOException {
        String url = HOST_URL; // http://localhost:8070
        String query = "{\"jsonrpc\": \"2.0\", \"method\": \"stop\", \"id\": 2}";
        String output = callHttp(url, HTTP_POST, query); // POST should be used always since the query is a JSON object.
        System.out.println(output);
    }
"""


def stop():
    # Make sure this change can be popped out in this scope
    global isWaiting
    isWaiting = False


# Set up logging
logging_file_name = None
if len(sys.argv) > 2:
    logging_file_name = sys.argv[2]
# Enable basic logging for the time being
logger.basicConfig(format='%(asctime)s %(message)s',
                   datefmt='%m/%d/%Y %I:%M:%S %p',
                   filename=logging_file_name,
                   level=logger.INFO)

logger.info("Scanpy information...")
# Info and hint
sc.settings.verbosity = 3
sc.settings.logfile = logging_file_name
sc.logging.print_versions()
sc.logging.print_version_and_date()
logger.info('Random state, {}, is used for all functions needed a randomState setting.'.format(analyzer.random_state))
logger.info('Default top_rank_genes is {}.'.format(analyzer.n_rank_genes))

# Define two global level variables so that we can control the server's behaviors
# Use 0 if we dont need to specify the port number
# Get get the port from sys params
port = 8999
if len(sys.argv) > 1:
    port = int(sys.argv[1])  # cast to int
logger.info("Port: {}".format(port))
server = SimpleJSONRPCServer(('localhost', port))
logger.info("Server initialized at {}".format(server.server_address))
# Server address is a tupe. The first element is the host and the second is the port number
isWaiting = True
# Start the server
# main()