from jsonrpclib.SimpleJSONRPCServer import SimpleJSONRPCServer
import logging as logger
from . import ScanpyWrapper as analyzer
import scvelo as scv
import scanpy as sc


def echo(text):
    return "You sent: " + text


def scv_open(file_name):
    adata = analyzer.scv_open(file_name)
    analyzer.cache_data(adata)
    return str(adata)


def scv_preprocess():
    adata = analyzer.get_loaded_data()
    analyzer.scv_preprocess(adata)
    # This is the same as the loaded data for scv
    analyzer.cache_processed_data(adata)
    return str(adata)


def scv_velocity(mode):
    adata = analyzer.get_processed_data()
    analyzer.scv_velocity(adata, mode=mode)
    return str(adata)


def scv_velocity_plot(gene: str):
    adata = analyzer.get_processed_data()
    if gene not in adata.var_names:
        return "error: " + gene + " cannot be found."
    file_name = gene + '_velocity.pdf'
    scv.pl.velocity(adata, gene, color='leiden', show=False, save=file_name)
    return "scvelo_" + file_name


def scv_rank_velocity_genes():
    adata = analyzer.get_processed_data()
    scv.tl.rank_velocity_genes(adata, groupby='leiden', n_genes=analyzer.n_rank_genes)
    return adata.uns['rank_velocity_genes']['names'].tolist()


def scv_rank_dynamic_genes():
    adata = analyzer.get_processed_data()
    # Have to make sure dynamic mode is used for RNA velocity analysis
    if adata.uns['velocity_params']['mode'] != 'dynamical':
        return "Error: The dynamical mode for RNA velocity analysis must be used to rank dynamic genes."
    scv.tl.rank_dynamical_genes(adata, groupby='leiden', n_genes=analyzer.n_rank_genes)
    return adata.uns['rank_dynamical_genes']['names'].tolist()


def scv_embedding(color_key=None):
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


def open_data(dir_name):
    adata = analyzer.open_10_genomics_data(dir_name)
    analyzer.cache_data(adata)
    # Just return a string for the client
    return str(adata)


def open_analyzed_data(file_name: str) -> str:
    """
    Open a processed adata writted by function write_data below.
    :param file_name:
    :return:
    """
    adata = sc.read(file_name)
    analyzer.cache_processed_data(adata)
    return str(adata)


def write_data(file_name: str) -> str:
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
            scv=False):
    adata = analyzer.get_processed_data()
    if adata is None:
        return "error: no pre-processed reference data is available."
    merged_data = analyzer.project(dir_name, adata, scv)
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
    print("imputation: ", imputation)
    print("regress_out_keys: ", regress_out_keys)
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
    # We should use pre-processed data for clustering analysis
    adata = analyzer.get_processed_data()
    if adata is None:
        return "error: no preprocessed data. Call open_data and preproces_data first."
    analyzer.cluster(adata, plot=False)  # Plot should be turn off
    # Expect to see more variables after clustering
    return str(adata)


def get_umap():
    adata = analyzer.get_processed_data()
    key = 'X_umap'
    if key not in adata.obsm.keys():
        return "error: no clustering data. Call open_data, preprocess_data, and cluster_data first."
    return adata.obsm[key].tolist()


def get_connectivites():
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
    adata = analyzer.get_processed_data()
    key = "cytotrace"
    if key not in adata.obs.keys():
        analyzer.cytotrace(adata)
    return adata.obs[key].tolist()


def get_cluster():
    return get_obs('leiden')


def get_obs(obs_name: str):
    """
    Get the values in the obs data frame for individual cells.
    :param obs_name:
    :return:
    """
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
    adata = analyzer.get_processed_data()
    if adata is None:
        return "error: no preprocessed data. Call open_data and preproces_data first."
    rtn = adata.obs.keys().to_list()
    # # cluster should be handled elsewhere
    # if 'leiden' in rtn :
    #     rtn.remove('leiden')
    return rtn


def get_cell_ids():
    # Have to use the processed data. Otherwise, cell ids may be too many
    adata = analyzer.get_processed_data()
    if adata is None:
        return "error: no preprocessed data. Call open_data and preproces_data first."
    return adata.obs.index.to_list()


def rank_genes_groups(groups='all',
                      reference='rest',
                      groupby='leiden') -> dict:
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
        if key1 is 'params':
            continue  # Don't need to expose this
        values = adata.uns[key][key1]
        values_converted = list()
        for value in values:
            values_converted.append(value.tolist())
        rtn[key1] = values_converted
    return rtn


def get_paga():
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
    adata = analyzer.get_processed_data()
    if adata is None:
        return "error: no preprocessed data. Call open_data, preprocess first."
    # Generate a list
    target_clusters = None
    if len(args) > 0:
        target_clusters = args
    return analyzer.infer_cell_root(adata, target_clusters)


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
    server.register_function(echo)
    server.register_function(stop)
    logger.info("Start server...")
    # server.serve_forever()
    start()


def start():
    while (isWaiting):
        server.handle_request()
        logger.info("isWaiting", isWaiting)
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


# Enable basic logging for the time being
logger.basicConfig()
# Define two global level variables so that we can control the server's behaviors
# Use 0 if we dont need to specify the port number
# TODO: Make sure this be configured or automatically assigned
port = 8999
server = SimpleJSONRPCServer(('localhost', port))
logger.info("Server initialized at ", server.server_address)
# Server address is a tupe. The first element is the host and the second is the port number
isWaiting = True

if __name__ == '__main__':
    main()
