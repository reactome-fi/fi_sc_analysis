from jsonrpclib.SimpleJSONRPCServer import SimpleJSONRPCServer
import logging as logger
import ScanpyWrapper as analyzer

def echo(text) :
    return "You sent: " + text

def open_data(dir_name) :
    adata = analyzer.open_10_genomics_data(dir_name)
    analyzer.cache_data(adata)
    # Just return a string for the client
    return str(adata)

def preprocess_data() :
    adata = analyzer.get_loaded_data()
    if adata is None :
        return "error: no data loaded. Call open_data first."
    processed = analyzer.preprocess(adata, copy=True, need_scale=True)
    analyzer.cache_processed_data(processed)
    return str(processed)

def cluster_data() :
    # We should use pre-processed data for clustering analysis
    adata = analyzer.get_processed_data()
    if adata is None :
        return "error: no preprocessed data. Call open_data and preproces_data first."
    analyzer.cluster(adata, plot=False) # Plot should be turn off
    # Expect to see more variables after clustering
    return str(adata)

def get_umap() :
    adata = analyzer.get_processed_data()
    key = 'X_umap'
    if key not in adata.obsm.keys() :
        return "error: no clustering data. Call open_data, preprocess_data, and cluster_data first."
    return adata.obsm[key].tolist()

def get_connectivites() :
    adata = analyzer.get_processed_data()
    key = 'connectivities'
    if key not in adata.obsp.keys() :
        return "error: no connectivities data. Call open_data, preprocess_data, and cluster_data first."
    # Use the network structure for output the connectivities
    import networkx as nx
    network = nx.Graph(adata.obsp[key])
    rtn = list()
    for edge in network.edges :
        rtn.append((str(edge[0]), str(edge[1]), str(network[edge[0]][edge[1]]['weight'])))
    return rtn

def get_cluster() :
    adata = analyzer.get_processed_data()
    key = 'leiden'
    if key not in adata.obs.keys() :
        return "error: no clustering data. Call open_data, preprocess_data, and cluster_data first."
    # Apparently we have to convert pandas.core.series.Series in to a list of string.
    # Otherwise, only class type is returned.
    return adata.obs[key].tolist()

def get_cell_ids() :
    # Have to use the processed data. Otherwise, cell ids may be too many
    adata = analyzer.get_processed_data();
    if adata is None :
        return "error: no preprocessed data. Call open_data and preproces_data first."
    return adata.obs.index.to_list()

def get_paga() :
    adata = analyzer.get_processed_data();
    key = 'paga'
    if key not in adata.uns.keys() :
        return "error: no clustering data. Call open_data, preprocess_data, and cluster_data first."
    # Need some process for json converting
    rtn = dict()
    # A list of list of double
    rtn['pos'] = adata.uns[key]['pos'].tolist()
    # Since this is a graph for clusters and the adjency matrix is not that spase,
    # using this should be fine. This should be a list of list of double for a n x n
    # matrix (n is the number of clusters)
    rtn['connectivities'] = adata.uns[key]['connectivities'].toarray().tolist()
    return rtn

def main() :
    # server = SimpleJSONRPCServer(('localhost', 8085))
    server.register_function(open_data)
    server.register_function(preprocess_data)
    server.register_function(cluster_data)
    server.register_function(get_umap)
    server.register_function(get_cluster)
    server.register_function(get_cell_ids)
    server.register_function(get_connectivites)
    server.register_function(get_paga)
    server.register_function(echo)
    server.register_function(stop)
    logger.info("Start server...")
    # server.serve_forever()
    start()

def start() :
    while (isWaiting) :
        server.handle_request()
        logger.info("isWaiting", isWaiting)
    server.server_close() # Don't use shutdown(). It will block the call.
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

def stop() :
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

if __name__=='__main__' :
    main()
