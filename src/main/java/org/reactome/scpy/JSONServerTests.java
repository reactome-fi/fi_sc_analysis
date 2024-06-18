package org.reactome.scpy;

import org.junit.Test;
import org.reactome.cytoscape.sc.server.JSONServerCaller;

/**
 * Some methods included here to test scpy2reactome outside of ReactomeFIViz.
 * @author wug
 *
 */
public class JSONServerTests {
    private JSONServerCaller caller;
    private boolean isDataLoaded = false;
    
    public JSONServerTests() {
        caller = new JSONServerCaller();
        caller.setIsStarted(true); // We want to control the server here.
    }
    
    private void loadData() throws Exception {
        if (isDataLoaded)
            return;
        String dir_17_5 = "/Users/wug/git/reactome-fi/fi_sc_analysis/cache/Users-wug-Documents-missy_single_cell-seq_data_v2-17_5_gfp-filtered_feature_bc_matrix-matrix.h5ad";
        String text = caller.openData(dir_17_5, "read_h5ad");
        
        System.out.println("Open data: " + text);
        text = caller.preprocessData(null, null);
        System.out.println("Preprocess data: " + text);
    }
    
    /**
     * Before run this method, make sure run testVegaTrain() first.
     * @throws Exception
     */
    @Test
    public void testGetPathwayScore() throws Exception {
        String pathway = "GABA synthesis";
        Object rtn = caller.callJSONServer("vega_pathway_score", pathway);
        System.out.println("Score for " + pathway + ":\n" + rtn);
    }
    
    /**
     * Before run this method, make sure testVegaTrain has been called first.
     * @throws Exception
     */
    @Test
    public void testVegaPathwayAnova() throws Exception {
        System.out.println("Clustering cells...");
        Object rtn = caller.clusterData();
        System.out.println("Done clsutering:\n" + rtn);
        System.out.println("Perform vega pathway anova...");
        rtn = caller.callJSONServer("vega_pathway_anova");
        System.out.println("Done with the following output:\n" + rtn);
    }
    
    @Test
    public void testVegaTrain() throws Exception {
        System.out.println("Loading data...");
        loadData();
        System.out.println("Test vega training...");
        String gmtFileName = "python/data/vega/MouseReactomePathways_Rel_75_122220.gmt";
        Object rtn = caller.callJSONServer("train_vega", 
                                           gmtFileName, 
                                           2 + "");
        System.out.println("Finished training with the following retruend results:\n" + rtn);
    }

}
