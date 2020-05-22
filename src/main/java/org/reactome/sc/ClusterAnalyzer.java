package org.reactome.sc;

import smile.clustering.GMeans;
import smile.clustering.HierarchicalClustering;
import smile.clustering.KMeans;
import smile.clustering.XMeans;
import smile.clustering.linkage.CompleteLinkage;
import smile.clustering.linkage.Linkage;
import smile.data.DataFrame;
import smile.data.vector.BaseVector;
import smile.data.vector.IntVector;

/**
 * This class uses the smile clustering features. The default algorithm is X-means
 * with maximum of clusters 50.
 * @author wug
 *
 */
@SuppressWarnings("rawtypes")
public class ClusterAnalyzer {
    
    /**
     * For the 6 group-simulated data UMAP plot, the XMeans and GMeans have the best
     * match with the original groups assignment and the UMAP plot results with almost
     * similar results. KMeans and Hierarchical clustering messed up the final results.
     * This analysis is based on clusterNumber = 6 with other default values.
     * @author wug
     *
     */
    static enum Type {
        XMeans,
        KMeans,
        GMeans,
        HierarchicalClustering;
    }
    
    private final int DEFAULT_MAX_CLUSTERS = 20;
    
    public ClusterAnalyzer() {
    }
    
    /**
     * The passed data should be in the cell x gene format.
     * @param data
     * @return
     */
    public BaseVector cluster(double[][] data) {
        XMeans xmeans = XMeans.fit(data, DEFAULT_MAX_CLUSTERS);
        return IntVector.of("Cluster", xmeans.y);
    }
    
    /**
     * Some algorithms don't need a pre-defined cluster name. For these algorithms,
     * the passed cluster will be used as the maximum cluster if provided.
     * @param data
     * @param type
     * @param clusterNumber
     * @return
     */
    public BaseVector cluster(double[][] data, 
                              Type type,
                              Integer clusterNumber) {
        if (clusterNumber == null)
            clusterNumber = DEFAULT_MAX_CLUSTERS;
        int[] clusters = null;
        switch (type) {
            case GMeans :  
                clusters = GMeans.fit(data, clusterNumber).y;
            case XMeans :
                clusters = XMeans.fit(data, clusterNumber).y;
            case HierarchicalClustering : 
                // Use CompleteLinkage as the default
                Linkage linkage = CompleteLinkage.of(data);
                HierarchicalClustering clustering = HierarchicalClustering.fit(linkage);
                // Need to check the total clusters we can get
                clusters = clustering.partition(clusterNumber);
            case KMeans : default :
                clusters = KMeans.fit(data, clusterNumber).y;
        }
        return IntVector.of("Cluster", clusters);
    }
    
    /**
     * Cluster the passed data frame and add a new column at the end for the clustering results.
     * @param df
     * @return
     */
    public void clusterDF(ScDataFrame df) {
        double[][] data = SmileUtilities.getCounts(df);
        BaseVector clusterCol = cluster(data);
        DataFrame newWrapped = df.getWrapped().merge(clusterCol);
        df.setWrapped(newWrapped);
    }

}
