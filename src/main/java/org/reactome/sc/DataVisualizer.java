package org.reactome.sc;

import java.awt.Color;
import java.lang.reflect.Field;
import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import javax.swing.JFrame;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import smile.data.DataFrame;
import smile.data.vector.BaseVector;
import smile.data.vector.DoubleVector;
import smile.data.vector.StringVector;
import smile.manifold.TSNE;
import smile.manifold.UMAP;
import smile.math.MathEx;
import smile.math.matrix.Matrix;
import smile.plot.swing.Canvas;
import smile.plot.swing.Legend;
import smile.plot.swing.ScatterPlot;
import smile.projection.PCA;

/**
 * Visualize the data.
 * @author wug
 *
 */
@SuppressWarnings("rawtypes")
public class DataVisualizer {
    private static final Logger logger = LoggerFactory.getLogger(DataVisualizer.class);

    public DataVisualizer() {
    }
    
    public static void main(String[] args) throws Exception {
        
        logger.debug("This is a test!");
        if (true)
            return;
        
//        String dir = "/Users/wug/Documents/missy_single_cell/seq_data_v2/17_5_gfp/filtered_feature_bc_matrix/";
//        String fileName = "matrix.mtx";
//        String dir = "data/simulated/";
//        String cellFileName = dir + "6_cell_info.csv";
//        String countFileName = dir + "6_group_true.csv";
//        countFileName = dir + "6_group_drop.csv";
//        countFileName = dir + "6_dca_results/mean.tsv";
//        
//        DataReader reader = new DataReader();
//        ScDataFrame df = reader.readCountsFiles(countFileName, 
//                                                CountMatrixLayout.GENE_TIMES_CELL,
//                                                cellFileName);
        
        String dir = "/Users/wug/Documents/missy_single_cell/seq_data_v2/17_5_gfp/filtered_feature_bc_matrix/";
        String countFileName = dir + "matrix.mtx";
        String barcodeFile = dir + "barcodes.tsv";
        String featureFile = dir + "features.tsv";
        String group = null;
        Integer top = 1000;
        boolean needPCA = true;
        Integer topPCs = 100;
        
        DataReader reader = new DataReader();
        ScDataFrame df = reader.readMatrixMarketFile(countFileName,
                                                     barcodeFile, 
                                                     featureFile,
                                                     1);
        
        DataVisualizer visualize = new DataVisualizer();
//        JFrame frame = visualize.visualizeInPCA(df, group, top);
        JFrame frame = visualize.visualizeInUMAP(df, 
                                                 group,
                                                 top, 
                                                 needPCA,
                                                 topPCs);
//        JFrame frame = visualize.visualizeInTSNE(df, 
//                                                 group, 
//                                                 top,
//                                                 needPCA,
//                                                 topPCs);
        frame.setTitle(countFileName);
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    }
    
    /**
     * The layout of the passed ScDataFrame should be cell x gene with cell meta.
     * The client can also specify if a PCA and how many PCs are needed first before tSNE.
     * To use all PCs, specify null for topPCs.
     * @param df
     * @param colorColName
     * @return
     * @throws Exception
     */
    public JFrame visualizeInTSNE(ScDataFrame df, // layout cell x gene
                                  String colorColName,
                                  Integer topVarFeatures,
                                  boolean needPCA,
                                  Integer topPCs) throws Exception {
        // Keep the color column first
        BaseVector colorCol = getColorCol(df, colorColName);
        Matrix normalizedMatrix = normalizeCounts(df);
        // Do tSNE
        double[][] values = getValues(normalizedMatrix, topVarFeatures);
        TSNE tsne = null;
        if (needPCA) {
            Matrix projection = prePCA(values, topPCs);
            tsne = new TSNE(convertMatrixToDoubleArray(projection), 2);
        }
        else
            tsne = new TSNE(values, 2);
        // Transpose so that we can use the ScatterPlot
        double[][] coordinates = MathEx.transpose(tsne.coordinates);
        return plot(coordinates,
                    "tSNE1",
                    "tSNE2",
                    colorCol);
    }

    /**
     * Get the values from the matrix as a two dimentional double array.
     * @param normalizedMatrix
     * @param topVarFeatures
     * @return
     */
    private double[][] getValues(Matrix normalizedMatrix, // cell x gene
                                 Integer topVarFeatures) {
        double[][] values = null;
        if (topVarFeatures == null)
            values = convertMatrixToDoubleArray(normalizedMatrix);
        else
            values = selectTopVariantGenes(normalizedMatrix, topVarFeatures);
        return values; // cell x gene
    }

    /**
     * Did a PCA for tSNE or UMAP.
     * @param values
     * @param topPCs
     * @return
     */
    private Matrix prePCA(double[][] values, Integer topPCs) {
        values = MathEx.transpose(values); // transpose from cell x gene to gene x cell
        // This layout of values is contradictory to the API doc:
        // The features (aka genes here) should be in rows and the samples should be columns
        PCA pca = PCA.fit(values).setProjection(topPCs == null ? values[1].length : topPCs);
        Matrix projection = pca.getProjection().transpose();
        return projection;
    }

    private JFrame plot(double[][] coordinates,
                        String xLabel,
                        String yLabel,
                        BaseVector colorCol) throws InterruptedException, InvocationTargetException {
        BaseVector[] cols = new BaseVector[2];
        // The order seems different from the conventional one
        cols[0] = DoubleVector.of(xLabel,
                                   coordinates[0]);
        cols[1] = DoubleVector.of(yLabel,
                                  coordinates[1]);
        DataFrame dataDF = DataFrame.of(cols);
        if (colorCol != null)
            dataDF = dataDF.merge(colorCol);
        ScatterPlot plot = null;
        if (colorCol == null) 
            plot = ScatterPlot.of(dataDF,
                                  xLabel, 
                                  yLabel, 
                                  '.',
                                  Color.black);
        else {
            plot = ScatterPlot.of(dataDF,
                                  xLabel, 
                                  yLabel, 
                                  colorCol == null ? null : colorCol.name(),
                    '.');
            // Sort the legends
            Legend[] legends = plot.legends().get();
            if (legends != null)
                Arrays.sort(legends, (l1, l2) -> getLegendText(l1).compareTo(getLegendText(l2)));
        }
        Canvas canvas = plot.canvas();
        canvas.setAxisLabels(xLabel, yLabel);
        return canvas.window();
    }

    public Matrix normalizeCounts(ScDataFrame df) {
        // Normalize the data
        DataNormalizer normalizer = new DataNormalizer();
        Matrix normalizedMatrix = normalizer.normalizePerCell(df.getWrapped(),
                                                              df.getCountStartCol(), 
                                                              df.getCountEndCol());
        return normalizedMatrix;
    }

    public BaseVector getColorCol(ScDataFrame df, String colorColName) {
        BaseVector colorCol = null;
        if (colorColName != null)
            colorCol = df.getWrapped().column(colorColName);
        return colorCol;
    }
    
    public JFrame visualizeInUMAP(ScDataFrame df,
                                  String colorColName,
                                  Integer topVarFeatures,
                                  boolean needPCA,
                                  Integer topPCs) throws Exception {
        BaseVector colorCol = getColorCol(df, colorColName);
        // Do PCA
        Matrix countMatrix = normalizeCounts(df);
        // df is in the cell x gene layout, which can be used by UMAP.
        double[][] values = getValues(countMatrix, topVarFeatures);
        UMAP umap = null;
        if (needPCA) {
            Matrix projection = prePCA(values, topPCs);
            umap = UMAP.of(convertMatrixToDoubleArray(projection));
        }
        else
            umap = UMAP.of(values);
        double[][] coordinates = MathEx.transpose(umap.coordinates);
        // From UMAP, we may lost some data points. Therefore, we need
        // to do this when the color column is assigned:
        if (colorCol != null) {
            String[] groups = new String[umap.index.length];
            for (int i = 0; i < groups.length; i++) 
                groups[i] = colorCol.get(umap.index[i]).toString();
            colorCol = StringVector.of(colorCol.name(), groups);
        }
        return plot(coordinates, "UMAP1", "UMAP2", colorCol);
    }
    
    /**
     * Select top variant genes. The passed matrix should be in cell x gene layout.
     * @param countMatrix
     * @param top
     * @return
     */
    private double[][] selectTopVariantGenes(Matrix countMatrix,
                                             int top) {
        double[][] values = SmileUtilities.convertToArrays(countMatrix);
        double[] sds = MathEx.colSds(values);
        // Sort the sds so that we can pick up the top
        List<Integer> indices = IntStream.range(0, sds.length).boxed().collect(Collectors.toList());
        Comparator<Integer> sorter = (i1, i2) -> new Double(sds[i2]).compareTo(sds[i1]);
        Collections.sort(indices, sorter);
        double[][] selected = new double[values.length][];
        for (int i = 0; i < selected.length; i++) {
            selected[i] = new double[top];
            double[] row = values[i];
            for (int j = 0; j < top; j++) {
                selected[i][j] = row[indices.get(j)];
            }
        }
        return selected;
    }
    
    /**
     * The layout of the passed ScDataFrame should be cell x gene with cell meta. This object
     * should have also countStart and countEnd columns set.
     * @param df
     * @param colorName
     * @return
     * @throws Exception
     */
    public JFrame visualizeInPCA(ScDataFrame df,
                                 String colorColName,
                                 Integer topVarGenes) throws Exception {
        BaseVector colorCol = getColorCol(df, colorColName);
        // Do PCA
        Matrix countMatrix = normalizeCounts(df);
        // df is in the cell x gene layout. However, for PCA we need to 
        // use the gene x cell layout.
        double[][] values = null;
        if (topVarGenes == null)
            values = convertMatrixToDoubleArray(countMatrix.transpose());
        else
            values = selectTopVariantGenes(countMatrix, topVarGenes);
        values = MathEx.transpose(values);
        PCA pca = PCA.fit(values).setProjection(2);
        // The project is a 2 x cell_number layout, which is used by the plot
        Matrix projection = pca.getProjection();
        double[][] coordinates = SmileUtilities.convertToArrays(projection);
        return plot(coordinates, "PC1", "PC2", colorCol);
    }
    
    private double[][] convertMatrixToDoubleArray(Matrix matrix) {
        logger.debug(String.format("Convert a matrix with %d rows, %d cols to a double[][]",
                                   matrix.nrows(),
                                   matrix.ncols()));
        return SmileUtilities.convertToArrays(matrix);
    }
    
    /**
     * Use Java reflection to get the text since it is not exposed.
     * @param legend
     * @return
     */
    private String getLegendText(Legend legend)  {
        try {
            Field field = legend.getClass().getDeclaredField("text");
            field.setAccessible(true); // By pass the access control 
            return field.get(legend).toString();
        }
        catch(Exception e) {
            logger.error(e.getMessage(), e);
        }
        return ""; // Do nothing
    }
    
}
