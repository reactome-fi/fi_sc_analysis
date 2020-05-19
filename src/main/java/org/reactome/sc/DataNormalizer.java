package org.reactome.sc;

import smile.data.DataFrame;
import smile.math.MathEx;
import smile.math.matrix.DenseMatrix;
import smile.math.matrix.Matrix;

/**
 * This class is used to normalize the data.
 * @author wug
 *
 */
public class DataNormalizer {
    
    public DataNormalizer() {
    }
    
    /**
     * This method is ported from the python code:
     *  https://scanpy.readthedocs.io/en/stable/api/scanpy.pp.normalize_per_cell.html.
     * It is assumed that the layout of the matrix is gene x cell
     */
    public DenseMatrix normalizePerCell(DenseMatrix matrix,
                                        CountMatrixLayout layout) {
        // Clone the matrix to avoid override the original matrix
        DenseMatrix dmatrix = matrix.clone();
        if (layout == CountMatrixLayout.CELL_TIMES_GENE)
            dmatrix = dmatrix.transpose();
        // The following code is implemented for gene x cell layout
        double[] colSums = dmatrix.colSums();
        // The following call will rearrange the input. We need to copy
        double[] copy = new double[colSums.length];
        MathEx.copy(colSums, copy);
        double median = MathEx.median(copy);
        for (int i = 0; i < colSums.length; i++)
            colSums[i] = colSums[i] / median;
        dmatrix = dmatrix.scale(null, colSums);
        log1p(dmatrix);
        // Convert it back 
        if (layout == CountMatrixLayout.CELL_TIMES_GENE)
            dmatrix = dmatrix.transpose();
        return dmatrix;
    }
    
    /**
     * Perform normalizion for a DataFrame. The returned DataFrame will have the columns having
     * counts only.
     * @param df
     * @param countsColStart
     * @param countsColEnd inclusive
     * @return
     */
    public Matrix normalizePerCell(DataFrame df,
                                   int countsColStart,
                                   int countsColEnd) {
        DataFrame countsDF = SmileUtilities.select(df, 
                                                   countsColStart, 
                                                   countsColEnd);
        DenseMatrix matrix = countsDF.toMatrix();
        matrix = normalizePerCell(matrix, 
                                  CountMatrixLayout.CELL_TIMES_GENE);
        return matrix;
    }
    
    /**
     * A simple log (natural) transformation after adding 1.0 for each element
     * in the matrix.
     * @param matrix
     */
    private void log1p(DenseMatrix matrix) {
        // Need to  log1 transformation
        for (int i = 0; i < matrix.nrows(); i ++)
            for (int j = 0; j < matrix.ncols(); j++)
                matrix.set(i, j, MathEx.log(matrix.get(i, j) + 1.0d));
    }

}
