package org.reactome.sc;

import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.IntStream;

import smile.data.DataFrame;
import smile.data.type.DataType;
import smile.data.type.DataTypes;
import smile.data.vector.BaseVector;
import smile.data.vector.DoubleVector;
import smile.math.matrix.Matrix;

/**
 * Utility methods to make data structure conversion easy.
 * @author wug
 *
 */
public class SmileUtilities {

    /**
     * Convert a matrix into a DataFrame object. If it is needed, some of columns
     * may be copied from the passed original DataFrame.
     * @param source
     * @param data
     * @param keptCols
     * @return
     */
    public static DataFrame convertMatrixToDataFrame(Matrix data,
                                                     String[] colNames,
                                                     DataFrame original,
                                                     Integer... keptCols) {
        if (colNames.length != data.ncols())
            throw new IllegalArgumentException("The length of colNames is not the same as the number of the matrix columns!");
        if (original != null && data.nrows() != original.nrows())
            throw new IllegalArgumentException("The number of the matrix rows are not the same as the original rows!");
        DoubleVector[] dataCols = new DoubleVector[colNames.length];
        for (int i = 0; i < dataCols.length; i++) {
            DoubleVector dv = DoubleVector.of(colNames[i],
                                              getColumn(data, i));
            dataCols[i] = dv;
        }
        DataFrame rtn = DataFrame.of(dataCols);
        // Need to get the more columns
        if (original != null && keptCols != null) {
            @SuppressWarnings("rawtypes")
            BaseVector[] otherCols = new BaseVector[keptCols.length];
            for (int i = 0; i < otherCols.length; i++)
                otherCols[i] = original.column(keptCols[i]);
            rtn.merge(otherCols);
        }
        return rtn;
    }
    
    /**
     * Get a column vector in a double array from a Matrix object.
     * @param matrix
     * @param col
     * @return
     */
    public static double[] getColumn(Matrix matrix, int col) {
        double[] rtn = new double[matrix.nrows()];
        for (int i = 0; i < rtn.length; i++) 
            rtn[i] = matrix.get(i, col);
        return rtn;
    }
    
    /**
     * Convert a Matrix object into an two-dimensional double array.
     * @param matrix
     * @return
     */
    public static double[][] convertToArrays(Matrix matrix) {
        // Put the matrix into a double array
        double[][] values = new double[matrix.nrows()][];
        for (int i = 0; i < values.length; i++) {
            values[i] = new double[matrix.ncols()];
            for (int j = 0; j < values[i].length; j++)
                values[i][j] = matrix.get(i, j);
        }
        return values;
    }

    /**
     * Transpose the matrix with the colnames as the rownames and the rownames as the colnames. The data matris is 
     * transposed. Since there is no actual row names in the DataFrame, we will use the first column as the row names.
     * @param df
     * @return
     */
    public static DataFrame transpose(DataFrame df) {
        MatrixDataFrame mdf = new MatrixDataFrame().fromDataFrame(df, 0);
        mdf = mdf.transpose();
        return mdf.toDataFrame();
    }
    
    /**
     * Select a subset of a DataFrame object.
     * @param df
     * @param startCol inclusive
     * @param endCol exclusive
     * @return
     */
    public static DataFrame select(DataFrame df, 
                                   int startCol,
                                   int endCol) {
        int[] selectedCols = IntStream.range(startCol, endCol + 1).toArray();
        DataFrame countsDF = df.select(selectedCols);
        return countsDF;
    }
    
}
