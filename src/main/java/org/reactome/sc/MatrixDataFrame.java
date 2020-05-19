package org.reactome.sc;

import java.util.Arrays;

import smile.data.DataFrame;
import smile.data.vector.BaseVector;
import smile.data.vector.DoubleVector;
import smile.math.matrix.Matrix;

/**
 * This is a wrapper for a smile Matrix object with row names and col names.
 * @author wug
 *
 */
public class MatrixDataFrame {
    private String[] colNames;
    private String[] rowNames;
    private Matrix matrix;
    
    public MatrixDataFrame() {
    }
    
    public MatrixDataFrame transpose() {
        MatrixDataFrame rtn = new MatrixDataFrame();
        rtn.matrix = matrix.clone().transpose();
        rtn.colNames = Arrays.copyOf(rowNames, rowNames.length);
        rtn.rowNames = Arrays.copyOf(colNames, colNames.length);
        return rtn;
    }
    
    /**
     * Convert it into a smile DataFrame object.
     * @return
     */
    @SuppressWarnings("rawtypes")
    public DataFrame toDataFrame() {
        BaseVector[] cols = new BaseVector[matrix.ncols()];
        for (int i = 0; i < cols.length; i++) {
            cols[i] = DoubleVector.of(colNames[i],
                                      SmileUtilities.getColumn(matrix, i));
        }
        return DataFrame.of(cols);
    }
    
    /**
     * Convert a DataFrame object into a MatrixDataFrame object.
     * @param frame
     * @return
     */
    @SuppressWarnings("rawtypes")
    public MatrixDataFrame fromDataFrame(DataFrame frame,
                                         int rowNameCol) {
        MatrixDataFrame rtn = new MatrixDataFrame();
        BaseVector rowNames = frame.column(rowNameCol);
        frame = frame.drop(rowNameCol);
        Matrix matrix = frame.toMatrix();
        rtn.matrix = matrix;
        rtn.colNames = frame.names();
        rtn.rowNames = new String[rowNames.size()];
        for (int i = 0; i < rtn.rowNames.length; i++)
            rtn.rowNames[i] = rowNames.get(i).toString();
        return rtn;
    }

    public String[] getColNames() {
        return colNames;
    }

    public void setColNames(String[] colNames) {
        this.colNames = colNames;
    }

    public String[] getRowNames() {
        return rowNames;
    }

    public void setRowNames(String[] rowNames) {
        this.rowNames = rowNames;
    }

    public Matrix getMatrix() {
        return matrix;
    }

    public void setMatrix(Matrix matrix) {
        this.matrix = matrix;
    }
    
}
