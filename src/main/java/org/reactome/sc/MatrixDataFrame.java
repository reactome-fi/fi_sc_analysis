package org.reactome.sc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import smile.data.DataFrame;
import smile.data.vector.BaseVector;
import smile.data.vector.DoubleVector;
import smile.math.MathEx;
import smile.math.matrix.Matrix;

/**
 * This is a wrapper for a smile Matrix object with row names and col names.
 * @author wug
 *
 */
@SuppressWarnings("rawtypes")
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
    
    public void setColNames(BaseVector vector) {
        this.colNames = new String[vector.size()];
        for (int i = 0; i < colNames.length; i++)
            colNames[i] = vector.get(i).toString();
    }

    public String[] getRowNames() {
        return rowNames;
    }

    public void setRowNames(String[] rowNames) {
        this.rowNames = rowNames;
    }
    
    public void setRowNames(BaseVector vector) {
        rowNames = new String[vector.size()];
        for (int i = 0; i < rowNames.length; i++)
            rowNames[i] = vector.get(i).toString();
    }

    public Matrix getMatrix() {
        return matrix;
    }

    public void setMatrix(Matrix matrix) {
        this.matrix = matrix;
    }
    
    public void filterRows(Double minValue) {
        double[][] values = SmileUtilities.convertToArrays(matrix);
        double[] rowSums = MathEx.rowSums(values);
        List<Integer> dropped = new ArrayList<>();
        for (int i = 0; i < rowSums.length; i++) {
            if (rowSums[i] < minValue) {
                dropped.add(i);
            }
        }
        if (dropped.size() == 0)
            return; // Nothing to do
        double[][] copy = new double[rowSums.length - dropped.size()][];
        String[] newRowNames = new String[rowNames.length - dropped.size()];
        int index = 0;
        for (int i = 0; i < values.length; i++) {
            if (dropped.contains(i))
                continue;
            copy[index] = values[i];
            newRowNames[index] = rowNames[i];
            index ++;
        }
        this.matrix = Matrix.of(copy);
        this.rowNames = newRowNames;
    }
    
    public void filterCols(Double minValue) {
        double[][] values = SmileUtilities.convertToArrays(matrix);
        double[] colSums = MathEx.colSums(values);
        List<Integer> dropped = new ArrayList<>();
        for (int i = 0; i < colSums.length; i++) {
            if (colSums[i] < minValue) {
                dropped.add(i);
            }
        }
        if (dropped.size() == 0)
            return; // Do nothing
        double[][] copy = new double[matrix.nrows()][];
        int index = 0;
        Set<Integer> dropSet = new HashSet<>(dropped);
        for (int i = 0; i < copy.length; i++) {
            copy[i] = new double[colSums.length - dropped.size()];
            index = 0;
            for (int j = 0; j < matrix.ncols(); j++) {
                if (dropSet.contains(j))
                    continue;
                copy[i][index] = matrix.get(i, j);
                index ++;
            }
        }
        String[] newColNames = new String[colNames.length - dropped.size()];
        index = 0;
        for (int i = 0; i < colNames.length; i++) {
            if (dropSet.contains(i))
                continue;
            newColNames[index ++] = colNames[i];
        }
        this.matrix = Matrix.of(copy);
        this.colNames = newColNames;
    }
    
}
