package org.reactome.sc;

import java.awt.Color;
import java.lang.reflect.Field;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import smile.data.DataFrame;
import smile.data.vector.BaseVector;
import smile.data.vector.DoubleVector;
import smile.math.matrix.DenseMatrix;
import smile.math.matrix.Matrix;
import smile.plot.swing.Legend;
import smile.plot.swing.Palette;
import smile.plot.swing.Point;
import smile.plot.swing.ScatterPlot;

/**
 * Utility methods to make data structure conversion easy.
 * @author wug
 *
 */
public class SmileUtilities {
    private static final Logger logger = LoggerFactory.getLogger(SmileUtilities.class);
    
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
    
    public static DenseMatrix getCountsMatrix(ScDataFrame df) {
        DataFrame countsDF = SmileUtilities.select(df.getWrapped(), 
                                                   df.getCountStartCol(), 
                                                   df.getCountEndCol());
        DenseMatrix matrix = countsDF.toMatrix();
        return matrix;
    }
    
    public static double[][] getCounts(ScDataFrame df) {
        DenseMatrix matrix = getCountsMatrix(df);
        return convertToArrays(matrix);
    }
    
    /**
     * This method is a merge of two static methods in the smile's ScatterPlot so that we can enhance the original
     * ScatterPlot by subclassing it.
     * @param data
     * @param x
     * @param y
     * @param category
     * @param mark
     * @return
     */
    private static ScatterPlot getScatterPlot(DataFrame data, 
                                              String x, 
                                              String y, 
                                              String category, 
                                              char mark,
                                              Color color) {
        int ix = data.columnIndex(x);
        int iy = data.columnIndex(y);
        double[][] xy = data.stream().map(row -> new double[]{row.getDouble(ix), row.getDouble(iy)}).toArray(double[][]::new);
        if (category == null)
            return new ScScatterPlot(new Point(xy, mark, color), data, x, y);
        String[] label = data.column(category).toStringArray();
        Map<String, List<Integer>> groups = IntStream.range(0, xy.length).boxed().collect(Collectors.groupingBy(i -> label[i]));
        Point[] points = new Point[groups.size()];
        Legend[] legends = new Legend[groups.size()];
        int k = 0;
        for (Map.Entry<String, List<Integer>> group : groups.entrySet()) {
            color = Palette.COLORS[k % Palette.COLORS.length];
            points[k] = new Point(
                    group.getValue().stream().map(i -> xy[i]).toArray(double[][]::new),
                    mark,
                    color
            );
            legends[k] = new Legend(group.getKey(), color);
            k++;
        }
        ScScatterPlot plot = new ScScatterPlot(points, legends, data, x, y, category);
        return plot;
    }
    
    public static ScatterPlot getScatterPlot(DataFrame data,
                                            String x,
                                            String y,
                                            String category,
                                            char mark) {
        return getScatterPlot(data, x, y, category, mark, null);
    }
    
    public static ScatterPlot getScatterPlot(DataFrame data, 
                                             String x,
                                             String y, 
                                             char mark,
                                             Color color) {
        return getScatterPlot(data, x, y, null, mark, color);
    }
    
    public static void setFieldViaReflection(Object obj, 
                                             String fieldName,
                                             Object newValue) {
        try {
            Class<?> cls = obj.getClass();
            // Need to use declared field 
            Field field = cls.getDeclaredField(fieldName);
            field.setAccessible(true);
            field.set(obj, newValue);
//            logger.debug("Reset mouse dragging coordinates.");
        }
        catch(Exception e) {
            logger.error(e.getMessage(), e);
        }
    }
}
