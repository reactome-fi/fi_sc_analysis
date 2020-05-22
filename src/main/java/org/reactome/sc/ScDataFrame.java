package org.reactome.sc;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;

import smile.data.DataFrame;

/**
 * A simple wrapper around the smile's DataFrame object to provide sc specific information. It is assumed that
 * this is a DataFrame object for counts with other meta information.
 * @author wug
 *
 */
public class ScDataFrame {
    // The wrapped DataFrame
    private DataFrame wrapped;
    // The first column for counts
    private Integer countStartCol;
    // The last column for counts
    private Integer countEndCol;
    // The layout of the data frame.
    private CountMatrixLayout layout = CountMatrixLayout.CELL_TIMES_GENE;
    // A label for plot
    private String name;
    
    public ScDataFrame() {
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public DataFrame getWrapped() {
        return wrapped;
    }

    public void setWrapped(DataFrame wrapped) {
        this.wrapped = wrapped;
    }

    public Integer getCountStartCol() {
        return countStartCol;
    }

    public void setCountStartCol(Integer countStartCol) {
        this.countStartCol = countStartCol;
    }

    public Integer getCountEndCol() {
        return countEndCol;
    }

    public void setCountEndCol(Integer countEndCol) {
        this.countEndCol = countEndCol;
    }

    public CountMatrixLayout getLayout() {
        return layout;
    }

    public void setLayout(CountMatrixLayout layout) {
        this.layout = layout;
    }
    
    /**
     * Write the wrapped ScDataFrame into a file, delimited by tab.
     * The layout should be cells x genes
     * @param fileName
     * @throws IOException
     */
    public void write(String fileName) throws IOException {
        PrintWriter writer = new PrintWriter(fileName);
        StringBuilder builder = new StringBuilder();
        Arrays.asList(wrapped.names()).forEach(name -> builder.append(name).append("\t"));
        builder.deleteCharAt(builder.length() - 1);
        writer.println(builder.toString());
        builder.setLength(0);
        for (int i = 0; i < wrapped.nrows(); i++) {
            for (int j = 0; j < wrapped.ncols();j ++)
                builder.append(wrapped.get(i, j)).append("\t");
            builder.deleteCharAt(builder.length() - 1);
            writer.println(builder.toString());
            builder.setLength(0);
        }
        writer.close();
    }

}
