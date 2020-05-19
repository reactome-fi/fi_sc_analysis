package org.reactome.sc;

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
    
    public ScDataFrame() {
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

}
