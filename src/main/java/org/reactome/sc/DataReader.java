package org.reactome.sc;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.log4j.Logger;

import smile.data.DataFrame;
import smile.data.type.DataTypes;
import smile.data.type.StructField;
import smile.data.type.StructType;
import smile.data.vector.BaseVector;
import smile.data.vector.DoubleVector;
import smile.data.vector.StringVector;
import smile.io.Read;
import smile.math.matrix.Matrix;

/**
 * This class is used to read and write matrix storing the expression counts.
 * @author wug
 *
 */
public class DataReader {
    private Logger logger = Logger.getLogger(DataReader.class);
    
    public DataReader() {
    }
    
    /**
     * Read an matrix in the MatrixMarket format.
     * @param file
     * @return
     */
    public Matrix readMatrixMarketFile(String file) throws IOException, ParseException {
        Matrix matrix = Matrix.market(Paths.get(file));
        logger.debug(String.format("Loaded matrix with %d rows, %d colmns:",
                     matrix.nrows(),
                     matrix.ncols()));
        return matrix; 
    }
    
    /**
     * This matrix file should be generated from R with rows for genes and columns for cells.
     * It is expected to have the first row for cell ids and the first column for gene ids.
     * The first cell (0, 0) should be empty. It is also assumed that there is no quotation marks used.
     * @param file
     * @return
     * @throws IOException
     * @throws ParseException
     */
    private DataFrame readMatrixInCSV(String file) throws IOException, ParseException {
        // Get the first line for structure type
        BufferedReader br = Files.newBufferedReader(Paths.get(file));
        String line = br.readLine();
        br.close();
        // This is the first line
        CSVParser parser = CSVParser.parse(line, CSVFormat.DEFAULT);
        List<StructField> fields = new ArrayList<>();
        for (CSVRecord record : parser) {
            for (int i = 0; i < record.size(); i++) {
                String text = record.get(i);
                if (i == 0)
                    fields.add(new StructField(text, DataTypes.StringType));
                else
                    fields.add(new StructField(text, DataTypes.IntegerType));
            }
        }
        CSVFormat format = CSVFormat.DEFAULT.withAllowMissingColumnNames().withFirstRecordAsHeader();
        StructType type = new StructType(fields);
        DataFrame df = Read.csv(Paths.get(file),
                                format,
                                type);
        return df;
    }
    
    /**
     * Read a matrix as a DataFrame object. The file type is defined by the extension. The first row
     * is the column names and the first column is for the row names.
     * @param fileName
     * @return
     * @throws IOException
     * @throws ParseException
     */
    private ScDataFrame readMatrixDataFrame(String fileName) throws IOException, ParseException {
        DataFrame df = null;
        if (isCSVFile(fileName))
            df = readMatrixInCSV(fileName);
        else
            df = readMatrixInTSV(fileName); // Default using the tab
        ScDataFrame rtn = new ScDataFrame();
        rtn.setWrapped(df);
        return rtn;
    }
    
    private boolean isCSVFile(String fileName) {
        // Check the extension
        File file = new File(fileName);
        String name = file.getName();
        int index = name.lastIndexOf(".");
        String ext = name.substring(index + 1).toLowerCase();
        return ext.equals("csv");
    }
    
    private DataFrame readMatrixInTSV(String file) throws IOException, ParseException {
        // Get the first line for structure type
        BufferedReader br = Files.newBufferedReader(Paths.get(file));
        String line = br.readLine();
        String[] colNames = line.split("\t");
        @SuppressWarnings("rawtypes")
        BaseVector[] cols = new BaseVector[colNames.length];
        List<String> names = new ArrayList<>(); // For the first column
        List<double[]> counts = new ArrayList<>();
        while ((line = br.readLine()) != null) {
            String[] tokens = line.split("\t");
            names.add(tokens[0]);
            double[] rowCount = new double[tokens.length - 1];
            for (int i = 1; i < tokens.length; i++)
                rowCount[i - 1] = new Double(tokens[i]);
            counts.add(rowCount);
        }
        cols[0] = StringVector.of(colNames[0], names.toArray(new String[] {}));
        for (int i = 1; i < colNames.length; i++) {
            double[] colCount = new double[counts.size()];
            for (int j = 0; j < colCount.length; j++)
                colCount[j] = counts.get(j)[i - 1];
            cols[i] = DoubleVector.of(colNames[i], colCount);
        }
        br.close();
        return DataFrame.of(cols);
    }
    
    /**
     * The first file should be the count file with col as cells and rows as genes. The first row is
     * for for column names and the first column for row names. The rows for cells and the columns for
     * column features in the second file with the first column and the first row are used for annotation.
     * @param countFile
     * @param cellFile
     * @return The layout is cell x gene with meta
     * @throws IOException
     * @throws URISyntaxException
     * @throws ParseException
     */
    public ScDataFrame readCountsFiles(String countFile,
                                       CountMatrixLayout layout,
                                       String cellFile) throws IOException, URISyntaxException, ParseException {
        ScDataFrame countDF = readMatrixDataFrame(countFile);
        // We are referring the actual wrapped DataFrame object.
        DataFrame df = countDF.getWrapped();
        if (layout == CountMatrixLayout.GENE_TIMES_CELL) {
            df = SmileUtilities.transpose(df); 
            // Need to re-assign
            countDF.setWrapped(df);
        }
        countDF.setLayout(CountMatrixLayout.CELL_TIMES_GENE);
        // Now we can set counts start and end columns
        countDF.setCountStartCol(1); // Since we have the first row names
        countDF.setCountEndCol(countDF.getWrapped().ncols() - 1); // Inclusive
        if (cellFile != null) {
            DataFrame cellDF = readDataFrame(cellFile);
            cellDF = cellDF.drop(0); // We don't need the first column
            DataFrame merged = df.merge(cellDF);
            countDF.setWrapped(merged);
        }
        return countDF;
    }
    
    /**
     * Read a generic DataFrame file. Currently two formats are supported:
     * csv with the file extension ends with .csv or tsv. For tsv files, the
     * file names may not be ended with .tsv.
     * @param fileName
     * @return
     * @throws IOException
     * @throws ParseException
     * @throws URISyntaxException
     */
    public DataFrame readDataFrame(String fileName) throws IOException, ParseException, URISyntaxException{
        DataFrame df = null;
        if (isCSVFile(fileName)) {
            // Need to read the cell annotation file
            CSVFormat format = CSVFormat.DEFAULT.withAllowMissingColumnNames().withFirstRecordAsHeader();
            df = Read.csv(fileName, format);
        }
        else
            df = readDataFrameInTSV(fileName);
        return df;
    }
    
    private DataFrame readDataFrameInTSV(String fileName) throws IOException {
        // Get the first line for structure type
        BufferedReader br = Files.newBufferedReader(Paths.get(fileName));
        String line = br.readLine();
        String[] colNames = line.split("\t");
        @SuppressWarnings("rawtypes")
        BaseVector[] cols = new BaseVector[colNames.length];
        List<String[]> rows = new ArrayList<>();
        while ((line = br.readLine()) != null) {
            String[] tokens = line.split("\t");
            rows.add(tokens);
        }
        for (int i = 0; i < colNames.length; i++) {
            String[] col = new String[rows.size()];
            for (int j = 0; j < col.length; j++)
                col[j] = rows.get(j)[i];
            cols[i] = StringVector.of(colNames[i], col);
        }
        br.close();
        return DataFrame.of(cols);
    }

}
