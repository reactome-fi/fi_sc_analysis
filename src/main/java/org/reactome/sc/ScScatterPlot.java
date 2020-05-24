package org.reactome.sc;

import java.awt.Component;
import java.awt.Rectangle;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.util.ArrayList;
import java.util.List;
import java.util.Optional;
import java.util.stream.Stream;

import javax.swing.JComponent;

import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import smile.data.DataFrame;
import smile.math.MathEx;
import smile.plot.swing.Canvas;
import smile.plot.swing.Legend;
import smile.plot.swing.PlotPanel;
import smile.plot.swing.Point;
import smile.plot.swing.ScatterPlot;

/**
 * To enhance the functions of the original ScatterPlot for scRNA-seq data.
 * @author wug
 *
 */
public class ScScatterPlot extends ScatterPlot {
    private static final Logger logger = LoggerFactory.getLogger(ScScatterPlot.class);
    private final int SENSING_DIST = 6; // Radius
    // Need to keep track its own PlotPane
    private Canvas canvas;
    private PlotPanel plotPanel;
    // For hit detect
    private DotQuadTree quadTree;
    
    public ScScatterPlot(Point[] points, 
                         Legend[] legends,
                         DataFrame df,
                         String x,
                         String y,
                         String category) {
        super(points, legends);
        setMetaInfo(df, x, y, category);
    }
    
    public ScScatterPlot(Point points,
                         DataFrame df,
                         String x,
                         String y) {
        super(points);
        setMetaInfo(df, x, y, null);
    }
    
    public void setCanvas(Canvas canvas, PlotPanel plotPane) {
        this.canvas = canvas;
        this.plotPanel = plotPane;
        customizeMouseListeners();
    }

    private void customizeMouseListeners() {
        JComponent canvasPane = getCanvasPane();
        if (canvasPane == null)
            return;
        MouseController controller = new MouseController();
        MouseListener[] canvasMouseListeners = canvasPane.getMouseListeners();
        if (canvasMouseListeners != null) {
            Stream.of(canvasMouseListeners).forEach(l -> canvasPane.removeMouseListener(l));
            controller.mouseListeners = canvasMouseListeners;
        }
        MouseMotionListener[] motionListeners = canvasPane.getMouseMotionListeners();
        if (motionListeners != null) {
            Stream.of(motionListeners).forEach(l -> canvasPane.removeMouseMotionListener(l));
            controller.motionListeners = motionListeners;
        }
        canvasPane.addMouseListener(controller);
        canvasPane.addMouseMotionListener(controller);
    }
    
    private JComponent getCanvasPane() {
        for (int i = 0; i < plotPanel.getComponentCount(); i ++) {
            Component comp = plotPanel.getComponent(i);
            if (comp instanceof JComponent &&
                comp.getClass().getName().endsWith("JCanvas")) {
                return (JComponent) comp;
            }
        }
        return null;
    }
    
    private void setMetaInfo(DataFrame df,
                             String x,
                             String y,
                             String category) {
        double[] xCoords = df.column(x).toDoubleArray();
        double[] yCoods = df.column(y).toDoubleArray();
        double minX = MathEx.min(xCoords);
        double minY = MathEx.min(yCoods);
        double maxX = MathEx.max(xCoords);
        double maxY = MathEx.max(yCoods);
        // Make sure we have some buffers so that all dots can be covered
        double bufferX = (maxX - minX) / xCoords.length;
        minX -= bufferX;
        maxX += bufferX;
        double bufferY = (maxY - minY) / yCoods.length;
        minY -= bufferY;
        maxY += bufferY;
        quadTree = new DotQuadTree(minX, minY, maxX, maxY);
        for (int i = 0; i < xCoords.length; i++) {
            double[] dot = {xCoords[i], yCoods[i]};
            NodeObject value = new NodeObject();
            value.dot = dot;
            quadTree.add(value);
        }
    }

    @Override
    public Optional<String> tooltip(double[] coord) {
//        return Optional.of(coord[0] + ", " + coord[1]);
        double[] mvRatios = getModelToViewRatios();
        if (mvRatios == null)
            return super.tooltip(coord);
        double dx = SENSING_DIST * mvRatios[0];
        double dy = SENSING_DIST * mvRatios[1];
        NodeObject value = quadTree.get(coord, 
                                        canvas.getLowerBounds(),
                                        canvas.getUpperBounds(),
                                        dx, 
                                        dy);
//        System.out.println("Found value: " + (value == null ? value : value.dot[0] + ", " + value.dot[1]));
        if (value == null)
            return super.tooltip(coord);
        return Optional.of(value.dot[0] + ", " + value.dot[1]);
    }
    
    private double[] getModelToViewRatios() {
        if (canvas == null || plotPanel == null)
            return null;
        double[] ratios = new double[2];
        int w = plotPanel.getWidth();
        int h = plotPanel.getHeight();
        ratios[0] = (canvas.getUpperBounds()[0] - canvas.getLowerBounds()[0]) / w;
        ratios[1] = (canvas.getUpperBounds()[1] - canvas.getLowerBounds()[1]) / h;
        return ratios;
    }
    
    private class MouseController implements MouseListener, MouseMotionListener {
        private MouseListener[] mouseListeners;
        private MouseMotionListener[] motionListeners;
        private int pressX;
        private int pressY;
        private int releaseX;
        private int releaseY;
        
        public MouseController() {
        }
        
        @Override
        public void mouseClicked(MouseEvent e) {
            Stream.of(mouseListeners).forEach(l -> l.mouseClicked(e));
        }
        
        public Rectangle getSelectionBox() {
            Rectangle rect = new Rectangle();
            rect.x = Math.min(pressX, releaseX);
            rect.y = Math.min(pressY, releaseY);
            rect.width = Math.abs(releaseX - pressX);
            rect.height = Math.abs(releaseY - pressY);
            return rect;
        }

        @Override
        public void mousePressed(MouseEvent e) {
            pressX = e.getX();
            pressY = e.getY();
            // Need to reset these two values
            SmileUtilities.setFieldViaReflection(e.getSource(),
                                                 "mouseDraggingX",
                                                 -1);
            SmileUtilities.setFieldViaReflection(e.getSource(),
                                                 "mouseDraggingY",
                                                 -1);
            Stream.of(mouseListeners).forEach(l -> l.mousePressed(e));
        }

        @Override
        public void mouseReleased(MouseEvent e) {
            if (e.isAltDown()) {
                releaseX = e.getX();
                releaseY = e.getY();
                logger.debug("Selection box: " + getSelectionBox());
                return; // We just want to draw the rectangle as a selection
            }
            Stream.of(mouseListeners).forEach(l -> l.mouseReleased(e));
        }

        @Override
        public void mouseEntered(MouseEvent e) {
            Stream.of(mouseListeners).forEach(l -> l.mouseEntered(e));
        }

        @Override
        public void mouseExited(MouseEvent e) {
            Stream.of(mouseListeners).forEach(l -> l.mouseExited(e));
        }
        
        @Override
        public void mouseDragged(MouseEvent e) {
            // Make sure the mouse dragged is not too sensitive to change the selection
            // box
            if (Math.abs(e.getX() - pressX) > 4 &&
                Math.abs(e.getY() - pressY) > 4) {
//                System.out.println("Mouse dragged fired: " + e.getX() + ", " + e.getY());
                Stream.of(motionListeners).forEach(l -> l.mouseDragged(e));
            }
        }

        @Override
        public void mouseMoved(MouseEvent e) {
            Stream.of(motionListeners).forEach(l -> l.mouseMoved(e));
        }
        
    }
    
    /**
     * A QuadTree used to fetch cells based on coordinates.
     * @author wug
     *
     */
    private class DotQuadTree {

        DotQuadTreeNode root = null;

        DotQuadTree(double minX,
                    double minY, 
                    double maxX, 
                    double maxY) {
            root = new DotQuadTreeNode(minX, minY, maxX, maxY);
        }

        void add(NodeObject value) {
            root.add(value);
        }
        
        NodeObject get(double[] dot,
                       double[] lowerBounds,
                       double[] upperBounds,
                       double dx,
                       double dy) {
            return root.get(dot, lowerBounds, upperBounds, dx, dy);
        }

    }
    
    private class DotQuadTreeNode {
        static final int maxValues = 32; // This is arbitrary now
        static final int maxLayer = 8;
        DotQuadTreeNode ne;
        DotQuadTreeNode nw;
        DotQuadTreeNode sw;
        DotQuadTreeNode se;
        boolean isLeaf = true;
        List<NodeObject> values;
        double minX;
        double minY;
        double maxX;
        double maxY;
        int layer = 0;

        DotQuadTreeNode(double minX, 
                        double minY, 
                        double maxX, 
                        double maxY) {
            this.minX = minX;
            this.minY = minY;
            this.maxX = maxX;
            this.maxY = maxY;
        }

        boolean contain(double x, double y) {
            // Min inclusive and max exclusive
            if (x >= minX && x < maxX && y < maxY && y >= minY)
                return true;
            return false;
        }

        boolean contain(double[] dot) {
            return contain(dot[0], dot[1]);
        }

        // Split this node into four sections
        void split() {
            double centerX = (minX + maxX) / 2.0d;
            double centerY = (minY + maxY) / 2.0d;
            ne = new DotQuadTreeNode(centerX, centerY, maxX, maxY);
            ne.layer = layer + 1;
            nw = new DotQuadTreeNode(minX, centerY, centerX, maxY);
            nw.layer = layer + 1;
            sw = new DotQuadTreeNode(minX, minY, centerX, centerY);
            sw.layer = layer + 1;
            se = new DotQuadTreeNode(centerX, minY, maxX, centerY);
            se.layer = layer + 1;
            isLeaf = false;
            values.forEach(value -> add(value));
            values.clear();
        }

        void add(NodeObject value) {
            if (isLeaf) {
                if (values == null)
                    values = new ArrayList<>();
                // If reach to the max layer, don't split again.
                if (values.size() < maxValues || layer == maxLayer) {
                    values.add(value);
                    return;
                }
                split();
            }
            DotQuadTreeNode chosen = choose(value.dot);
            chosen.add(value);
        }
        
        NodeObject get(double[] dot,
                       double[] lowerBounds,
                       double[] upperBounds,
                       double dx,
                       double dy) {
            if (isLeaf) {
                if (values == null)
                    return null;
                for (NodeObject value : values) {
                    if (value.dot[0] < lowerBounds[0] ||
                        value.dot[0] > upperBounds[0] ||
                        value.dot[1] < lowerBounds[1] ||
                        value.dot[1] > upperBounds[1])
                        continue;
                    if (Math.abs(dot[0] - value.dot[0]) <= dx &&
                        Math.abs(dot[1] - value.dot[1]) <= dy) {
                            return value;
                        }
                }
                return null; // Nothing we can do.
            }
            DotQuadTreeNode child = choose(dot);
            if (child == null)
                return null; // This should never occur
            return child.get(dot, lowerBounds, upperBounds, dx, dy);
        }

        DotQuadTreeNode choose(double[] dot) {
            if (ne.contain(dot))
                return ne;
            if (nw.contain(dot))
                return nw;
            if (sw.contain(dot))
                return sw;
            return se;
        }

    }
    
    private class NodeObject {
        double[] dot;
    }
    
}
