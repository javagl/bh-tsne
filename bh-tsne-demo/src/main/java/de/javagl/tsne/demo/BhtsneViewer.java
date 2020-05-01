/*
 * www.javagl.de - bh-tsne
 *
 * Copyright (c) 2020 Marco Hutter - http://www.javagl.de
 */ 
package de.javagl.tsne.demo;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.GridLayout;
import java.awt.Point;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseWheelEvent;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.Point2D;
import java.awt.image.BufferedImage;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.stream.DoubleStream;

import javax.swing.JPanel;

import de.javagl.colors.brewer.ColorBrewer;
import de.javagl.colors.brewer.ColorSchemeType;
import de.javagl.colors.maps.indexed.IndexedColorMap1D;
import de.javagl.colors.maps.indexed.IndexedColorMaps;
import de.javagl.viewer.Painter;
import de.javagl.viewer.Viewer;

class BhtsneViewer extends JPanel
{
    private static final long serialVersionUID = 1L;

    private final Viewer viewer;
    private final IndexedColorMap1D colorMap;
    private final Map<String, Integer> colorIndices;

    private double[][] input;
    private double[][] result;
    private String[] labels;

    private Point popupPosition = new Point();
    private BufferedImage popupImage;
    
    public BhtsneViewer()
    {
        super(new GridLayout(1, 1));
        this.viewer = new Viewer();
        this.colorMap = IndexedColorMaps.wrapping(
            IndexedColorMaps.create(ColorBrewer.get(
                ColorSchemeType.QUALITATIVE, "Paired", 10)));
        
        this.colorIndices = new LinkedHashMap<String, Integer>();
        
        Painter painter = new Painter()
        {
            @Override
            public void paint(Graphics2D g, AffineTransform worldToScreen,
                double w, double h)
            {
                paintResults(g, worldToScreen);
            }
        };
        viewer.addPainter(painter);
        
        Painter popupPainter = new Painter()
        {
            @Override
            public void paint(Graphics2D g, AffineTransform worldToScreen,
                double w, double h)
            {
                if (popupImage != null)
                {
                    g.drawImage(popupImage, 
                        popupPosition.x - popupImage.getWidth(),
                        popupPosition.y - popupImage.getHeight(), null);
                }
            }
        };
        viewer.addPainter(popupPainter);
        
        viewer.setMaintainAspectRatio(false);
        add(viewer);
        
        MouseAdapter mouseAdapter = new MouseAdapter()
        {
            @Override
            public void mouseWheelMoved(MouseWheelEvent e)
            {
                popupImage = null;
                viewer.repaint();
            }

            @Override
            public void mouseDragged(MouseEvent e)
            {
                popupImage = null;
                viewer.repaint();
            }
            
            @Override
            public void mouseMoved(MouseEvent e)
            {
                popupImage = null;
                viewer.repaint();
                if (input == null)
                {
                    return;
                }
                if (input[0].length != 784)
                {
                    // This works only for MNIST...
                    return;
                }
                AffineTransform wts = viewer.getWorldToScreen();
                Point2D p = new Point2D.Double(0, 0);
                double r = 6.0;
                double rr = r * r;
                for (int i = 0; i < result.length; i++)
                {
                    double x = result[i][0];
                    double y = result[i][1];
                    p.setLocation(x, y);
                    wts.transform(p, p);
                    if (p.distanceSq(e.getX(), e.getY()) < rr)
                    {
                        popupImage = createMnistImage(input[i]);
                        popupPosition.setLocation(e.getPoint());
                    }
                }
            }
        };
        viewer.addMouseWheelListener(mouseAdapter);
        viewer.addMouseMotionListener(mouseAdapter);
    }
    
    private static BufferedImage createMnistImage(double data[])
    {
        double min = DoubleStream.of(data).min().orElse(0.0);
        double max = DoubleStream.of(data).max().orElse(1.0);
        int w = 28;
        int h = 28;
        BufferedImage image = new BufferedImage(
            w, h, BufferedImage.TYPE_INT_ARGB);
        for (int y = 0; y < h; y++)
        {
            for (int x = 0; x < w; x++)
            {
                double d = data[y * w + x];
                double n = (d - min) / (max - min);
                int b = (int)(n * 255);
                int argb = 0xFF000000 | (b << 16) | (b << 8) | b;
                image.setRGB(x, y, argb);
            }
        }
        return image;
    }

    private void paintResults(Graphics2D g, AffineTransform worldToScreen)
    {
        if (result == null)
        {
            return;
        }
        g.setColor(Color.BLUE);
        Point2D p = new Point2D.Double();
        Ellipse2D e = new Ellipse2D.Double(0, 0, 0, 0);
        double r = 3.0;
        for (int i = 0; i < result.length; i++)
        {
            if (labels != null)
            {
                String label = labels[i];
                Integer colorIndex = colorIndices.get(label);
                if (colorIndex == null)
                {
                    colorIndex = colorIndices.size();
                    colorIndices.put(label, colorIndex);
                }
                Color color = colorMap.getColor(colorIndex);
                g.setColor(color);
            }
            p.setLocation(result[i][0], result[i][1]);
            worldToScreen.transform(p, p);
            e.setFrame(p.getX() - r, p.getY() - r, r + r, r + r);
            g.fill(e);
        }
    }
    
    void setResult(double input[][], double result[][], String[] labels)
    {
        this.input = input;
        this.result = result;
        this.labels = labels;
        this.colorIndices.clear();
        
        if (result == null)
        {
            viewer.setDisplayedWorldArea(0, 0, 1, 1);
        }
        else
        {
            double minX = Double.POSITIVE_INFINITY;
            double minY = Double.POSITIVE_INFINITY;
            double maxX = Double.NEGATIVE_INFINITY;
            double maxY = Double.NEGATIVE_INFINITY;
            for (int i = 0; i < result.length; i++)
            {
                double x = result[i][0];
                double y = result[i][1];
                minX = Math.min(minX, x);
                minY = Math.min(minY, y);
                maxX = Math.max(maxX, x);
                maxY = Math.max(maxY, y);
            }
            double w = maxX - minX;
            double h = maxY - minY;
            double x0 = minX - 0.05 * w;
            double y0 = minY - 0.05 * h;
            double x1 = maxX + 0.05 * w;
            double y1 = maxY + 0.05 * h;
            //System.out.println(x0 + " " + y0 + " " + x1 + " " + y1 + " size "
            //    + (x1 - x0) + ", " + (y1 - y0));
            viewer.setDisplayedWorldArea(x0, y0, x1 - x0, y1 - y0);
        }
        repaint();
    }
}
