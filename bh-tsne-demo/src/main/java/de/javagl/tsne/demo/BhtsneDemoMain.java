/*
 * www.javagl.de - bh-tsne
 *
 * Copyright (c) 2020 Marco Hutter - http://www.javagl.de
 */
package de.javagl.tsne.demo;

import java.awt.BorderLayout;

import javax.swing.JFrame;
import javax.swing.SwingUtilities;

public class BhtsneDemoMain
{
    public static void main(String[] args)
    {
        SwingUtilities.invokeLater(() -> createAndShowGui());
    }

    private static void createAndShowGui()
    {
        JFrame f = new JFrame("bh-tsne demo");
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.getContentPane().setLayout(new BorderLayout());

        BhtsneDemo bhtsneDemo = new BhtsneDemo();
        f.getContentPane().add(bhtsneDemo.getMainPanel(), BorderLayout.CENTER);

        f.setSize(1024, 768);
        f.setLocationRelativeTo(null);
        f.setVisible(true);
    }
}
