/*
 * www.javagl.de - bh-tsne
 *
 * Copyright (c) 2020 Marco Hutter - http://www.javagl.de
 */ 
package de.javagl.tsne.demo;

/**
 * An interface for this package that represents input data for t-SNE
 */
interface TsneData
{
    /**
     * Returns the actual data
     * 
     * @return The data
     */
    double[][] getData();
    
    /**
     * Returns the labels for the data. If there are no labels, then this
     * array contains only <code>null</code> entries
     * 
     * @return The labels
     */
    String[] getLabels();
}
