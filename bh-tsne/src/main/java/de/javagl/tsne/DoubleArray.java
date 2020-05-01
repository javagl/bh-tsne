/*
 * www.javagl.de - bh-tsne
 *
 * Copyright (c) 2020 Marco Hutter - http://www.javagl.de
 */ 
package de.javagl.tsne;

// (bh-tsne) : Introduced to reduce the number of double[] allocations
interface DoubleArray
{
    double get(int index);
    int length();
}