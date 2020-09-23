/*
 * www.javagl.de - bh-tsne
 *
 * Copyright (c) 2020 Marco Hutter - http://www.javagl.de
 */ 
package de.javagl.tsne.demo;

import java.io.IOException;
import java.util.Arrays;

import de.javagl.tsne.Tsne;

public class BhtsneTest
{
    private static String fileName = 
        "src/main/resources/data/"
        + "training-2500-with-250-from-each-normalized.arff.gz";
    
    public static void main(String[] args) throws IOException
    {
        testSimple();
    }

    private static void testSimple() throws IOException
    {
        TsneData tsneData = WekaUtils.readFile(fileName);
        Tsne tsne = new Tsne();
        double[][] Y = tsne.run(tsneData.getData());

        System.out.println(Arrays.toString(Y[0]));
    	System.out.println(Arrays.toString(Y[1]));
    	System.out.println(Arrays.toString(Y[2]));
    }
    

}
