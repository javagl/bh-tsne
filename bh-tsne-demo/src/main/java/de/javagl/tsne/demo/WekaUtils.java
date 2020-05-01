/*
 * www.javagl.de - bh-tsne
 *
 * Copyright (c) 2020 Marco Hutter - http://www.javagl.de
 */ 
package de.javagl.tsne.demo;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.List;

import weka.core.Attribute;
import weka.core.Instance;
import weka.core.Instances;
import weka.core.converters.ArffLoader;

/**
 * Weka utilities for the demo package
 */
class WekaUtils
{
    /**
     * Read the {@link TsneData} from the given file
     * 
     * @param inputFile The input file
     * @return The data
     * @throws IOException If an IO error occurs
     */
    static TsneData read(Path inputFile) throws IOException
    {
        InputStream fis = null;
        try
        {
            ArffLoader arffLoader = new ArffLoader();
            arffLoader.setFile(inputFile.toFile());
            Instances instances = arffLoader.getDataSet();
            return extractData(instances);
        }
        finally
        {
            if (fis != null)
            {
                fis.close();
            }
        }
    }
    
    /**
     * Extract the {@link TsneData} from the given Weka instances
     * 
     * @param instances The instances
     * @return The {@link TsneData}
     */
    private static TsneData extractData(Instances instances)
    {
        List<Integer> numericIndices = new ArrayList<Integer>();
        List<Integer> nominalIndices = new ArrayList<Integer>();
        for (int i = 0; i < instances.numAttributes(); i++)
        {
            Attribute attribute = instances.attribute(i);
            if (attribute.isNominal())
            {
                nominalIndices.add(i);
            }
            if (attribute.isNumeric())
            {
                numericIndices.add(i);
            }
        }
        
        double data[][] = new double[instances.size()][numericIndices.size()];
        String labels[] = new String[instances.size()];
        for (int i = 0; i < instances.size(); i++)
        {
            Instance instance = instances.get(i);
            for (int a = 0; a < numericIndices.size(); a++)
            {
                int index = numericIndices.get(a);
                data[i][a] = instance.value(index);
            }
            if (!nominalIndices.isEmpty())
            {
                int index = nominalIndices.get(0);
                int label = (int)instance.value(index);
                labels[i] = String.valueOf(label);
            }
        }
        
        return new TsneData()
        {
            @Override
            public double[][] getData()
            {
                return data;
            }

            @Override
            public String[] getLabels()
            {
                return labels;
            }
            
        };
    }
    
    /**
     * Private constructor to prevent instantiation
     */
    private WekaUtils()
    {
        // Private constructor to prevent instantiation
    }
}
