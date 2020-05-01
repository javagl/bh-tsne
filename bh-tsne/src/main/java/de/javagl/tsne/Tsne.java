/*
 * www.javagl.de - bh-tsne
 *
 * Copyright (c) 2020 Marco Hutter - http://www.javagl.de
 */ 
package de.javagl.tsne;

import java.util.function.Consumer;
import java.util.function.DoubleConsumer;

/**
 * A class offering the interface to the underlying t-SNE implementation.
 */
public final class Tsne
{
    /**
     * The number of output dimensions
     */
    private int outputDims;
    
    /**
     * The number of initial dimensions
     */
    private int initialDims;
    
    /**
     * The perplexity
     */
    private double perplexity;

    /**
     * Theta. Well. Look at the paper. Whatever.
     */
    private double theta;
    
    /**
     * The maximum number of iterations
     */
    private int maxIterations;
    
    /**
     * The random seed
     */
    private long randomSeed;
    
    /**
     * A consumer for the iteration progress
     */
    private DoubleConsumer progressConsumer;
    
    /**
     * A consumer for the log messages
     */
    private Consumer<? super String> messageConsumer;
    
    /**
     * Creates a new instance.<br>
     * <br>
     * The {@link #setOutputDims(int) output dimensions} will be 2.<br>
     * The {@link #setMaxIterations(int) maximum iterations} will be 1000.<br>
     * The {@link #setPerplexity(double) perplexity} will be 20.0.<br>
     * The {@link #setTheta(double) theta} will be 0.5.<br>
     * The {@link #setInitialDims(int) initial dimensions} will be 55.<br>
     * The {@link #setRandomSeed(int) random seed} will be 0.<br>
     * <br>
     * Messages will be written to the standard output. To disable the output, 
     * call {@link #setMessageConsumer(Consumer) setMessageConsumer(null)}.
     */
    public Tsne()
    {
        this.outputDims = 2;
        this.initialDims = 55;
        this.perplexity = 20.0;
        this.theta = 0.5;
        this.maxIterations = 1000;
        this.randomSeed = 0;
        
        this.messageConsumer = message -> System.out.println(message);
    }
    
    /**
     * Set the output dimensions. This should (or will) be 2 in most cases.
     * 
     * @param outputDims The output dimensions
     * @throws IllegalArgumentException If the given number is not positive
     */
    public void setOutputDims(int outputDims)
    {
        if (outputDims <= 0)
        {
            throw new IllegalArgumentException(
                "The outputDims must be positive, but is " + outputDims);
        }
        this.outputDims = outputDims;
    }
    
    /**
     * Returns the number of output dimensions
     * 
     * @return The number of output dimensions
     */
    public int getOutputDims()
    {
        return outputDims;
    }
    
    /**
     * Set the maximum number of iterations that will be performed.<br>
     * <br>
     * Note that due to some learning parameters that are used in the
     * underlying implementation (and that are not accessible from the
     * outside), the number should be ~"considerably larger" than 250.
     * 
     * @param maxIterations The maximum number of iterations
     * @throws IllegalArgumentException If the given number is not positive
     */
    public void setMaxIterations(int maxIterations)
    {
        if (maxIterations <= 0)
        {
            throw new IllegalArgumentException(
                "The maxIterations must be positive, but is " + maxIterations);
        }
        this.maxIterations = maxIterations;
    }
    
    /**
     * Returns the maximum number of iterations that will be performed.
     * See {@link #setMaxIterations(int)} for details.
     * 
     * @return The maximum number of iterations.
     */
    public int getMaxIterations()
    {
        return maxIterations;
    }
    
    /**
     * Set the number of initial dimensions. <br>
     * <br>
     * This is the number of dimensions that the input data should be reduced 
     * to (using a PCA). If this value is not positive, or if it is larger
     * than the dimensionality of the actual input data, then no PCA will 
     * be performed and the dimensionality of the input data will be used.  
     * 
     * @param initialDims The initial dimensions
     */
    public void setInitialDims(int initialDims)
    {
        this.initialDims = initialDims;
    }
    
    /**
     * Returns the initial dimensions. See {@link #setInitialDims(int)} for
     * details.
     * 
     * @return The initial dimensions.
     */
    public int getInitialDims()
    {
        return initialDims;
    }

    /**
     * Set the perplexity. <br>
     * <br>
     * According to the t-SNE FAQ, this can be seen as a measure for the 
     * number of effective nearest neighbors. It should be larger for 
     * larger/denser data sets. The typical value range is between 
     * 5.0 and 50.0.<br>
     * <br>
     * Note that <code>numberOfInputPoints - 1 &gt; 3 * perplexity</code>
     * must hold, otherwise, the underlying implementation will throw
     * an exception.
     * 
     * @param perplexity The perplexity
     */
    public void setPerplexity(double perplexity)
    {
        this.perplexity = perplexity;
    }
    
    /**
     * Returns the perplexity. See {@link #setPerplexity(double)} for details.
     * 
     * @return The perplexity.
     */
    public double getPerplexity()
    {
        return perplexity;
    }
    
    /**
     * Set this theta parameter. Don't do this if you don't know what it does.
     * 
     * @param theta The theta parameter.
     * @throws IllegalArgumentException If the given number is not positive
     */
    public void setTheta(double theta)
    {
        if (theta <= 0.0)
        {
            throw new IllegalArgumentException(
                "The theta must be positive, but is " + theta);
        }
        this.theta = theta;
    }
    
    /**
     * Returns the theta parameter. See {@link #setTheta(double)} for details.
     * 
     * @return The theta parameter
     */
    public double getTheta()
    {
        return theta;
    }

    /**
     * Set the random seed. By default, this will be <code>0</code>, causing
     * the same results to be generated each time. In order to generate
     * different results, the random seed may be set to different values.
     * 
     * @param randomSeed The random seed
     */
    public void setRandomSeed(long randomSeed)
    {
        this.randomSeed = randomSeed;
    }
    
    /**
     * Returns the random seed
     * 
     * @return The random seed
     */
    public long getRandomSeed()
    {
        return randomSeed;
    }
    
    /**
     * Set an optional consumer for the progress of the iterations. This 
     * will receive a value in [0.0, 1.0] indicating the progress of the
     * iterations.
     * 
     * @param progressConsumer The progress consumer
     */
    public void setProgressConsumer(DoubleConsumer progressConsumer)
    {
        this.progressConsumer = progressConsumer;
    }
    
    /**
     * Returns the progress consumer
     * 
     * @return The progress consumer
     */
    public DoubleConsumer getProgressConsumer()
    {
        return progressConsumer;
    }
    
    /**
     * Set an optional consumer for the messages that are emitted by the 
     * underlying implementation.
     * 
     * @param messageConsumer The message consumer
     */
    public void setMessageConsumer(Consumer<? super String> messageConsumer)
    {
        this.messageConsumer = messageConsumer;
    }
    
    /**
     * Returns the consumer for the messages
     * 
     * @return The message consumer
     */
    public Consumer<? super String> getMessageConsumer()
    {
        return messageConsumer;
    }

    /**
     * Performs the t-SNE on the given input data, and returns the result.
     * 
     * @param input The input data
     * @return The result, or <code>null</code> if the calling thread was
     * interrupted
     */
    public double[][] run(double input[][])
    {
        boolean usePca = false;
        int actualDims = input[0].length;
        int pcaDims = actualDims; 
        if (initialDims > 0 && initialDims < actualDims)
        {
            usePca = true;
            pcaDims = initialDims;
        }
        
        boolean silent = true;
        boolean printErrors = false;
        
        if (progressConsumer != null || messageConsumer != null)
        {
            silent = false;
            printErrors = true;
        }
        
        TSneConfiguration config =
            new TSneConfig(input, outputDims, pcaDims, perplexity,
                maxIterations, usePca, theta, silent, printErrors);
        double[][] Y = BHTSne.tsne(config, randomSeed,
            progressConsumer, messageConsumer);
        return Y;
    }

}
