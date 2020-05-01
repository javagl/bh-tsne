/*
 * See BHTSne.java for copyright and license information!  
 */
package de.javagl.tsne;

class TSneConfig implements TSneConfiguration
{
    private double[][] xin;
    private int outputDims;

    private int initial_dims;
    private double perplexity;
    private int max_iter;
    private boolean use_pca;
    private double theta;
    private boolean silent;
    private boolean print_error;

    TSneConfig(double[][] xin, int outputDims, int initial_dims,
        double perplexity, int max_iter, boolean use_pca, double theta,
        boolean silent, boolean print_error)
    {
        this.xin = xin;
        this.outputDims = outputDims;
        this.initial_dims = initial_dims;
        this.perplexity = perplexity;
        this.max_iter = max_iter;
        this.use_pca = use_pca;
        this.theta = theta;
        this.silent = silent;
        this.print_error = print_error;
    }

    @Override
    public double[][] getXin()
    {
        return xin;
    }

    @Override
    public void setXin(double[][] xin)
    {
        this.xin = xin;
    }

    @Override
    public int getOutputDims()
    {
        return outputDims;
    }

    @Override
    public void setOutputDims(int n)
    {
        this.outputDims = n;
    }

    @Override
    public int getInitialDims()
    {
        return initial_dims;
    }

    @Override
    public void setInitialDims(int initial_dims)
    {
        this.initial_dims = initial_dims;
    }

    @Override
    public double getPerplexity()
    {
        return perplexity;
    }

    @Override
    public void setPerplexity(double perplexity)
    {
        this.perplexity = perplexity;
    }

    @Override
    public int getMaxIter()
    {
        return max_iter;
    }

    @Override
    public void setMaxIter(int max_iter)
    {
        this.max_iter = max_iter;
    }

    @Override
    public boolean usePca()
    {
        return use_pca;
    }

    @Override
    public void setUsePca(boolean use_pca)
    {
        this.use_pca = use_pca;
    }

    @Override
    public double getTheta()
    {
        return theta;
    }

    @Override
    public void setTheta(double theta)
    {
        this.theta = theta;
    }

    @Override
    public boolean silent()
    {
        return silent;
    }

    @Override
    public void setSilent(boolean silent)
    {
        this.silent = silent;
    }

    @Override
    public boolean printError()
    {
        return print_error;
    }

    @Override
    public void setPrintError(boolean print_error)
    {
        this.print_error = print_error;
    }

    @Override
    public int getXStartDim()
    {
        return xin[0].length;
    }

    @Override
    public int getNrRows()
    {
        return xin.length;
    }

}