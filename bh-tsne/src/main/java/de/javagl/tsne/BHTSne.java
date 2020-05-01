/*
 * www.javagl.de - bh-tsne 
 *
 * Original implementation: Copyright (c) Leif Jonsson 2016 and 
 * Copyright (c) 2014, Laurens van der Maaten
 * 
 * Modified implementation: Copyright (c) 2020 Marco Hutter 
 * http://www.javagl.de
 * 
 * All files in this package have originally been part of the 
 * Barnes-Hut-TSNE implementation by Leif Jonsson 
 * https://github.com/lejon/T-SNE-Java 
 * (see the copyright and license information below!) 
 * 
 * The result is published under any license that is compatible with
 * the BSD license that is quoted below: 
 */

/*
 *
 * This Java port of Barnes Hut t-SNE is Copyright (c) Leif Jonsson 2016 and 
 * Copyright (c) 2014, Laurens van der Maaten (Delft University of Technology)
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. All advertising materials mentioning features or use of this software
 *    must display the following acknowledgement:
 *    This product includes software developed by the Delft University of Technology.
 * 4. Neither the name of the Delft University of Technology nor the names of
 *    its contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ''AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO
 * EVENT SHALL LAURENS VAN DER MAATEN BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
 * IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY
 * OF SUCH DAMAGE.
 *
 */
package de.javagl.tsne;

import static java.lang.Math.exp;
import static java.lang.Math.log;

import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Random;
import java.util.function.Consumer;
import java.util.function.DoubleConsumer;
import java.util.logging.Logger;
import java.util.stream.DoubleStream;
import java.util.stream.IntStream;

class BHTSne
{
    /**
     * The logger used in this class
     */
    private static final Logger logger =
        Logger.getLogger(BHTSne.class.getName());

    // (bh-tsne) : A trivial parallelization. This is not related
    // to the parallel implementation of the original library.
    private static final boolean PARALLEL = true;

    private static final Distance distance = new EuclideanDistance();

    // (bh-tsne) This is the only remaining non-private method, used
    // to launch the actual computation
    static double[][] tsne(TSneConfiguration config, long randomSeed,
        DoubleConsumer progressConsumer,
        Consumer<? super String> messageConsumer)
    {
        return run(config, randomSeed, progressConsumer, messageConsumer);
    }

    private static double[] flatten(double[][] x)
    {
        int noCols = x[0].length;
        double[] flat = new double[x.length * x[0].length];
        for (int i = 0; i < x.length; i++)
        {
            for (int j = 0; j < x[i].length; j++)
            {
                flat[i * noCols + j] = x[i][j];
            }
        }
        return flat;
    }

    private static double[][] expand(double[] x, int N, int D)
    {
        double[][] expanded = new double[N][D];
        for (int row = 0; row < N; row++)
        {
            for (int col = 0; col < D; col++)
            {
                expanded[row][col] = x[row * D + col];
            }
        }
        return expanded;
    }

    private static double sign_tsne(double x)
    {
        return (x == .0 ? .0 : (x < .0 ? -1.0 : 1.0));
    }

    // Perform t-SNE
    private static double[][] run(TSneConfiguration parameterObject,
        long randomSeed, DoubleConsumer progressConsumer,
        Consumer<? super String> messageConsumer)
    {
        int D = parameterObject.getXStartDim();
        double[][] Xin = parameterObject.getXin();
        boolean exactError = (parameterObject.getTheta() == .0);
        if (exactError)
        {
            throw new IllegalArgumentException(
                "The Barnes Hut implementation does not support "
                + "exact inference yet (theta==0.0), if you want exact "
                + "t-SNE please use one of the standard t-SNE "
                + "implementations (FastTSne for instance)");
        }

        if (parameterObject.usePca() && D > parameterObject.getInitialDims()
            && parameterObject.getInitialDims() > 0)
        {
            PrincipalComponentAnalysis pca = new PrincipalComponentAnalysis();
            Xin = pca.pca(Xin, parameterObject.getInitialDims());
            D = parameterObject.getInitialDims();
            printLog(messageConsumer, "X:Shape after PCA is = " + Xin.length
                + " x " + Xin[0].length);
        }

        double[] X = flatten(Xin);
        int N = parameterObject.getNrRows();
        int no_dims = parameterObject.getOutputDims();

        double[] Y = new double[N * no_dims];
        printLog(messageConsumer, "X:Shape is = " + N + " x " + D);
        
        // Determine whether we are using an exact algorithm
        double perplexity = parameterObject.getPerplexity();
        if (N - 1 < 3 * perplexity)
        {
            throw new IllegalArgumentException(
                "Perplexity too large for the number of data points! " 
                + "The number of data points should be more than three "
                + "times the perplexity. There are " + N 
                + " data points, and the perplexity is " + perplexity + ".");
        }
        printLog(messageConsumer,
            "Using no_dims = %d, perplexity = %f, and theta = %f", no_dims,
            perplexity, parameterObject.getTheta());

        Random random = new Random(randomSeed);
        
        // Set learning parameters
        double total_time = 0;
        int stop_lying_iter = 250, mom_switch_iter = 250;
        double momentum = .5, final_momentum = .8;
        double eta = 200.0;

        // Allocate some memory
        double[] dY = new double[N * no_dims];
        double[] uY = new double[N * no_dims];
        double[] gains = new double[N * no_dims];
        for (int i = 0; i < N * no_dims; i++)
            gains[i] = 1.0;

        // Normalize input data (to prevent numerical problems)
        printLog(messageConsumer, "Computing input similarities...");
        long start = System.currentTimeMillis();
        // zeroMean(X, N, D);
        double max_X = .0;
        for (int i = 0; i < N * D; i++)
        {
            if (X[i] > max_X)
                max_X = X[i];
        }

        for (int i = 0; i < N * D; i++)
            X[i] /= max_X;

        double[] P = null;
        int K = (int) (3 * perplexity);
        int[] row_P = new int[N + 1];
        int[] col_P = new int[N * K];
        double[] val_P = new double[N * K];

        // Compute input similarities for approximate t-SNE
        // Compute asymmetric pairwise input similarities
        computeGaussianPerplexity(X, N, D, row_P, col_P, val_P, perplexity, K,
            random, messageConsumer);
        if (Thread.currentThread().isInterrupted())
        {
            return null;
        }

        // Verified that val_P,col_P,row_P is the same at this point

        // Symmetrize input similarities
        SymResult res = symmetrizeMatrix(row_P, col_P, val_P, N);
        row_P = res.sym_row_P;
        col_P = res.sym_col_P;
        val_P = res.sym_val_P;

        double sum_P = .0;
        for (int i = 0; i < row_P[N]; i++)
            sum_P += val_P[i];
        for (int i = 0; i < row_P[N]; i++)
            val_P[i] /= sum_P;
        long end = System.currentTimeMillis();

        // Lie about the P-values
        for (int i = 0; i < row_P[N]; i++)
            val_P[i] *= 12.0;

        // Initialize solution (randomly)
        for (int i = 0; i < N * no_dims; i++)
            Y[i] = random.nextDouble() * 0.0001;

        // Perform main training loop
        printLog(messageConsumer,
            "Done in %4.2f seconds (sparsity = %f)!\nLearning embedding...",
            (end - start) / 1000.0,
            (double) row_P[N] / ((double) N * (double) N));
        start = System.currentTimeMillis();
        for (int iter = 0; iter < parameterObject.getMaxIter(); iter++)
        {
            if (Thread.currentThread().isInterrupted())
            {
                return null;
            }

            // Compute (approximate) gradient
            computeGradient(P, row_P, col_P, val_P, Y, N, no_dims, dY,
                parameterObject.getTheta());

            updateGradient(N, no_dims, Y, momentum, eta, dY, uY, gains);

            // Make solution zero-mean
            zeroMean(Y, N, no_dims);
            
            // Stop lying about the P-values after a while, and switch momentum
            if (iter == stop_lying_iter)
            {
                for (int i = 0; i < row_P[N]; i++)
                    val_P[i] /= 12.0;
            }
            if (iter == mom_switch_iter)
                momentum = final_momentum;

            if (progressConsumer != null)
            {
                double progress =
                    (double) iter / (parameterObject.getMaxIter() - 1);
                progressConsumer.accept(progress);
            }

            // Print out progress
            if (((iter > 0 && iter % 50 == 0)
                || iter == parameterObject.getMaxIter() - 1)
                && !parameterObject.silent())
            {
                end = System.currentTimeMillis();
                String err_string = "not_calculated";
                if (parameterObject.printError())
                {
                    double C = .0;
                    // doing approximate computation here!
                    C = evaluateError(row_P, col_P, val_P, Y, N, no_dims,
                        parameterObject.getTheta());
                    err_string = "" + C;
                }
                if (iter == 0)
                    printLog(messageConsumer, "Iteration %d: error is %s",
                        iter + 1, err_string);
                else
                {
                    total_time += (end - start) / 1000.0;
                    printLog(messageConsumer,
                        "Iteration %d: error is %s (50 iterations in %4.2f seconds)",
                        iter, err_string, (end - start) / 1000.0);
                }
                start = System.currentTimeMillis();
            }
        }
        end = System.currentTimeMillis();
        total_time += (end - start) / 1000.0;

        printLog(messageConsumer, "Fitting performed in %4.2f seconds.",
            total_time);
        return expand(Y, N, no_dims);
    }

    private static void updateGradient(int N, int no_dims, double[] Y,
        double momentum, double eta, double[] dY, double[] uY, double[] gains)
    {
        for (int i = 0; i < N * no_dims; i++)
        {
            // Update gains
            gains[i] = (sign_tsne(dY[i]) != sign_tsne(uY[i])) ? (gains[i] + .2)
                : (gains[i] * .8);
            if (gains[i] < .01)
                gains[i] = .01;

            // Perform gradient update (with momentum and gains)
            Y[i] = Y[i] + uY[i];
            uY[i] = momentum * uY[i] - eta * gains[i] * dY[i];
        }
    }

    // Compute gradient of the t-SNE cost function (using Barnes-Hut algorithm)
    private static void computeGradient(double[] P, int[] inp_row_P,
        int[] inp_col_P, double[] inp_val_P, double[] Y, int N, int D,
        double[] dC, double theta)
    {
        // Construct space-partitioning tree on current map
        SPTree tree = new SPTree(D, Y, N);
        
        double totalSum_Q = 0.0;
        double[] sum_Q = new double[N];
        double[] pos_f = new double[N * D];
        double[][] neg_f = new double[N][D];
        double[][] buff = new double[N][D];
        tree.computeEdgeForces(inp_row_P, inp_col_P, inp_val_P, N, pos_f);

        if (!PARALLEL)
        {
            // Compute all terms required for t-SNE gradient
            for (int n = 0; n < N; n++)
                tree.computeNonEdgeForces(n, theta, neg_f[n], buff[n], sum_Q);
        }
        else
        {
            // Compute all terms required for t-SNE gradient
            IntStream.range(0, N).parallel().forEach(i -> {
                tree.computeNonEdgeForces(i, theta, neg_f[i], buff[i], sum_Q);
            });
        }
        totalSum_Q = DoubleStream.of(sum_Q).sum();

        // Compute final t-SNE gradient
        for (int n = 0; n < N; n++)
        {
            for (int d = 0; d < D; d++)
            {
                dC[n * D + d] = pos_f[n * D + d] - (neg_f[n][d] / totalSum_Q);
            }
        }
    }

    // Evaluate t-SNE cost function (approximately)
    private static double evaluateError(int[] row_P, int[] col_P,
        double[] val_P, double[] Y, int N, int D, double theta)
    {
        // Get estimate of normalization term
        SPTree tree = new SPTree(D, Y, N);
        double[] buff = new double[D];
        double[][] buffs = new double[N][D];
        double totalSum_Q = 0.0;
        double[] sum_Q = new double[N];

        if (!PARALLEL)
        {
            // Compute all terms required for t-SNE gradient
            for (int n = 0; n < N; n++)
                tree.computeNonEdgeForces(n, theta, buff, buffs[n], sum_Q);
        }
        else
        {
            // Compute all terms required for t-SNE gradient
            IntStream.range(0, N).parallel().forEach(i -> {
                tree.computeNonEdgeForces(i, theta, buff, buffs[i], sum_Q);
            });
        }
        totalSum_Q = DoubleStream.of(sum_Q).sum();

        // Loop over all edges to compute t-SNE error
        int ind1, ind2;
        double C = .0, Q;
        for (int n = 0; n < N; n++)
        {
            ind1 = n * D;
            for (int i = row_P[n]; i < row_P[n + 1]; i++)
            {
                Q = .0;
                ind2 = col_P[i] * D;
                for (int d = 0; d < D; d++)
                    buff[d] = Y[ind1 + d];
                for (int d = 0; d < D; d++)
                    buff[d] -= Y[ind2 + d];
                for (int d = 0; d < D; d++)
                    Q += buff[d] * buff[d];
                Q = (1.0 / (1.0 + Q)) / totalSum_Q;
                C += val_P[i] * log(
                    (val_P[i] + Double.MIN_VALUE) / (Q + Double.MIN_VALUE));
            }
        }

        return C;
    }

    // Compute input similarities with a fixed perplexity using ball trees
    private static void computeGaussianPerplexity(double[] X, int N, int D,
        int[] _row_P, int[] _col_P, double[] _val_P, double perplexity, int K,
        Random random, Consumer<? super String> messageConsumer)
    {
        if (perplexity > K)
        {
            logger.warning("Perplexity should be lower than K!");
            // Come on, tell the user what's wrong, it's not so hard:
            logger.warning("perplexity=" + perplexity + ", K=" + K);
        }

        // Allocate the memory we need
        int[] row_P = _row_P;
        int[] col_P = _col_P;
        double[] val_P = _val_P;
        double[] cur_P = new double[N - 1];

        row_P[0] = 0;
        for (int n = 0; n < N; n++)
            row_P[n + 1] = row_P[n] + K;

        // Build ball tree on data set
        VpTree<DataPoint> tree = new VpTree<DataPoint>(distance, random);
        final DataPoint[] obj_X = new DataPoint[N];
        for (int n = 0; n < N; n++)
        {
            DoubleArray row = MatrixOps.extractRowViewFromFlatMatrix(X, n, D);
            obj_X[n] = new DataPoint(D, n, row);
        }
        tree.create(obj_X);

        // Loop over all points to find nearest neighbors
        printLog(messageConsumer, "Building tree...");
        long beforeNs = System.nanoTime();
        List<DataPoint> indices = new ArrayList<>();
        List<Double> distances = new ArrayList<>();
        for (int n = 0; n < N; n++)
        {
            if (Thread.currentThread().isInterrupted())
            {
                return;
            }
            if (n % 10000 == 0)
                printLog(messageConsumer, " - point %d of %d", n, N);

            // Find nearest neighbors
            indices.clear();
            distances.clear();
            tree.search(obj_X[n], K + 1, indices, distances);
            
            // Initialize some variables for binary search
            boolean found = false;
            double beta = 1.0;
            double min_beta = -Double.MAX_VALUE;
            double max_beta = Double.MAX_VALUE;
            double tol = 1e-5;

            // Iterate until we found a good perplexity
            int iter = 0;
            double sum_P = 0.;
            while (!found && iter < 200)
            {

                // Compute Gaussian kernel row and entropy of current row
                sum_P = Double.MIN_VALUE;
                double H = .0;
                for (int m = 0; m < K; m++)
                {
                    cur_P[m] = exp(-beta * distances.get(m + 1));
                    sum_P += cur_P[m];
                    H += beta * (distances.get(m + 1) * cur_P[m]);
                }
                H = (H / sum_P) + log(sum_P);

                // Evaluate whether the entropy is within the tolerance level
                double Hdiff = H - log(perplexity);
                if (Hdiff < tol && -Hdiff < tol)
                {
                    found = true;
                }
                else
                {
                    if (Hdiff > 0)
                    {
                        min_beta = beta;
                        if (max_beta == Double.MAX_VALUE
                            || max_beta == -Double.MAX_VALUE)
                            beta *= 2.0;
                        else
                            beta = (beta + max_beta) / 2.0;
                    }
                    else
                    {
                        max_beta = beta;
                        if (min_beta == -Double.MAX_VALUE
                            || min_beta == Double.MAX_VALUE)
                            beta /= 2.0;
                        else
                            beta = (beta + min_beta) / 2.0;
                    }
                }

                // Update iteration counter. Oh. Really?
                iter++;
            }

            // Row-normalize current row of P and store in matrix
            for (int m = 0; m < K; m++)
            {
                cur_P[m] /= sum_P;
                col_P[row_P[n] + m] = indices.get(m + 1).index();
                val_P[row_P[n] + m] = cur_P[m];
            }
        }
        long afterNs = System.nanoTime();
        printLog(messageConsumer, "Building tree took %f ms",
            (afterNs - beforeNs) / 1e6);
    }

    private static class SymResult
    {

        int[] sym_row_P;

        int[] sym_col_P;

        double[] sym_val_P;

        private SymResult(int[] sym_row_P, int[] sym_col_P, double[] sym_val_P)
        {
            super();
            this.sym_row_P = sym_row_P;
            this.sym_col_P = sym_col_P;
            this.sym_val_P = sym_val_P;
        }
    }

    private static SymResult symmetrizeMatrix(int[] _row_P, int[] _col_P,
        double[] _val_P, int N)
    {
        // Get sparse matrix
        int[] row_P = _row_P;
        int[] col_P = _col_P;
        double[] val_P = _val_P;

        // Count number of elements and row counts of symmetric matrix
        int[] row_counts = new int[N];
        for (int n = 0; n < N; n++)
        {
            for (int i = row_P[n]; i < row_P[n + 1]; i++)
            {

                // Check whether element (col_P[i], n) is present
                boolean present = false;
                for (int m = row_P[col_P[i]]; m < row_P[col_P[i] + 1]; m++)
                {
                    if (col_P[m] == n)
                        present = true;
                }
                if (present)
                    row_counts[n]++;
                else
                {
                    row_counts[n]++;
                    row_counts[col_P[i]]++;
                }
            }
        }
        int no_elem = 0;
        for (int n = 0; n < N; n++)
            no_elem += row_counts[n];

        // Allocate memory for symmetrized matrix
        int[] sym_row_P = new int[N + 1];
        int[] sym_col_P = new int[no_elem];
        double[] sym_val_P = new double[no_elem];

        // Construct new row indices for symmetric matrix
        sym_row_P[0] = 0;
        for (int n = 0; n < N; n++)
            sym_row_P[n + 1] = sym_row_P[n] + row_counts[n];

        // Fill the result matrix
        int[] offset = new int[N];
        for (int n = 0; n < N; n++)
        {
            for (int i = row_P[n]; i < row_P[n + 1]; i++)
            { // considering element(n, col_P[i])

                // Check whether element (col_P[i], n) is present
                boolean present = false;
                for (int m = row_P[col_P[i]]; m < row_P[col_P[i] + 1]; m++)
                {
                    if (col_P[m] == n)
                    {
                        present = true;
                        if (n <= col_P[i])
                        { // make sure we do not add elements
                          // twice
                            sym_col_P[sym_row_P[n] + offset[n]] = col_P[i];
                            sym_col_P[sym_row_P[col_P[i]] + offset[col_P[i]]] =
                                n;
                            sym_val_P[sym_row_P[n] + offset[n]] =
                                val_P[i] + val_P[m];
                            sym_val_P[sym_row_P[col_P[i]] + offset[col_P[i]]] =
                                val_P[i] + val_P[m];
                        }
                    }
                }

                // If (col_P[i], n) is not present, there is no addition
                // involved
                if (!present)
                {
                    sym_col_P[sym_row_P[n] + offset[n]] = col_P[i];
                    sym_col_P[sym_row_P[col_P[i]] + offset[col_P[i]]] = n;
                    sym_val_P[sym_row_P[n] + offset[n]] = val_P[i];
                    sym_val_P[sym_row_P[col_P[i]] + offset[col_P[i]]] =
                        val_P[i];
                }

                // Update offsets
                if (!present || (present && n <= col_P[i]))
                {
                    offset[n]++;
                    if (col_P[i] != n)
                        offset[col_P[i]]++;
                }
            }
        }

        // Divide the result by two
        for (int i = 0; i < no_elem; i++)
            sym_val_P[i] /= 2.0;

        return new SymResult(sym_row_P, sym_col_P, sym_val_P);
    }

    // Makes data zero-mean
    private static void zeroMean(double[] X, int N, int D)
    {

        // Compute data mean
        double[] mean = new double[D];
        for (int n = 0; n < N; n++)
        {
            for (int d = 0; d < D; d++)
            {
                mean[d] += X[n * D + d];
            }
        }
        for (int d = 0; d < D; d++)
        {
            mean[d] /= (double) N;
        }

        // Subtract data mean
        for (int n = 0; n < N; n++)
        {
            for (int d = 0; d < D; d++)
            {
                X[n * D + d] -= mean[d];
            }
        }
    }

    private static void printLog(Consumer<? super String> messageConsumer,
        String format, Object... args)
    {
        if (messageConsumer != null)
        {
            String message = String.format(Locale.ENGLISH, format, args);
            messageConsumer.accept(message);
        }
    }
}
