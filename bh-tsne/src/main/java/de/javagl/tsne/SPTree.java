/*
 * See BHTSne.java for copyright and license information!  
 */
package de.javagl.tsne;

import static java.lang.Math.max;
import static java.lang.Math.sqrt;

class SPTree
{

    // Fixed constants

    private final static int QT_NODE_CAPACITY = 1;

    private int dimension;

    private boolean is_leaf;

    private int size;

    private int cum_size;

    // Axis-aligned bounding box stored as a center with half-dimensions to
    // represent the boundaries of this quad tree

    private Cell boundary;

    // Indices in this space-partitioning tree node, corresponding
    // center-of-mass, and list of all children

    private double[] data;

    private double[] center_of_mass;

    private int[] index = new int[QT_NODE_CAPACITY];

    // Children

    private SPTree[] children;

    private int no_children;

    private double[] buff;
    
    SPTree(int D, double[] inp_data, int N)
    {
        // Compute mean, width, and height of current map (boundaries of SPTree)
        int nD = 0;
        double[] mean_Y = new double[D];
        double[] min_Y = new double[D];
        double[] max_Y = new double[D];
        for (int d = 0; d < D; d++)
        {
            min_Y[d] = Double.POSITIVE_INFINITY;
            max_Y[d] = Double.NEGATIVE_INFINITY;
        }
        for (int n = 0; n < N; n++)
        {
            for (int d = 0; d < D; d++)
            {
                mean_Y[d] += inp_data[n * D + d];
                if (inp_data[nD + d] < min_Y[d])
                    min_Y[d] = inp_data[nD + d];
                if (inp_data[nD + d] > max_Y[d])
                    max_Y[d] = inp_data[nD + d];
            }
            nD += D;
        }
        for (int d = 0; d < D; d++)
            mean_Y[d] /= (double) N;

        // Construct SPTree
        double[] width = new double[D];
        for (int d = 0; d < D; d++)
            width[d] = max(max_Y[d] - mean_Y[d], mean_Y[d] - min_Y[d]) + 1e-5;
        init(null, D, inp_data, mean_Y, width);
        fill(N);
    }

    // Main initialization function
    private void init(SPTree inp_parent, int D, double[] inp_data,
        double[] inp_corner, double[] inp_width)
    {
        dimension = D;
        no_children = 2;
        for (int d = 1; d < D; d++)
            no_children *= 2;
        data = inp_data;
        is_leaf = true;
        size = 0;
        cum_size = 0;

        center_of_mass = new double[D];
        boundary = new Cell(dimension, inp_corner, inp_width);

        children = getTreeArray(no_children);
        for (int i = 0; i < no_children; i++)
            children[i] = null;
        
        buff = new double[dimension];
    }

    // Constructor for SPTree with particular size and parent (do not fill tree)
    private SPTree(SPTree inp_parent, int D, double[] inp_data,
        double[] inp_corner, double[] inp_width)
    {
        init(inp_parent, D, inp_data, inp_corner, inp_width);
    }

    private SPTree[] getTreeArray(int no_children)
    {
        return new SPTree[no_children];
    }

    // Insert a point into the SPTree
    private boolean insert(int new_index)
    {
        // Ignore objects which do not belong in this quad tree
        DoubleArray point =
            MatrixOps.extractRowViewFromFlatMatrix(data, new_index, dimension);

        if (!boundary.containsPoint(point))
            return false;

        // Online update of cumulative size and center-of-mass
        cum_size++;
        double mult1 = (double) (cum_size - 1) / (double) cum_size;
        double mult2 = 1.0 / (double) cum_size;
        for (int d = 0; d < dimension; d++)
        {
            center_of_mass[d] *= mult1;
            center_of_mass[d] += mult2 * point.get(d);
        }

        // If there is space in this quad tree and it is a leaf, add the object
        // here
        if (is_leaf && size < QT_NODE_CAPACITY)
        {
            index[size] = new_index;
            size++;
            return true;
        }

        // Don't add duplicates for now (this is not very nice)
        boolean any_duplicate = false;
        for (int n = 0; n < size; n++)
        {
            boolean duplicate = true;
            for (int d = 0; d < dimension; d++)
            {
                if (point.get(d) != data[index[n] * dimension + d])
                {
                    duplicate = false;
                    break;
                }
            }
            any_duplicate = any_duplicate || duplicate;
        }
        if (any_duplicate)
            return true;

        // Otherwise, we need to subdivide the current cell
        if (is_leaf)
            subdivide();

        // Find out where the point can be inserted
        for (int i = 0; i < no_children; i++)
        {
            if (children[i].insert(new_index))
                return true;
        }

        // Otherwise, the point cannot be inserted (this should never happen)
        assert false;
        return false;
    }

    // Create four children which fully divide this cell into four quads of
    // equal area
    private void subdivide()
    {

        // Create new children
        for (int i = 0; i < no_children; i++)
        {
            double[] new_corner = new double[dimension];
            double[] new_width = new double[dimension];
            int div = 1;
            for (int d = 0; d < dimension; d++)
            {
                new_width[d] = .5 * boundary.getWidth(d);
                if ((i / div) % 2 == 1)
                    new_corner[d] =
                        boundary.getCorner(d) - .5 * boundary.getWidth(d);
                else
                    new_corner[d] =
                        boundary.getCorner(d) + .5 * boundary.getWidth(d);
                div *= 2;
            }
            children[i] = getNewTree(this, new_corner, new_width);
        }

        // Move existing points to correct children
        for (int i = 0; i < size; i++)
        {
            boolean success = false;
            for (int j = 0; j < no_children; j++)
            {
                if (!success)
                    success = children[j].insert(index[i]);
            }
            index[i] = -1;
        }

        // Empty parent node
        size = 0;
        is_leaf = false;
    }

    private SPTree getNewTree(SPTree root, double[] new_corner,
        double[] new_width)
    {
        return new SPTree(root, dimension, data, new_corner, new_width);
    }

    // Build SPTree on dataset
    private void fill(int N)
    {
        for (int i = 0; i < N; i++)
            insert(i);
    }

    
    // Compute non-edge forces using Barnes-Hut algorithm
    double computeNonEdgeForces(int point_index, double theta, double[] neg_f, 
        double buff[], double[] sum_Q)
    {
        //double[] buff = new double[dimension];

        // Make sure that we spend no time on empty nodes or self-interactions
        if (cum_size == 0 || (is_leaf && size == 1 && index[0] == point_index))
            return 0.0;

        // Compute distance between point and center-of-mass
        double D = .0;
        int ind = point_index * dimension;
        // Check whether we can use this node as a "summary"
        double max_width = 0.0;
        double cur_width;
        for (int d = 0; d < dimension; d++)
        {
            buff[d] = data[ind + d] - center_of_mass[d];
            D += buff[d] * buff[d];
            cur_width = boundary.getWidth(d);
            max_width = (max_width > cur_width) ? max_width : cur_width;
        }

        if (is_leaf || max_width / sqrt(D) < theta)
        {
            // Compute and add t-SNE force between point and current node
            D = 1.0 / (1.0 + D);
            double mult = cum_size * D;
            sum_Q[point_index] += mult;
            mult *= D;
            for (int d = 0; d < dimension; d++)
                neg_f[d] += mult * buff[d];
        }
        else
        {

            // Recursively apply Barnes-Hut to children
            for (int i = 0; i < no_children; i++)
                children[i].computeNonEdgeForces(point_index, theta, neg_f,
                    buff, sum_Q);
        }
        return sum_Q[point_index];
    }

    // Computes edge forces
    void computeEdgeForces(int[] row_P, int[] col_P, double[] val_P, int N,
        double[] pos_f)
    {
        // Loop over all edges in the graph
        //double[] buff = new double[dimension];
        int ind1 = 0;
        int ind2 = 0;
        double D;
        for (int n = 0; n < N; n++)
        {
            for (int i = row_P[n]; i < row_P[n + 1]; i++)
            {

                // Compute pairwise distance and Q-value
                D = 1.0;
                ind2 = col_P[i] * dimension;
                for (int d = 0; d < dimension; d++)
                {
                    buff[d] = data[ind1 + d] - data[ind2 + d];
                    D += buff[d] * buff[d];
                }
                D = val_P[i] / D;

                // Sum positive force
                for (int d = 0; d < dimension; d++)
                    pos_f[ind1 + d] += D * buff[d];
            }
            ind1 += dimension;
        }
    }

    private class Cell
    {

        int dimension;

        double[] corner;

        double[] width;

        // Constructs cell
        private Cell(int inp_dimension, double corner[], double width[])
        {
            dimension = inp_dimension;
            this.corner = corner;
            this.width = width;
        }

        private double getCorner(int d)
        {
            return corner[d];
        }

        private double getWidth(int d)
        {
            return width[d];
        }

        // Checks whether a point lies in a cell
        private boolean containsPoint(DoubleArray point)
        {
            for (int d = 0; d < dimension; d++)
            {
                if (corner[d] - width[d] > point.get(d))
                    return false;
                if (corner[d] + width[d] < point.get(d))
                    return false;
            }
            return true;
        }
        
    }
}
