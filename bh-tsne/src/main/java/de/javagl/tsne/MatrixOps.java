/*
 * See BHTSne.java for copyright and license information!  
 */
package de.javagl.tsne;

class MatrixOps
{
    static DoubleArray extractRowViewFromFlatMatrix(double[] flatMatrix,
        int rowIdx, int dimension)
    {
        DoubleArray row = new DoubleArray()
        {
            @Override
            public double get(int index)
            {
                return flatMatrix[rowIdx * dimension + index];
            }

            @Override
            public int length()
            {
                return dimension;
            }

        };
        return row;
    }

}
