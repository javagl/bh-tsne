/*
 * See BHTSne.java for copyright and license information!  
 */
package de.javagl.tsne;

import static java.lang.Math.sqrt;

class EuclideanDistance implements Distance
{

    public EuclideanDistance()
    {
    }

    @Override
    public double distance(DataPoint d1, DataPoint d2)
    {
        double dd = .0;
        DoubleArray x1 = d1._x;
        DoubleArray x2 = d2._x;
        double diff;
        for (int d = 0; d < d1._D; d++)
        {
            diff = (x1.get(d) - x2.get(d));
            dd += diff * diff;
        }
        return sqrt(dd);

    }

}
