/*
 * See BHTSne.java for copyright and license information!  
 */
package de.javagl.tsne;

import static java.lang.Math.min;

class DataPoint
{

    private int _ind;
    DoubleArray _x;
    int _D;

    DataPoint()
    {
        _D = 1;
        _ind = -1;
    }

    DataPoint(int D, int ind, DoubleArray x)
    {
        _D = D;
        _ind = ind;
        _x = x;
    }

    @Override
    public String toString()
    {
        String xStr = "";
        for (int i = 0; i < min(20, _x.length()); i++)
        {
            xStr += _x.get(i) + ", ";
        }
        return "DataPoint (index=" + _ind + ", Dim=" + _D + ", point=" + xStr
            + ")";
    }

    int index()
    {
        return _ind;
    }

}

