/*
 * See BHTSne.java for copyright and license information!  
 */
package de.javagl.tsne;

import java.util.Comparator;

class DistanceComparator implements Comparator<DataPoint>
{

    private DataPoint refItem;

    private Distance dist;

    DistanceComparator(DataPoint refItem, Distance dist)
    {
        this.refItem = refItem;
        this.dist = dist;
    }

    @Override
    public int compare(DataPoint o1, DataPoint o2)
    {
        return dist.distance(o1, refItem) < dist.distance(o2, refItem) ? -1
            : (dist.distance(o1, refItem) > dist.distance(o2, refItem) ? 1 : 0);
    }
}