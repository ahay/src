/*
 * Copyright 1987 the Board of Trustees of the Leland Stanford Junior
 * University. Official permission to use this software is included in
 * the documentation. It authorizes you to use this file for any
 * non-commercial purpose, provided that this copyright notice is not
 * removed and that any modifications made to this file are commented
 * and dated in the style of my example below.
 */

/*
 *
 *  source file:   ./filters/utilities/sort.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

void sort (int *vec, int n)
/*< sort the elements of vec into ascending order >*/
{
register int   *above, *below, *last;
register int    temp;
    last = vec + (n - 1);
    for (above = vec; above < last; above++)
    {
	for (below = above + 1; below <= last; below++)
	{
	    if (*above > *below)
	    {
		temp = *above;
		*above = *below;
		*below = temp;
	    }
	}
    }
}
