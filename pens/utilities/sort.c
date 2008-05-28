/*
  Copyright (C) 1987 The Board of Trustees of Stanford University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
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
