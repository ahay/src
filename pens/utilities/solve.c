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
 *  source file:   ./filters/utilities/solve.c
 *
 * Joe Dellinger (SEP), June 11 1987
 *	Inserted this sample edit history entry.
 *	Please log any further modifications made to this file:
 */

int solve (int pnot, int p1, int q1, int p2, int q2)
/*< compute intersection - floating point version >*/
{
double          invslope;
int    qnot;
float  tempf;
    if (pnot == p1)
	return (q1);
    if (pnot == p2)
	return (q2);
    if (q1 == q2)
	return (q1);
    invslope = (q1 - q2) / ((double) (p1 - p2));
    tempf = (pnot - p1) * invslope + (double) q1;
    qnot = (tempf > 0) ? tempf + 0.5 : tempf - 0.5;
    return (qnot);
}
