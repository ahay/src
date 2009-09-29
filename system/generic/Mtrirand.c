/* Edit points for triangulation by removing similar and randomizing. */
/*
  Copyright (C) 2004 University of Texas at Austin

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

#include <stdlib.h>
#include <rsf.h>

static int RandomInteger(int low, int high);
static void Shuffle(float **xyz, int n);

int main(int argc, char* argv[])
{
    int nd, three, i, j, id;
    float **xyz=NULL, **xyz1=NULL, xi, yi, xj, yj;
    sf_file in=NULL, out=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&three) || 3 != three) 
	sf_error("Need n1=3 in input");
    if (!sf_histint(in,"n2",&nd)) nd=1;

    xyz = sf_floatalloc2(3,nd);
    xyz1 = sf_floatalloc2(3,nd);

    sf_floatread(xyz[0],3*nd,in);

    sf_pqueue_init (3*nd);
    sf_pqueue_start();
    for (i=0; i < nd; i++) {
	sf_pqueue_insert (xyz[i]);
    }

    id = sf_pqueue_extract()-xyz[0];
    xyz1[0][0] = xj = xyz[0][id];
    xyz1[0][1] = yj = xyz[0][id+1];
    xyz1[0][2] = xyz[0][id+2];
    for (j=1,i=1; i < nd; i++) {
	id = sf_pqueue_extract()-xyz[0];
	xi = xyz[0][id];
	yi = xyz[0][id+1];
	if (xi != xj || yi != yj) {
	    xj = xi; yj = xi;
	    xyz1[j][0] = xi;
	    xyz1[j][1] = yi;
	    xyz1[j][2] = xyz[0][id+2];
	    j++;
	}
    }
    sf_pqueue_close ();

    nd = j;
    sf_putint(out,"n2",nd);

    Shuffle (xyz1,nd);

    /* after shuffle, the array is no longer contiguous */

    for (i=0; i < nd; i++) {
	sf_floatwrite(xyz1[i],3,out);
    }

    exit(0);
}


/*
 * Function: RandomInteger
 * -----------------------
 * This function first obtains a random integer in
 * the range [0..RAND_MAX] by applying four steps:
 * (1) Generate a real number between 0 and 1.
 * (2) Scale it to the appropriate range size.
 * (3) Truncate the value to an integer.
 * (4) Translate it to the appropriate starting point.
 */

static int RandomInteger(int low, int high)
{
    int k;
    double d;

    d = (double) rand() / ((double) RAND_MAX + 1);
    k = (int) (d * (high - low + 1));
    return (low + k);
}

/*
 * Implementation notes: Shuffle
 * -----------------------------
 * This function shuffles the array xyz of size m by n.
 * This is done by slightly modifing the selection sort algorithm:
 * instead of choosing the smallest element in each cycle, pick a random
 * element.
 */

static void Shuffle(float **xyz, int n)
{
    int lh, rh;
    float *temp;

    for (lh = 0; lh < n; lh++) {
	rh = RandomInteger(lh, n - 1);
	temp = xyz[lh]; 
	xyz[lh] = xyz[rh]; 
	xyz[rh] = temp;
    }
}
