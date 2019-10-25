/* Finding the easy path */
/*
  Copyright (C) 2018 University of Texas at Austin
   
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

#include <rsf.h>

#include "easypath.h"

int main(int argc, char *argv[])
{
    int n1, n2, n3, i;
    int num_samples;                // number of data samples
    int *indexes, *path_array;
    float **data;
    sf_file in, out;

    sf_init(argc, argv);

    in      = sf_input("in");    
    /* input 3-D data */
    out     = sf_output("out");
    /* output path vector */

    /* check dimension */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");
    num_samples = n1*n2*n3;

    /* set up the output path file */
    sf_settype(out, SF_INT);
    sf_putint(out, "n1", n2*n3);
    sf_putint(out, "n2", 1);
    sf_putint(out, "n3", 1);

    data    = sf_floatalloc2(n1, n2*n3);
    indexes = sf_intalloc(n2*n3);

    /* set indexes for each traces */
    for (i = 0; i < n2*n3; i++)
        indexes[i] = i+1;

    /* read 3-D data */
    sf_floatread(data[0], num_samples, in);

    /* get path vector */
    path_array = pathway(data, indexes, n1, n2, n3);

    sf_intwrite(path_array, n2*n3, out);

    exit(0);
}
