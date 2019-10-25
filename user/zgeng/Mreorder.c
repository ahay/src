/* Reorder the data according to the path */
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
    int n1, n2, n3;
    int num_samples;            // number of data samples
    int *path_array = NULL;     // input path vector
    float **data, **newdata;
    bool inv;
    sf_file in, out, path;

    sf_init(argc, argv);

    in   = sf_input("in");    
    /* input data */
    path = sf_input("path");
    /* input path vector */
    out  = sf_output("out");
    /* output reordered data */

    /* check dimension */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");
    num_samples = n1*n2*n3;

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, do inverse transform */

    data       = sf_floatalloc2(n1, n2*n3);
    newdata    = sf_floatalloc2(n1, n2*n3);
    path_array = sf_intalloc(n2*n3);

    /* read data and path */
    sf_floatread(data[0], num_samples, in);
    sf_intread(path_array, n2*n3, path);

    /* reorder data */
    newdata = reorder(data, path_array, n1, n2, n3, inv);

    sf_floatwrite(newdata[0], num_samples, out);

    exit(0);
}
