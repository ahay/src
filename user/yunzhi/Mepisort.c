/* Sort microseismic surface array recording traces by their distances to a given epicenter. */
/*
  Copyright (C) 2013 University of Texas at Austin

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
#include <math.h>
#include <stdlib.h>
#define RAD2DEG 1.0/3.1415926534*180

// Sorting algorithm adopted from system/main/headersort.c
struct skey {
    int ikey; // sorting index
    float value; // sorting value
};

static int float_key_compare (const void *k1, const void *k2) {
    float f1 = ((struct skey*) k1)->value;
    float f2 = ((struct skey*) k2)->value;
    return (f1<f2)?-1 : (f1>f2)?1 : 0;
    // by default ascending; change signs if need descending
}

int main(int argc, char *argv[])
{
    int n1,n2,n3,n12; // n1 is trace time length, n2 is trace number
    float **origin_data; // n1*n2 origin dataset
    float **sorted_data; // n1*n2 sorted dataset
    float *abs_x,*abs_y,epi_x,epi_y; // absolute and center coordinates
    float *origin_d,*origin_a; // unsorted distances and angles
    float *sorted_d,*sorted_a; // sorted distances and angles
    int nx,ny,i1,i2,i3; // auxiliary indexes
    float x,y;
    struct skey *sorted;

    sf_file in,out,x_file,y_file,dist_file,theta_file;
    sf_init(argc,argv);
    in = sf_input("in"); // input unsorted traces
    /* DEFAULT: input */
    x_file = sf_input("x");
    /* input: x coordinates extraced from seg-y header */
    y_file = sf_input("y");
    /* input: y coordinates extraced from seg-y header */
    out = sf_output("out"); // output sorted traces
    /* DEFAULT: output */
    dist_file = sf_output("dist");
    /* output: distances of traces after sorting */
    theta_file = sf_output("theta");
    /* output: angles of traces after sorting */
    if (!sf_getfloat("epi_x",&epi_x)) sf_error("Need epicenter x.");
    /* parameter: referenced center point x. */
    if (!sf_getfloat("epi_y",&epi_y)) sf_error("Need epicenter y.");
    /* parameter: referenced center point y. */

    // Check input traces dimensions
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input.");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input.");
    n12 = n1 * n2;
    n3 = sf_leftsize(in,2);

    // Check x and y dimensions
    if (!sf_histint(x_file,"n1",&nx) || !sf_histint(y_file,"n1",&ny)) {
        sf_error("Check x or y dimension.");
    } else if ((nx!=1) || (ny!=1)) {
        sf_error("x and y dimension should n1=1.");
    } else if (!sf_histint(x_file,"n2",&nx) ||
            !sf_histint(y_file,"n2",&ny)) {
        sf_error("No n2= in x or y.");
    } else if ((nx!=n2) || (ny!=n2)) {
        sf_error("x or y dimension n2 mismatches.");
    } else if ((sf_leftsize(x_file,2)!=n3) ||
            (sf_leftsize(y_file,2)!=n3)) {
        sf_error("x or y leftsize dimension mismatches.");
    }

    // Set output traces dimensions
    // As default, same as input, no change

    // Set output distances and angles dimensions
    sf_putint(dist_file,"n1",1);
    sf_putint(dist_file,"n2",n2);
    sf_putint(dist_file,"n3",n3);
    sf_putint(theta_file,"n1",1);
    sf_putint(theta_file,"n2",n2);
    sf_putint(theta_file,"n3",n3);

    // Memory allocation for all arrays
    origin_data = sf_floatalloc2(n1,n2);
    abs_x = sf_floatalloc(n2);
    abs_y = sf_floatalloc(n2);
    sorted_data = sf_floatalloc2(n1,n2);
    origin_d = sf_floatalloc(n2);
    origin_a = sf_floatalloc(n2);
    sorted_d = sf_floatalloc(n2);
    sorted_a = sf_floatalloc(n2);
    sorted = (struct skey*) sf_alloc(n2,sizeof(struct skey));

    // Most outside loop on n3
    for (i3=0; i3<n3; i3++) {

        // Read in input files
        sf_floatread(origin_data[0],n12,in);
        sf_floatread(abs_x,n2,x_file);
        sf_floatread(abs_y,n2,y_file);

        // Initialize other arrays
        for (i2=0; i2<n2; i2++) {
            origin_d[i2] = 0.0;
            origin_a[i2] = 0.0;
            sorted_d[i2] = 0.0;
            sorted_a[i2] = 0.0;
            for (i1=0; i1<n1; i1++) {
                sorted_data[i2][i1] = 0.0;
            }
        }

        // Essential computations
        for (i2=0; i2<n2; i2++) {
            // Calculate relative coordinates
            x = abs_x[i2]-epi_x;
            y = abs_y[i2]-epi_y;
            // Convert to polar coordinates
            origin_d[i2] = hypotf(x,y);
                // hypotf: distance with x and y
            origin_a[i2] = RAD2DEG*atan2f(y,x);
                // atan2f: arctan(y/x) with correct quadrant
            sorted[i2].value = origin_d[i2];
            sorted[i2].ikey = i2;
        }

        qsort(sorted,n2,sizeof(struct skey),float_key_compare);
            // qsort: C library function for d&c quick sort

        for (i2=0; i2<n2; i2++) {
            // Sort according to ikey
            for (i1=0; i1<n1; i1++) {
                sorted_data[i2][i1] = origin_data[sorted[i2].ikey][i1];
            }
            sorted_d[i2] = origin_d[sorted[i2].ikey];
            sorted_a[i2] = origin_a[sorted[i2].ikey];
        }

        // Write into output files
        sf_floatwrite(sorted_data[0],n12,out);
        sf_floatwrite(sorted_d,n2,dist_file);
        sf_floatwrite(sorted_a,n2,theta_file);
    }

    exit(0);
}
