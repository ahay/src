/* Difference profile of two data */
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

#include <rsf.h>
#include <math.h>

int main(int argc, char *argv[])
{
    int n1, n2, sn1, sn2, i, j;
    float *data1, *data2, *difference;
    float maxdata1;
    sf_file in, out, sub;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");
    sub = sf_input("subtracter");

    /* get data size */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_histint(sub,"n1",&sn1)) sf_error("No n1= in subtracter");
    sn2 = sf_leftsize(sub,1);

    if (n1!=sn1 || n2!=sn2) sf_error("Different size in two data");

    data1 = sf_floatalloc(n1*n2);
    data2 = sf_floatalloc(n1*n2);
    difference = sf_floatalloc(n1*n2);

    sf_floatread(data1,n1*n2,in);
    sf_floatread(data2,n1*n2,sub);

       /* loop over traces */
    for (i=0; i < n2; i++) {
        for (j=0; j < n1; j++) {
            difference[n1*i+j]= data1[n1*i+j]-data2[n1*i+j];
        }
    }

    maxdata1=0.;
    for (i=0; i<n1*n2; i++) {
    if (maxdata1 < fabs(data1[i])) {
        maxdata1 = fabs(data1[i]);
        }
    }

    difference[n1*n2-1]=maxdata1;
    sf_floatwrite(difference,n1*n2,out);


    exit(0);
}
/* 	$Id: Mdifference.c 4796 2009-09-29 19:39:07Z ivlad $	 */
