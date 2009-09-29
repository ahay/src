/* Initialize a data */
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
    int n1, n2, i, j;
    bool sign;
    float *data;
    sf_file in, out;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");

    /* get data size */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(in,1);

    if (!sf_getbool("sign",&sign)) sign=true;
    /* if y, initialize data with "1", or with "0" */

    data = sf_floatalloc(n1*n2);
    sf_floatread(data,n1*n2,in);

    if(sign){

      /* loop over traces */
      for (i=0; i < n2; i++) {
        for (j=0; j < n1; j++) {
	  data[n1*i+j]= 1.;
        }
      }
    } else {
      /* loop over traces */
      for (i=0; i < n2; i++) {
        for (j=0; j < n1; j++) {
	  data[n1*i+j]= 0.;
        }
      }
    }
    sf_floatwrite(data,n1*n2,out);


    exit(0);
}
/* 	$Id$	 */
