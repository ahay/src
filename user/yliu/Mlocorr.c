/* Local correlation measure between two datasets. */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

int main(int argc, char* argv[])
{
    int dim, dim1, i, n1, i1, i2, n2, n[SF_MAX_DIM]; 
    int nw, iw, ii, jj, kk, ind;
    int rect[SF_MAX_DIM], half[SF_MAX_DIM];
    char key[6];	
    float done, dtwo;
    float *one, *two, *win1, *win2, *cor;
    bool verb;
    sf_file in, other, out;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    other = sf_input("other");

    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity */

    if (SF_FLOAT != sf_gettype(in) ||
        SF_FLOAT != sf_gettype(other)) sf_error("Need float input");

    for (i=0; i < 3; i++) {
	rect[i] = 1;
	n[i] = 1;
    }

    nw = 1;

    dim = sf_filedims(in,n);

    dim1 = -1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i) || i > 2) rect[i]=1;
	/*( rect#=(1,1,1) smoothing radius on #-th axis )*/ 
	if (rect[i] > 1) dim1 = i+1;
	if (0 == rect[i]%2) {
	    rect[i] ++;
	    sf_warning("rect%d has been changed to be %d", i+1, rect[i]);
	}
    }

    for (i=0; i < 3; i++) {
	half[i] = (rect[i]-1)/2;
    }

    n1 = n2 = 1;
    for (i=0; i < dim; i++) {
        if (i < dim1) {
            n1 *= n[i];
	    nw *= rect[i];
        } else {
            n2 *= n[i];
        }
    }

    one = sf_floatalloc(n1);
    two = sf_floatalloc(n1);

    cor = sf_floatalloc(n1);

    win1 = sf_floatalloc(nw);
    win2 = sf_floatalloc(nw);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(one,n1,in);
        sf_floatread(two,n1,other);

	for (i1=0; i1 < n1; i1++) {
	    if (verb) sf_warning("record %d of %d;",i1+1,n1);

	    for (ii=0; ii < rect[2]; ii++) {
		for (jj=0; jj < rect[1]; jj++) {
		    for (kk=0; kk < rect[0]; kk++) {
			ind = (ii-half[2])*n[1]*n[0]+(jj-half[1])*n[0]
			    +i1+(kk-half[0]);
			if (ind >= 0 && ind < n1) {
			    win1[ii*rect[1]*rect[0]+jj*rect[0]+kk]=one[ind];
			    win2[ii*rect[1]*rect[0]+jj*rect[0]+kk]=two[ind];
			} else {
			    win1[ii*rect[1]*rect[0]+jj*rect[0]+kk]=0.;
			    win2[ii*rect[1]*rect[0]+jj*rect[0]+kk]=0.;
			}
		    }
		}
	    }

	    done = 0.;
	    dtwo = 0.;
	    cor[i1] = 0;
	    for (iw=0; iw < nw; iw++) {
		done += win1[iw]*win1[iw];
		dtwo += win2[iw]*win2[iw];
		cor[i1] += win1[iw]*win2[iw];
	    }		
	    cor[i1] = fabsf(cor[i1])/(sqrtf(done*dtwo)+FLT_EPSILON);	
	}
        sf_floatwrite(cor,n1,out);
    }
    if (verb) sf_warning(".");

    exit(0);
}
/* 	$Id$	 */
