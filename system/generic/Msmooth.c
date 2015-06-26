/* Multi-dimensional triangle smoothing. 

January 2012 program of the month:
http://ahay.org/rsflog/index.php?/archives/280-Program-of-the-month-sfsmooth.html
*/
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

int main (int argc, char* argv[]) 
{
    int dim, dim1, i, j, n[SF_MAX_DIM], rect[SF_MAX_DIM], s[SF_MAX_DIM];
    int nrep, irep, n1, n2, i2, i0;
    bool adj, diff[SF_MAX_DIM], box[SF_MAX_DIM];
    char key[6];
    float* data;
    sf_triangle tr;
    sf_file in, out;

    sf_init (argc, argv);
    in  = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    dim = sf_filedims (in,n);
    dim1 = -1;
    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (!sf_getint(key,rect+i)) rect[i]=1;
	/*( rect#=(1,1,...) smoothing radius on #-th axis )*/ 
	if (rect[i] > 1) dim1 = i;
	snprintf(key,6,"diff%d",i+1);
	if (!sf_getbool(key,diff+i)) diff[i]=false;
	/*( diff#=(n,n,...) differentiation on #-th axis )*/
	snprintf(key,5,"box%d",i+1);
	if (!sf_getbool(key,box+i)) box[i]=false;
	/*( box#=(n,n,...) box (rather than triangle) on #-th axis )*/
    }

    n1 = n2 = 1;
    for (i=0; i < dim; i++) {
	if (i <= dim1) {
	    s[i] = n1;
	    n1 *= n[i];
	} else {
	    n2 *= n[i];
	}
    }

    data = sf_floatalloc (n1);

    if (!sf_getint("repeat",&nrep)) nrep=1;
    /* repeat filtering several times */

    if (!sf_getbool("adj",&adj)) adj=false;
    /* run in the adjoint mode */

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data,n1,in);

	for (i=0; i <= dim1; i++) {
	    if (rect[i] <= 1) continue;
	    tr = sf_triangle_init (rect[i],n[i],box[i]);
	    for (j=0; j < n1/n[i]; j++) {
		i0 = sf_first_index (i,j,dim1+1,n,s);
		for (irep=0; irep < nrep; irep++) {
		    if (adj) {
			sf_smooth (tr,i0,s[i],diff[i],data);
		    } else {
			sf_smooth2 (tr,i0,s[i],diff[i],data);
		    }
		}
	    }
	    sf_triangle_close(tr);
	}
	
	sf_floatwrite(data,n1,out);
    }    


    exit (0);
}

/* 	$Id$	 */
