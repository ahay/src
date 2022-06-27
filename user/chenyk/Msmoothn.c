/* N-D non-stationary smoothing (fixed version of sfnsmooth.) 
Note Fomel's Mnsmooth.c may contain bugs in user/fomel/ntriangle.c
t[i] -> t[o+i*d]; s[i]->s[o+i*d] in triple 
Example: chenyk/ortho_non/nsmooth
*/
/*
  Copyright (C) 2020 Zhejiang University & University of Texas at Austin

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
    int *sft[SF_MAX_DIM];
    int box[SF_MAX_DIM], n[SF_MAX_DIM];
    int dim, dim1, i, n1, n2, i1, i2, b, nrep;
    float *data, *smoo, *rct[SF_MAX_DIM];
    char key[8];
    sf_file in, out, rect[SF_MAX_DIM], shift[SF_MAX_DIM];

    sf_init (argc, argv);
    in = sf_input ("in");
    out = sf_output ("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
 
    if (!sf_getint("repeat",&nrep)) nrep=1;
    /* repeat filtering several times */

    dim = sf_filedims (in,n);

    dim1 = -1;

    for (i=0; i < dim; i++) {
	snprintf(key,6,"rect%d",i+1);
	if (NULL != sf_getstring(key)) {
	    /*( rect# size of the smoothing stencil in #-th dimension /auxiliary input file/ )*/
	    rect[i] = sf_input(key);
	    if (SF_FLOAT != sf_gettype(rect[i])) sf_error("Need float %s",key);
	    dim1 = i;
	    snprintf(key,8,"shift%d",i+1);
	    if (NULL != sf_getstring(key)) {
		/*( shift# shifting of the smoothing stencil in #-th dimension /auxiliary input file/ )*/
		shift[i] = sf_input(key);
		if (SF_INT != sf_gettype(shift[i])) sf_error("Need int %s",key);
	    } else {
		shift[i] = NULL;
	    }
	} else {
	    rect[i] = NULL;
	    shift[i] = NULL;
	}
    }

    n1 = n2 = 1;
    for (i=0; i < dim; i++) {
	if (i <= dim1) {
	    n1 *= n[i];
	} else {
	    n2 *= n[i];
	}
    }
	
    data = sf_floatalloc (n1);
    smoo = sf_floatalloc (n1);

    for (i=0; i <= dim1; i++) {
	box[i] = 1;
	if (NULL != rect[i]) {
	    rct[i] = sf_floatalloc (n1);
	    sft[i] = sf_intalloc (n1);

	    sf_floatread(rct[i],n1,rect[i]);
	    sf_fileclose(rect[i]);

	    if (NULL != shift[i]) {
		sf_intread(sft[i],n1,shift[i]);
		sf_fileclose(shift[i]);
	    } else {
		for (i1=0; i1 < n1; i1++) {
		    sft[i][i1] = 0;
		}
	    }

		
	    for (i1=0; i1 < n1; i1++) {
		b = ceilf(rct[i][i1])+SF_ABS(sft[i][i1]);
		if (b > box[i]) box[i] = b;
	    }	    
	} else {
	    rct[i] = NULL;
	    sft[i] = NULL;
	}
    }

    sf_ntrianglen_init(dim1+1,box,n,rct,sft,nrep);

    for (i2=0; i2 < n2; i2++) {
	sf_floatread(data,n1,in);

	sf_ntrianglen_lop(false,false,n1,n1,data,smoo);
	
	sf_floatwrite(smoo,n1,out);
    }    


    exit (0);
}


