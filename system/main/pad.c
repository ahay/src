/* Pad a dataset with zeros.

   n#out is equivalent to n#, both of them overwrite end#.

   Other parameters from the command line are passed to the output (similar to sfput).
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

#include <stdio.h>
#include <string.h>

#include <rsf.h>

int main (int argc, char* argv[])
{
    int i, j, nj, dim, ntr, itr, esize;
    int beg[SF_MAX_DIM], end[SF_MAX_DIM];
    off_t n[SF_MAX_DIM], n2[SF_MAX_DIM];
    size_t n0, n20, beg0, end0;
    sf_file in, out;
    float o, d;
    char key[10], key2[10], *tr, *zero;
    bool inside;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    dim = sf_largefiledims(in,n);

    sf_expandpars(out);
    
    ntr=1;
    for (j=0; j < SF_MAX_DIM; j++) {
	i=j+1;
	if (j>=dim) {
	    n[j]=1;
	}
	sprintf(key,"beg%d",i);
	if(!sf_getint(key,beg+j)) {
	    /*( beg#=0 the number of zeros to add before the beginning of #-th axis )*/
	    beg[j]=0;
	} else if (beg[j]<0) {
	    sf_error("negative beg=%d",beg[j]);
	}
	sprintf(key,"end%d",i);
	sprintf(key2,"n%d",i);
	if(!sf_getint(key,end+j)) {
	    /*( end#=0 the number of zeros to add after the end of #-th axis )*/
	    sprintf(key,"n%dout",i);
	    if (sf_getint(key,&nj) || sf_getint(key2,&nj)) {
		/*( n# the output length of #-th axis - padding at the end )*/
		if (0==nj) 
		    for (nj++; nj < n[j]; nj *= 2) ;
			
		end[j]=nj-n[j]-beg[j];
		if (end[j]<0)
		    sf_error("negative end=%d",end[j]);
	    } else {
		end[j]=0;
	    }
	}
	n2[j]=n[j]+beg[j]+end[j];
	if (j>0) ntr *= n2[j];

	if (n2[j] != n[j]) sf_putint(out,key2,n2[j]);

	if (beg[j] > 0) {
	    sprintf(key,"o%d",i);
	    if (sf_histfloat(in,key,&o)) {
		sprintf(key2,"d%d",i);
		if (sf_histfloat(in,key2,&d))
		    sf_putfloat(out,key,o-d*beg[j]);
	    }
	}
	if (n2[j] > 1 && j >= dim) dim=i;
    }

    if (!sf_histint(in,"esize",&esize)) {
	esize=sf_esize(in);
    } else if (0>=esize) {
	sf_error("cannot handle esize=%d",esize);
    }

    sf_fileflush(out,in);
    sf_setform(in,SF_NATIVE);
    sf_setform(out,SF_NATIVE);

    n[0]   *= esize; n0   = (size_t) n[0];
    n2[0]  *= esize; n20  = (size_t) n2[0];
    beg[0] *= esize; beg0 = (size_t) beg[0];
    end[0] *= esize; end0 = (size_t) end[0];

    tr   = sf_charalloc(n20);
    zero = sf_charalloc(n20);
    memset(zero,0,n20);

    if (beg0>0) memcpy(tr,zero,beg0);
    if (end0>0) memcpy(tr+beg0+n0,zero,end0);

    for (itr=0; itr < ntr; itr++) {
    	inside = true;
    	for (nj=j=1; j < dim; nj *= n2[j], j++) {
  	    i = (itr/nj)%n2[j];
   	    if (i < beg[j] || i >= beg[j]+n[j]) {
       		inside = false;
      		break;
        }
    	}
	    if (inside) {
	      sf_charread (tr+beg0,n0,in);
	      sf_charwrite(tr,n20,out);
	    } else {
	      sf_charwrite(zero,n20,out);
	    }
    }
    exit (0);
}

/* 	$Id$	 */
