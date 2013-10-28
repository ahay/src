/* Apply Spitz filter operator in patches */
/*
  Copyright (C) 2010 Politecnico di Milano
  
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

#include <rsf.h>
/*^*/

#include <rsfgee.h>

#include "patch1.h"
#include "spitzbl.h"
#include "patching1bl.h"

void patching1bl(float* modl /* input */, 
		 float* data      /* output */, 
		 int dim          /* number of dimensions */, 
		 int* npatch      /* number of patches [dim] */, 
		 int* nwall       /* data size in[dim] */, 
		 int* nwind       /* patch size in[dim] */, 
		 int* nwall_out   /* data size out [dim] */, 
		 int* nwind_out   /* patch size out [dim] */, 
		 float* windwt    /* window weight */,
		 int order        /* linear PEF order */,
		 int ntraces      /* number of traces to be interpolated*/,
		 float *f		   /* frequency bandwitch normalized [0 0.5]*/,
		 bool verb		   /* verbosity flag*/,
		 bool norm        /* normalization flag*/	)
/*< patch spitz filter >*/
{

    float *winmodl, *windata, *wallwt;
    int i, j, iw, ip, np, n, nw, n_out, nw_out;
    //char WB[]="|/-\\|/-\\";

    np = n = nw = n_out = nw_out = 1; 
    for (j=0; j < dim; j++) {
	np *= npatch[j];
	n *= nwall[j];
	nw *= nwind[j];
	n_out *= nwall_out[j];
	nw_out *= nwind_out[j];
    }
  
    winmodl = sf_floatalloc(nw);
    windata = sf_floatalloc(nw_out);
    wallwt  = sf_floatalloc(n_out);

    for (i=0; i < n_out; i++) data[i] = 0.;

    patch_init(dim, npatch, nwall, nwind);
    patch1_init(dim, npatch, nwall_out, nwind_out);    

    for (ip = 0; ip < np; ip++) {
	/* modl -> winmodl */
	patch_lop(false, false, n, nw, modl, winmodl);	
	/* winmodl -> windata */

	spitzbl(winmodl,windata, nwind, order, ntraces, f, norm);

	/* apply window weighting */
	for (iw=0; iw < nw_out; iw++) windata[iw] *= windwt[iw];
	/* data <- windata */
	patch1_lop(true, true, n_out, nw_out, data, windata);
	patch_close();
	patch1_close();
	if (verb==true) sf_warning("##### Spitz on Patch %d/%d #####",ip+1,np); 
    	//if (verb==true) sf_warning("\r %3.2f%%      %c",(float)100*(ip+1)/np,WB[ip%8]);
    	//fflush(stdout);
    }

	

    /* windwt -> wallwt */
    mkwallwt(dim, npatch, nwall_out, nwind_out, windwt, wallwt);

    /* apply wall weighting */
    for (i=0; i < n_out; i++) data[i] *= wallwt[i];

    free (winmodl);
    free (windata);
    free (wallwt);
}

