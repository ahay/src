/* Interface for image-domain waveform tomography (L-BFGS). */
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

#include "iwidip.h"
#include "iwimodl.h"
#include "iwigrad.h"

#include "iwilbfgs.h"

static int nn[3], dorder;
static sf_fslice sfile, rfile;
static float *image, *pimage, *pipz, *piph;

void iwilbfgs_init(char *order,
		   int npml, float vpml, 
		   int n1, int n2, 
		   float d1, float d2,
		   int nh, int ns, 
		   float ow, float dw, int nw,
		   sf_file source, sf_file data,
		   bool load, char *datapath,
		   int uts,
		   int prect1, int prect2, int prect3,
		   int porder, int pniter, int pliter,
		   int dorder0)
/*< initialization >*/
{
    nn[0] = n1; nn[1] = n2; nn[2] = 2*nh+1;
    dorder = dorder0;

    /* open temporary file */
    sfile = sf_fslice_init(n1*n2*ns,nw,sizeof(sf_complex));
    rfile = sf_fslice_init(n1*n2*ns,nw,sizeof(sf_complex));

    /* allocate temporary memory */
    image  = sf_floatalloc(n1*n2*(2*nh+1));
    pimage = sf_floatalloc(n1*n2*(2*nh+1));

    pipz = sf_floatalloc(n1*n2*(2*nh+1));
    piph = sf_floatalloc(n1*n2*(2*nh+1));

    /* forward modeling */
    iwimodl_init(order,npml,vpml,
		 n1,n2, d1,d2,
		 nh,ns, ow,dw,nw,
		 source,data, sfile,rfile,
		 load,datapath, uts);

    /* PWD dip estimator */
    iwidip_init(n1,n2,nh, d1,d2,
		prect1,prect2,prect3,
		porder,pniter,pliter);

    /* tomography operator */
    iwigrad_init(order,npml,vpml,
		 n1,n2, d1,d2,
		 nh,ns, ow,dw,nw,
		 sfile,rfile,
		 load,datapath, uts);
}

void iwilbfgs_free()
/*< free allocated memory >*/
{
    sf_fslice_close(sfile);
    sf_fslice_close(rfile);

    free(image);
    free(pimage);
    free(pipz);
    free(piph);
}

void iwilbfgs_eval(float **vel,
		   float ***mask, float ***wght)
/*< forward modeling and evaluate objective function >*/
{
    int i;

    /* clean temporary image */
    iwimodl_clean();
    
    /* forward modeling */
    iwimodl_modl(vel,image);
    
    if (mask != NULL) {
	for (i=0; i < nn[0]*nn[1]*nn[2]; i++) {
	    image[i] *= mask[0][0][i];
	}
    }

    /* evaluate objective function */
}

void iwilbfgs_grad(float **vel,
		   float ***wght, float **prec,
		   int ng, float *grad,
		   int ni, float *image)
/*< prepare image perturbation and compute gradient >*/
{
    int i1, i2, i3;
    float *din, *dout;

    /* set-up linear operator */
    iwigrad_set(vel, wght,prec);

    /* estimate slope */
    iwidip_both(image, pimage);

    /* partial i partial z */
    din  = sf_floatalloc(nn[0]);
    dout = sf_floatalloc(nn[0]);

    sf_deriv_init(nn[0], dorder, 0.);

    for (i3=0; i3 < nn[2]; i3++) {
	for (i2=0; i2 < nn[1]; i2++) {
	    for (i1=0; i1 < nn[0]; i1++) {
		din[i1] = image[i3*nn[0]*nn[1]+i2*nn[0]+i1];
	    }
	    sf_deriv(din,dout);
	    for (i1=0; i1 < nn[0]; i1++) {
		pipz[i3*nn[0]*nn[1]+i2*nn[0]+i1] = dout[i1];
	    }
	}
    }

    sf_deriv_free();

    free(din);
    free(dout);

    /* partial i partial h */
    din  = sf_floatalloc(nn[2]);
    dout = sf_floatalloc(nn[2]);

    sf_deriv_init(nn[2], dorder, 0.);

    for (i2=0; i2 < nn[1]; i2++) {
	for (i1=0; i1 < nn[0]; i1++) {
	    for (i3=0; i3 < nn[2]; i3++) {
		din[i3] = image[i3*nn[0]*nn[1]+i2*nn[0]+i1];
	    }
	    sf_deriv(din,dout);
	    for (i3=0; i3 < nn[2]; i3++) {
		pipz[i3*nn[0]*nn[1]+i2*nn[0]+i1] = dout[i3];
	    }
	}
    }

    sf_deriv_free();

    free(din);
    free(dout);

    /* assemble image derivative */
    for (i3=0; i3 < nn[2]; i3++) {
	for (i2=0; i2 < nn[1]; i2++) {
	    for (i1=0; i1 < nn[0]; i1++) {
		
	    }
	}
    }

    /* conjugate-gradient */
    iwigrad_oper(true, false, ng, ni, grad, pimage);
}
