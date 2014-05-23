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

#include "lbfgs.h"

#include "iwimodl.h"
#include "iwigrad.h"

#include "iwilbfgs.h"

static bool verb, deriv, bound;
static int nn[3], ss[2], grect[2];
static float dd[3], **vel;
static sf_fslice sfile, rfile;
static float *image, *pimage, *pipz, *piph;
static lbfgsfloatval_t gscale, lower, upper;

void didz(bool adj, float *img)
/* apply derivative */
{
    int i1, i2, i3;
    float *din, *dout;

    din  = sf_floatalloc(nn[0]);
    dout = sf_floatalloc(nn[0]);

    for (i3=0; i3 < nn[2]; i3++) {
	for (i2=0; i2 < nn[1]; i2++) {
	    for (i1=0; i1 < nn[0]; i1++) {
		din[i1] = img[i3*nn[0]*nn[1]+i2*nn[0]+i1];
	    }
	    
	    if (adj) {
		sf_igrad1_lop (true, false,nn[0],nn[0],dout,din);
	    } else {
		sf_igrad1_lop (false,false,nn[0],nn[0],din,dout);
	    }

	    for (i1=0; i1 < nn[0]; i1++) {
		img[i3*nn[0]*nn[1]+i2*nn[0]+i1] = dout[i1]/dd[0];
	    }
	}
    }
    
    free(din); free(dout);
}

void scale(const lbfgsfloatval_t *x, lbfgsfloatval_t *g)
/* re-scale gradient */
{
    int i1, i2;
    lbfgsfloatval_t scale=0.;

    /* scan for maximum absolute ratio */
    for (i2=0; i2 < nn[1]; i2++) {
	for (i1=0; i1 < nn[0]; i1++) {
	    scale = (fabs(g[i2*nn[0]+i1])/x[i2*nn[0]+i1])>scale? 
		fabs(g[i2*nn[0]+i1])/x[i2*nn[0]+i1]: scale;
	}
    }

    if (scale <= gscale) return;

    /* re-scale to pre-defined ratio */
    scale /= gscale;

    for (i2=0; i2 < nn[1]; i2++) {
	for (i1=0; i1 < nn[0]; i1++) {
	    g[i2*nn[0]+i1] /= scale;
	}
    }
}

void iwilbfgs_init(bool verb0,
		   char *order,
		   int npml,
		   int n1, int n2,
		   float d1, float d2,
		   int nh, int ns, 
		   float ow, float dw, int nw,
		   sf_file source, sf_file data,
		   bool deriv0,
		   bool load, char *datapath,
		   int uts,
		   int grect1, int grect2,
		   float gscale0,
		   float lower0, float upper0)
/*< initialization >*/
{
    verb = verb0;
    deriv = deriv0;

    /* model */
    nn[0] = n1; nn[1] = n2; nn[2] = 2*nh+1;
    ss[0] = 1;  ss[1] = n1;
    dd[0] = d1; dd[1] = d2; dd[2] = d2;

    /* control */
    grect[0] = grect1; grect[1] = grect2;
    gscale = (lbfgsfloatval_t) gscale0;
    lower = (lbfgsfloatval_t) lower0;
    upper = (lbfgsfloatval_t) upper0;

    /* open temporary file */
    sfile = sf_fslice_init(n1*n2*ns,nw,sizeof(sf_complex));
    rfile = sf_fslice_init(n1*n2*ns,nw,sizeof(sf_complex));

    /* allocate temporary memory */
    vel = sf_floatalloc2(n1,n2);

    image  = sf_floatalloc(n1*n2*(2*nh+1));
    pimage = sf_floatalloc(n1*n2*(2*nh+1));

    pipz = sf_floatalloc(n1*n2*(2*nh+1));
    piph = sf_floatalloc(n1*n2*(2*nh+1));

    /* forward modeling */
    iwimodl_init(order,npml,
		 n1,n2, d1,d2,
		 nh,ns, ow,dw,nw,
		 source,data, sfile,rfile,
		 load,datapath, uts);

    /* tomography operator */
    iwigrad_init(order,npml,
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

lbfgsfloatval_t iwilbfgs_eval(const lbfgsfloatval_t *x,
			      float ***wdso, float ***wstk)
/*< forward modeling and evaluate objective function >*/
{
    int i1, i2, i3, ii;
    lbfgsfloatval_t fx=0.;

    bound = false;
    for (i2=0; i2 < nn[1]; i2++) {
	for (i1=0; i1 < nn[0]; i1++) {
	    if (x[i2*nn[0]+i1] < lower) {
		bound = true;
		return DBL_MAX;
	    }
	    if (x[i2*nn[0]+i1] > upper) {
		bound = true;
		return DBL_MAX;
	    }

	    vel[i2][i1] = (float) x[i2*nn[0]+i1];
	}
    }

    if (verb) sf_warning("Forward modeling...");

    /* clean temporary image */
    iwimodl_clean();

    /* forward modeling */
    iwimodl_modl(vel,image);

    /* apply image derivative */
    if (deriv) didz(false,image);
    
    /* evaluate objective function */
    for (i3=0; i3 < nn[2]; i3++) {
	for (i2=0; i2 < nn[1]; i2++) {
	    for (i1=0; i1 < nn[0]; i1++) {
		ii = i3*nn[0]*nn[1]+i2*nn[0]+i1;

		/* DSO */
		fx += (lbfgsfloatval_t) 0.5*
		    image[ii]*wdso[i3][i2][i1]*image[ii]*wdso[i3][i2][i1];

		/* STK */
		if (wstk != NULL) {
		    fx -= (lbfgsfloatval_t) 0.5*
			image[ii]*wstk[i3][i2][i1]*image[ii]*wstk[i3][i2][i1];
		}
	    }
	}
    }

    return fx;
}

void iwilbfgs_grad(const lbfgsfloatval_t *x,
		   float ***wdso, float ***wstk, float **prec,
		   lbfgsfloatval_t *g)
/*< prepare image perturbation and compute gradient >*/
{
    int i1, i2, i3, ii;
    sf_triangle tr;

    if (bound) {
	for (i2=0; i2 < nn[1]; i2++) {
	    for (i1=0; i1 < nn[0]; i1++) {
		g[i2*nn[0]+i1] = DBL_MAX;
	    }
	}

	return;
    }

    if (verb) sf_warning("Computing gradient...");

    /* set-up linear operator */
    iwigrad_set(vel, NULL,prec);    

    /* prepare image perturbation */
    if (deriv) didz(true,image);

    for (i3=0; i3 < nn[2]; i3++) {
	for (i2=0; i2 < nn[1]; i2++) {
	    for (i1=0; i1 < nn[0]; i1++) {
		ii = i3*nn[0]*nn[1]+i2*nn[0]+i1;

		/* DSO */
		pimage[ii] = image[ii]*wdso[i3][i2][i1]*wdso[i3][i2][i1];

		/* STK */
		if (wstk != NULL) {
		    pimage[ii] -= image[ii]*wstk[i3][i2][i1]*wstk[i3][i2][i1];
		}
	    }
	}
    }

    /* steepest descent */
    iwigrad_oper(true,false, nn[0]*nn[1],nn[0]*nn[1]*nn[2], image,pimage);

    /* smooth gradient */
    tr = sf_triangle_init(grect[0],nn[0]);
    for (i2=0; i2 < nn[1]; i2++) {
	i1 = sf_first_index(0,i2,2,nn,ss);
	
	sf_smooth (tr,i1,ss[0],false,false,image);
	sf_smooth2(tr,i1,ss[0],false,false,image);
    }
    sf_triangle_close(tr);
    
    tr = sf_triangle_init(grect[1],nn[1]);
    for (i1=0; i1 < nn[0]; i1++) {
	i2 = sf_first_index(1,i1,2,nn,ss);
	
	sf_smooth (tr,i2,ss[1],false,false,image);
	sf_smooth2(tr,i2,ss[1],false,false,image);
    }
    sf_triangle_close(tr);    

    /* convert to velocity */
    for (i2=0; i2 < nn[1]; i2++) {
	for (i1=0; i1 < nn[0]; i1++) {
	    g[i2*nn[0]+i1] = (lbfgsfloatval_t) 
		(-0.5*powf(x[i2*nn[0]+i1],3.)*image[i2*nn[0]+i1]);
	}
    }

    /* re-scale */
    if (0. < gscale && gscale < 1.) scale(x,g);
}
