/* Interface for image-domain waveform tomography (Non-linear CG). */
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

#include "iwinlcg.h"

static bool verb;
static int nn[3], ss[2], dorder, grect[2], gliter;
static float dd[3], **vel, plower, pupper, geps, gscale;
static sf_fslice sfile, rfile;
static float *image, *pimage, *pipz, *piph;
static float lower, upper;

void scale(const float *x, float *g)
/* re-scale gradient */
{
    int i1, i2;
    float scale=0.;

    for (i2=0; i2 < nn[1]; i2++) {
	for (i1=0; i1 < nn[0]; i1++) {
	    scale = (fabsf(g[i2*nn[0]+i1])/x[i2*nn[0]+i1])>scale? 
		fabsf(g[i2*nn[0]+i1])/x[i2*nn[0]+i1]: scale;
	}
    }

    if (scale <= gscale) return;

    scale /= gscale;

    for (i2=0; i2 < nn[1]; i2++) {
	for (i1=0; i1 < nn[0]; i1++) {
	    g[i2*nn[0]+i1] /= scale;
	}
    }
}

void iwinlcg_init(bool verb0,
		  char *order,
		  int npml, float vpml,
		  int n1, int n2,
		  float d1, float d2,
		  int nh, int ns, 
		  float ow, float dw, int nw,
		  sf_file source, sf_file data,
		  bool load, char *datapath,
		  int uts,
		  int prect1, int prect2, int prect3,
		  int pliter,
		  float plower0, float pupper0, 
		  int dorder0,
		  int grect1, int grect2,
		  int gliter0, float geps0, float gscale0,
		  float lower0, float upper0)
/*< initialization >*/
{
    verb = verb0;

    /* model */
    nn[0] = n1; nn[1] = n2; nn[2] = 2*nh+1;
    ss[0] = 1;  ss[1] = n1;
    dd[0] = d1; dd[1] = d2; dd[2] = d2;

    /* control */
    plower = plower0;
    pupper = pupper0;
    dorder = dorder0;
    grect[0] = grect1; grect[1] = grect2;
    gliter = gliter0; geps = geps0;
    gscale = gscale0;    
    lower = lower0; upper = upper0;

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
    iwimodl_init(order,npml,vpml,
		 n1,n2, d1,d2,
		 nh,ns, ow,dw,nw,
		 source,data, sfile,rfile,
		 load,datapath, uts);

    /* dip estimator */
    iwidip_init(n1,n2,nh, d1,d2,
		prect1,prect2,prect3,
		pliter);

    /* tomography operator */
    iwigrad_init(order,npml,vpml,
		 n1,n2, d1,d2,
		 nh,ns, ow,dw,nw,
		 sfile,rfile,
		 load,datapath, uts);
}

void iwinlcg_free()
/*< free allocated memory >*/
{
    sf_fslice_close(sfile);
    sf_fslice_close(rfile);

    free(image);
    free(pimage);
    free(pipz);
    free(piph);
}

float iwinlcg_eval(const float *x,
		   float ***mask, float ***wght)
/*< forward modeling and evaluate objective function >*/
{
    int i, j;
    float fx=0.;

    for (j=0; j < nn[1]; j++) {
	for (i=0; i < nn[0]; i++) {
	    if (x[j*nn[0]+i] < lower) {
		return SF_HUGE;
	    }
	    if (x[j*nn[0]+i] > upper) {
		return SF_HUGE;
	    }

	    vel[j][i] = x[j*nn[0]+i];
	}
    }

    if (verb) sf_warning("Forward modeling...");

    /* clean temporary image */
    iwimodl_clean();

    /* forward modeling */
    iwimodl_modl(vel,image);
    
    if (mask != NULL) {
	for (i=0; i < nn[0]*nn[1]*nn[2]; i++) {
	    image[i] *= mask[0][0][i];
	}
    }

    /* estimate slope */
    iwidip_fdip(image, pimage);

    /* evaluate objective function */
    for (i=0; i < nn[0]*nn[1]*nn[2]; i++) {
	/* thresholding */
	if (fabsf(pimage[i]) >= plower && 
	    fabsf(pimage[i]) <= pupper) {
	    fx += image[i]*wght[0][0][i]*image[i]*wght[0][0][i];
	}
    }

    return fx;
}

void iwinlcg_grad(const float *x,
		  float ***wght, float **prec,
		  float *g)
/*< prepare image perturbation and compute gradient >*/
{
    int i1, i2, i3, ii;
    float *din, *dout, *p;
    sf_triangle tr;

    for (i2=0; i2 < nn[1]; i2++) {
	for (i1=0; i1 < nn[0]; i1++) {
	    vel[i2][i1] = x[i2*nn[0]+i1];
	}
    }

    if (verb) sf_warning("Computing gradient...");

    /* set-up linear operator */
    iwigrad_set(vel, wght,prec);

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
		pipz[i3*nn[0]*nn[1]+i2*nn[0]+i1] = dout[i1]/dd[0];
	    }
	}
    }

    sf_deriv_free(); free(din); free(dout);

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
		piph[i3*nn[0]*nn[1]+i2*nn[0]+i1] = dout[i3]/dd[2];
	    }
	}
    }

    sf_deriv_free(); free(din); free(dout);

    /* prepare image perturbation */
    for (i3=0; i3 < nn[2]; i3++) {
	for (i2=0; i2 < nn[1]; i2++) {
	    for (i1=0; i1 < nn[0]; i1++) {
		ii = i3*nn[0]*nn[1]+i2*nn[0]+i1;

		/* thresholding */
		if (fabsf(pimage[ii]) < plower || 
		    fabsf(pimage[ii]) > pupper) {
		    pimage[ii] = 0.;
		    continue;
		}

		/* non-stationary focusing */
		if (i3 < (nn[2]-1)/2) {
		    pimage[ii] = (pimage[ii]<0.? 1.: -1.)
			*(pipz[ii]-pimage[ii]*piph[ii])
			/sqrtf(1.+pimage[ii]*pimage[ii])
			*wght[0][0][ii];
		} else if (i3 == (nn[2]-1)/2) {
		    pimage[ii] = 0.;
		} else {
		    pimage[ii] = (pimage[ii]>0.? 1.: -1.)
			*(pipz[ii]-pimage[ii]*piph[ii])
			/sqrtf(1.+pimage[ii]*pimage[ii])
			*wght[0][0][ii];
		}
	    }
	}
    }

    if (gliter == 0) {
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
    } else {
	/* shaping regularization */
	sf_trianglen_init(2,grect,nn);
	sf_repeat_init(nn[0]*nn[1],1,sf_trianglen_lop);

	sf_conjgrad_init(nn[0]*nn[1],nn[0]*nn[1],
			 nn[0]*nn[1]*nn[2],nn[0]*nn[1]*nn[2],
			 geps,1.e-6,false,false);
	p = sf_floatalloc(nn[0]*nn[1]);
	
	/* Gauss-Newton */
	sf_conjgrad(NULL, iwigrad_oper,sf_repeat_lop, p,image,pimage,gliter);
	
	free(p); sf_conjgrad_close(); sf_trianglen_close();
    }

    /* convert to velocity */
    for (i2=0; i2 < nn[1]; i2++) {
	for (i1=0; i1 < nn[0]; i1++) {
	    g[i2*nn[0]+i1] = (-0.5*powf(x[i2*nn[0]+i1],3.)*image[i2*nn[0]+i1]);
	}
    }

    /* re-scale */
    scale(x,g);
}

void iwinlcg_image(sf_file fimage)
/*< output image >*/
{
    sf_floatwrite(image,nn[0]*nn[1]*nn[2],fimage);
}
