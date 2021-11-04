/* Convert RMS to interval velocity using LS and shaping regularization with non-stationary smoothing */
/*
  Copyright (C) 2004 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that ixwt will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <rsf.h>
#include <math.h>
#include "nsmoothder.h"

int main(int argc, char* argv[])
{
    int i, niter, nd, dim, n1, n2, i1, i2;

    int n[SF_MAX_DIM];
    float **vr, **vi, **wt, **v0, wti;
    char key[6];
    sf_file vrms, vint, weight, vout;
    float *sft[SF_MAX_DIM]; /* storing non-stationary shifting size */
    float *rct[SF_MAX_DIM]; /* storing non-stationary smoothing radii */
    int box[SF_MAX_DIM];    /*box means maximum (edge) padding for triangle smoothing, box[i]=max(rect[i])*/
    sf_file rect[SF_MAX_DIM], shift[SF_MAX_DIM];
    int dim1, b; 
    int n1_1, n2_1 ; 


    sf_init(argc,argv);
    vrms = sf_input("in");
    vint = sf_output("out");
    
    if (NULL != sf_getstring("weight")) {
  weight = sf_input("weight");
    } else {
  weight = NULL;
    }

    dim = sf_filedims (vrms,n);

    nd = 1;
    for (i=0; i < dim; i++) {
  nd *= n[i];
    }

    //n1 = n[0];
    //n2 = nd/n1;
    n1_1 = n[0];
    n2_1 = nd/n1_1;

      /*Calculate dim1*/
    dim1 = -1;
    for (i=0; i < dim; i++) {
  snprintf(key,6,"rect%d",i+1);
  if (NULL != sf_getstring(key)) {
    /*( rect# size of the smoothing stencil in #-th dimension /auxiliary input file/ )*/
      rect[i] = sf_input(key);
      if (SF_FLOAT != sf_gettype(rect[i])) sf_error("Need float %s",key);
      dim1 = i;
      snprintf(key,6,"shift%d",i+1);
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


      /*Calculate n1*/
    n1 = n2 = 1;
    for (i=0; i < dim; i++) {
  if (i <= dim1) {
      n1 *= n[i];
  } else {
      n2 *= n[i];
  }
    }

      /*reading the non-stationary smoothing radii*/
    for (i=0; i <= dim1; i++) {
  box[i] = 1;
  if (NULL != rect[i]) {
      rct[i] = sf_floatalloc (n1);
      sft[i] = sf_floatalloc (n1);

      sf_floatread(rct[i],n1,rect[i]);
      sf_fileclose(rect[i]);

      if (NULL != shift[i]) {
    sf_floatread(sft[i],n1,shift[i]);
    sf_fileclose(shift[i]);
      } else {
    for (i1=0; i1 < n1; i1++) {
        sft[i][i1] = 0;
    }
      }

      for (i1=0; i1 < n1; i1++) {
	  b = ceilf(rct[i][i1]+SF_ABS(sft[i][i1]));
	  if (b > box[i]) box[i] = b;
      }     
  } else {
      rct[i] = NULL;
      sft[i] = NULL;
  }
    }

    nsmoothder_init(nd, dim, box, rct, sft, n,n1,n2);


    vr = sf_floatalloc2(n1_1,n2_1);
    vi = sf_floatalloc2(n1_1,n2_1);
    wt = sf_floatalloc2(n1_1,n2_1);
    v0 = sf_floatalloc2(n1_1,n2_1);

    sf_floatread(vr[0],nd,vrms);

    if (NULL != weight) {
  sf_floatread(wt[0],nd,weight);
    } else {
  for (i1=0; i1 < nd; i1++) {
      wt[0][i1] = 1.0f;
  }
    }

    if (!sf_getint("niter",&niter)) niter=100;
    /* maximum number of iterations */

    wti = 0.;
    for (i2=0; i2 < n2_1; i2++) {
  for (i1=0; i1 < n1_1; i1++) {
      wti += wt[i2][i1]*wt[i2][i1];
  }
    }
    if (wti > 0.) wti = sqrtf(n1_1*n2_1/wti);

    for (i2=0; i2 < n2_1; i2++) {
  for (i1=0; i1 < n1_1; i1++) {
      vr[i2][i1] *= vr[i2][i1]*(i1+1.0f); /* vrms^2*t - data */
      wt[i2][i1] *= wti/(i1+1.0f); /* decrease weight with time */   
      v0[i2][i1] = -vr[i2][0];
  }
    }
    
    sf_repeat_lop(false,true,nd,nd,v0[0],vr[0]);

    nsmoothder(niter, wt[0], vr[0], vi[0]);
 
    for (i2=0; i2 < n2_1; i2++) {
  for (i1=0; i1 < n1_1; i1++) {
      vi[i2][i1] -= v0[i2][i1];
  }
    }

    sf_repeat_lop(false,false,nd,nd,vi[0],vr[0]);

    for (i2=0; i2 < n2_1; i2++) {
  for (i1=0; i1 < n1_1; i1++) {
      vr[i2][i1] = sqrtf(fabsf(vr[i2][i1]/(i1+1.0f)));
      vi[i2][i1] = sqrtf(fabsf(vi[i2][i1]));
  }
    }

    sf_floatwrite(vi[0],nd,vint);

    if (NULL != sf_getstring("vrmsout")) {
  /* optionally, output predicted vrms */
  vout = sf_output("vrmsout");

  sf_floatwrite(vr[0],nd,vout);
    }

    exit(0);
}

/*  $Id$   */
