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

#include "iwimodl.h"
#include "iwigrad.h"
#include "iwilbfgs.h"

static sf_fslice sfile, rfile;

void iwilbfgs_init(char *order,
		   int npml, float vpml, 
		   int n1, int n2, 
		   float d1, float d2,
		   int nh, int ns, 
		   float ow, float dw, int nw,
		   sf_file source, sf_file data,
		   bool load, char *datapath,
		   int uts)
/*< initialization >*/
{
    /* open temporary file */
    sfile = sf_fslice_init(n1*n2*ns,nw,sizeof(sf_complex));
    rfile = sf_fslice_init(n1*n2*ns,nw,sizeof(sf_complex));

    /* forward modeling */
    iwimodl_init(order,npml,vpml,
		 n1,n2, d1,d2,
		 nh,ns, ow,dw,nw,
		 source,data, sfile,rfile,
		 load,datapath, uts);

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
}

void iwilbfgs_eval(float **vel,
		   float ****image)
/*< forward modeling and evaluate objective function >*/
{
    iwimodl_modl(vel,image);
}

void iwilbfgs_grad(float **vel,
		   float ***wght, float **prec,
		   int ng, float **grad,
		   int ni, float ****image)
/*< prepare image perturbation and compute gradient >*/
{
    iwigrad_set(vel,wght,prec);
    iwigrad_oper(true, false, ng, ni, grad[0], image[0][0][0]);
}
