/* Least-squares interpolation for 2D validition
 */
/*
  Copyright (C) 2014  Xi'an Jiaotong University, UT Austin (Pengliang Yang)

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
#include <math.h>

#include "masklop.h"
#include "fftn.h"

int main(int argc, char* argv[])
{
    bool verb;
    int i, iter, niter, nouter, n[2]; 
    float eps;
    sf_file Fin, Fout, Fmask;/* mask and I/O files*/ 

    /* define temporary variables */
    int n1,n2;
    float *w, *mask;
    sf_complex *dobs, *mm, *dd;

    sf_init(argc,argv);	/* Madagascar initialization */

    /* setup I/O files */
    Fin=sf_input("in");	/* read the data to be interpolated */
    Fmask=sf_input("mask");  	/* read the 2-D mask for 3-D data */
    Fout=sf_output("out"); 	/* output the reconstructed data */
 
    if(!sf_getbool("verb",&verb))    	verb=false;
    /* verbosity */
    if (!sf_getint("niter",&niter)) 	niter=100;
    /* inner iterations */
    if (!sf_getint("nouter",&nouter)) 	nouter=5;
    /* outer iterations */
    if (!sf_getfloat("eps",&eps)) 	eps=1.e-2;
    /* regularization parameter */

    /* Read the data size */
    if (!sf_histint(Fin,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(Fin,"n2",&n2)) sf_error("No n2= in input");

    /* allocate data and mask arrays */
    w=sf_floatalloc(n1*n2);
    dobs=sf_complexalloc(n1*n2);
    mm=sf_complexalloc(n1*n2);
    dd=sf_complexalloc(n1*n2);

    sf_floatread(w,n1*n2,Fin);
    for(i=0; i<n1*n2; i++) dobs[i]=w[i];
    if (NULL != sf_getstring("mask")){
	mask=sf_floatalloc(n2);
	sf_floatread(mask,n2,Fmask);
    } else {
	mask=NULL;
    }

    n[0]=n1; n[1]=n2;
    fftn_init(2, n);
    mask_init(n1, n2, mask);
    for(i=0; i<n1*n2; i++) w[i]=1.0;

    for(iter=0; iter<nouter; iter++){
	sf_csolver_prec(mask_lop, sf_ccgstep, fftn_lop, n1*n2, n1*n2, n1*n2, dd, dobs, niter, eps, "mwt",w,"xp",mm, "verb",verb,"end");
    	sf_ccgstep_close();
	for(i=0; i<n1*n2; i++) w[i]=cabsf(mm[i]); 
    }
    fftn_close();

    for(i=0; i<n1*n2; i++) w[i]=crealf(dd[i]);
    sf_floatwrite(w,n1*n2,Fout); /* output reconstructed seismograms */

    free(w);
    free(dobs);
    free(mm);
    free(dd);

    exit(0);
}
