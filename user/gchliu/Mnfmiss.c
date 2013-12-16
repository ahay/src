/* Missing data interpolation in freq domain. */
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

#include "nfmis.h"

int main(int argc, char* argv[])
{
    int n1, na2,na12, niter, xniter, nf,i,f1,iy,ny;
    sf_complex *xx, *aa;
    float *kk;
    bool *known,exact, verb;
    sf_file in, out, filt,mask=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    filt = sf_input("filt");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nf)) sf_error("No n2= in input");
    ny = sf_leftsize(in,2);

	if (!sf_histint(filt,"n2",&na2)) sf_error("No n2= in filtin");
	na12=n1*na2;


    if (!sf_getbool("exact",&exact)) exact=true;
    /* If y, preserve the known data values */

    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

    
    if (!sf_getint("niter",&niter)) niter=n1;
    /* number of iterations */
    
    if (!sf_getint("xniter",&xniter)) xniter=1;
    /* number of iterations */
    
    if (NULL != sf_getstring("mask")) {
	/* optional input mask file for known data */
	mask = sf_input("mask");
    }
   
    xx    = sf_complexalloc(n1);
    known = sf_boolalloc(n1);
    aa    = sf_complexalloc(na12);
    kk    = sf_floatalloc(n1);


    for (iy=0;iy<ny;iy++){ // y
	           sf_warning("y slice %d of %d",iy+1,ny);
         if (NULL != mask) {
              sf_floatread(kk,n1,mask);
              for (i=0; i < n1; i++) {
      	          known[i] = (bool) (cabsf(kk[i]) != 0.);
	          }
         }
	
        
         for (f1=0;f1 < nf; f1++) { // freq
	           sf_warning("Frequency slice %d of %d",f1+1,nf);
               sf_complexread(aa,na12,filt);
               sf_complexread(xx,n1,in);

               if (NULL == mask) { 

	                for (i=0; i < n1; i++) {
	                    known[i] = (bool) (cabsf(xx[i]) != 0.);
	                }
               }


               nfmis(niter, n1, na2, aa, xx, known, verb);           
               sf_complexwrite (xx,n1,out);
	     } /* freq end*/
    } /*y end*/
	
	
    exit(0);
}

/* 	$Id: Mmiss1.c 7107 2011-04-10 02:04:14Z ivlad $	 */
