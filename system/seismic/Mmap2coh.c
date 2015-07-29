/* From parameter's attribute map (veltran) to coherency-like plots. 
   (eventually masked) */
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

#include <math.h>
#include <float.h>
#include <rsf.h>

int main (int argc, char* argv[])
{
    int it,ix,i2,iv, nt,nx, n2, nv, ntv, nw, imin2, imax2;
    float dt, t0, t,o2, d2, v0, dv, min2, max2;
    float *M=NULL, *ord=NULL, *v=NULL, *v2=NULL;
    float **coord=NULL;
    sf_file cmp=NULL, map=NULL, coh=NULL;
    
    sf_init (argc,argv);
    cmp = sf_input("in");
    
    /* Auxiliary inputs */

    if (NULL != sf_getstring("map")) map = sf_input("map"); 
    /*parameters map */
    else sf_error("Need parameter map");
    
    coh = sf_output("out");
    
    if (SF_FLOAT != sf_gettype(cmp)) sf_error("Need float input");
    
    if (!sf_histint(cmp,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(cmp,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(cmp,"o1",&t0)) sf_error("No o1= in input");
    
    if (!sf_histint(cmp,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histfloat(cmp,"d2",&d2)) sf_error("No d2= in input");
    if (!sf_histfloat(cmp,"o2",&o2)) sf_error("No o2= in input");
    
    /* parameter sampling */
    if (!sf_getint("nv",&nv)) sf_error("Need nv=");       /* number of velocities */
    if (!sf_getfloat("v0",&v0)) sf_error("Need v0=");     /* velocity origin */
    if (!sf_getfloat("dv",&dv)) sf_error("Need dv=");     /* velocity sampling */
    
    /* slope/offset window */
    if (!sf_getfloat("min2",&min2)) min2=o2;     	   /* min2 */
    if (min2<o2) min2=o2; 
	
    if (!sf_getfloat("max2",&max2)) max2=o2+d2*(n2-1);     /* max2 */    
    if (max2>o2+d2*(n2-1)) max2=o2+d2*(n2-1);     

    imin2=(int) floorf((min2-o2)/d2);
    imax2=(int) ceilf((max2-o2)/d2);
    
    if (min2>max2) sf_error("min2=%g must be less than max2=%f",min2,max2);
    sf_warning("imin2=%d imax2=%d",imin2,imax2);
    /* adding the second dimension to the ouput files*/
    sf_putint(coh,"n2",nv);
    sf_putfloat(coh,"o2",v0);
    sf_putfloat(coh,"d2",dv);

    ntv = nt*nv;
    
    v  = sf_floatalloc(ntv);
    v2 = sf_floatalloc(ntv);
    ord = sf_floatalloc(nt);	    
    M = sf_floatalloc(nt);

    coord = sf_floatalloc2(2,nt);

    /* reading the number of cmp in the data */
    nx = sf_leftsize(cmp,2);
    
    if (!sf_getint("nw",&nw)) nw=4;
    /* interpolator size (2,3,4,6,8) */
    
    for (ix = 0; ix < nx; ix++) { /* CMP loop*/
	
	for (iv=0; iv < ntv; iv++) {
	    v[iv]=0.;
	}

	sf_seek(cmp,   imin2 * nt * sizeof(float),  SEEK_SET);
	sf_seek(map,   imin2 * nt * sizeof(float),  SEEK_SET);
        
	for (i2 = imin2; i2 <= imax2; i2++) { /* slope p loop*/
	    
	    

	    sf_floatread (ord, nt ,cmp);
	    sf_floatread (M  , nt ,map);	    
	    for (it=0; it < nt; it++) { /* time tau loop*/
                t = t0 + it*dt;
		coord[it][0] = t; /* T0[it]; */
		coord[it][1] = M[it];
				
            } /* END t loop */
	    
	    
	    sf_int2_init (coord, t0,v0, dt,dv, nt,nv, sf_spline_int, nw, nt);
	    sf_int2_lop (true,true,ntv,nt,v,ord);
	    
        } /* END slope p loop */
        
	/* sf_warning("QUI"); */
	/* from spline coefficients to model */
	if (nw > 2) { /* loop spline */
	    for (iv=0; iv < nv; iv++) {
		sf_spline_post (nw, iv*nt, 1, nt, v, v2);
	    }
	    
	    for (it=0; it < nt; it++) {
		sf_spline_post (nw, it, nt, nv, v2, v);
		
	    }
	} /* END loop spline */
        
	sf_floatwrite (v,ntv,coh);
	    
    } /* END CMP loop */
 
    exit (0);
} /* END of main */
