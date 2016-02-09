/* Slope-based tau-p 3D velocity transform for elliptical anisotropy.
 * The program returns the squared velocity vx,vy,vxy spectra or maps
 */
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
    bool map,interval;
    sf_map4 nmo;
    int it,ix,ip1,ip2,iv, nt,nx, np1,np2, nvx=0, ntvx=0, nw, nvy=0, ntvy=0, nvxy=0, ntvxy=0;
    float dt, t0, p1, p10, p2, p20, t, dp1,dp2, vx0, dvx, vy0, dvy , vxy0 , dvxy;
    float N,D;
    float *v1=NULL, *v2=NULL, *v3=NULL; /* to compute velocity attribute panels */
    float *Rx=NULL, *Ry=NULL, *Rxy=NULL, **coord1=NULL, **coord2=NULL, **coord3=NULL, *ord=NULL;
    float *vx=NULL, *vx2=NULL, *vy=NULL, *vy2=NULL, *vxy=NULL, *vxy2=NULL, *TAU0=NULL; /* *TAU0t=NULL; */

    sf_file input=NULL, cmp=NULL, dipx=NULL, dipy=NULL, dipxy=NULL;
    /* Input files */
    sf_file velx=NULL, vely=NULL, velxy=NULL;
    /* Output Files */

    float eps = 100*FLT_EPSILON;

    sf_init (argc,argv);
    input = sf_input("in");

    velx = sf_output("out");      /*vx velocity*/
    vely = sf_output("vely");     /*vy velocity*/
    velxy = sf_output("velxy");   /*vxy velocity*/

    if (!sf_getbool("map",&map)) map=false;
    /* output maps instead of coherency panels */

    if (!map) {
    	if (NULL != sf_getstring("cmp")) cmp = sf_input("cmp");
    	else sf_error("Need cmp input when map=y");
    }

    if (!sf_getbool("interval",&interval)) interval=false;
    /* interval values by 3D stripping equations */
    
    if (SF_FLOAT != sf_gettype(input)) sf_error("Need float input");

    if (!sf_histint(input,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histfloat(input,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(input,"o1",&t0)) sf_error("No o1= in input");
    
    if (!sf_histint(input,"n2",&np1)) sf_error("No n2= in input");
    if (!sf_histfloat(input,"d2",&dp1)) sf_error("No d2= in input");
    if (!sf_histfloat(input,"o2",&p10)) sf_error("No o2= in input");

    if (!sf_histint(input,"n3",&np2)) sf_error("No n3= in input");
    if (!sf_histfloat(input,"d3",&dp2)) sf_error("No d3= in input");
    if (!sf_histfloat(input,"o3",&p20)) sf_error("No o3= in input");

    p10 /= dp1; /* slope axis is normalized: no need for slope normalization dp1=1 and dp2=1*/
    p20 /= dp2;

    v1  = sf_floatalloc(nt);
    v2  = sf_floatalloc(nt);
    v3  = sf_floatalloc(nt);
    
    if (map) {    
	nmo = sf_stretch4_init (nt, t0, dt, nt, 0.01);

    } else {
	nmo = NULL;

	ord = sf_floatalloc(nt);	
	coord1 = sf_floatalloc2(2,nt);
	coord2 = sf_floatalloc2(2,nt);
	coord3 = sf_floatalloc2(2,nt);    

	/* first velocity... vx */
	if (!sf_getint("nvx",&nvx)) sf_error("Need nvx=");        /* number of vx squared velocities */
	if (!sf_getfloat("vx0",&vx0)) sf_error("Need vx0=");     /*vx squared velocity origin */
	if (!sf_getfloat("dvx",&dvx)) sf_error("Need dvx=");     /*vx squared velocity sampling */
	/* second velocity... vy */
	if (!sf_getint("nvy",&nvy)) nvy=nvx;     	/* number of vy squared velocities  */
	if (!sf_getfloat("vy0",&vy0)) vy0=vx0;     /* vy squared  velocity origin */
	if (!sf_getfloat("dvy",&dvy)) dvy=dvx;     /* vy squared  velocity sampling */
	/* third velocity... vxy */
	if (!sf_getint("nvxy",&nvxy))   nvxy = 101;   /* number of vxy velocities */
	if (!sf_getfloat("vxy0",&vxy0)) vxy0 =-0.1;   /* vxy   velocity origin */
	if (!sf_getfloat("dvxy",&dvxy)) dvxy = 0.1;   /* vxy   velocity sampling */
    
    
	/* adding the second dimension to the ouput files*/
	sf_putint(velx,"n2",nvx);
	sf_putfloat(velx,"o2",vx0);
	sf_putfloat(velx,"d2",dvx);
	
	sf_putint(vely,"n2",nvy);
	sf_putfloat(vely,"o2",vy0);
	sf_putfloat(vely,"d2",dvy);
	
	sf_putint(velxy,"n2",nvxy);
	sf_putfloat(velxy,"o2",vxy0);
	sf_putfloat(velxy,"d2",dvxy);

	if(!map) {
	    sf_putint(velx,"n3",1);
	    sf_putint(vely,"n3",1);
	    sf_putint(velxy,"n3",1);
	}


	ntvx  = nt * nvx;
	ntvy  = nt * nvy;
	ntvxy = nt * nvxy;

	vx  = sf_floatalloc(ntvx);
	vx2 = sf_floatalloc(ntvx);
	vy  = sf_floatalloc(ntvy);
	vy2 = sf_floatalloc(ntvy);
	vxy  = sf_floatalloc(ntvxy);
	vxy2 = sf_floatalloc(ntvxy);
    }
    
    /* reading the number of cmp in the data */
    nx = sf_leftsize(input,3);
    
    if (!sf_getint("nw",&nw)) nw=4;
    /* interpolator size (2,3,4,6,8) */

    /* Auxiliary inputs */
    TAU0 = sf_floatalloc(nt);

    if (NULL != sf_getstring("dipx")) dipx = sf_input("dipx");
    else sf_error("Need dipx input");

    if (NULL != sf_getstring("dipy")) dipy = sf_input("dipy");
    else sf_error("Need dipy input");

    if (NULL != sf_getstring("dipxy")) dipxy = sf_input("dipxy");
    else sf_error("Need dipxy input");

    Rx = sf_floatalloc(nt);
    Ry = sf_floatalloc(nt);
    Rxy = sf_floatalloc(nt);

    
    for (ix = 0; ix < nx; ix++) { /* CMP loop*/
	
	if (!map) {
	    for (iv=0; iv < ntvx; iv++) {
		vx[iv]=0.;
	    }
	    for (iv=0; iv < ntvy; iv++) {
		vx[iv]=0.;
	    }
	    for (iv=0; iv< ntvxy; iv++) {
		vxy[iv]=0.;
	    }
	}
	for (ip2 = 0; ip2 < np2; ip2++) { /* slope 2 (dimension 3)*/
	    p2 = p20+ip2;

	    for (ip1 = 0; ip1 < np1; ip1++) { /* slope 1 (dimension 2)*/
		p1 = p10+ip1;

		/* sf_warning("ip1=%d ip2=%d np2=%d",ip1,ip2,np2); */

		sf_floatread (TAU0, nt,input);
		sf_floatread (Rx, nt, dipx);
		sf_floatread (Ry, nt, dipy);
		sf_floatread (Rxy,nt, dipxy);

		if (!map) sf_floatread (ord, nt, cmp);


		for (it=0; it < nt; it++) { /* time tau loop*/

		    t = t0 + it*dt;


		    if (TAU0[it] <= 0.) {
                	v1[it] = 0.0;
                	v2[it] = 0.0;
                	v3[it] = 0.0;
		    } else {
                	if (interval) {
			    N = (Rxy[it]  +   Rx[it] * Ry[it] ) / dp1 /dp2;
			    D = (Rx[it] * p1 +  Ry[it]  * p2 - 1 );

			    v1[it] = ( Rx[it] / dp1 - N * p2*dp2 ) / (D*p1*dp1 + eps);
			    v2[it] = ( Ry[it] / dp2 - N * p1*dp1 ) / (D*p2*dp2 + eps);
			    v3[it] = N/(D+eps);

/*
  tau_ = (Rx[it] * p1  +  Ry[it] * p2  - 1.0);
  v3[it] =  ( (Rxy[it]  +  Rx[it] * Ry[it])  ) / dp1 /dp2 ;

  v1[it] = -1.0 * ( (Rx[it] * Ry[it] * 1.0/dp1 * 1.0/dp2 + Rxy[it] * 1.0/dp1/dp2)*p2*dp2 - Rx[it] * 1.0/dp1 )/(p1*dp1*tau_+eps);
  v2[it] = -1.0 * ( (Rx[it] * Ry[it] * 1.0/dp1 * 1.0/dp2 + Rxy[it] * 1.0/dp1/dp2)*p1*dp1 - Ry[it] * 1.0/dp2 )/(p2*dp2*tau_+eps);
*/
			}
			else {
			    N = (Rxy[it] *dt +  (1.0/t+eps) * Rx[it] * Ry[it] *dt *dt) / dp1 /dp2;
			    D = (Rx[it]* dt * p1 +  Ry[it] * dt * p2 - t );

			    v1[it] = ( Rx[it]* dt/dp1 - N * p2*dp2 ) / (D*p1*dp1 + eps);
			    v2[it] = ( Ry[it]* dt/dp2 - N * p1*dp1 ) / (D*p2*dp2 + eps);
			    v3[it] = N/(D+eps);



			    /*tau_ = (-1.0) * (t / (TAU0[it] * TAU0[it] + eps) );*/
/*
  v3[it] =  1.0/(tau_+eps) * (Rxy[it] * dt/dp1/dp2 + 1/(t+eps)* Rx[it] * Ry[it] * dt/dp1 * dt/dp2);
  v1[it] = ( 1.0/(tau_+eps) * Rx[it] *(dt/dp1) - v3[it] * p2*dp2) / ( p1*dp1 + eps);
  v2[it] = ( 1.0/(tau_+eps) * Ry[it] *(dt/dp2) - v3[it] * p1*dp1) / ( p2*dp2 + eps);
*/
			}
		    } /* end of if tau0 >= 0 */
		    if (!map) {
			coord1[it][0] = TAU0[it];
			coord2[it][0] = TAU0[it];
			coord3[it][0] = TAU0[it];

			coord1[it][1] = v1[it];
			coord2[it][1] = v2[it];
			coord3[it][1] = v3[it];
		    }
		} /* END tau t loop */

		if (map) {
		    sf_stretch4_define (nmo,TAU0);

		    sf_stretch4_apply (false,nmo,v1,v1); sf_floatwrite (v1,nt,velx);
		    sf_stretch4_apply (false,nmo,v2,v2); sf_floatwrite (v2,nt,vely);
		    sf_stretch4_apply (false,nmo,v3,v3); sf_floatwrite (v3,nt,velxy);
		} else {
		    sf_int2_init (coord1, t0,vx0, dt,dvx, nt,nvx, sf_spline_int, nw, nt);
		    sf_int2_lop (true,true,ntvx,nt,vx,ord);

		    sf_int2_init (coord2, t0,vy0, dt,dvy, nt,nvy, sf_spline_int, nw, nt);
		    sf_int2_lop (true,true,ntvy,nt,vy,ord);

		    sf_int2_init (coord3, t0,vxy0, dt,dvxy, nt,nvxy, sf_spline_int, nw, nt);
		    sf_int2_lop (true,true,ntvy,nt,vxy,ord);
		}
	    } /* END slope p1 loop */
        } /* END slope p2 loop */
        
	if (!map) {
	    /* from spline coefficients to model */
	    if (nw > 2) { /* loop spline */
		for (iv=0; iv < nvx; iv++) {
		    sf_spline_post (nw, iv*nt, 1, nt, vx, vx2);
		}
		for (iv=0; iv < nvy; iv++) {
		    sf_spline_post (nw, iv*nt, 1, nt, vy, vy2);
		}
		
		for (iv=0; iv < nvxy; iv++) {
		    sf_spline_post (nw, iv*nt, 1, nt, vxy, vxy2);
		}
		
		for (it=0; it < nt; it++) {
		    sf_spline_post (nw, it, nt, nvx, vx2, vx);
		    sf_spline_post (nw, it, nt, nvy, vy2, vy);
		    sf_spline_post (nw, it, nt, nvxy, vxy2, vxy);
		}
	    } /* END loop spline */
        
	    sf_floatwrite (vx,ntvx,velx);
	    sf_floatwrite (vy,ntvy,vely);
	    sf_floatwrite (vxy,ntvxy,velxy);
	}
    } /* END CMP loop */
 
    exit (0);
} /* END of main */
