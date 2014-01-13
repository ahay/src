/* 2-D simplest-form post-stack Kirchhoff time modeling and migration. 
Suppose the input_image & output_data or input_data & output_image have the same dimensions, samplings.
The dottest has been past. */

/*
  Copyright (C) 2014 University of Texas at Austin
  
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

static int nt, nx, sw;
static float t0, dt, dx, *vrms;

void aa_kirch_init (float *vrms_in /* RMS velocity */, 
		    float t0_in    /* time origin */, 
		    float dt_in    /* time sampling */, 
		    float dx_in    /* midpoint sampling */, 
		    int nt_in      /* time samples */, 
		    int nx_in      /* midpoint samples */, 
		    int sw_in      /* branch to compute */)
/*< Initialize >*/
{
    vrms = vrms_in;
    t0 = t0_in;
    dt = dt_in;
    dx = dx_in;
    nt = nt_in;
    nx = nx_in;
	sw = sw_in;
}

void kirch_init (float *vrms_in /* RMS velocity */, 
		    float t0_in    /* time origin */, 
		    float dt_in    /* time sampling */, 
		    float dx_in    /* midpoint sampling */, 
		    int nt_in      /* time samples */, 
		    int nx_in      /* midpoint samples */)
/*< Initialize >*/
{
    vrms = vrms_in;
    t0 = t0_in;
    dt = dt_in;
    dx = dx_in;
    nt = nt_in;
    nx = nx_in;
}

void aa_kirch_lop (bool adj, bool add, int nm, int nd, 
		   float *modl, float *data)
/*< Antialiaising kirchhoff operator by reciprocal parameterization >*/
{
    int im, id1, id2, ix,iz,it,ib,iy, minx[2],maxx[2], is,i;
    float amp,t,z,b,db,f,g;

    sf_adjnull(adj,add,nm,nd,modl,data);

    maxx[0] = nx;	
    minx[1] = 1;
    for (iz=0; iz < nt-1; iz++) {     
	z = t0 + dt * iz;		/* vertical traveltime */
	for (it=nt-1; it >= iz+1; it--) { 
	    t = t0 + dt * it;		/* time shift */
	    b = sqrtf(t*t - z*z); 
	    db = dx*b*2./(vrms[iz]*t);
	    if(db < dt || sw == 1) break;

	    f = 0.5*vrms[iz]*b/dx; 
	    iy = f; f = f-iy; 
	    i = iy+1; g = 1.-f;

	    if(i >= nx)	continue;

	    amp = dt / db;

	    minx[0] = i;  
	    maxx[1] = nx-i;
	    for (is=0; is < 2; is++) {  
		iy = -iy; 
		i = -i;		/* two branches of hyperbola */
		for (ix=minx[is]; ix < maxx[is]; ix++) {
		    im = ix*nt+iz;
		    id1 = (ix+iy)*nt+it;
		    id2 = (ix+i)*nt+it;
		    if( adj) {	
			    modl[im] += amp*(data[id1]*g + data[id2]*f);
		    } else {
			    data[id1] += modl[im]*amp*g;
			    data[id2] += modl[im]*amp*f;
		    }
		}
	    }
	}

	for (ib=0; ib < nx; ib++) {	   
	    b = dx*ib*2./vrms[iz]; /* space shift */ 
	    iy = ib; 
	    t = hypotf(z,b); 
	    db = dx*b*2./(vrms[iz]*t);
	    if(db > dt || sw == 2) break;

	    f = (t-t0)/dt; 
	    it = f; f = f-it; 
	    i = it+1; g = 1.-f;
	    if(it >= nt) break;

	    amp = 1.;
	    if(ib == 0) amp *= 0.5;

	    minx[0] = iy; 
	    maxx[1] = nx-iy;
	    for (is=0; is < 2; is++) {	
		iy = -iy; /* two branches of hyperbola */
		for (ix=minx[is]; ix <	maxx[is]; ix++) {
		    im = ix*nt+iz;
		    id1 = (ix+iy)*nt+it;
		    id2 = (ix+iy)*nt+i;

		    if( adj) {
			    modl[im] += amp*(data[id1]*g + data[id2]*f);
		    } else {
			    data[id1] += modl[im]*amp*g;
			    data[id2] += modl[im]*amp*f;
		    }
		}
	    }
	}
    }
}


void kirch_lop (bool adj, bool add, int nm, int nd, 
		   float *modl, float *data)
/*< Simplest-form Kirchhoff operator >*/
{
    int im, id1, id2, ix,iz,it,ib,i;
    float t,z,b,f,g;

    sf_adjnull(adj,add,nm,nd,modl,data);


    for (iz=0; iz < nt-1; iz++) {     
		z = t0 + dt * iz;		/* vertical traveltime (t_0) */

		for (ib=0; ib < nx; ib++) {	   /* Loop over input trace */
	    	b = dx*ib*2./vrms[iz]; 
	    	t = hypotf(z,b); 

	    	f = (t-t0)/dt; 
	    	it = f; f = f-it; 
	    	i = it+1; g = 1.-f;
	    	if(it >= nt) break;


			for (ix=0; ix <	nx; ix++) { /* Loop over output trace */
		    im = ix*nt+iz;
		    id1 = ix*nt+it;
		    id2 = ix*nt+i;

		    	if( adj) {
			    	modl[im] += data[id1]*g + data[id2]*f;
		    	} else {
			    	data[id1] += modl[im]*g;
			    	data[id2] += modl[im]*f;
		    	}/* if */
	    	}/* ix */
		}/* ib */
    }/* iz */
}



int main(int argc, char* argv[])
{
    int n12, n1, n2, n3, i1, i3, sw;
    bool adj,aa;
    float **data, **modl, *vrms, o1,d1,o2,d2, v0;
    char *test;
    sf_file in, out, vel;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in,2);
    
    if (!sf_getbool("adj",&adj)) adj=true;
    /* yes: migration, no: modeling */


	if (!sf_getbool("aa",&aa)) aa=false;	
	/* yes: apply reciprocal antialiaising operator */
	if(aa)
	{    
		if (!sf_getint("sw",&sw)) sw=0;
   		/* if > 0, select a branch of the antialiasing operation */
	}
    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"o2",&o2)) sf_error("No o2= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input");

    vrms = sf_floatalloc(n1);

    if (NULL != (test = sf_getstring("velocity"))) { 
	/* velocity file */
	free(test);
	vel = sf_input("velocity");
	sf_floatread(vrms,n1,vel);
	sf_fileclose(vel);
    } else {
	if (!sf_getfloat("v0",&v0)) sf_error("Need velocity= or v0=");
	/* constant velocity (if no velocity=) */
	for (i1=0; i1 < n1; i1++) {
	    vrms[i1] = v0;
	}
    }

    n12 = n1*n2;
    data = sf_floatalloc2(n1,n2);
    modl = sf_floatalloc2(n1,n2);
	
	if(aa)
	{
    	aa_kirch_init (vrms, o1, d1, d2, n1, n2, sw);
	}else{
    	kirch_init (vrms, o1, d1, d2, n1, n2);
	}

    for (i3=0; i3 < n3; i3++) {
	if (adj) {
	    sf_floatread (data[0],n12,in);
	} else {
	    sf_floatread (modl[0],n12,in);
	}
	
	if(aa)
	{
		aa_kirch_lop (adj,false,n12,n12,modl[0],data[0]);
	}else{
		kirch_lop (adj,false,n12,n12,modl[0],data[0]);
	}

	if (adj) {
	    sf_floatwrite (modl[0],n12,out);
	} else {
	    sf_floatwrite (data[0],n12,out);
	}
    }


    exit(0);
}
