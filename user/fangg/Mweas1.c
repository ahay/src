/* 1-D analytic solution for acoustic wave equation */
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
#include <math.h>


float trape( float *f,
	     int n,
	     float d
             )
{
    int i;
    float sum = 0.;
    float w = .5*d;
    for(i=0; i<n-1; i++) {
	sum = sum + w*(f[i]+f[i+1]);
    }
 
    return sum;
}


float simpson(float *f, int n, float d) 
{
    int i;
    float we, wo;
    float sum0, sum1, sum2;
    we = 4./3.;
    wo = 1./3.;
    sum0 = 0.;
    sum1 = 0.;
    sum2 = 0.;
    
    sum1 = sum0 + (3./8*f[0]
		   + (9./8-1./3)*f[1]
		   + (9./8-4./3)*f[2]
		   + (3./8-1./3)*f[3])*d;
    for (i=2; i<n-1; i++) {
	sum2 = sum0 + we*f[i-1] + wo*(f[i] + f[i-2]);
	sum0 = sum1;
	sum1 = sum2;
    }
    
    return sum2;
    
}

int main(int argc, char* argv[])
{
    int nt, nx, ix, it;
    float dt, dx, *f, *fker, *u;
    float vel;
    float t, tau, dist;
    float x, ox, x0;
    int spx; //source point in x
    int kt; 
    char *rule;
    sf_file Fin, Fout;
    sf_axis at;

    
    sf_init(argc, argv);
    Fin = sf_input("in");
    
    if (SF_FLOAT != sf_gettype(Fin) ) sf_error("Need type=float in input");
    
    rule = sf_getstring("rule");
    /* t, s : quadrature rules */ 
    if (!sf_getfloat("vel", &vel)) sf_error("Need vel");
    if (!sf_getint("nx", &nx)) sf_error("Need nx");
    if (!sf_getfloat("dx", &dx)) sf_error("Need dx");
    if (!sf_getint("spx", &spx)) sf_error("Need spx");
    /*source point in x*/
    if (!sf_getint("kt", &kt)) sf_error("Need time");
    /*selected time*/
    if (!sf_getfloat("ox", &ox)) ox=0.0;

    Fout= sf_output("out");
    sf_putint(Fout,"n1",nx);
    sf_putfloat(Fout,"d1",dx);

    at = sf_iaxa(Fin, 1); nt = sf_n(at); dt = sf_d(at);
    
    f = sf_floatalloc(nt);
    fker = sf_floatalloc(nt);
    u = sf_floatalloc(nx);
    sf_floatread(f,nt,Fin);
    
    x0 = spx*dx;
    t   = kt*dt;
    sf_warning("=====================");
    sf_warning("nx=%d nt=%d dx=%f dt=%f kt=%d", nx, nt, dx, dt, kt);
    sf_warning("vel=%f ", vel);
    sf_warning("x0=%f t=%f", x0, t);
    
    
    /* Analytic solution*/
    for(ix=0; ix<nx; ix++) {
	x = ix*dx;
	dist = fabs(x-x0);
	if (dist == vel*t) {
	    u[ix] = f[0]*dt;
	} 
	else if(dist > vel*t) {
	    u[ix] = 0.0;
	} 
	else {
	    // integral kernel
	    for (it=0; it<nt; it++) {
		tau = it*dt;
		if (tau>=0 && tau<= t && dist<vel*(t-tau) ) {
		    fker[it] = f[it];
		} else {
		    fker[it] = 0.0;
		}
	    }
	}
	
	if (rule[0] == 't') {
	    u[ix] = trape(fker, nt, dt);
	} else {
	    u[ix] = simpson(fker, nt, dt);
	}
	//u[ix]=ix*dx;
    }
    
    sf_floatwrite(u,nx,Fout);
    
    return 0;
}    


    
    

    
