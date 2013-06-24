/* 1D Source for the method of manufactured solution */
/*
  Copyright (C) 2008 University of Texas at Austin
  
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

int main(int argc, char* argv[])
{
    /*grid parameters*/
    int nx, nt;
    float dx, dt;
    int ix, it;
    
    sf_axis at, ax;

    float *vel;          /*input velocity*/
    float **src, **slt;  /*output MMS source and solution*/
    
    /*parameters*/
    float slx;
    float alpha, beta, beta2;
    float dist2, xx, tt, tmp;
    float cosbdt2;
    
    sf_file Fvel, Fsrc, Fslt;
    sf_init(argc, argv);
    Fvel = sf_input("in");
    Fsrc = sf_output("out");
    Fslt = sf_output("mslt"); /*the manufactured solution*/
    
    if (SF_FLOAT != sf_gettype(Fvel)) sf_error("Need float input");
    ax = sf_iaxa(Fvel, 1); nx = sf_n(ax); dx = sf_d(ax); 
    

    if (!sf_getint("nt",&nt)) sf_error("Need nt");
    /*number of time step*/
    if (!sf_getfloat("dt", &dt)) sf_error("Need dt");
    /*time step*/
    if (!sf_getfloat("slx", &slx)) slx = nx*dx*0.5;
    /*center of source location: x*/
    if (!sf_getfloat("alpha", &alpha)) alpha = 1.0e-2;
    /*source parameter*/
    if (!sf_getfloat("beta", &beta)) beta = 1.0;
    /*source parameter*/
    
    at = sf_maxa(nt, 0.0, dt);
        
    vel = sf_floatalloc(nx);
    src = sf_floatalloc2(nx, nt);
    slt = sf_floatalloc2(nx, nt);
    
    sf_floatread(vel, nx, Fvel);

    /*Set output axis*/
    sf_oaxa(Fsrc, ax, 1);
    sf_oaxa(Fsrc, at, 2);
    
    sf_oaxa(Fslt, ax, 1);
    sf_oaxa(Fslt, at, 2);

    /*Manufactured Solution Source*/
    beta2 = beta*beta;
    cosbdt2 = cosf(beta*dt*0.5);
    for (it=0; it<nt; it++) {
	//tt = sinf(it*dt*beta);
	tt = (1.0*it+0.5)*dt;
	tt = cosf(tt*beta);
	//tt = expf(-1*(it*dt-0.15)*beta);
	for (ix=0; ix<nx; ix++) {
	    xx = ix*dx;
	    dist2 = (xx-slx)*(xx-slx);
	    tmp = expf(-1*alpha*dist2)/beta;
	    //src[it][ix] = -1*tmp*(beta2*tt+2*alpha*vel[ix]*vel[ix]*(2*alpha*dist2-1)*(tt)); 
	    src[it][ix] = tmp*(beta2*tt+2*alpha*vel[ix]*vel[ix]*(2*alpha*dist2-1)*(tt-cosbdt2)); 
	    //src[it][ix] = tmp*(-1*beta2+4*alpha*vel[ix]*vel[ix]*(alpha*dist2-1)); 
	}
	sf_floatwrite(src[it], nx, Fsrc);
    }

    /*solution*/
    for (it=0; it<nt; it++) {
	//tt = cosf(it*dt*beta);
	tt =sinf(it*dt*beta);
	//tt = expf(-1*(it*dt-0.15)*beta);
	for (ix=0; ix<nx; ix++) {
	    xx = ix*dx;
	    dist2 = (xx-slx)*(xx-slx);
	    slt[it][ix] = expf(-1*alpha*dist2)*tt;
	}
	sf_floatwrite(slt[it], nx, Fslt);
    }
    
    
}
					  
