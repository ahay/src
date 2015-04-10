/* 1D method of manufactured solution using Gaussian pulsa*/
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


static float *vel, *dvel, *velhf, *dvelhf, *den, *denhf;
static float alpha;
static float dt, dx, slx;

float gauss(int ix, float xx, float tt, bool hf)
{
    float tmp;
    if (hf) {
	tmp = xx-slx-velhf[ix]*tt;
    } else {
	tmp = xx-slx-vel[ix]*tt;
    }
    return expf(-1*alpha*tmp*tmp);
}


int main(int argc, char* argv[])
{

    /*I/O*/
    sf_file Fvel, Fdvel, Fden;
    sf_file Fvelhf, Fdvelhf, Fdenhf;
    sf_file Fpsrc, Fvsrc, Fpint, Fvint, Fmms;
    
    /*I/O array*/
    float **psrc, **vsrc, *pint, *vint, **mms; 

    /*grid parameters*/
    int nx, nt;
    int ix, it;
    
    sf_axis at, ax;

    /*parameters*/
    
    float gg, dist;
    float xx, tt;
    
    bool hf = true;
    bool mg = false;
  
    sf_init(argc, argv);
    Fvel  = sf_input("in");
    Fdvel = sf_input("dvel");
    Fden  = sf_input("den");
    Fvelhf= sf_input("velhf");
    Fdvelhf=sf_input("dvelhf");
    Fdenhf =sf_input("denhf");
    
    Fpsrc = sf_output("presrc");
    Fvsrc = sf_output("velsrc");
    Fpint = sf_output("preinit");
    Fvint = sf_output("velinit");
    Fmms  = sf_output("out");
    
    if (SF_FLOAT != sf_gettype(Fvel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(Fdvel)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(Fden)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(Fvelhf)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(Fdvelhf)) sf_error("Need float input");
    if (SF_FLOAT != sf_gettype(Fdenhf)) sf_error("Need float input");
    
    ax = sf_iaxa(Fvel, 1); nx = sf_n(ax); dx = sf_d(ax); 

    /*parameters*/
    if (!sf_getint("nt",&nt)) sf_error("Need nt");
    /*number of time step*/
    if (!sf_getfloat("dt", &dt)) sf_error("Need dt");
    /*time step*/
    if (!sf_getfloat("slx", &slx)) slx = nx*dx*0.5;
    /*center of source location: x*/
    if (!sf_getfloat("alpha", &alpha)) alpha = 1.0e-2;
    /*source parameter*/
    
    
    /*set axis*/
    at = sf_maxa(nt, 0.0, dt);

    /*set output axis*/
    sf_oaxa(Fpsrc, ax, 1);
    sf_oaxa(Fpsrc, at, 2);
    sf_oaxa(Fvsrc, ax, 1);
    sf_oaxa(Fvsrc, at, 2);

    sf_oaxa(Fpint, ax, 1);
    sf_oaxa(Fvint, ax, 1);

    sf_oaxa(Fmms, ax, 1);
    sf_oaxa(Fmms, at, 2);

    /*allocate memory*/
    vel  = sf_floatalloc(nx);
    velhf= sf_floatalloc(nx);
    dvel = sf_floatalloc(nx);
    dvelhf=sf_floatalloc(nx);
    den  = sf_floatalloc(nx);
    denhf =sf_floatalloc(nx);
    
    psrc = sf_floatalloc2(nx, nt);
    vsrc = sf_floatalloc2(nx, nt);
    mms  = sf_floatalloc2(nx, nt);
    
    pint = sf_floatalloc(nx);
    vint = sf_floatalloc(nx);

    /* read */
    sf_floatread(vel, nx, Fvel);
    sf_floatread(dvel, nx, Fdvel);
    sf_floatread(den, nx, Fden);
    sf_floatread(velhf, nx, Fvelhf);
    sf_floatread(dvelhf, nx, Fdvelhf);
    sf_floatread(denhf, nx, Fdenhf);

    /* velocity source */
    for (it=0; it<nt; it++) {
	tt = it*dt; /*time derivative is caculted on main grid */
	for (ix=0; ix<nx; ix++) {
	    xx = (ix+0.5)*dx;
	    dist = xx - slx;
	    gg = gauss(ix, xx, tt, hf);    
	    vsrc[it][ix] = 2.0*alpha*(dist-velhf[ix]*tt)*(velhf[ix]+(dvelhf[ix]*tt-1)/denhf[ix])*gg;
	}
	sf_floatwrite(vsrc[it], nx, Fvsrc);
    }

    /* pressure source */
    for (it=0; it<nt; it++) {
	tt = (it+0.5)*dt; /*time derivative is caculted on half grid */
	for (ix=0; ix<nx; ix++) {
	    xx = ix*dx;
	    dist = xx - slx;
	    gg = gauss(ix, xx, tt, mg);
	    psrc[it][ix] = 2.0*alpha*(dist-vel[ix]*tt)*(vel[ix]+den[ix]*vel[ix]*vel[ix]* \
							(dvel[ix]*tt-1))*gg;
	}
	sf_floatwrite(psrc[it], nx, Fpsrc);
    }
    
    /* initial condtion*/
    /* velocity:  U(x, -dt/2)*/
    tt = -0.5*dt;
    for (ix=0; ix<nx; ix++) {
	xx = (ix+0.5)*dx;
	vint[ix] = gauss(ix, xx, tt, hf);
    }
    sf_floatwrite(vint, nx, Fvint);
    
    /*pressure: P(x, 0)*/
    for (ix=0; ix<nx; ix++) {
	xx = (ix)*dx;
	dist = xx - slx;
	dist = dist*dist;
	pint[ix] = expf(-1*alpha*dist);
    }
    sf_floatwrite(pint, nx, Fpint);

    /* Manufactured solution: P(x, t), pressure on main grid*/
    for (it=0; it<nt; it++) {
	tt = it*dt;
	for (ix=0; ix<nx; ix++) {
	    xx = ix*dx;
	    mms[it][ix] = gauss(ix, xx, tt, mg);
	}
	sf_floatwrite(mms[it], nx, Fmms);
    }
    exit(0);

    
}
					  
