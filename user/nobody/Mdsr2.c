/* 2-D prestack modeling/migration with split-step DSR. */
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

#include "dsr2.h"

int main (int argc, char *argv[])
{
    int nz;		/* depth samples */
    int ny;             /* lateral samples */
    int nx;		/* number of midpoints 	*/
    int nh;             /* number of offsets */
    int nw;		/* number of frequencies */	
    int nt, ntx, nth;   /* boundary taper size */
    int nr;             /* number of reference velocities */
    int npad;           /* padding on offset wavenumber */

    float z0, dz;	/* depth origin, sampling interval */
    float w0, dw;	/* frequency origin, sampling interval */
    float x0, dx;	/* midpoint origin, sampling interval	*/
    float y0, dy;       /* spatial origin, sampling interval	*/
    float h0, dh;       /* offset origin, sampling interval */
    float dt;           /* time error */

    float          **slow;

    bool inv;             /* modeling or migration        */
    bool verb;            /* verbosity */
    float eps;            /* dip filter constant          */  
    sf_file in, out, vel;
    slice imag;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    vel = sf_input("slowness");

    if (!sf_getbool("inv",&inv)) inv = false;
    /* If y, modeling; if n, migration */
    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    if (!sf_getfloat("eps",&eps)) eps = 0.01;
    /* stability parameter */

    if (!sf_histint  (in,"n1",&nh)) nh = 1;
    if (!sf_histfloat(in,"d1",&dh)) sf_error ("No d1= in input");
    if (!sf_histfloat(in,"o1",&h0)) h0=0.;

    if (!sf_histint  (in,"n2",&nx)) nx = 1;
    if (!sf_histfloat(in,"d2",&dx)) sf_error ("No d2= in input");
    if (!sf_histfloat(in,"o2",&x0)) x0=0.;

    if (!sf_getint("nt",&nt)) nt = 1;
    /* taper size */
    ntx = SF_MIN(nt,nx-1);
    nth = SF_MIN(nt,nh-1);
    
    if (!sf_getint("nr",&nr)) nr = 1;
    /* maximum number of references */

    if (!sf_getfloat("dt",&dt)) dt=0.004;
    /* time error */

    if (!sf_getint("npad",&npad)) npad = 0;
    /* padding on offset wavenumber */

    if (inv) { /* modeling */
	if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
	sf_settype(out,SF_COMPLEX);

	if (!sf_histint  (in,"n3",&nz)) sf_error ("No n3= in input");
	if (!sf_histfloat(in,"d3",&dz)) sf_error ("No d3= in input");
	
	if (!sf_histint  (in,"n2",&ny)) sf_error ("No n2= in input");
	if (!sf_histfloat(in,"d2",&dy)) sf_error ("No d2= in input");
	if (!sf_histfloat(in,"o2",&y0)) sf_error ("No o2= in input");

	if (!sf_getint  ("nw",&nw)) sf_error ("Need nw=");
	/* Length of frequency axis (for modeling) */ 
	if (!sf_getfloat("dw",&dw)) sf_error ("Need dw=");
	/* Frequency sampling (for modeling) */
	if (!sf_getfloat("w0",&w0)) w0=0.;
	/* Frequency origin (for modeling) */
	sf_putint  (out,"n3",nw);
	sf_putfloat(out,"d3",dw);
	sf_putfloat(out,"o3",w0);
    } else { /* migration */
	if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
	sf_settype(out,SF_FLOAT);

	if (!sf_histint  (in,"n3",&nw)) sf_error ("No n3= in input");
	if (!sf_histfloat(in,"d3",&dw)) sf_error ("No d3= in input");
	if (!sf_histfloat(in,"o3",&w0)) sf_error ("No o3= in input");

	if (!sf_histint  (vel,"n2",&nz)) sf_error ("No n2= in slowness");
	if (!sf_histfloat(vel,"d2",&dz)) sf_error ("No d2= in slowness");
	if (!sf_histfloat(vel,"o2",&z0)) z0=0.; 

	if (!sf_histint  (vel,"n1",&ny)) sf_error ("No n1= in slowness");
	if (!sf_histfloat(vel,"d1",&dy)) sf_error ("No d1= in slowness");
	if (!sf_histfloat(vel,"o1",&y0)) y0=0.; 

	sf_putint  (out,"n3",nz);
	sf_putfloat(out,"d3",dz);
	sf_putfloat(out,"o3",z0);
    }
    /* from hertz to radian */
    dw *= 2.*SF_PI; 
    w0 *= 2.*SF_PI;

    slow = sf_floatalloc2(ny,nz);
    sf_floatread(slow[0],ny*nz,vel);
    sf_fileclose(vel);

    imag = slice_init( inv ? in:out,nh,ny,nz);

    dsr2_init (nz,dz, 
	       nh,dh,h0, 
	       nx,dx,x0, 
	       ny,dy,y0, 
	       ntx,nth,
	       nr,
	       npad);
    dsr2      (verb, inv, eps,  
	       nw, dw, w0,
	       inv ? out:in, 
	       imag, 
	       slow, 
	       dt);
    dsr2_close();

    exit (0);
}
