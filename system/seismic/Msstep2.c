/* 3-D post-stack modeling/migration with extended split step. */
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
#include "split2.h"

int main (int argc, char *argv[])
{
    int nw;		/* number of frequencies */
    int nz;		/* number of migrated time samples */
    int nx, ny;		/* number of midpoints 	*/
    int nt, ntx, nty;   /* boundary taper size */
    int padx, pady;     /* padding */
    int nr;             /* number of reference velocities */

    float w0;		/* frequency origin 	*/
    float z0, dz;	/* migrated time sampling interval */
    float dw;	        /* frequency sampling interval */
    float dx,dy;	/* spatial sampling interval	*/
    float dt;           /* time error */

    bool inv;             /* modeling or migration        */
    bool verb;            /* verbosity */
    float eps;            /* dip filter constant          */
    sf_file in=NULL, out=NULL, vel=NULL;
    sf_slice imag, slow;

    sf_init(argc,argv);
    in  = sf_input ("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv = false;
    /* If y, modeling; if n, migration */
    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */
    if (!sf_getfloat("eps",&eps)) eps = 0.01;
    /* stability parameter */

    if (!sf_histint  (in,"n1",&ny)) ny = 1;
    if (!sf_histfloat(in,"d1",&dy)) sf_error ("No d1= in input");

    if (!sf_histint  (in,"n2",&nx)) nx = 1;
    if (!sf_histfloat(in,"d2",&dx)) dx=dy;

    if (!sf_getint("nt",&nt)) nt = 1;
    /* taper size */
    ntx = SF_MIN(nt,nx-1);
    nty = SF_MIN(nt,ny-1);

    if (!sf_getint("nr",&nr)) nr = 1;
    /* maximum number of references */

    if (!sf_getint("padx",&padx)) padx = 0;
    /* cross-line padding */
    if (!sf_getint("pady",&pady)) pady = 0;
    /* in-line padding */

    if (!sf_getfloat("dt",&dt)) dt=0.004;
    /* time error */

    vel = sf_input("slowness");

    if (inv) { /* modeling */
	if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
	sf_settype(out,SF_COMPLEX);

	if (!sf_histint  (in,"n3",&nz)) sf_error ("No n3= in input");
	if (!sf_histfloat(in,"d3",&dz)) sf_error ("No d3= in input");
	if (!sf_getint   (   "nw",&nw)) sf_error ("Need nw=");
	/* Length of frequency axis (for modeling) */ 
	if (!sf_getfloat (   "dw",&dw)) sf_error ("Need dw=");
	/* Frequency sampling (for modeling) */
	if (!sf_getfloat (   "w0",&w0)) w0=0.;
	/* Frequency origin (for modeling) */
	sf_putint  (out,"n3",nw);
	sf_putfloat(out,"d3",dw);
	sf_putfloat(out,"o3",w0);
    } else { /* migration */
	if (SF_COMPLEX != sf_gettype(in)) sf_error("Need complex input");
	sf_settype(out,SF_FLOAT);

	if (!sf_histint  (in, "n3",&nw)) sf_error ("No n3= in input");
	if (!sf_histfloat(in, "d3",&dw)) sf_error ("No d3= in input");
	if (!sf_histfloat(in, "o3",&w0)) sf_error ("No o3= in input");
	if (!sf_histint  (vel,"n3",&nz)) sf_error ("No n3= in slowness");
	if (!sf_histfloat(vel,"d3",&dz)) sf_error ("No d3= in slowness");

	if (!sf_histfloat(vel,"o3",&z0)) z0=0.; 
	sf_putint  (out,"n3",nz);
	sf_putfloat(out,"d3",dz);
	sf_putfloat(out,"o3",z0);
    }
    /* from hertz to radian */
    dw *= 2.*SF_PI; 
    w0 *= 2.*SF_PI;

    /* allocate space for slowness and image */
    slow = sf_slice_init(vel,ny,nx,nz);
    imag = sf_slice_init(inv? in:out,ny,nx,nz);

    /* initialize split-step */
    split2_init(nz,dz,ny,dy,nx,dx,ntx,nty,padx,pady,nr,dt);
    split2 (verb, inv, eps,  nw, dw, w0, inv? out:in, imag, slow);

    exit (0);
}
