/* 3-D common-azimuth prestack modeling/migration with split-step DSR. */
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

#include "cam.h"

int main (int argc, char *argv[])
{
    bool inv;             /* modeling or migration */
    bool verb;            /* verbosity */
    float eps;            /* dip filter constant   */  

    int nt,ntx,nty,nth;   /* boundary taper size */
    int nr;               /* number of reference velocities */
    int npad;             /* padding on offset wavenumber */
    float dt;             /* time error */

    slice imag;
    slice slow;
    slice data;

    axa az,ay,ax,aw,au,av,ah;
    sf_file Fi; /* image    file */
    sf_file Fd; /* data     file */
    sf_file Fs; /* slowness file */

    /*------------------------------------------------------------*/
    sf_init(argc,argv);

    if (!sf_getbool( "inv",&inv ))  inv = false; /* y=modeling; n=migration */
    if (!sf_getbool("verb",&verb)) verb = false; /* verbosity flag */
    if (!sf_getfloat("eps",&eps ))  eps =  0.01; /* stability parameter */
    if (!sf_getint  ( "nr",&nr  ))   nr =     1; /* maximum number of references */
    if (!sf_getfloat( "dt",&dt  ))   dt = 0.004; /* time error */
    if (!sf_getint( "npad",&npad)) npad =     0; /* padding on offset wavenumber */
    if (!sf_getint(   "nt",&nt  ))   nt =     1; /* taper size */
    
    Fi = inv ? sf_input ( "in"): sf_output("out");
    Fd = inv ? sf_output("out"): sf_input ( "in"); 
    Fs = sf_input ("slowness");
    
    /*     data[nh][ny][nx][nw] */
    /*    image[nh][nu][nv][nz] */
    /* slowness    [nu][nv][nz] */
    if (inv) { /* modeling */
	if (SF_FLOAT != sf_gettype(Fi)) sf_error("Need float image");
	sf_settype(Fd,SF_COMPLEX);

	if (!sf_getint  ("nw",&aw.n)) sf_error ("Need nw=");
	if (!sf_getfloat("dw",&aw.d)) sf_error ("Need dw=");
	if (!sf_getfloat("w0",&aw.o)) aw.o=0.;

	iaxa(Fi,&ah,1);        oaxa(Fd,&ah,1);
	iaxa(Fi,&au,2); ay=au; oaxa(Fd,&ay,2);
	iaxa(Fi,&av,3); ax=av; oaxa(Fd,&ax,3);
	iaxa(Fi,&az,4);        oaxa(Fd,&aw,4);

    } else { /* migration */
	if (SF_COMPLEX != sf_gettype(Fd)) sf_error("Need complex data");
	sf_settype(Fi,SF_FLOAT);

	iaxa(Fd,&ah,1); oaxa(Fi,&ah,1);
	iaxa(Fd,&ay,2);
	iaxa(Fd,&ax,3);
	iaxa(Fd,&aw,4);

	iaxa(Fs,&au,1); oaxa(Fi,&au,2);
	iaxa(Fs,&av,2); oaxa(Fi,&av,3);
	iaxa(Fs,&az,3); oaxa(Fi,&az,4);
    }

    /* taper */
    ntx = SF_MIN(nt,ax.n-1);
    nty = SF_MIN(nt,ay.n-1);
    nth = SF_MIN(nt,ah.n-1);

    /* from hertz to radian */
    aw.d *= 2.*SF_PI; 
    aw.o *= 2.*SF_PI;

    slow = slice_init(Fs,     au.n,av.n,az.n);
    imag = slice_init(Fi,ah.n,au.n*av.n,az.n);
    data = slice_init(Fd,ah.n,ay.n*ax.n,aw.n);

    cam_init(verb,eps,dt,
	     az,aw,ah,ay,ax,au,av,
	     ntx,nty,nth,nr,npad);
    cam      (inv,data,imag,slow);
    cam_close();
    
    exit (0);
}
