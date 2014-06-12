/* 2-D zero-offset reverse-time migration
 */
/*
  Copyright (C) 2014  Xi'an Jiaotong University, UT Austin (Pengliang Yang)

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

#include "rtm2d.h"

int main(int argc, char* argv[])
{
    	bool adj;    
    	int n1, n2, nb, nt, n0, nx, nz;
    	float dt, dx, dz, o1, o2;
    	float *mod, *dat, **v0;      

    	sf_file data, imag, modl;/* I/O files */

    	/* initialize Madagascar */
    	sf_init(argc,argv);

    	modl = sf_input ("vel");/* velocity model */
    
    	if (!sf_histint(modl,"n1",&n1)) sf_error("n1");
	/* 1st dimension size */
    	if (!sf_histint(modl,"n2",&n2)) sf_error("n2");
	/* 2nd dimension size */
    	if (!sf_histfloat(modl,"d1",&dz)) sf_error("d1");
	/* d1 */
    	if (!sf_histfloat(modl,"d2",&dx)) sf_error("d2");
	/* d2 */
    	if (!sf_histfloat(modl,"o1",&o1)) sf_error("o1");
	/* o1 */
    	if (!sf_histfloat(modl,"o2",&o2)) sf_error("o2");
	/* o2 */
    	if (!sf_getbool("adj",&adj)) adj=false;
	/* if y, migration; else, modeling */
    	if (!sf_getint("nb",&nb)) nb=20;
	/* number (thickness) of ABC boundary grid on each side */
    	if (!sf_getint("n0",&n0)) n0=0;
	/* shot depth in the grid */

	if(adj){/* migration */
	    data = sf_input ("in"); /* seismic data */
	    imag = sf_output("out");  /* output image */
	    if (!sf_histint(data,"n1",&nt)) sf_error("n1");
	    /* number of time steps */
	    if (!sf_histfloat(data,"d1",&dt)) sf_error("d1");
	    /* time sampling interval: dt */
	    if (!sf_histint(data,"n2",&nx) || nx != n2) 
		sf_error("Need n2=%d in data",n2);
	    sf_putint(imag,"n1",n1);
	    sf_putint(imag,"n2",n2);
	    sf_putfloat(imag,"d1",dz);
	    sf_putfloat(imag,"d2",dx);
	    sf_putfloat(imag,"o1",o1);
	    sf_putfloat(imag,"o2",o2);
	    sf_putstring(imag,"label1","Depth");
	    sf_putstring(imag,"label2","Distance");
	}else{/* modeling */
	    imag = sf_input ("in"); /* input image */
	    data = sf_output("out");  /* output seismic data */
	    if (!sf_histint(imag,"n1",&nz) || nz != n1)
		    sf_error("Need n1=%d in imag",n1);
	    if (!sf_histint(imag,"n2",&nx) || nx != n2) 
		    sf_error("Need n2=%d in imag",n2);
	    if (!sf_getint("nt",&nt)) sf_error("nt");
	    /* number of time steps */
	    if (!sf_getfloat("dt",&dt)) sf_error("dt");
	    /* time sampling interval: dt */
	    sf_putint(data, "n1", nt);
	    sf_putint(data, "n2", n2);
	    sf_putfloat(data,"d1",dt);
	    sf_putfloat(data,"d2",dx);
	    sf_putfloat(data,"o1",0);
	    sf_putfloat(data,"o2",o2);
	    sf_putstring(data,"label1","Time");
	    sf_putstring(data,"label2","Distance");
	}

	/* In rtm, v0 is the velocity model [modl], which is input parameter; 
	   mod is the image/reflectivity [imag]; dat is seismogram [data]! */
    	v0 = sf_floatalloc2(n1,n2);
    	mod = sf_floatalloc(n1*n2);
    	dat = sf_floatalloc(nt*n2);

    	sf_floatread(v0[0],n1*n2,modl);
    	if(adj){/* migration */
		sf_floatread(dat,nt*n2,data);
    	}else{ /* modeling */
		sf_floatread(mod,n1*n2,imag);
    	}

	rtm2d_init(dz, dx, dt, n0, n1, n2, nb, nt, v0, mod, dat);
	rtm2d_lop(adj, false, n1*n2, nt*n2, mod, dat);
	rtm2d_close();

    	if(adj) sf_floatwrite(mod, n1*n2, imag);/*output image */
    	else	sf_floatwrite(dat, nt*n2, data);/* output data */

	free(*v0); free(v0);
	free(mod);
	free(dat);    
    
	exit (0);
}
