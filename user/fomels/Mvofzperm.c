/* V(z) prestack exploditing reflector. */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

#include "cosft3.h"
#include "lowrank1.h"

int main(int argc, char* argv[])
{
    bool mig, stat;
    int it, nt, ix, nx=0, iz, nz=0, ih, nh, it1, it2, its, ik, n1, n2, nk, nr, **siz;
    float dt, dx, dz, dh, c;
    float ***curr, **img, **dat, ***lft, ***mid, ***rht;
    sf_file data, image, size, left, middle, right;

    sf_init(argc,argv);

    if (!sf_getbool("mig",&mig)) mig=false;
    /* if n, modeling; if y, migration */

    if (mig) { /* migration */
	data = sf_input("in");
	image = sf_output("out");

	stat = sf_histint(data,"n1",&nh);
	stat = stat && sf_histfloat(data,"d1",&dh);

	stat = stat && sf_histint(data,"n2",&nx);
	stat = stat && sf_histfloat(data,"d2",&dx);

	stat = stat && sf_histint(data,"n3",&nt);
	stat = stat && sf_histfloat(data,"d3",&dt);

	if (!stat) sf_error("Need n1,n2,etc. in input");

	if (!sf_getint("nz",&nz)) sf_error("Need nz=");
	/* time samples (if migration) */
	if (!sf_getfloat("dz",&dz)) sf_error("Need dz=");
	/* time sampling (if migration) */

	sf_putint(image,"n1",nz);
	sf_putfloat(image,"d1",dz);
	sf_putstring(image,"label1","Depth");

	sf_putint(image,"n3",1); /* stack for now */
    } else { /* modeling */
	image = sf_input("in");
	data = sf_output("out");

	stat = sf_histint(image,"n1",&nz);
	stat = stat && sf_histfloat(image,"d1",&dz);

	stat = stat && sf_histint(image,"n2",&nx);
	stat = stat && sf_histfloat(image,"d2",&dx);

	if (!stat) sf_error("Need n1,n2,etc. in input");

	if (!sf_getint("nt",&nt)) sf_error("Need nt=");
	/* time samples (if modeling) */
	if (!sf_getfloat("dt",&dt)) sf_error("Need dt=");
	/* time sampling (if modeling) */

	if (!sf_getint("nh",&nh)) sf_error("Need nh=");
        /* offset samples (if modeling) */
	if (!sf_getfloat("dh",&dh)) sf_error("Need dh=");
	/* offset sampling (if modeling) */

	sf_putint(data,"n1",nh);
	sf_putfloat(data,"d1",dh);
	sf_putstring(data,"label1","Half-Offset");

	sf_putint(data,"n3",nt);
	sf_putfloat(data,"d3",dt);
	sf_putstring(data,"label3","Time");
	sf_putstring(data,"unit3","s");
    }

    img = sf_floatalloc2(nz,nx);
    dat = sf_floatalloc2(nh,nx);

    curr = sf_floatalloc3(nz,nh,nx);

    for (ix=0; ix < nx; ix++) {
	for (ih=0; ih < nh; ih++) {
	    for (iz=0; iz < nz; iz++) {
		curr[ix][ih][iz] = 0.;
	    }
	}
    }

    nk = nx*nh;

    size = sf_input("size");
    if (SF_INT != sf_gettype(size) ||
	!sf_histint(size,"n1",&n1) || 2  != n1 ||
	!sf_histint(size,"n2",&n2) || nk != n2) sf_error("Wrong size in size=");
    siz = sf_intalloc2(2,nk);
    sf_intread(siz[0],2*nk,size);
    sf_fileclose(size);

    left =  sf_input("left");
    if (SF_FLOAT != sf_gettype(left) ||
	!sf_histint(left,"n1",&n1) || nz != n1 ||
	!sf_histint(left,"n3",&n2) || nk != n2) sf_error("Wrong size in left=");
    if (!sf_histint(left,"n2",&nr)) sf_error("No n2= in left");
    lft = sf_floatalloc3(nz,nr,nk);
    sf_floatread(lft[0][0],nz*nr*nk,left);
    sf_fileclose(left);

    middle =  sf_input("middle");
    if (SF_FLOAT != sf_gettype(middle) ||
	!sf_histint(middle,"n1",&n1) || nr != n1 ||
	!sf_histint(middle,"n3",&n2) || nk != n2) sf_error("Wrong size in middle=");
    mid = sf_floatalloc3(nr,nr,nk);
    sf_floatread(mid[0][0],nr*nr*nk,middle);
    sf_fileclose(middle);

    right =  sf_input("right");
    if (SF_FLOAT != sf_gettype(right) ||
	!sf_histint(right,"n2",&n1) || nz != n1 ||
	!sf_histint(right,"n3",&n2) || nk != n2) sf_error("Wrong size in right=");
    rht = sf_floatalloc3(nr,nz,nk);
    sf_floatread(rht[0][0],nr*nz*nk,right);
    sf_fileclose(right);

    if (mig) { /* migration */
	/* initialize image */
	for (iz=0; iz < nz; iz++) {
	    for (ix=0; ix < nx; ix++) {
		img[ix][iz] = 0.;
	    }
	}

	/* step backward in time */
	it1 = nt-1;
	it2 = -1;
	its = -1;

    } else { /* modeling */
	sf_floatread(img[0],nz*nx,image);

	/* Initialize model */
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		curr[ix][0][iz] = img[ix][iz];
	    }
	}
	cosft12(false,nz,nh,nx,curr);
	
	/* step forward in time */
	it1 = 0;
	it2 = nt;
	its = +1;
    }

    lowrank1_init(nz,nr);

    /* time stepping */
    for (it=it1; it != it2; it += its) {
	sf_warning("it=%d;",it);

	if (mig) { /* migration <- read data */
	    sf_floatread(dat[0],nh*nx,data);
	    cosft2(false,nh,nx,dat);	    
	} else {
	    for (ix=0; ix < nx; ix++) {
		for (ih=0; ih < nh; ih++) {
		    dat[ix][ih] = 0.;
		}
	    }
	}

	sf_cosft_init(nz);
	
	ik=0;
	for (ix=0; ix < nx; ix++) {
	    for (ih=0; ih < nh; ih++) {		
		if (mig) {
		    curr[ix][ih][0] += dat[ix][ih];
		} else {
		    dat[ix][ih] += curr[ix][ih][0];
		}
		
		lowrank1_step(siz[ik][0],siz[ik][1],lft[ik],mid[ik],rht[ik],curr[ix][ih]);
		ik++;
	    }
	}

	sf_cosft_close();

	if (!mig) { /* modeling -> write out data */
	    cosft2(true,nh,nx,dat);
	    sf_floatwrite(dat[0],nh*nx,data);
	}
    }
    sf_warning(".");

    if (mig) {
	for (ix=0; ix < nx; ix++) {
	    for (ih=0; ih < nh; ih++) {
		for (iz=1; iz < nz; iz++) {
		    c = curr[ix][ih][iz];
		    img[ix][iz] += (iz==nz-1)? c*0.5: c;
		}
	    }
	}

	cosft2(true,nz,nx,img);
	sf_floatwrite(img[0],nz*nx,image);
    }

    exit(0);
}
