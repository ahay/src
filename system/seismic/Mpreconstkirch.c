/* Prestack Kirchhoff modeling/migration in constant velocity. 

Requires the input to be in (time,cmp x,cmp y,offset)
*/
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

#include <math.h>
#include <rsf.h>

int main(int argc, char* argv[])
{
    bool inv, zero, aal;
    int nt,nx,ny, nh, ix,iy,it,ih, ixin,iyin, n123, ix1, ix2, iy1, iy2;
    float dt,dx,dy, t0,x,y, vel0, rx,ry, dh, h0, h, hx, t, sq, ti, t1,t2;
    float *time=NULL, *str=NULL, *add=NULL, *tx=NULL, *amp=NULL, ***cinp=NULL, ***cout=NULL, ***stack=NULL;
    sf_file in=NULL, out=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (!sf_getbool("inv",&inv)) inv=false;
    /* if y, modeling; if n, migration */
    if (!sf_getbool("zero",&zero)) zero=false;
    /* if y, stack in migration */
    if (!sf_getbool("aal",&aal)) aal=true;
    /* if y, apply antialiasing */

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(in,"n3",&ny)) ny=1;
    n123 = nt*nx*ny;

    if (inv && zero) {
	if (!sf_getint ("nh",&nh)) sf_error("Need nh=");
	/* number of offsets */
	if (!sf_getfloat ("dh",&dh)) sf_error("Need dh=");
	/* offset sampling */
	if (!sf_getfloat ("h0",&h0)) sf_error("Need h0="); 
	/* offset origin */

	sf_putint(out,"n4",nh);
	sf_putfloat(out,"d4",dh);
	sf_putfloat(out,"o4",h0);
    } else {
	if (!sf_histint(in,"n4",&nh)) sf_error("No n4= in input");
	if (!sf_histfloat(in,"d4",&dh)) sf_error("No d4= in input");
	if (!sf_histfloat(in,"o4",&h0)) sf_error("No o4= in input");
    }

    if (!inv && zero) sf_putint(out,"n4",1);

    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"d3",&dy)) dy=dx;

    if (!sf_getfloat ("vel",&vel0)) sf_error("Need vel=");
    /* velocity */

    vel0 *= 0.5;
    dx /= vel0;
    dy /= vel0;
    dh /= vel0;
    h0 /= vel0;

    time = sf_floatalloc(nt);
    str = sf_floatalloc(nt);
    add = sf_floatalloc(nt);
    tx = sf_floatalloc(nt);
    amp = sf_floatalloc(nt);

    cinp = sf_floatalloc3(nt,nx,ny);
    cout = sf_floatalloc3(nt,nx,ny);

    if (!inv && zero) {
	stack = sf_floatalloc3(nt,nx,ny);
	for (it=0; it<n123; it++) {
	    stack[0][0][it] = 0.;
	}
    }

    for (it=0; it<nt; it++) {
	t = t0+it*dt;
	time[it] = t*t;
	tx[it] = 0.;
    }

    sf_aastretch_init (false, nt, t0, dt, nt);

    if (inv && zero) sf_floatread (cinp[0][0],n123,in);

    for (ih=0; ih < nh; ih++) {
	h = h0 + ih*dh;
	h *= h;
	sf_warning("offset %d of %d;",ih+1, nh);

	if (!inv || !zero) sf_floatread (cinp[0][0],n123,in);
	
	for (it=0; it<n123; it++) {
	    cout[0][0][it] = 0.;
	}

	for (iy = 1-ny; iy <= ny-1; iy++) {
	    if (iy > 0) {
		iy1 = 0;
		iy2 = ny-iy;
	    } else {
		iy1 = -iy;
		iy2 = ny;
	    }

	    y = iy*dy;
	    ry = fabsf(y*dy);
	    y *= y;
	    for (ix = 1-nx; ix <= nx-1; ix++) {
		if (ix > 0) {
		    ix1 = 0;
		    ix2 = nx-ix;
		} else {
		    ix1 = -ix;
		    ix2 = nx;
		}

		x = ix*dx;
		rx = fabsf(x*dx);
		x *= x;
		hx = h*x;
		x = x + y + h;
		for (it=0; it < nt; it++) {
		    if (inv) { /* modeling */
			t = time[it] + x;
			sq = t*t - 4.*hx;
			if (sq > 0. && t > 0.) {
			    sq = sqrtf(sq);
			    ti = sqrtf(0.5*(t + sq));
			    str[it] = ti;
			    t1 = rx / ti * 0.5 * (1 + (t - 2.*h)/sq);
			    t2 = ry * ti / sq;
			    if (aal) tx[it] = SF_MAX(t1,t2);
			    amp[it] = 
				sqrtf(time[it] * nt * dt) * ti * sqrtf(ti)/
				(sq * sqrtf(sq));
			} else {
			    str[it] = t0 - 2.*dt;
			    tx[it] = 0.;
			    amp[it] = 0.;
			}
		    } else { /* migration */
			t = time[it];
			if (t > h) {
			    sq = t - x + hx/t;
			    if (sq > 0.) {
				sq = sqrtf(sq);
				str[it] = sq;
				t1 = rx*(1.-h/t);
				t2 = ry;
				if (aal) tx[it] = SF_MAX(t1,t2)/sq;
				amp[it] = 1.;
			    } else {
				str[it] = t0 - 2.*dt;
				tx[it] = 0.;
				amp[it] = 0.;
			    }
			} else {
			    str[it] = t0 - 2.*dt;
			    tx[it] = 0.;
			    amp[it] = 0.;
			}
		    } /* modeling - migration */ 
		}  /* it */

		sf_aastretch_define (str, tx, amp);
		
		for (iyin=iy1; iyin < iy2; iyin++) {
		    for (ixin=ix1; ixin < ix2; ixin++) {
			sf_aastretch_lop (false,false,nt,nt,
				       cinp[iyin][ixin], add);
			for (it=0; it < nt; it++) {
			    cout[iyin+iy][ixin+ix][it] += add[it];
			}
		    }
		}

	    } /* x */
	} /* y */

	if (!inv && zero) {
	    for (it=0; it < n123; it++) {
		stack[0][0][it] += cout[0][0][it];
	    }
	} else {
	    sf_floatwrite (cout[0][0],n123,out);
	}
    } /* h */
    sf_warning(".");

    if (!inv && zero)
	sf_floatwrite (stack[0][0],n123,out);

    exit(0);
}
