/* 3-D Kirchhoff time migration with antialiasing. */
/*
  Copyright (C) 2007 University of Texas at Austin
  
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

int main(int argc, char* argv[])
{
    int nt, nx, ny, ix, iy, n1, it, iz, itp, itm, ntr, itr, i;
    char *antialias;
    float hdr[3], *trace, tt, ***out;
    float h_in, x_in, y_in, hx, ftm, ftp, imp;
    float x,y, z, t, dx, dy, ox, oy, rx, ry, ft;
    float h, sq, ti, tx, ot, dt, tm, tp, vel;
    sf_file in, mig, head;

    sf_init (argc,argv);
    in = sf_input("in");
    head = sf_input("hdr");
    mig = sf_output("out");

    if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
    ntr = sf_leftsize(in,1);

    if (!sf_histfloat(in,"o1",&ot)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
    
    if (!sf_getint("n2",&nx))   sf_error("Need n2="); sf_putint(mig,"n2",nx);
    if (!sf_getfloat("d2",&dx)) sf_error("Need d2="); sf_putfloat(mig,"d2",dx);
    if (!sf_getfloat("o2",&ox)) sf_error("Need o2="); sf_putfloat(mig,"o2",ox);

    if (!sf_getint("n3",&ny))   sf_error("Need n3="); sf_putint(mig,"n3",ny);
    if (!sf_getfloat("d3",&dy)) sf_error("Need d3="); sf_putfloat(mig,"d3",dy);
    if (!sf_getfloat("o3",&oy)) sf_error("Need o3="); sf_putfloat(mig,"o3",oy);

    if (!sf_getint("n1",&n1))   sf_error("Need n1="); sf_putint(mig,"n1",n1);
    
    if (NULL == (antialias = sf_getstring("antialias"))) antialias="triangle";
    /* antialiasing type [triangle,flat,steep,none] */

    if (!sf_getfloat("vel",&vel)) sf_error("Need vel=");
    /* migration velocity */
    vel = 2./vel; 
    vel *= vel;

    trace = sf_floatalloc(nt);
    out = sf_floatalloc3(n1,nx,ny);

    for (i=0; i < n1*nx*ny; i++) {
	out[0][0][i] = 0.;
    }

    for (itr=0; itr < ntr; itr++) {
	sf_floatread (trace,nt,in);
	sf_floatread (hdr,3,head);

	x_in = hdr[0];
	y_in = hdr[1];
	h_in = hdr[2];
	
	h = h_in * h_in * vel;

	if ('t'==antialias[0]) { 
            /* double integration */
	    tt = trace[0];
	    for (it=1; it < nt; it++) {
		tt += trace[it];
		trace[it] = tt;
	    }
	    tt = trace[nt-1];
	    for (it=nt-2; it >=0; it--) {
		tt += trace[it];
		trace[it] = tt;
	    }
	}

	for (iy=0; iy < ny; iy++) {
	    y = oy + iy*dy - y_in;
	    ry = fabsf(y*dy) * vel;
	    y *= y * vel;
        
	    for (ix=0; ix < nx; ix++) {
		x = ox + ix*dx - x_in;
		rx = fabsf(x*dx) * vel;
		x *= x * vel;
		hx = h*x;
		x += y + h;

		for (iz=0; iz < n1; iz++) {
		    z = ot + iz*dt;              

		    t = z*z + x; 
		    if (t <= 0.) continue;
              
		    sq = t*t - 4.*hx; 
		    if (sq <= 0.) continue;

		    sq = sqrtf(sq);
		    ti = sqrtf(0.5*(t + sq));

		    tx = SF_MAX(
			rx / ti * 0.5 * (1 + (t - 2.*h)/sq), 
			ry * ti / sq);

		    if ('f'==antialias[0] && tx > dt) continue;
		    
		    ft = (ti-ot)/dt; it = ft; ft -= it; 

		    if ( it < 0) continue;
		    if ( it >= nt-1) break; 
              
		    if ('t'==antialias[0]) {
			tm = ti-tx-dt;
			ftm = (tm-ot)/dt; itm = ftm; ftm = ftm - itm; 
			if (itm < 0) continue;
                 
			tp = ti+tx+dt;
			ftp = (tp-ot)/dt; itp = ftp; ftp = ftp - itp; 
			if (itp >= nt-1) continue;

			imp = dt/(dt+tp-tm);
			imp *= imp;
			
			out[iy][ix][iz] += imp*(
			    2.*(1.-ft)*trace[it] + 2.*ft*trace[it+1] -
			    (1.-ftm)*trace[itm] - ftm*trace[itm+1]    - 
			    (1.-ftp)*trace[itp] - ftp*trace[itp+1]);
		    } else {
			out[iy][ix][iz] +=
			    (1.-ft)*trace[it] + ft*trace[it+1];
		    }
		} /* iz */
	    } /* ix */
	} /* iy */
    } /* itr */

    sf_floatwrite(out[0][0],n1*nx*ny,mig);        

    exit(0);
}
