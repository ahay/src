/* Shot-recorder Oriented prestack migration data is (shots x recs x time). */
/*
  Copyright (C) 2009 University of Texas at Austin
  
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
#include <assert.h>

int main(int argc, char* argv[])
{
    int iz, nz, ix, nx, ia, ja, na, it, nt, ih, nh, is, ns;
    float t, h, s, dt, dh, ds, t0, s0, h0, dz, z0, tolz;
    float zi, zj, xi, ti, tj;
    float *tim, *dis, *dep, ***dat, **img;
    sf_file time, dist, dept, data, imag;
    bool is_offset;

    sf_init(argc,argv);

    dist = sf_input("in");
    time = sf_input("time");
    dept = sf_input("dept");
    imag = sf_output("out");

    /* get image dimensions */
    if (!sf_histint(time,"n1",&na)) sf_error("No n1= in input");
    if (!sf_histint(time,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(time,"n3",&nz)) sf_error("No n3= in input");

    if (!sf_histfloat(time,"d3",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(time,"o3",&z0)) sf_error("No d1= in input");

    /* tolerance depth */
    if (!sf_getfloat("tolz",&tolz)) tolz=2.0;
    tolz *= dz;

    /* tolerance depth */
    if (!sf_getbool("is_offset",&is_offset)) is_offset=false;

    tim = sf_floatalloc(na);
    dis = sf_floatalloc(na);
    dep = sf_floatalloc(na);

    data = sf_input("data");

    /* get data dimensions */
    if (!sf_histint(data,"n1",&nt)) sf_error("No n1= in data");
    if (!sf_histint(data,"n2",&nh)) sf_error("No n2= in data");
    if (!sf_histint(data,"n3",&ns)) sf_error("No n3= in data");

    if (!sf_histfloat(data,"d1",&dt)) sf_error("No d1= in data");
    if (!sf_histfloat(data,"d2",&dh)) sf_error("No d2= in data");
    if (!sf_histfloat(data,"d3",&ds)) sf_error("No d3= in data");

    if (!sf_histfloat(data,"o1",&t0)) sf_error("No o1= in data");
    if (!sf_histfloat(data,"o2",&h0)) sf_error("No o2= in data");
    if (!sf_histfloat(data,"o3",&s0)) sf_error("No o3= in data");

    /* read the whole dataset in memory */
    dat = sf_floatalloc3(nt,nh,ns);
    sf_floatread(dat[0][0],nt*nh*ns,data);
    sf_fileclose(data);

    /* The output is (a, a, x, z) */
    sf_shiftdim(dist, imag, 1);
    sf_putint(imag,"n1",na);

    img = sf_floatalloc2(na,na);

    for (iz=0; iz < nz; iz++) {

	sf_warning("depth %d of %d",iz+1,nz);

	for (ja=0; ja < na*na; ja++) 
	    img[0][ja] = 0.0;    

	for (ix=0; ix < nx; ix++) {

	    sf_floatread(tim,na,time);
	    sf_floatread(dis,na,dist);
	    sf_floatread(dep,na,dept);

	    for (ia=0; ia < na; ia++) {

		/* check that depth is on the surface */

		zi = dep[ia];
		
		if (zi > z0 + tolz) {
		    for (ja=0/*ia*/; ja < na; ja++) {
			img[ia][ja] = 0.0;
		    }
		    continue;
		}

		/* get source location, interpolate */

		xi = dis[ia];

		s = (xi-s0)/ds;
		is = floorf(s);
		if (is < 0 || is > ns-1) {
		    for (ja=0/*ia*/; ja < na; ja++) {
			img[ia][ja] = 0.0;
		    }
		    continue;
		}
		s -= is;

		ti = tim[ia];		

		for (ja=0; ja < na; ja++) {
		    zj = dep[ja];
		    
		    if (zj > z0 + tolz) {
			img[ia][ja] = 0.0;
			continue;
		    }
		    
		    /* get receiver location, interpolate */

		    if (is_offset)
			h = (dis[ja]-xi-h0)/dh; 
		    else
			h = (dis[ja] - h0)/dh;

		    ih = floorf(h);
		    if (ih < 0 || ih > nh-1) {
			img[ia][ja] = 0.0;
			continue;
		    }
		    h -= ih;

		    /* get time, interpolate */

		    tj = tim[ja];

		    t = (ti+tj-t0)/dt; 
		    it = floorf(t);
		    if (it < 0 || it >= nt-1) {
			img[ia][ja] = 0.0;
			continue;
		    }
		    t -= it;

		    if (is+1 > ns-1 && ih+1 > nh-1) {

			img[ia][ja] = 
			    dat[is][ih][it  ] * (1.-t) +
			    dat[is][ih][it+1] *     t;
		    }
		    else {
			if (is+1 > ns-1 && ih+1 <= nh-1) {

			    img[ia][ja] = 
				
				dat[is][ih][it]*(1.-h)*(1.-t) +
				dat[is][ih][it+1]*(1.-h)*t +
				dat[is][ih+1][it]*h*(1.-t) +
				dat[is][ih+1][it+1]*h*t;
			}
			else {
			    if (is + 1 <= ns -1) {
				/* trilinear interpolation from the data */
				
				img[ia][ja] = 
				
				    dat[is][ih][it]*(1.-s)*(1.-h)*(1.-t) +
				    dat[is][ih][it+1]*(1.-s)*(1.-h)*t +
				    dat[is][ih+1][it]*(1.-s)*h*(1.-t) +
				    dat[is+1][ih][it]*s*(1.-h)*(1.-t) +
				    dat[is+1][ih][it+1]*s*(1.-h)*t +
				    dat[is][ih+1][it+1]*(1.-s)*h*t +
				    dat[is+1][ih+1][it]*s*h*(1.-t) +
				    dat[is+1][ih+1][it+1]*s*h*t;
			    }
			}
		    }
		}
	    }

	    sf_floatwrite(img[0],na*na,imag);
	}
    }

    exit(0);
}
