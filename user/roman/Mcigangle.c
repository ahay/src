/* src-receiver to angle gathers */
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

int main(int argc, char* argv[])
{
    int iz, nz, ix, nx, ia, ja, na, it, nt, ih, nh, is, ns;
    float t, h, s, dt, dh, ds, t0, s0, h0, dz, z0;
    float zi, zj, xi, ti, tj;
    float *tim, *dis, *dep, ***dat, ***img;
    sf_file time, dist, dept, data, imag;

    int nalpha, ngamma, ialpha, igamma;
    float tolz, da, a0, alpha, gamma, alpha0, gamma0, dalpha, dgamma;

    sf_init(argc,argv);

    dist = sf_input("in");
    time = sf_input("time");
    dept = sf_input("dept");
    imag = sf_output("out");

    /* get image dimensions */
    if (!sf_histint(time,"n1",&na)) sf_error("No n1= in input");
    if (!sf_histint(time,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(time,"n3",&nz)) sf_error("No n3= in input");

    if (!sf_histfloat(time,"d3",&dz)) sf_error("No d3= in time");
    if (!sf_histfloat(time,"d1",&da)) sf_error("No d1= in time");

    if (!sf_histfloat(time,"o3",&z0)) sf_error("No d1= in input");
    if (!sf_histfloat(time,"o1",&a0)) sf_error("No o1= in time");

    if (!sf_getint("nalpha",&nalpha)) nalpha=90;

    if (!sf_getfloat("tolz",&tolz)) tolz=1.f;
    /* surface depth */
    tolz *= dz;

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

    /* The output is (alpha, gamma, t) */
    //nalpha = 360;
    ngamma = nalpha;
    alpha0 = 0.;
    gamma0 = 0.;
    dalpha = 360/nalpha;
    dgamma = 360/ngamma;

    //sf_shiftdim(dist, imag, 1);
    sf_putint(imag,"n3",nalpha);
    sf_putint(imag,"n2",ngamma);
    sf_putint(imag,"n1",nt);

    sf_putfloat(imag,"o3",alpha0);
    sf_putfloat(imag,"o2",gamma0);
    sf_putfloat(imag,"o1",t0);
    //sf_putint(imag,"n1",nt);

    sf_putfloat(imag,"d3",dalpha);
    sf_putfloat(imag,"d2",dgamma);
    sf_putfloat(imag,"d1",dt);

    img = sf_floatalloc3(nt, ngamma, nalpha);

    for (iz=0; iz < nalpha*ngamma*nt; iz++) 
	img[0][0][iz] = 0.f;

    for (iz=0; iz < nz; iz++) {

	sf_warning("depth %d of %d",iz+1,nz);

	for (ix=0; ix < nx; ix++) {
	    
	    sf_floatread(tim,na,time);
	    sf_floatread(dis,na,dist);
	    sf_floatread(dep,na,dept);


	    for (ia=0; ia < na; ia++) {

		/* check that depth is on the surface */

		zi = dep[ia];
		
		if (zi > z0 + tolz) {
		    continue;
		}

		/* get source location, interpolate */

		xi = dis[ia];

		s = (xi-s0)/ds;
		is = floorf(s);

		if (is < 0 || is >= ns-1) {
		    continue;
		}
		s -= is;

		ti = tim[ia];		

		for (ja=0; ja < na; ja++) {
		    zj = dep[ja];
		    
		    if (zj > z0 + tolz) {
			continue;
		    }
		    
		    /* get receiver location, interpolate */

		    //h = (dis[ja]-xi-h0)/dh;
		    h = ( dis[ja] - h0 )/dh;

		    ih = floorf(h);
		    if (ih < 0 || ih >= nh-1) {
			continue;
		    }
		    h -= ih;

		    /* get time, interpolate */
		    tj = tim[ja];

		    t = (ti+tj-t0)/dt; 
		    it = floorf(t);
		    if (it < 0 || it >= nt-1) {
			continue;
		    }
		    t -= it;

		    /* trilinear interpolation from the data */
		    alpha = a0 + 0.5f*da*(ia + ja);
		    gamma = 0.5f * da * (ia - ja);

		    ialpha = floor( (alpha - alpha0)/dalpha);
		    igamma = floor( (gamma - gamma0)/dgamma);

		    if (ialpha < 0) ialpha += nalpha;
		    if (igamma < 0) igamma += ngamma;

		    if (ialpha >= nalpha) ialpha -= nalpha;
		    if (igamma >= ngamma) igamma -= ngamma;
		    
		    if (ialpha >= 0 && ialpha < nalpha && igamma >=0 && igamma < ngamma)

			img[ialpha][igamma][it] = 

			dat[is][ih][it]*(1.-s)*(1.-h)*(1.-t) +
			dat[is][ih][it+1]*(1.-s)*(1.-h)*t +
			dat[is][ih+1][it]*(1.-s)*h*(1.-t) +
			dat[is+1][ih][it]*s*(1.-h)*(1.-t) +
			dat[is+1][ih][it+1]*s*(1.-h)*t +
			dat[is][ih+1][it+1]*(1.-s)*h*t +
			dat[is+1][ih+1][it]*s*h*(1.-t) +
			dat[is+1][ih+1][it+1]*s*h*t;
		}
	    } // ia
	} // ix
    } // iz

    sf_floatwrite(img[0][0],nalpha*ngamma*nt,imag);


    exit(0);
}
