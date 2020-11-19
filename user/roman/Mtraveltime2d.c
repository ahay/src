/* Oriented zero-offset migration. */
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

void select_t_l_x(float l, float t, float src_dist, float dx, float *p_src, float *p_len_min, float *p_tr_time)
{
    /*
	if ( (*p_src_neg > src_dist && src_dist > dx) ||
	     (*p_src_neg > dx && src_dist < dx) ||
	     (*p_src_neg < dx && src_dist < dx && *p_len_min_neg > l)) 
    */
    if ( (*p_src > src_dist + dx) || (fabs(*p_src-src_dist)<dx && *p_len_min > l))
        {
                *p_src = src_dist;
		*p_len_min = l;
		*p_tr_time = t;
	}
}
int iz, nz, ix, nx, ia, na, interpolate;
float dz, dx, z0, x0;

void set_image(sf_file *minpath)
{
    sf_putint(*minpath,"n1",nx);
    sf_putint(*minpath,"n2",nz);
    sf_putint(*minpath,"n3",1);

    sf_putfloat(*minpath,"d1",dx);
    sf_putfloat(*minpath,"d2",dz);

    sf_putfloat(*minpath,"o1",x0);
    sf_putfloat(*minpath,"o2",z0);

    sf_putstring(*minpath,"label1","x dist");
    sf_putstring(*minpath,"label2","z depth");
    sf_putstring(*minpath,"unit1","km");	
    sf_putstring(*minpath,"unit2","km");	
}
int main(int argc, char* argv[])
{
    float t, z, x, xsrc, l, tolz_dz,tolz, xprev, lprev, tprev, zprev, eps;
    float *tim, *dis, *dep, *len;
/*
    float **len_min_neg, **tr_time_neg, **src_neg;
    float **len_min_pos, **tr_time_pos, **src_pos;

    float *p_len_min_neg, *p_tr_time_neg, *p_src_neg;
    float *p_len_min_pos, *p_tr_time_pos, *p_src_pos;
*/
    int is_prev_sign, is_sign=0;
    float **len_minpath, **tr_time_minpath, **tr_time_pos, ***tr_time_z0, **t_eik, **front, 
	front_t0, front_dt, front_eps;
    int k, front_nt;

    sf_file time, dist, dept, lent, imagt,timez0, minpath, eiksf, frontsf;
    char *timez0_name, *minpath_name, *eik_name, *front_name;
    

    sf_init(argc,argv);

    dist = sf_input("in");
    time = sf_input("time");
    dept = sf_input("dept");
    lent = sf_input("len");

    if (!sf_histint(time,"n3",&nz)) sf_error("No n1= in input");
    if (!sf_histint(time,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(time,"n1",&na)) sf_error("No n3= in input");
    if (!sf_histfloat(time,"d3",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(time,"o3",&z0)) sf_error("No d1= in input");
    if (!sf_histfloat(time,"d2",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(time,"o2",&x0)) sf_error("No d1= in input");

    if (!sf_getfloat("xsrc",&xsrc)) sf_error("No xsrc= ");
    if (!sf_getfloat("tolz",&tolz)) sf_error("No tolz= (float)");
    //if (!sf_getfloat("tolx",&tolx)) sf_error("No tolx= (float)");
    if (!sf_getfloat("frontt0",&front_t0)) sf_error("No frontt0= (float)");
    if (!sf_getfloat("frontdt",&front_dt)) sf_error("No frontdt= (float)");
    if (!sf_getint("frontnt",&front_nt)) sf_error("No frontnt= (int)");
    if (!sf_getfloat("fronteps",&front_eps)) sf_error("No fronteps= (float)");

    if (!sf_getint("interpolate",&interpolate)) 
	sf_error("<'2' works BEST!> No interpolate= (int) 0-direct rays 1-escape dist to source 2-escape ray-length");

    /* first arrivals */
    imagt = sf_output("out");

    tim = sf_floatalloc(na);
    dis = sf_floatalloc(na);
    dep = sf_floatalloc(na);
    len = sf_floatalloc(na);

    /* Output: */
    timez0_name = sf_getstring ("timez0");    
    assert (timez0_name);
    timez0 = sf_output(timez0_name);   

    minpath_name = sf_getstring ("minpath");
    assert (minpath_name);
    minpath = sf_output(minpath_name);

    front_name = sf_getstring ("front");
    assert (front_name);
    frontsf = sf_output(front_name);

    /* For drawing I tested this stuff  */
    t_eik = (float**)0;
    eik_name = sf_getstring ("eik");    
    if (eik_name) {
	eiksf = sf_input(eik_name);
	t_eik = sf_floatalloc2(nx,nz);
	sf_floatread(t_eik[0],nx*nz,eiksf);
    }

    /* min-path rays arrivals */
    set_image(&minpath);
    /* fronts of all arrivals */
    set_image(&frontsf);
    /* first arrivals */
    set_image(&imagt);

    /* allocate */
    tr_time_z0 = sf_floatalloc3(na, nx, nz); 
    tr_time_pos = sf_floatalloc2(nx, nz); 
    tr_time_minpath = sf_floatalloc2(nx, nz); 
    len_minpath = sf_floatalloc2(nx, nz); 
    front = sf_floatalloc2(nx, nz); 

    tolz_dz = tolz * dz;
    //tolx_dx = tolx * dx;

    assert(tr_time_z0);

    for (ix=0; ix < nx*nz; ix++) 
	tr_time_minpath[0][ix] = tr_time_pos[0][ix] = front[0][ix]= -1e-2f;

    for (ix=0; ix < nx*nz*na; ix++) 
	tr_time_z0[0][0][ix]=-0.01;

    for (iz=0; iz < nz; iz++) {
	
	for (ix=0; ix < nx; ix++) {					

	    sf_floatread(tim,na,time);
	    sf_floatread(dis,na,dist);
	    sf_floatread(dep,na,dept);
	    sf_floatread(len,na,lent);


	    /* take care of angle ia==0: */

	    is_prev_sign = 0; /* 0-uninit,'-1'-less, '1'-greater then x-src */

	    z = dep[na-1];
	    xprev = dis[na-1];
	    lprev = len[na-1];
	    tprev = tim[na-1];

	    if (z > (z0+ tolz_dz)) {
		is_prev_sign = 0; 
		zprev=1e10;
	    }
	    else {
		if (xprev < xsrc - 1e-5f) {
		    is_prev_sign = -1;
		    zprev=z;
		}
		else {
		    if (xprev > xsrc + 1e-5f) {
			is_prev_sign = 1;
			zprev = z;
		    }
		    else {
			/* rare case- hit src exactly */
			is_prev_sign = 2;
			tr_time_z0[iz][ix][na-1]=tprev;
			zprev=1e10;
		    }
		}
	    }

	    for (ia=0; ia < na; ia++) {

		z = dep[ia];
		/* dep[iz] = 0.0; */
		x = dis[ia];

		l = len[ia];

		t = tim[ia];

		if (z > (z0+ 10*tolz_dz)) { /* || fabs(xsrc-x)>2*dx)  */
		    
		    if (t_eik)
			tr_time_z0[iz][ix][ia]=t_eik[iz][ix];//-0.01f;
		    else
			tr_time_z0[iz][ix][ia]=-0.01f;

		    is_sign = 0; 
		}
		else {
		    if (x < xsrc - 1e-5f) {
			is_sign = -1;
		    }
		    else {
			if (x > xsrc + 1e-5f) {
			    is_sign = 1;
			}
			else {
			    /* rare case- hit src exactly */
			    if (z<z0+tolz_dz) {
				is_sign = 2;
				tr_time_z0[iz][ix][ia]=t;
				
				if (tr_time_pos[iz][ix] < 0 || tr_time_pos[iz][ix]  > t)
				    tr_time_pos[iz][ix]  = t;
				
				if (tr_time_minpath[iz][ix] < 0 || len_minpath[iz][ix]  > l) {
				    tr_time_minpath[iz][ix]  = t;
				    len_minpath[iz][ix]=l;
				}
			    }
			}
		    }

		    if (((is_prev_sign == -1 && is_sign == 1) || 
			 (is_prev_sign == 1 && is_sign == -1)) &&
			(z < z0+tolz_dz || zprev < z0+tolz_dz))
		    {
			/* xprev = dis[ia-1]; */
			
			eps = fabs(x - xsrc)/fabs(x - xprev);
			//d1=(x-xsrc)*(x-xsrc)+(z-zsrc)*(z-zsrc);
			//d2=(xprev-xsrc)*(xprev-xsrc)+(zprev-zsrc)*(zprev-zsrc);

			//assert (eps + 1e-6f > 0 && eps < 1.f + 1e-6);
			
			tr_time_z0[iz][ix][ia] = tprev * eps + t*(1.f-eps);

			/* lprev = len[ia-1]; */
			l = lprev*eps + l*(1.f - eps);
			
			if (tr_time_pos[iz][ix] < 0 || tr_time_pos[iz][ix]  > tr_time_z0[iz][ix][ia])
			    tr_time_pos[iz][ix]  = tr_time_z0[iz][ix][ia];

			if (tr_time_minpath[iz][ix] < 0 || len_minpath[iz][ix]  > l) {
			    tr_time_minpath[iz][ix]  = tr_time_z0[iz][ix][ia];
			    len_minpath[iz][ix]=l;
			}
		    }
		} /* hit surface */
 
		is_prev_sign = is_sign;
		tprev = t;
		xprev = x;
		zprev = z;
		lprev = len[ia];
		    
	    } /* ia */
	} /* iz */
    } /* ix */
	
    /* DID not WORK! */
    for (k = 0; k < front_nt; k++) {

	l = front_t0 + k* front_dt;

	for (iz=0; iz < nz; iz++) 
	    for (ix=0; ix < nx; ix++) 
	    for (ia=0; ia < na; ia++) {
		/* if (tr_time_z0[iz][ix][ia] < 0.f)
		   tr_time_z0[iz][ix][ia]=tr_time_minpath[iz][ix];*/
		if (fabs(l-tr_time_z0[iz][ix][ia]) < front_eps)
		    front[iz][ix]=k;//1.0; front_val;
                //else
		//  front[iz][ix]=0.0;
	    }
    }

    sf_floatwrite(tr_time_z0[0][0],nz*nx*na,timez0);
    sf_floatwrite(tr_time_pos[0],nz*nx,imagt);
    sf_floatwrite(tr_time_minpath[0],nz*nx,minpath);
    sf_floatwrite(front[0],nz*nx,frontsf);

    sf_close();
    exit(0);
}
