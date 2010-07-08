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
int main(int argc, char* argv[])
{
    int iz, nz, ix, nx, ia, na, interpolate;
    float t, dz, z, x, dx, z0, x0, xsrc, l, src_dist, l1, l2,tolz_dz,tolz, tolx, tolx_dx;
    float *tim, *dis, *dep, *len;
    float **len_min_neg, **tr_time_neg, **src_neg;
    float **len_min_pos, **tr_time_pos, **src_pos;

    float *p_len_min_neg, *p_tr_time_neg, *p_src_neg;
    float *p_len_min_pos, *p_tr_time_pos, *p_src_pos;
    float **tr_time_z0;

    sf_file time, dist, dept, lent, imagt,timez0;
    char *timez0_name;
    

    sf_init(argc,argv);

    dist = sf_input("in");
    time = sf_input("time");
    dept = sf_input("dept");
    lent = sf_input("len");

    if (!sf_histint(time,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histint(time,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(time,"n3",&na)) sf_error("No n3= in input");
    if (!sf_histfloat(time,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(time,"o1",&z0)) sf_error("No d1= in input");
    if (!sf_histfloat(time,"d2",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(time,"o2",&x0)) sf_error("No d1= in input");

    if (!sf_getfloat("xsrc",&xsrc)) sf_error("No xsrc= ");
    if (!sf_getfloat("tolz",&tolz)) sf_error("No tolz= (float)");
    if (!sf_getfloat("tolx",&tolx)) sf_error("No tolx= (float)");

    if (!sf_getint("interpolate",&interpolate)) 
	sf_error("<'2' works BEST!> No interpolate= (int) 0-direct rays 1-escape dist to source 2-escape ray-length");

    imagt = sf_output("out");
    sf_putint(imagt,"n1",nz);
    sf_putint(imagt,"n2",nx);
    sf_putint(imagt,"n3",1);
    sf_putfloat(imagt,"d1",dz);
    sf_putfloat(imagt,"d2",dx);
    sf_putfloat(imagt,"o1",z0);
    sf_putfloat(imagt,"o2",x0);
    sf_putstring(imagt,"label1","z depth");
    sf_putstring(imagt,"label2","x dist");
    sf_putstring(imagt,"unit1","km");	
    sf_putstring(imagt,"unit2","km");	

    tim = sf_floatalloc(nz);
    dis = sf_floatalloc(nz);
    dep = sf_floatalloc(nz);
    len = sf_floatalloc(nz);

    timez0_name = sf_getstring ("timez0");    
    if (timez0_name) {
	tr_time_z0 = sf_floatalloc2(nz, nx); 
	timez0 = sf_output(timez0_name);
    }

    len_min_neg = sf_floatalloc2(nz, nx);
    src_neg = sf_floatalloc2(nz, nx);
    tr_time_neg = sf_floatalloc2(nz, nx);

    len_min_pos = sf_floatalloc2(nz, nx);
    src_pos = sf_floatalloc2(nz, nx);
    tr_time_pos = sf_floatalloc2(nz, nx);

    tolz_dz = tolz * dz;
    tolx_dx = tolx * dx;

    for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		len_min_neg[ix][iz] = 1e10; 
		src_neg[ix][iz] = 1e10; 
		tr_time_neg[ix][iz] = -1.f;		

		len_min_pos[ix][iz] = 1e10; 
		src_pos[ix][iz] = 1e10; 
		tr_time_pos[ix][iz] = -1.f;		
	    }
    }

    for (ia=0; ia < na; ia++) {

	p_len_min_neg = &(len_min_neg[0][0]);
	p_src_neg = &(src_neg[0][0]);
	p_tr_time_neg = &(tr_time_neg[0][0]);

	p_len_min_pos = &(len_min_pos[0][0]);
	p_src_pos = &(src_pos[0][0]);
	p_tr_time_pos = &(tr_time_pos[0][0]);

	for (ix=0; ix < nx; ix++) {					

	    sf_floatread(tim,nz,time);
	    sf_floatread(dis,nz,dist);
	    sf_floatread(dep,nz,dept);
	    sf_floatread(len,nz,lent);

	/* shooting 360+- angles may not be enough:
	best way: if xsrc \in [a, a+-1] => add t=(t_a + t_a+-1)/2,
	here now: keep the closest to the src */
	    for (iz=0; iz < nz; iz++, p_len_min_neg++, p_len_min_pos++, p_tr_time_neg++, p_tr_time_pos++, p_src_neg++, p_src_pos++) {
		z = dep[iz];
		/* dep[iz] = 0.0; */
		x = dis[iz];

		if (z > (z0+ tolz_dz)) { /* || fabs(xsrc-x)>2*dx)  */
		    
		    if (timez0_name) 
			tr_time_z0[ix][iz]=-0.01f;

		    continue;
		}

		if (timez0_name) {
		    if (fabs(x - xsrc) < tolx_dx) {
			tr_time_z0[ix][iz]=tim[iz];
		    }
		    else {
			tr_time_z0[ix][iz]=-0.01f;
		    }
		}

		l = len[iz];

		t = tim[iz];

		src_dist = fabs(xsrc - x);

		if (x < xsrc) {
			select_t_l_x(l, t, src_dist, tolx_dx, p_src_neg, p_len_min_neg, p_tr_time_neg);
		} 
		else {
			select_t_l_x(l, t, src_dist, tolx_dx, p_src_pos, p_len_min_pos, p_tr_time_pos);
		}
	    } /* iz */
	} /* ix */
	
	if (timez0_name) 
	    sf_floatwrite(tr_time_z0[0],nz*nx,timez0);

    } /* ia */

	p_src_neg = src_neg[0];
	p_src_pos = src_pos[0];
	p_tr_time_neg = tr_time_neg[0];
	p_tr_time_pos = tr_time_pos[0];

	p_len_min_neg = len_min_neg[0];
	p_len_min_pos = len_min_pos[0];

	/* linear interpolate and return time_pos */
	for (ix=0; ix < nx; ix++) {

	    for (iz=0; iz < nz; iz++, 
		     p_src_neg++, p_src_pos++, 
		     p_tr_time_neg++, p_tr_time_pos++,
		     p_len_min_neg++, p_len_min_pos++) 
	    {
		if (*p_tr_time_neg < 0.f) {
			; /* *p_tr_time_pos = pos */
		}
		else { 
			if (*p_tr_time_pos < 0.f) {
				*p_tr_time_pos = *p_tr_time_neg;
			}
			else {
				/* if (*p_src_neg < 2*dx && *p_src_pos < 2*dx) { */
				/* linear interpolate by straight rays i.e. triangle */
				switch (interpolate) {

				    case 0: /* interpolate by straight lines from the source */
					x = x0 + ix*dx;
					z = z0 + iz*dz;

					l1 = (xsrc - *p_src_neg) - x;
					l1 = sqrt(l1*l1 + z*z);
					
					l2 = (xsrc + *p_src_pos) - x;
					l2 = sqrt(l2*l2 + z*z);
					
					l = xsrc - x;
					l = sqrt(l*l + z*z);				

					if (fabs(l2 - l1) > 1e-5) { // 1e-2f) { 
					    l = (l - l1) / (l2 - l1);
					    /* bad: l = (*p_src_neg) / (*p_src_neg + *p_src_pos); */
					}
					else 
					    l = 0;					
					break;
				    case 1: /* by travel ray-length from the source */
					l = (*p_len_min_neg) / (*p_len_min_neg + *p_len_min_neg); 
					break;
				    case 2: /* by escape distance to the  source */
					l = (*p_src_neg) / (*p_src_neg + *p_src_pos); 
					break;
				    case 3: /* by escape distance to the  source */
					l = 0.5; 
					break;
				    case 4: /* by escape distance to the  source */
					l = 0.0; 
					break;					
				    default: 

					assert(0);
				}

				*p_tr_time_pos = *p_tr_time_neg * (1.f - l) + *p_tr_time_pos * l;				
			}
		}
	    }
	}

    sf_floatwrite(tr_time_pos[0],nz*nx,imagt);

    sf_close();
    exit(0);
}
