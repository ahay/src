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

void select_t_l_x(float l, float t, float src_dist, float dx, float *p_src_neg, float *p_len_min_neg, float *p_tr_time_neg)
{
    if ( (*p_src_neg > src_dist && src_dist > dx) ||
	 (*p_src_neg > dx && src_dist < dx) ||
	 (*p_src_neg < dx && src_dist < dx && *p_len_min_neg > l)) 
    {
	*p_src_neg = src_dist;
	*p_len_min_neg = l;
	*p_tr_time_neg = t;
    }
}

int nz, nx, ny, na, nb;
float dz, z0, dx, x0, dy,  y00, xsrc, ysrc;

float ***len_min1, ***tr_time1, ***src1;
float ***len_min2, ***tr_time2, ***src2;
float ***len_min3, ***tr_time3, ***src3;
float ***len_min4, ***tr_time4, ***src4;

float *p_len_min1, *p_tr_time1, *p_src1, *p_srcY1;
float *p_len_min2, *p_tr_time2, *p_src2, *p_srcY2;
float *p_len_min3, *p_tr_time3, *p_src3,  *p_srcY3;
float *p_len_min4, *p_tr_time4, *p_src4, *p_srcY4;

void init_ptrs()
{
    p_len_min1 = len_min1[0][0];
    p_src1 = src1[0][0];
    p_tr_time1 = tr_time1[0][0];

    p_len_min2 = len_min2[0][0];
    p_src2 = src2[0][0];
    p_tr_time1 = tr_time1[0][0];

    p_len_min3 = len_min3[0][0];
    p_src3 = src3[0][0];
    p_tr_time3 = tr_time3[0][0];

    p_len_min4 = len_min4[0][0];
    p_src4 = src4[0][0];
    p_tr_time4 = tr_time4[0][0];
}

void init_vals()
{
    int ix, iy, iz;
    len_min1 = sf_floatalloc3(nz, nx, ny);
    src1 = sf_floatalloc3(nz, nx, ny);
    tr_time1 = sf_floatalloc3(nz, nx, ny);

    len_min2 = sf_floatalloc3(nz, nx, ny);
    src2 = sf_floatalloc3(nz, nx, ny);
    tr_time2 = sf_floatalloc3(nz, nx, ny);

    len_min3 = sf_floatalloc3(nz, nx, ny);
    src3 = sf_floatalloc3(nz, nx, ny);
    tr_time3 = sf_floatalloc3(nz, nx, ny);

    len_min4 = sf_floatalloc3(nz, nx, ny);
    src4 = sf_floatalloc3(nz, nx, ny);
    tr_time4 = sf_floatalloc3(nz, nx, ny);


    for (iy=0; iy < ny; iy++) {
	for (ix=0; ix < nx; ix++) {
	    for (iz=0; iz < nz; iz++) {
		len_min1[iy][ix][iz] = 1e10f; 
		tr_time1[iy][ix][iz] = -1.f;
		src1[iy][ix][iz] = 1e10f; 

		len_min2[iy][ix][iz] = 1e10f; 
		tr_time2[iy][ix][iz] = -1.f;
		src2[iy][ix][iz] = 1e10f; 

		len_min3[iy][ix][iz] = 1e10f; 
		tr_time3[iy][ix][iz] = -1.f;
		src3[iy][ix][iz] = 1e10f; 

		len_min4[iy][ix][iz] = 1e10f; 
		tr_time4[iy][ix][iz] = -1.f;
		src4[iy][ix][iz] = 1e10f; 
	    }
	}
    }
}
void move_ptrs()
{
/*		     p_len_min1++, p_tr_time1++, p_src1++,
		     p_len_min2++, p_tr_time2++, p_src2++,
		     p_len_min3++, p_tr_time3++, p_src3++,
		     p_len_min4++, p_tr_time4++, p_src4++) 
*/
    p_len_min1++;
    p_tr_time1++;
    p_src1++;
    p_len_min2++;
    p_tr_time2++;
    p_src2++;
    p_len_min3++;
    p_tr_time3++;
    p_src3++;
    p_len_min4++;
    p_tr_time4++;
    p_src4++;
}
float dist_to (float z, float *p_tr_time1, float *p_src1, float * p_srcY1)
{
    if (*p_tr_time1 < 0.f) 
    {
	*p_tr_time1 = 0.f;		    		    
	return 1e10f;		    
    }
    
    return sqrt(z*z + (xsrc - *p_src1)* (xsrc - *p_src1) +  (ysrc - *p_srcY1)* (ysrc - *p_srcY1));
}
int main(int argc, char* argv[])
{
    int ix, iy, iz, ia, ib;

    float t, z, x, y, l, src_dist2, dr2;
    float *tim, *dis, *disY, *dep, *len;

    sf_file time, dist, dept, lent, imagt, distY;
    
    sf_init(argc,argv);

    dist = sf_input("in");
    distY = sf_input("disty");
    time = sf_input("time");
    dept = sf_input("dept");
    lent = sf_input("len");

    if (!sf_histint(time,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histint(time,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(time,"n3",&ny)) sf_error("No n3= in input");
    if (!sf_histint(time,"n4",&na)) sf_error("No n4= in input");
    if (!sf_histint(time,"n5",&nb)) sf_error("No n5= in input");

    if (!sf_histfloat(time,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(time,"o1",&z0)) sf_error("No d1= in input");
    if (!sf_histfloat(time,"d2",&dx)) sf_error("No d2= in input");
    if (!sf_histfloat(time,"o2",&x0)) sf_error("No d2= in input");
    if (!sf_histfloat(time,"d3",&dy)) sf_error("No d3= in input");
    if (!sf_histfloat(time,"o3",&y00)) sf_error("No d3= in input");

    if (!sf_getfloat("xsrc",&xsrc)) sf_error("No xsrc= ");
    if (!sf_getfloat("ysrc",&ysrc)) sf_error("No ysrc= ");

    imagt = sf_output("out");

    sf_putint(imagt,"n1",nz);
    sf_putint(imagt,"n2",nx);
    sf_putint(imagt,"n3",ny);
    sf_putint(imagt,"n4",1);

    sf_putfloat(imagt,"d1",dz);
    sf_putfloat(imagt,"d2",dx);
    sf_putfloat(imagt,"d3",dy);

    sf_putfloat(imagt,"o1",z0);
    sf_putfloat(imagt,"o2",x0);
    sf_putfloat(imagt,"o3",y00);

    sf_putstring(imagt,"label1","z depth");
    sf_putstring(imagt,"label2","x dist");
    sf_putstring(imagt,"label3","y dist");
    sf_putstring(imagt,"unit1","km");	
    sf_putstring(imagt,"unit2","km");	
    sf_putstring(imagt,"unit3","km");	

    tim = sf_floatalloc(nz);
    dis = sf_floatalloc(nz);
    disY = sf_floatalloc(nz);
    dep = sf_floatalloc(nz);
    len = sf_floatalloc(nz);

    init_vals();

    dr2 =  dx*dx + dy*dy + dz*dz;

    for (ib=0; ib < nb; ib++) {
	
	for (ia=0; ia < na; ia++) {

	    init_ptrs();

	    for (iy=0; iy < ny; iy++) {					

		for (ix=0; ix < nx; ix++) {					

		    sf_floatread(tim,nz,time);
		    sf_floatread(dis,nz,dist);
		    sf_floatread(disY,nz,distY);
		    sf_floatread(dep,nz,dept);
		    sf_floatread(len,nz,lent);

		    /* shooting 360+- angles may not be dense enough:
		       best way: if xsrc \in [a, a+-1] => add t=(t_a + t_a+-1)/2,
		       here now: keep the closest to the src */
		    for (iz=0; iz < nz; iz++, move_ptrs())
		    {
			z = dep[iz];

			x = dis[iz];

			y = disY[iz];

			if (z > (z0+ dz)) /* || fabs(xsrc-x)>2*dx)  */
			    continue;

			l = len[iz];

			t = tim[iz];

			src_dist2 = (xsrc - x)*(xsrc - x) + (ysrc - y)*(ysrc - y);

			if (x < xsrc) {
			    if (y < ysrc) 
				select_t_l_x(l, t, src_dist2, dr2, p_src1, p_len_min1, p_tr_time1);
			    else
				select_t_l_x(l, t, src_dist2, dr2, p_src2, p_len_min2, p_tr_time2);
			}
			else {
			    if (y < ysrc) 
				select_t_l_x(l, t, src_dist2, dr2, p_src3, p_len_min3, p_tr_time3);
			    else
				select_t_l_x(l, t, src_dist2, dr2, p_src4, p_len_min4, p_tr_time4);
			}
		    } /* iz */
		} /* ix */
	    } /* iy */
	} /* ia */
    } /* ib */

    init_ptrs();

    /* linear interpolate and return time_pos */
    for (iy=0; iy < ny; iy++) 
    {
	y = y00 + iy * dy;
	
	for (ix=0; ix < nx; ix++) 
	{
	    x = x0 + ix*dx;
	    
	    for (iz=0; iz < nz; iz++, move_ptrs())
	    {
		z = z0 + iz*dz;

		*p_tr_time1 = fmaxf(fmaxf(*p_tr_time1,*p_tr_time2), fmaxf(*p_tr_time3, *p_tr_time4));
	    }
	}
    }

    sf_floatwrite(tr_time1[0][0],ny*nz*nx,imagt);

    sf_close();
    exit(0);
}
