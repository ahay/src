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
/*------------------------------------------*/

#define IMIN(X,Z) ( (X < Z) ? X : Z)
#define IMAX(X,Z) ( (X > Z) ? X : Z)

/* Phase space grid */
static int nz,nx,na;
static float oz,ox,oa;
static float dz,dx,da;


typedef enum {OUTSIDE=0, IZ_IX_UNINIT=0, IZ_IX_Z_N1=1, IZ_IX_Z_0=2, IZ_IX_N1_X=3, IZ_IX_0_X=4} colors;
const float INIT_VALUE = 1e10;

////////////////////////////////////////////////////
float ISIGN(float X) 
{
    if (X + 1e-6f > 0.f)
		return 1.f;
	return -1.f;
}
int bc_padd(int ix, int nx)
{
    if (ix < 0)
	return 0;
    if (ix >= nx)
	return nx-1;
    return ix;
}
int bc_periodic(int ia1, int na)
{
	    while (ia1 < 0) 
		ia1 += na; 
	    while (ia1 >= na) 
		ia1 -= na; 
	    return ia1;
}
float bc_periodicf(float ray_a1)
{
    ray_a1 =  fmodf(ray_a1, 2.*SF_PI);

    assert (fabs(ray_a1) < 2.0*SF_PI);

    if (ray_a1 < 0.f) 
	ray_a1 += 2.*SF_PI; 
    
    return ray_a1;
}
// float fx(float x, float z, float a)  {}
////////////////////////////////////////////////////
void step_z(float ray_z0, float DDz, float oz, float dz,
	    float *p_ray_z1, int *p_iz1, int *p_iz2, float *p_dist_z1)
{
	// step in z
        *p_ray_z1 = ray_z0 - DDz; // fz / fx * DDx;

	
	*p_iz1 = floorf((*p_ray_z1 - oz)/dz);
	assert(*p_ray_z1 +1e-4f > oz + (*p_iz1)*dz);

	*p_iz2 = *p_iz1 + 1;
	assert(*p_ray_z1 < 1e-4f + oz + (*p_iz2)*dz);

	*p_dist_z1 = (*p_ray_z1 - oz)/dz - *p_iz1;
	assert(*p_dist_z1 + 1e-4f > 0.f && *p_dist_z1 <= 1.0f);
	
	// merge two indices if dist==0 or ==1
	if (*p_dist_z1 < 1e-4) {
	    *p_dist_z1 = 0.f;
	    *p_iz2 = *p_iz1;
	    *p_ray_z1 = oz + (*p_iz2)*dz;
	}
	else  if (*p_dist_z1 + 1e-4 > 1.f) {
	    *p_dist_z1 = 1.f;
	    *p_iz1 = *p_iz2;
	    *p_ray_z1 = oz + (*p_iz2)*dz;
	}
}

float	propagate_x(int dir,
		    float fx,  float dx, 
		    float fz,  float dz, 
		    float fa,  float da,
		    float ox,
		    float oz,
		    float oa,
		    int nx,
		    int nz,
		    int na,
		    float ray_x0,  
		    float ray_z0, 
		    float ray_a0,
		    int * p_ix1, 
		    int * p_iz1, int *p_iz2,
		    int * p_ia1, int *p_ia2,
		    float * p_dist_z1,
		    float * p_dist_a1,
		    float * p_ray_x1, 
		    float * p_ray_z1,
		    float * p_ray_a1)
{
    assert (fabs(fx) > 0.1); // propagate in x

    const int ixs = dir * ISIGN(fx);

    //assert((ixs > 0 && ray_x0  > ox) || (ixs < 0 && ray_x0 < ox+(nx-1)*dx));

    if (ixs > 0) {						
	*p_ix1 = floorf ( (ray_x0-ox) / dx );
	if ((ray_x0 - (ox + (*p_ix1*dx)) < 0.01*dx) && *p_ix1 > 0) {
	    (*p_ix1) = (*p_ix1) - 1;
	}
    }
    else {
	*p_ix1 = 1 + floorf ((ray_x0-ox) / dx);
	if (((ox + (*p_ix1*dx)) - ray_x0 < 0.01*dx) && *p_ix1 < nx - 1) {
	    (*p_ix1) = (*p_ix1) + 1;
	}
    }
    //assert(dir < 0 || (*p_ix1 >= 0 && *p_ix1 < nx));

	*p_ray_x1 = ox + (*p_ix1) * dx;
	//assert(dir < 0 || (fabs(*p_ray_x1 - ray_x0) > 1e-6f));

	const float DDx = ray_x0 - *p_ray_x1;

	step_z(ray_z0, fz / fx * DDx, oz, dz, p_ray_z1, p_iz1, p_iz2, p_dist_z1);

	step_z(ray_a0, fa / fx * DDx, oa, da, p_ray_a1, p_ia1, p_ia2, p_dist_a1);

	return fmax( fabs(DDx)/dx, fmax( fabs(fz/fx*DDx)/dz, fabs(fa / fx * DDx)/da));
}


//
// propagate in x direction and in z direction the solution fx T_x dx + fz T_z dz = f
//
// here fx,fz > 0 but ixs, izs >< 0
//
// fx!=0 => propagate in x, 
// fz!=0 => propagate in z
// choose smallest value and smallest color (is someone was un-init => remain un-init)
float propagate_x_z (int dir,
                     int iq,
		     float *** t, int *** t_color, float **s, float ** sx, float **sz,
		     float dx, float dz, float da, 
		     int iz0, int ix0, int ia0, 
		     int * color,
                     float **imagt_rays, int imagt_rays_col) 
//int * is_short_ray)
{
    float 
	ray_x0 = ox + ix0 * dx,
        ray_z0 = oz + iz0 * dz,
	ray_a0 = oa + ia0 * da,
	accum_val = 0.f,
	ff2 = s[ix0][iz0] * s[ix0][iz0],
	ff1 = s[ix0][iz0],
	sn = sinf(ray_a0),
	cs = cosf(ray_a0),
	fx = sn * s[ix0][iz0],
	fz = cs * s[ix0][iz0],
	fa = (cs*sx[ix0][iz0] - sn*sz[ix0][iz0]);
    int
	ix1 = ix0,
	iz1 = iz0,
	ia1 = ia0;
   
    int ik = 0;

    while (ik < nx + nz + nz)
    {
	float  ss1, ss2, sx1, sx2, sz1, sz2, dist_x1z1, ray_x1 = ray_x0, ray_z1 = ray_z0, ray_a1 = ray_a0, dist_a1 = 0.f;
	int ix2, iz2, ia2 = ia0;
	ik++;
	
	if (fabs(fx) > fabs(fz)) {

	    const float max_delta_x = 

		propagate_x( dir,
			     fx, dx, fz, dz, fa, da, 
			     ox, oz, oa, nx, nz, na,
			     ray_x0, ray_z0, ray_a0,
			     &ix1, 
			     &iz1, &iz2, &ia1, &ia2,
			     &dist_x1z1, &dist_a1,
			     &ray_x1, &ray_z1, &ray_a1);
	    

	    assert(dir < 0 || max_delta_x > 1e-4f);
	    assert(dir < 0 || ISIGN(ray_x0 - ray_x1) == dir * ISIGN(fx));
	    assert( 1e-4 > fabs(ox + ix1*dx - ray_x1) );
	    assert(ray_z1 + 1e-4f > oz + iz1*dz && ray_z1 < 1e-4f + oz + iz2*dz);
	    //assert(ix1 >= 0 && ix1 < nx && iz1 >= 0 && iz1 < nz && iz2 >= 0 && iz2 < nz);
 
	     if (iq == 2) 
		 accum_val += ff2 / fx * (ray_x0 - ray_x1);
	     if (iq == 4) 
		 accum_val += ff1 / fx * (ray_x0 - ray_x1);
	     
	     ray_a1 = bc_periodicf(ray_a1);

		 ia1 = bc_periodic(ia1, na);
		 ia2 = bc_periodic(ia2, na);

		 ix1 = bc_padd(ix1, nx);

		 iz1 = bc_padd(iz1, nz);
		 iz2 = bc_padd(iz2, nz);


	     if (ix1<=0 || ix1>=nx-1 || ((iz1 <= 0 || iz1 >= nz-1) && (iz2<=0 || iz2 >= nz-1)))		 
	     {		
		 if (t_color)
		     *color = IMIN(IMIN(t_color[ia1][ix1][iz1], t_color[ia1][ix1][iz2]), IMIN(t_color[ia2][ix1][iz1], t_color[ia2][ix1][iz2]));
		 else
		     *color = 0;

		 if (t)
		     return  accum_val + 
			 (1.0f - dist_x1z1) * (t[ia1][ix1][iz1] * (1.0f - dist_a1) +  t[ia2][ix1][iz1] * dist_a1) +
			 dist_x1z1 * (t[ia1][ix1][iz2] * (1.0f - dist_a1) +  t[ia2][ix1][iz2] * dist_a1);
		 return -1;
	     }
	     else {
		
		 ss1 = s[ix1][iz1];ss2 = s[ix1][iz2];
		 
		 sx1 = sx[ix1][iz1];sx2 = sx[ix1][iz2];

		 sz1 = sz[ix1][iz1];sz2 = sz[ix1][iz2];
	     }
        }
	else { // propagate Z

	    const float max_delta_z = 

		propagate_x( dir,
			     fz, dz, fx, dx, fa, da, 
			     oz, ox, oa, nz, nx, na,
			     ray_z0, ray_x0, ray_a0,
			     &iz1, 
			     &ix1, &ix2, &ia1, &ia2,
			     &dist_x1z1, &dist_a1,
			     &ray_z1, &ray_x1, &ray_a1);
	    

	    assert(dir < 0 || max_delta_z > 1e-4f);
	    assert(dir < 0 || ISIGN(ray_z0 - ray_z1) == dir * ISIGN(fz));
	    assert(1e-4 > fabs(oz + iz1*dz - ray_z1));
	    assert(ray_x1  + 1e-4f > ox + ix1*dx && ray_x1 < 1e-4f +  ox + ix2*dx);
 
	     if (iq == 2) 
		 accum_val += ff2 / fz * (ray_z0 - ray_z1);
	     if (iq == 4) 
		 accum_val += ff1 / fz * (ray_z0 - ray_z1);
	     
	     ray_a1 = bc_periodicf(ray_a1);

		 ia1 = bc_periodic(ia1, na);
		 ia2 = bc_periodic(ia2, na);

		 iz1 = bc_padd(iz1, nz);

		 ix1 = bc_padd(ix1, nx);
		 ix2 = bc_padd(ix2, nx);

	     if (iz1 <= 0 || iz1 >= nz-1  || ((ix1<=0 || ix1>=nx-1) && (ix2<=0 || ix2 >= nx-1)))		 
	     {
		 if (t_color)
		     *color = IMIN(IMIN(t_color[ia1][ix1][iz1], t_color[ia1][ix2][iz1]), IMIN(t_color[ia2][ix1][iz1], t_color[ia2][ix2][iz1]));
		 else
		     *color = 0;
		 
		 if (t)
		     return  accum_val + 
			 (1.0f - dist_x1z1) * (t[ia1][ix1][iz1] * (1.0f - dist_a1) +  t[ia2][ix1][iz1] * dist_a1) +
			 dist_x1z1 * (t[ia1][ix2][iz1] * (1.0f - dist_a1) +  t[ia2][ix2][iz1] * dist_a1);
		 return -1;
	    }	    
	    else {

		ss1 = s[ix1][iz1]; ss2 = s[ix2][iz1];

		sx1 = sx[ix1][iz1];sx2 = sx[ix2][iz1];

		sz1 = sz[ix1][iz1];sz2 = sz[ix2][iz1];
		
	    }
	}
	assert(ray_x1 >= ox && ray_x1 <= ox + (nx-1)*dx);
	assert(ray_z1 >= oz && ray_z1 <= oz + (nz-1)*dz);

	ray_x0 = ray_x1;
	ray_z0 = ray_z1;
	ray_a0 = ray_a1;

	if (imagt_rays) {
	    imagt_rays[iz1][ix1] = imagt_rays_col;
	}

	const float 
	    ss = ss1 * (1.0f - dist_x1z1) + ss2 * dist_x1z1,
	    ssx = sx1 * (1.0f - dist_x1z1) + sx2 * dist_x1z1,
	    ssz = sz1 * (1.0f - dist_x1z1) + sz2 * dist_x1z1;

	sn = sinf(ray_a0);
	cs = cosf(ray_a0);

	fx = sn * ss;
	fz = cs * ss;

	ff2 = ss*ss;
	ff1 = ss;

	fa = (cs*ssx - sn*ssz);
    } // end while

    assert(0);
    //sf_warning(" ray tracing failed ");
    *color = IZ_IX_UNINIT;
    return 1e10f;
}

/*----------------------------------------*/
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
/*************************88

input: 
1. escape x, a, t, z 
2. xref, zref, aref
outpur:
image[ix x iz] == 1 at rays' points
=0 otherwise
****************************/
int main(int argc, char* argv[])
{
    int stride=1;
    int num_refs, i,ix, k;//, interpolate; iz, ia, 
    float **imagt_rays;//**len_minpath, **tr_time_minpath, **tr_time_pos, ***tr_time_z0;

    /*float escx1 = 0.25, escz1=0.0, esca1=0.5;
    float escx2 = 0.75, escz2=0.0, esca2=-0.5;
    float refx1=0.5, refz1=0.75, refa1=0.5;*/
    int ix0, iz0, ia0;

    float **s, **sx, **sz;
    float **m=NULL;//m[3][400];
    sf_file time;
    /* sf_file angl, time, dist, dept, timez0, minpath,lent, */
    sf_file imagt, slow, slowz, slowx, fref1, fref2, fref3;
    char *ref1_name, *ref2_name, *ref3_name;
    

    sf_init(argc,argv);

    time = sf_input("in"); // just to no the sizes
    /* time = sf_input("time");
    dept = sf_input("dept");
    angl = sf_input("ang");*/
    if (!sf_getint("stride",&stride)) stride=1;
    assert(stride>=1);

    if (!sf_histint(time,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histint(time,"n2",&nx)) sf_error("No n2= in input");
    if (!sf_histint(time,"n3",&na)) sf_error("No n3= in input");

    if (!sf_histfloat(time,"d1",&dz)) sf_error("No d1= in input");
    if (!sf_histfloat(time,"o1",&oz)) sf_error("No d1= in input");

    if (!sf_histfloat(time,"d2",&dx)) sf_error("No d1= in input");
    if (!sf_histfloat(time,"o2",&ox)) sf_error("No d1= in input");

    if (!sf_histfloat(time,"d3",&da)) sf_error("No d1= in input");
    if (!sf_histfloat(time,"o3",&oa)) sf_error("No d1= in input");
    /*
    if (!sf_getfloat("tolz",&tolz)) sf_error("No tolz= (float)");
    if (!sf_getfloat("tolx",&tolx)) sf_error("No tolx= (float)");
    if (!sf_getfloat("tola",&tola)) sf_error("No tola= (float)");
    tolz_dz = tolz * dz;
    tolx_dx = tolx * dx;
    tola_da = tola * da;
    */

    /* read auxiliary slowness file */
    slow = sf_input("slow");

    s = sf_floatalloc2(nz,nx);
    sf_floatread(s[0],nz*nx,slow);

    /* read auxiliary slowness z-gradient file */
    slowz = sf_input("slowz");
    sz = sf_floatalloc2(nz,nx);
    sf_floatread(sz[0],nz*nx,slowz);

    /* read auxiliary slowness x-gradient file */
    slowx = sf_input("slowx");
    sx = sf_floatalloc2(nz,nx);
    sf_floatread(sx[0],nz*nx,slowx);

    /* read ref1 = {x, z, a1, a2} */
    ref1_name = sf_getstring ("ref1");
    if (ref1_name) {
	fref1 = sf_input(ref1_name);

	if (!sf_histint(fref1,"n1",&num_refs)) sf_error("No num_refs= in ref1");    
	if (!m)
	    m=sf_floatalloc2(num_refs,3);

	sf_floatread(m[0], num_refs,fref1);
    }

    ref2_name = sf_getstring ("ref2");
    if (ref2_name) {
	fref2 = sf_input(ref2_name);

	if (!sf_histint(fref2,"n1",&num_refs)) sf_error("No num_refs= in ref2");    
	if (!m)
	    m=sf_floatalloc2(num_refs,3);

	sf_floatread(m[1], num_refs,fref2);
    }

    ref3_name = sf_getstring ("ref3");
    if (ref3_name) {
	fref3 = sf_input(ref3_name);

	if (!sf_histint(fref3,"n1",&num_refs)) sf_error("No num_refs= in ref3");
	if (!m)
	    m=sf_floatalloc2(num_refs,3);

	sf_floatread(m[2], num_refs,fref3);    
    }

    //fref3 = sf_input("ref3");
    //sf_floatread(m[2], num_refs,fref3);

    /* convert to radians */
    oa = oa*SF_PI/180.;
    da = da*SF_PI/180.;


    /* rays  */
    imagt = sf_output("out");
    imagt_rays = sf_floatalloc2(nx, nz); 

    sf_putint(imagt,"n2",nz);
    sf_putint(imagt,"n1",nx);
    sf_putint(imagt,"n3",1);

    sf_putfloat(imagt,"d2",dz);
    sf_putfloat(imagt,"d1",dx);
    sf_putfloat(imagt,"o2",oz);
    sf_putfloat(imagt,"o1",ox);

    sf_putstring(imagt,"label2","z depth");
    sf_putstring(imagt,"label1","x dist");
    sf_putstring(imagt,"unit2","km");	
    sf_putstring(imagt,"unit1","km");	


    for (ix=0; ix < nx*nz; ix++) 
	imagt_rays[0][ix] = 0.f;
    
	
    int ts_color;

    sf_warning("\n------------\n num-refs=%d\n",num_refs/4);

    for (i=0; i<3;i++) {

	if ((i== 0 && ref1_name)||
	    (i== 1 && ref2_name)||
	    (i== 2 && ref3_name))
	{

	    for (k=0; k<num_refs/4; k++) {

		ix0 = floorf(0.5+(m[i][4*k  ]-ox)/dx);
		iz0 = floorf(0.5+(m[i][4*k+1]-oz)/dz);

		if (stride > num_refs/4) {
		    imagt_rays[iz0][ix0]=i+3;
		    continue;
		}

		if (k%stride!=0) {
		    continue;
		}

	    ia0 = floorf(0.5+(m[i][4*k+2]-oa)/da);

	    sf_warning("ref=%d(%d): x=%g z=%g as=%g ar=%g\n",k,num_refs/4,m[i][4*k],m[i][4*k+1],m[i][4*k+2],m[i][4*k+3]);
	    if (ix0<=0 || iz0<=0 || ix0>=nx || iz0>= nz)
		continue;

	    if (ia0 < 0)
		ia0 += na;
	    if (ia0 > na-1)
		ia0 -= na;

	    imagt_rays[iz0][ix0]=i+3;

	    (void) propagate_x_z (1, 2/*iq*/, (float***)0/*t*/, (int***)0/*t_colors*/, s, sx, sz,  dx, dz, da, iz0, ix0, ia0, &ts_color, imagt_rays, i+3);

	    ix0 = floorf(0.5+(m[i][4*k  ]-ox)/dx);
	    iz0 = floorf(0.5+(m[i][4*k+1]-oz)/dz);
	    ia0 = floorf(0.5+(m[i][4*k+3]-oa)/da);

	    if (ia0 < 0)
		ia0 += na;
	    if (ia0 > na-1)
		ia0 -= na;

	    (void) propagate_x_z (1, 2/*iq*/, (float***)0/*t*/, (int***)0/*t_colors*/, s, sx, sz,  dx, dz, da, iz0, ix0, ia0, &ts_color, imagt_rays, i+3);
	    }
	}
    }
    sf_floatwrite(imagt_rays[0],nz*nx,imagt);

    sf_close();
    exit(0);
}
