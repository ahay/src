/* Ray tracing interface for general anisotropic media. */
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

#include <rsf.h>

#include "araytracegenani.h"
#include "eno3.h"
#include "grid3genani.h"
#include "bond.h"

#ifndef _araytracegenani_h

typedef struct aRayTrace* araytrace;
/* abstract data type */
/*^*/

#endif

struct aRayTrace {
    int nt, nz, ny, nx, order, *rungecount;
    float dt, z0, dz, y0, dy, x0, dx;
    float **oldg;
    bool wantsine, wantpolar;
    eno3 p11, p12, p13, p14, p15, p16, p22, p23, p24, p25, p26, p33, p34, p35, p36, p44, p45, p46, p55, p56, p66, pthe, pphi;
    grid3genani grd3;
};
/* concrete data type */

static void aniso_p_rhs(void* par, float* y, float* f)
/* right-hand side for qP anisotropic raytracing */
{    
    araytrace rt;

    rt = (araytrace) par;

	grid3genani_p_rhs(rt->grd3,y,f,rt->rungecount,rt->oldg);

}

static void aniso_s1_rhs(void* par, float* y, float* f)
/* right-hand side for qS1 anisotropic raytracing */
{    
    araytrace rt;

    rt = (araytrace) par;

	grid3genani_s1_rhs(rt->grd3,y,f,rt->rungecount,rt->oldg);

}

static void aniso_s2_rhs(void* par, float* y, float* f)
/* right-hand side for qS2 anisotropic raytracing */
{    
    araytrace rt;

    rt = (araytrace) par;

	grid3genani_s2_rhs(rt->grd3,y,f,rt->rungecount,rt->oldg);

}

static void aniso_s1c_rhs(void* par, float* y, float* f)
/* right-hand side for most coupled qS1 anisotropic raytracing */
{    
    araytrace rt;

    rt = (araytrace) par;

	grid3genani_s1c_rhs(rt->grd3,y,f,rt->rungecount,rt->oldg);

}

static void aniso_s2c_rhs(void* par, float* y, float* f)
/* right-hand side for most coupled qS2 anisotropic raytracing */
{    
    araytrace rt;

    rt = (araytrace) par;

	grid3genani_s2c_rhs(rt->grd3,y,f,rt->rungecount,rt->oldg);

}

static int term(void* par, float* y)
/* grid termination */
{
    araytrace rt;

    rt = (araytrace) par;

	return grid3genani_term(rt->grd3,y);

}

araytrace araytracegenani_init(int dim            /* dimensionality (2 or 3) */, 
			 int nt             /* number of ray tracing steps */, 
			 float dt           /* ray tracing step (in time) */,
			 int* n             /* velocity dimensions [dim] */, 
			 float* o, float* d /* velocity grid [dim] */,
			 float* c11         /* [n3*n2*n1] */, 
			 float* c12         /* [n3*n2*n1] */,
			 float* c13         /* [n3*n2*n1] */,
			 float* c14         /* [n3*n2*n1] */, 
			 float* c15         /* [n3*n2*n1] */,
			 float* c16         /* [n3*n2*n1] */,
			 float* c22         /* [n3*n2*n1] */, 
			 float* c23         /* [n3*n2*n1] */,
			 float* c24         /* [n3*n2*n1] */,
			 float* c25         /* [n3*n2*n1] */, 
			 float* c26         /* [n3*n2*n1] */,
			 float* c33         /* [n3*n2*n1] */,
			 float* c34         /* [n3*n2*n1] */, 
			 float* c35         /* [n3*n2*n1] */,
			 float* c36         /* [n3*n2*n1] */,
			 float* c44         /* [n3*n2*n1] */, 
			 float* c45         /* [n3*n2*n1] */,
			 float* c46         /* [n3*n2*n1] */,
			 float* c55         /* [n3*n2*n1] */, 
			 float* c56         /* [n3*n2*n1] */,
			 float* c66         /* [n3*n2*n1] */,
			 float* medthe      /* [n3*n2*n1] */,
			 float* medphi      /* [n3*n2*n1] */,
			 int order          /* interpolation order */, 
			 bool wantsine      /* want sine computation*/,
			 bool wantpolar      /* want polarization computation*/)
/*< Initialize ray tracing object. 
 * Increasing order increases accuracy but
 decreases efficiency. Recommended values: 3 or 4.
 * vz2, vx2, and q can be changed or deallocated after
 araytrace_init.
 >*/
{
    araytrace rt;
    
    rt = (araytrace) sf_alloc (1,sizeof(*rt));
    
    rt->wantsine = wantsine;
    rt->wantpolar = wantpolar;
    
    rt->nt = nt;
    rt->dt = dt;
    rt->z0 = o[0]; rt->nz = n[0]; rt->dz = d[0];
    rt->y0 = o[1]; rt->ny = n[1]; rt->dy = d[1];
    rt->x0 = o[2]; rt->nx = n[2]; rt->dx = d[2];
    rt->order = order;
    
    rt->rungecount=sf_intalloc(1);
    rt->oldg = sf_floatalloc2(3,nt);

	rt->grd3 = grid3genani_init (n[0], o[0], d[0], // z
				n[1], o[1], d[1], // y
				n[2], o[2], d[2], // x
				c11, c12, c13,
				c14, c15, c16,
				c22, c23, c24,
				c25, c26, c33,
				c34, c35, c36,
				c44, c45, c46,
				c55, c56, c66,
				medthe,medphi,order);
				
    rt->p11 = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->p12 = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->p13 = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->p14 = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->p15 = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->p16 = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->p22 = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->p23 = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->p24 = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->p25 = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->p26 = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->p33 = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->p34 = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->p35 = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->p36 = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->p44 = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->p45 = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->p46 = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->p55 = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->p56 = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->p66 = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->pthe = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    rt->pphi = eno3_init (rt->order, rt->nz, rt->ny, rt->nx);
    
    eno3_set1 (rt->p11, c11);
    eno3_set1 (rt->p12, c12);
    eno3_set1 (rt->p13, c13);
    if (NULL!=c14) eno3_set1 (rt->p14, c14);
    if (NULL!=c15) eno3_set1 (rt->p15, c15);
    if (NULL!=c16) eno3_set1 (rt->p16, c16);
    eno3_set1 (rt->p22, c22);
    eno3_set1 (rt->p23, c23);
    if (NULL!=c24) eno3_set1 (rt->p24, c24);
    if (NULL!=c25) eno3_set1 (rt->p25, c25);
    if (NULL!=c26) eno3_set1 (rt->p26, c26);
    eno3_set1 (rt->p33, c33);
    if (NULL!=c34) eno3_set1 (rt->p34, c34);
    if (NULL!=c35) eno3_set1 (rt->p35, c35);
    if (NULL!=c36) eno3_set1 (rt->p36, c36);
    eno3_set1 (rt->p44, c44);
    if (NULL!=c45) eno3_set1 (rt->p45, c45);
    if (NULL!=c46) eno3_set1 (rt->p46, c46);
    eno3_set1 (rt->p55, c55);
    if (NULL!=c56) eno3_set1 (rt->p56, c56);
    eno3_set1 (rt->p66, c66);
    if (NULL!=medthe) eno3_set1 (rt->pthe, medthe);
    if (NULL!=medphi) eno3_set1 (rt->pphi, medphi);

    return rt;
}

void araytracegenani_close (araytrace rt)
/*< Free internal storage >*/
{
	grid3genani_close (rt->grd3);
    free (rt);
}

int trace_p_aray (araytrace rt  /* ray tracing object */, 
		float* x      /* point location {z,y,x} [dim] */, 
		float* p      /* ray parameter vector [dim] */, 
		float** traj  /* output ray trajectory [nt+1,2*dim] */,
		float** polar   /* polarization vectors [nt+1,dim] */)
/*< Trace a ray.
 * Values of x and p are changed inside the function.
 * The trajectory traj is stored as follows:
 {z0,y0,z1,y1,z2,y2,...} in 2-D
 {z0,y0,x0,z1,y1,x1,...} in 3-D
 * Vector p points in the direction of the ray. 
 The length of the vector is not important.
 Example initialization:
 p[0] = cos(a); p[1] = sin(a) in 2-D, a is between 0 and 2*pi radians
 p[0] = cos(b); p[1] = sin(b)*sin(a); p[2] = sin(b)*cos(a) in 3-D
 b is inclination between 0 and   pi radians
 a is azimuth     between 0 and 2*pi radians
 * The output code for it = trace_ray(...)
 it=0 - ray traced to the end without leaving the grid
 it>0 - ray exited at the top of the grid
 it<0 - ray exited at the side or bottom of the grid
 * The total traveltime along the ray is 
 nt*dt if (it = 0); abs(it)*dt otherwise 
 >*/
{
    int i, dim, it=0, nt;
    float y[6];

    dim = 3;
    nt = rt->nt;

    for (i=0; i < dim; i++) {
	y[i] = x[i];
	y[i+dim] = p[i];
    }

	rt->rungecount[0]=0; // Count the number of runge-kutta solve
    sf_runge_init(2*dim, nt, rt->dt);
    it = sf_ode23_step (y, rt, aniso_p_rhs, term, traj);
    sf_runge_close();
    
    for (i=0; i < dim; i++) {
	x[i] = y[i];
	p[i] = y[i+dim];
    }
    
    if (rt->wantpolar) {
    int time;
    
    for(time=0;time<nt+1;time++) {
        if (time < nt) { 
        	polar[time][0] = rt->oldg[time][2]; polar[time][1] = rt->oldg[time][1]; polar[time][2] = rt->oldg[time][0];
        } else {
        	polar[time][0] = rt->oldg[time-1][2]; polar[time][1] = rt->oldg[time-1][1]; polar[time][2] = rt->oldg[time-1][0]; // Duplicate the last two ray points
        }
    }
    }
    
    if (it > 0 && x[0] > rt->z0) {
	return (-it); /* exit through the side or bottom */
    } else {
	return it;
    }
}

int trace_s1_aray (araytrace rt  /* ray tracing object */, 
		float* x      /* point location {z,y,x} [dim] */, 
		float* p      /* ray parameter vector [dim] */, 
		float** traj  /* output ray trajectory [nt+1,2*dim] */,
		float* sine   /* output proximity to point singularity*/,
		float** polar   /* polarization vectors*/)
/*< Trace a ray.
 * Values of x and p are changed inside the function.
 * The trajectory traj is stored as follows:
 {z0,y0,z1,y1,z2,y2,...} in 2-D
 {z0,y0,x0,z1,y1,x1,...} in 3-D
 * Vector p points in the direction of the ray. 
 The length of the vector is not important.
 Example initialization:
 p[0] = cos(a); p[1] = sin(a) in 2-D, a is between 0 and 2*pi radians
 p[0] = cos(b); p[1] = sin(b)*sin(a); p[2] = sin(b)*cos(a) in 3-D
 b is inclination between 0 and   pi radians
 a is azimuth     between 0 and 2*pi radians
 * The output code for it = trace_ray(...)
 it=0 - ray traced to the end without leaving the grid
 it>0 - ray exited at the top of the grid
 it<0 - ray exited at the side or bottom of the grid
 * The total traveltime along the ray is 
 nt*dt if (it = 0); abs(it)*dt otherwise 
 >*/
{
    int i, dim, it=0, nt;
    float y[6];

    dim = 3;
    nt = rt->nt;

    for (i=0; i < dim; i++) {
	y[i] = x[i];
	y[i+dim] = p[i];
    }

	rt->rungecount[0]=0; // Count the number of runge-kutta solve
    sf_runge_init(2*dim, nt, rt->dt);
    it = sf_ode23_step (y, rt, aniso_s1_rhs, term, traj);
    sf_runge_close();
    
    for (i=0; i < dim; i++) {
	x[i] = y[i];
	p[i] = y[i+dim];
    }
    
    if (rt->wantpolar || rt-> wantsine) {
    /* Compute sin nu/3 for low-symmtey media only*/
    float c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,the,phi;
    float n1, n2, n3;
    int ii, jj, kk, time;
    float xx, yy, zz;
    
    for(time=0;time<nt+1;time++) {
        xx = (traj[time][2]-rt->x0)/rt->dx; ii = xx; xx -= ii;
        yy = (traj[time][1]-rt->y0)/rt->dy; jj = yy; yy -= jj;
        zz = (traj[time][0]-rt->z0)/rt->dz; kk = zz; zz -= kk;
       
        eno3_apply(rt->p11, kk, jj, ii, zz, yy, xx, &c11, NULL, FUNC);
        eno3_apply(rt->p12, kk, jj, ii, zz, yy, xx, &c12, NULL, FUNC);
        eno3_apply(rt->p13, kk, jj, ii, zz, yy, xx, &c13, NULL, FUNC);
        eno3_apply(rt->p14, kk, jj, ii, zz, yy, xx, &c14, NULL, FUNC);
        eno3_apply(rt->p15, kk, jj, ii, zz, yy, xx, &c15, NULL, FUNC);
        eno3_apply(rt->p16, kk, jj, ii, zz, yy, xx, &c16, NULL, FUNC);
        eno3_apply(rt->p22, kk, jj, ii, zz, yy, xx, &c22, NULL, FUNC);
        eno3_apply(rt->p23, kk, jj, ii, zz, yy, xx, &c23, NULL, FUNC);
        eno3_apply(rt->p24, kk, jj, ii, zz, yy, xx, &c24, NULL, FUNC);
        eno3_apply(rt->p25, kk, jj, ii, zz, yy, xx, &c25, NULL, FUNC);
        eno3_apply(rt->p26, kk, jj, ii, zz, yy, xx, &c26, NULL, FUNC);
        eno3_apply(rt->p33, kk, jj, ii, zz, yy, xx, &c33, NULL, FUNC);
        eno3_apply(rt->p34, kk, jj, ii, zz, yy, xx, &c34, NULL, FUNC);
        eno3_apply(rt->p35, kk, jj, ii, zz, yy, xx, &c35, NULL, FUNC);
        eno3_apply(rt->p36, kk, jj, ii, zz, yy, xx, &c36, NULL, FUNC);
        eno3_apply(rt->p44, kk, jj, ii, zz, yy, xx, &c44, NULL, FUNC);
        eno3_apply(rt->p45, kk, jj, ii, zz, yy, xx, &c45, NULL, FUNC);
        eno3_apply(rt->p46, kk, jj, ii, zz, yy, xx, &c46, NULL, FUNC);
        eno3_apply(rt->p55, kk, jj, ii, zz, yy, xx, &c55, NULL, FUNC);
        eno3_apply(rt->p56, kk, jj, ii, zz, yy, xx, &c56, NULL, FUNC);
        eno3_apply(rt->p66, kk, jj, ii, zz, yy, xx, &c66, NULL, FUNC);
        eno3_apply(rt->pthe, kk, jj, ii, zz, yy, xx, &the, NULL, FUNC);
        eno3_apply(rt->pphi, kk, jj, ii, zz, yy, xx, &phi, NULL, FUNC);
        
	    double  Chr[9];  // Lapack SVD array 
	    
        float one = sqrt(traj[time][3]*traj[time][3] + traj[time][4]*traj[time][4] + traj[time][5]*traj[time][5]);
        n1 = traj[time][5]/one; n2 = traj[time][4]/one; n3 = traj[time][3]/one;
        
        
        /*Bond transformation*/
        bond(&phi,&the,&c11,&c12,&c13,&c14,&c15,&c16,&c22,&c23,&c24,&c25,&c26,&c33,&c34,&c35,&c36,&c44,&c45,&c46,&c55,&c56,&c66);
        
        /* Christoffel matrix */
	    Chr[0] = c11*n1*n1 + c66*n2*n2 + c55*n3*n3 + 2*c16*n1*n2 + 2*c15*n1*n3 + 2*c56*n2*n3;
	    Chr[4] = c66*n1*n1 + c22*n2*n2 + c44*n3*n3 + 2*c26*n1*n2 + 2*c46*n1*n3 + 2*c24*n2*n3;
	    Chr[8] = c55*n1*n1 + c44*n2*n2 + c33*n3*n3 + 2*c45*n1*n2 + 2*c35*n1*n3 + 2*c34*n2*n3;
	    Chr[1] = Chr[3] = c16*n1*n1 + c26*n2*n2 + c45*n3*n3 + (c12+c66)*n1*n2 + (c14+c56)*n1*n3 + (c25+c46)*n2*n3; 
	    Chr[2] = Chr[6] = c15*n1*n1 + c46*n2*n2 + c35*n3*n3 + (c14+c56)*n1*n2 + (c13+c55)*n1*n3 + (c36+c45)*n2*n3; 
	    Chr[5] = Chr[7] = c56*n1*n1 + c24*n2*n2 + c34*n3*n3 + (c25+c46)*n1*n2 + (c36+c45)*n1*n3 + (c23+c44)*n2*n3; 

        double A = -(Chr[0]+Chr[4]+Chr[8]);
	    double B = Chr[0]*Chr[4]+Chr[0]*Chr[8]+Chr[4]*Chr[8]-Chr[1]*Chr[1]-Chr[2]*Chr[2]-Chr[5]*Chr[5];
	    double C = Chr[0]*Chr[5]*Chr[5]+Chr[4]*Chr[2]*Chr[2]+Chr[8]*Chr[1]*Chr[1]-Chr[0]*Chr[4]*Chr[8]-2*Chr[1]*Chr[2]*Chr[5];
	    double D = -A*A/3 + B;
	    double Q = 2*pow(A/3,3)-A*B/3+C;
	    double nu = acos(-Q/(2*sqrt(pow(-D/3,3))));
	    double coss = cos(nu/3);
	    sine[time] = sqrt(1-coss*coss);
	    
        if (time < nt) { 
        	polar[time][0] = rt->oldg[time][2]; polar[time][1] = rt->oldg[time][1]; polar[time][2] = rt->oldg[time][0];
        } else {
        	polar[time][0] = rt->oldg[time-1][2]; polar[time][1] = rt->oldg[time-1][1]; polar[time][2] = rt->oldg[time-1][0]; // Duplicate the last two ray points
        }
    }
    }
    
    if (it > 0 && x[0] > rt->z0) {
	return (-it); /* exit through the side or bottom */
    } else {
	return it;
    }
}

int trace_s2_aray (araytrace rt  /* ray tracing object */, 
		float* x      /* point location {z,y,x} [dim] */, 
		float* p      /* ray parameter vector [dim] */, 
		float** traj  /* output ray trajectory [nt+1,2*dim] */,
		float* sine   /* output proximity to point singularity*/,
		float** polar   /* polarization vectors*/)
/*< Trace a ray.
 * Values of x and p are changed inside the function.
 * The trajectory traj is stored as follows:
 {z0,y0,z1,y1,z2,y2,...} in 2-D
 {z0,y0,x0,z1,y1,x1,...} in 3-D
 * Vector p points in the direction of the ray. 
 The length of the vector is not important.
 Example initialization:
 p[0] = cos(a); p[1] = sin(a) in 2-D, a is between 0 and 2*pi radians
 p[0] = cos(b); p[1] = sin(b)*cos(a); p[2] = sin(b)*sin(a) in 3-D
 b is inclination between 0 and   pi radians
 a is azimuth     between 0 and 2*pi radians
 * The output code for it = trace_ray(...)
 it=0 - ray traced to the end without leaving the grid
 it>0 - ray exited at the top of the grid
 it<0 - ray exited at the side or bottom of the grid
 * The total traveltime along the ray is 
 nt*dt if (it = 0); abs(it)*dt otherwise 
 >*/
{
    int i, dim, it=0, nt;
    float y[6];

    dim = 3;
    nt = rt->nt;

    for (i=0; i < dim; i++) {
	y[i] = x[i];
	y[i+dim] = p[i];
    }

	rt->rungecount[0]=0; // Count the number of runge-kutta solve
    sf_runge_init(2*dim, nt, rt->dt);
    it = sf_ode23_step (y, rt, aniso_s2_rhs, term, traj);
    sf_runge_close();
    
    for (i=0; i < dim; i++) {
	x[i] = y[i];
	p[i] = y[i+dim];
    }
    
    if (rt->wantpolar || rt-> wantsine) {
	/* Compute sin nu/3 for low-symmtey media only*/
    float c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,the,phi;
    float n1, n2, n3;
    int ii, jj, kk, time;
    float xx, yy, zz;
    
    for(time=0;time<nt+1;time++) {
        xx = (traj[time][2]-rt->x0)/rt->dx; ii = xx; xx -= ii;
        yy = (traj[time][1]-rt->y0)/rt->dy; jj = yy; yy -= jj;
        zz = (traj[time][0]-rt->z0)/rt->dz; kk = zz; zz -= kk;
       
        eno3_apply(rt->p11, kk, jj, ii, zz, yy, xx, &c11, NULL, FUNC);
        eno3_apply(rt->p12, kk, jj, ii, zz, yy, xx, &c12, NULL, FUNC);
        eno3_apply(rt->p13, kk, jj, ii, zz, yy, xx, &c13, NULL, FUNC);
        eno3_apply(rt->p14, kk, jj, ii, zz, yy, xx, &c14, NULL, FUNC);
        eno3_apply(rt->p15, kk, jj, ii, zz, yy, xx, &c15, NULL, FUNC);
        eno3_apply(rt->p16, kk, jj, ii, zz, yy, xx, &c16, NULL, FUNC);
        eno3_apply(rt->p22, kk, jj, ii, zz, yy, xx, &c22, NULL, FUNC);
        eno3_apply(rt->p23, kk, jj, ii, zz, yy, xx, &c23, NULL, FUNC);
        eno3_apply(rt->p24, kk, jj, ii, zz, yy, xx, &c24, NULL, FUNC);
        eno3_apply(rt->p25, kk, jj, ii, zz, yy, xx, &c25, NULL, FUNC);
        eno3_apply(rt->p26, kk, jj, ii, zz, yy, xx, &c26, NULL, FUNC);
        eno3_apply(rt->p33, kk, jj, ii, zz, yy, xx, &c33, NULL, FUNC);
        eno3_apply(rt->p34, kk, jj, ii, zz, yy, xx, &c34, NULL, FUNC);
        eno3_apply(rt->p35, kk, jj, ii, zz, yy, xx, &c35, NULL, FUNC);
        eno3_apply(rt->p36, kk, jj, ii, zz, yy, xx, &c36, NULL, FUNC);
        eno3_apply(rt->p44, kk, jj, ii, zz, yy, xx, &c44, NULL, FUNC);
        eno3_apply(rt->p45, kk, jj, ii, zz, yy, xx, &c45, NULL, FUNC);
        eno3_apply(rt->p46, kk, jj, ii, zz, yy, xx, &c46, NULL, FUNC);
        eno3_apply(rt->p55, kk, jj, ii, zz, yy, xx, &c55, NULL, FUNC);
        eno3_apply(rt->p56, kk, jj, ii, zz, yy, xx, &c56, NULL, FUNC);
        eno3_apply(rt->p66, kk, jj, ii, zz, yy, xx, &c66, NULL, FUNC);
        eno3_apply(rt->pthe, kk, jj, ii, zz, yy, xx, &the, NULL, FUNC);
        eno3_apply(rt->pphi, kk, jj, ii, zz, yy, xx, &phi, NULL, FUNC);
        
        
	    double  Chr[9];  // Lapack SVD array 
	    
        float one = sqrt(traj[time][3]*traj[time][3] + traj[time][4]*traj[time][4] + traj[time][5]*traj[time][5]);
        n1 = traj[time][5]/one; n2 = traj[time][4]/one; n3 = traj[time][3]/one;
        
        /*Bond transformation*/
        bond(&phi,&the,&c11,&c12,&c13,&c14,&c15,&c16,&c22,&c23,&c24,&c25,&c26,&c33,&c34,&c35,&c36,&c44,&c45,&c46,&c55,&c56,&c66);
        
        /* Christoffel matrix */
	    Chr[0] = c11*n1*n1 + c66*n2*n2 + c55*n3*n3 + 2*c16*n1*n2 + 2*c15*n1*n3 + 2*c56*n2*n3;
	    Chr[4] = c66*n1*n1 + c22*n2*n2 + c44*n3*n3 + 2*c26*n1*n2 + 2*c46*n1*n3 + 2*c24*n2*n3;
	    Chr[8] = c55*n1*n1 + c44*n2*n2 + c33*n3*n3 + 2*c45*n1*n2 + 2*c35*n1*n3 + 2*c34*n2*n3;
	    Chr[1] = Chr[3] = c16*n1*n1 + c26*n2*n2 + c45*n3*n3 + (c12+c66)*n1*n2 + (c14+c56)*n1*n3 + (c25+c46)*n2*n3; 
	    Chr[2] = Chr[6] = c15*n1*n1 + c46*n2*n2 + c35*n3*n3 + (c14+c56)*n1*n2 + (c13+c55)*n1*n3 + (c36+c45)*n2*n3; 
	    Chr[5] = Chr[7] = c56*n1*n1 + c24*n2*n2 + c34*n3*n3 + (c25+c46)*n1*n2 + (c36+c45)*n1*n3 + (c23+c44)*n2*n3; 

        double A = -(Chr[0]+Chr[4]+Chr[8]);
	    double B = Chr[0]*Chr[4]+Chr[0]*Chr[8]+Chr[4]*Chr[8]-Chr[1]*Chr[1]-Chr[2]*Chr[2]-Chr[5]*Chr[5];
	    double C = Chr[0]*Chr[5]*Chr[5]+Chr[4]*Chr[2]*Chr[2]+Chr[8]*Chr[1]*Chr[1]-Chr[0]*Chr[4]*Chr[8]-2*Chr[1]*Chr[2]*Chr[5];
	    double D = -A*A/3 + B;
	    double Q = 2*pow(A/3,3)-A*B/3+C;
	    double nu = acos(-Q/(2*sqrt(pow(-D/3,3))));
	    double coss = cos(nu/3);
	    sine[time] = sqrt(1-coss*coss);
	    
        if (time < nt) { 
        	polar[time][0] = rt->oldg[time][2]; polar[time][1] = rt->oldg[time][1]; polar[time][2] = rt->oldg[time][0];
        } else {
        	polar[time][0] = rt->oldg[time-1][2]; polar[time][1] = rt->oldg[time-1][1]; polar[time][2] = rt->oldg[time-1][0]; // Duplicate the last two ray points
        }
    }
    }
    
    if (it > 0 && x[0] > rt->z0) {
	return (-it); /* exit through the side or bottom */
    } else {
	return it;
    }
}

/*Most coupled via dot product*/
int trace_s1c_aray (araytrace rt  /* ray tracing object */, 
		float* x      /* point location {z,y,x} [dim] */, 
		float* p      /* ray parameter vector [dim] */, 
		float** traj  /* output ray trajectory [nt+1,2*dim] */,
		float* sine   /* output proximity to point singularity*/,
		float** polar   /* polarization vectors*/)
/*< Trace a ray.
 * Values of x and p are changed inside the function.
 * The trajectory traj is stored as follows:
 {z0,y0,z1,y1,z2,y2,...} in 2-D
 {z0,y0,x0,z1,y1,x1,...} in 3-D
 * Vector p points in the direction of the ray. 
 The length of the vector is not important.
 Example initialization:
 p[0] = cos(a); p[1] = sin(a) in 2-D, a is between 0 and 2*pi radians
 p[0] = cos(b); p[1] = sin(b)*sin(a); p[2] = sin(b)*cos(a) in 3-D
 b is inclination between 0 and   pi radians
 a is azimuth     between 0 and 2*pi radians
 * The output code for it = trace_ray(...)
 it=0 - ray traced to the end without leaving the grid
 it>0 - ray exited at the top of the grid
 it<0 - ray exited at the side or bottom of the grid
 * The total traveltime along the ray is 
 nt*dt if (it = 0); abs(it)*dt otherwise 
 >*/
{
    int i, dim, it=0, nt;
    float y[6];

    dim = 3;
    nt = rt->nt;

    for (i=0; i < dim; i++) {
	y[i] = x[i];
	y[i+dim] = p[i];
    }
	
	rt->rungecount[0]=0; // Count the number of runge-kutta solve
    rt->oldg[0][0]=0; rt->oldg[0][1]=0; rt->oldg[0][2]=0;
    
    sf_runge_init(2*dim, nt, rt->dt);
    it = sf_ode23_step (y, rt, aniso_s1c_rhs, term, traj);
    sf_runge_close();
    
    for (i=0; i < dim; i++) {
	x[i] = y[i];
	p[i] = y[i+dim];
    }
    
    if (rt->wantpolar || rt-> wantsine) {
	/* Compute sin nu/3 for low-symmtey media only*/
    float c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,the,phi;
    float n1, n2, n3;
    int ii, jj, kk, time;
    float xx, yy, zz;
    
    for(time=0;time<nt+1;time++) {
        xx = (traj[time][2]-rt->x0)/rt->dx; ii = xx; xx -= ii;
        yy = (traj[time][1]-rt->y0)/rt->dy; jj = yy; yy -= jj;
        zz = (traj[time][0]-rt->z0)/rt->dz; kk = zz; zz -= kk;
       
        eno3_apply(rt->p11, kk, jj, ii, zz, yy, xx, &c11, NULL, FUNC);
        eno3_apply(rt->p12, kk, jj, ii, zz, yy, xx, &c12, NULL, FUNC);
        eno3_apply(rt->p13, kk, jj, ii, zz, yy, xx, &c13, NULL, FUNC);
        eno3_apply(rt->p14, kk, jj, ii, zz, yy, xx, &c14, NULL, FUNC);
        eno3_apply(rt->p15, kk, jj, ii, zz, yy, xx, &c15, NULL, FUNC);
        eno3_apply(rt->p16, kk, jj, ii, zz, yy, xx, &c16, NULL, FUNC);
        eno3_apply(rt->p22, kk, jj, ii, zz, yy, xx, &c22, NULL, FUNC);
        eno3_apply(rt->p23, kk, jj, ii, zz, yy, xx, &c23, NULL, FUNC);
        eno3_apply(rt->p24, kk, jj, ii, zz, yy, xx, &c24, NULL, FUNC);
        eno3_apply(rt->p25, kk, jj, ii, zz, yy, xx, &c25, NULL, FUNC);
        eno3_apply(rt->p26, kk, jj, ii, zz, yy, xx, &c26, NULL, FUNC);
        eno3_apply(rt->p33, kk, jj, ii, zz, yy, xx, &c33, NULL, FUNC);
        eno3_apply(rt->p34, kk, jj, ii, zz, yy, xx, &c34, NULL, FUNC);
        eno3_apply(rt->p35, kk, jj, ii, zz, yy, xx, &c35, NULL, FUNC);
        eno3_apply(rt->p36, kk, jj, ii, zz, yy, xx, &c36, NULL, FUNC);
        eno3_apply(rt->p44, kk, jj, ii, zz, yy, xx, &c44, NULL, FUNC);
        eno3_apply(rt->p45, kk, jj, ii, zz, yy, xx, &c45, NULL, FUNC);
        eno3_apply(rt->p46, kk, jj, ii, zz, yy, xx, &c46, NULL, FUNC);
        eno3_apply(rt->p55, kk, jj, ii, zz, yy, xx, &c55, NULL, FUNC);
        eno3_apply(rt->p56, kk, jj, ii, zz, yy, xx, &c56, NULL, FUNC);
        eno3_apply(rt->p66, kk, jj, ii, zz, yy, xx, &c66, NULL, FUNC);
        eno3_apply(rt->pthe, kk, jj, ii, zz, yy, xx, &the, NULL, FUNC);
        eno3_apply(rt->pphi, kk, jj, ii, zz, yy, xx, &phi, NULL, FUNC);
        
        
	    double  Chr[9];  // Lapack SVD array 
	    
        float one = sqrt(traj[time][3]*traj[time][3] + traj[time][4]*traj[time][4] + traj[time][5]*traj[time][5]);
        n1 = traj[time][5]/one; n2 = traj[time][4]/one; n3 = traj[time][3]/one;
        
        /*Bond transformation*/
        bond(&phi,&the,&c11,&c12,&c13,&c14,&c15,&c16,&c22,&c23,&c24,&c25,&c26,&c33,&c34,&c35,&c36,&c44,&c45,&c46,&c55,&c56,&c66);
        
        /* Christoffel matrix */
	    Chr[0] = c11*n1*n1 + c66*n2*n2 + c55*n3*n3 + 2*c16*n1*n2 + 2*c15*n1*n3 + 2*c56*n2*n3;
	    Chr[4] = c66*n1*n1 + c22*n2*n2 + c44*n3*n3 + 2*c26*n1*n2 + 2*c46*n1*n3 + 2*c24*n2*n3;
	    Chr[8] = c55*n1*n1 + c44*n2*n2 + c33*n3*n3 + 2*c45*n1*n2 + 2*c35*n1*n3 + 2*c34*n2*n3;
	    Chr[1] = Chr[3] = c16*n1*n1 + c26*n2*n2 + c45*n3*n3 + (c12+c66)*n1*n2 + (c14+c56)*n1*n3 + (c25+c46)*n2*n3; 
	    Chr[2] = Chr[6] = c15*n1*n1 + c46*n2*n2 + c35*n3*n3 + (c14+c56)*n1*n2 + (c13+c55)*n1*n3 + (c36+c45)*n2*n3; 
	    Chr[5] = Chr[7] = c56*n1*n1 + c24*n2*n2 + c34*n3*n3 + (c25+c46)*n1*n2 + (c36+c45)*n1*n3 + (c23+c44)*n2*n3; 

        double A = -(Chr[0]+Chr[4]+Chr[8]);
	    double B = Chr[0]*Chr[4]+Chr[0]*Chr[8]+Chr[4]*Chr[8]-Chr[1]*Chr[1]-Chr[2]*Chr[2]-Chr[5]*Chr[5];
	    double C = Chr[0]*Chr[5]*Chr[5]+Chr[4]*Chr[2]*Chr[2]+Chr[8]*Chr[1]*Chr[1]-Chr[0]*Chr[4]*Chr[8]-2*Chr[1]*Chr[2]*Chr[5];
	    double D = -A*A/3 + B;
	    double Q = 2*pow(A/3,3)-A*B/3+C;
	    double nu = acos(-Q/(2*sqrt(pow(-D/3,3))));
	    double coss = cos(nu/3);
	    sine[time] = sqrt(1-coss*coss);
	    
        if (time < nt) { 
        	polar[time][0] = rt->oldg[time][2]; polar[time][1] = rt->oldg[time][1]; polar[time][2] = rt->oldg[time][0];
        } else {
        	polar[time][0] = rt->oldg[time-1][2]; polar[time][1] = rt->oldg[time-1][1]; polar[time][2] = rt->oldg[time-1][0]; // Duplicate the last two ray points
        }
    }
    }
    
    if (it > 0 && x[0] > rt->z0) {
	return (-it); /* exit through the side or bottom */
    } else {
	return it;
    }
}

/*Most coupled via dot product*/
int trace_s2c_aray (araytrace rt  /* ray tracing object */, 
		float* x      /* point location {z,y,x} [dim] */, 
		float* p      /* ray parameter vector [dim] */, 
		float** traj  /* output ray trajectory [nt+1,2*dim] */,
		float* sine   /* output proximity to point singularity*/,
		float** polar   /* polarization vectors*/)
/*< Trace a ray.
 * Values of x and p are changed inside the function.
 * The trajectory traj is stored as follows:
 {z0,y0,z1,y1,z2,y2,...} in 2-D
 {z0,y0,x0,z1,y1,x1,...} in 3-D
 * Vector p points in the direction of the ray. 
 The length of the vector is not important.
 Example initialization:
 p[0] = cos(a); p[1] = sin(a) in 2-D, a is between 0 and 2*pi radians
 p[0] = cos(b); p[1] = sin(b)*cos(a); p[2] = sin(b)*sin(a) in 3-D
 b is inclination between 0 and   pi radians
 a is azimuth     between 0 and 2*pi radians
 * The output code for it = trace_ray(...)
 it=0 - ray traced to the end without leaving the grid
 it>0 - ray exited at the top of the grid
 it<0 - ray exited at the side or bottom of the grid
 * The total traveltime along the ray is 
 nt*dt if (it = 0); abs(it)*dt otherwise 
 >*/
{
    int i, dim, it=0, nt;
    float y[6];

    dim = 3;
    nt = rt->nt;

    for (i=0; i < dim; i++) {
	y[i] = x[i];
	y[i+dim] = p[i];
    }
    
	rt->rungecount[0]=0; // Count the number of runge-kutta solve
    rt->oldg[0][0]=0; rt->oldg[0][1]=0; rt->oldg[0][2]=0;
    
    sf_runge_init(2*dim, nt, rt->dt);
    it = sf_ode23_step (y, rt, aniso_s2c_rhs, term, traj);
    sf_runge_close();
    
    for (i=0; i < dim; i++) {
	x[i] = y[i];
	p[i] = y[i+dim];
    }
    
    if (rt->wantpolar || rt-> wantsine) {
	/* Compute sin nu/3 for low-symmtey media only*/
    float c11,c12,c13,c14,c15,c16,c22,c23,c24,c25,c26,c33,c34,c35,c36,c44,c45,c46,c55,c56,c66,the,phi;
    float n1, n2, n3;
    int ii, jj, kk, time;
    float xx, yy, zz;
    
    for(time=0;time<nt+1;time++) {
        xx = (traj[time][2]-rt->x0)/rt->dx; ii = xx; xx -= ii;
        yy = (traj[time][1]-rt->y0)/rt->dy; jj = yy; yy -= jj;
        zz = (traj[time][0]-rt->z0)/rt->dz; kk = zz; zz -= kk;
       
        eno3_apply(rt->p11, kk, jj, ii, zz, yy, xx, &c11, NULL, FUNC);
        eno3_apply(rt->p12, kk, jj, ii, zz, yy, xx, &c12, NULL, FUNC);
        eno3_apply(rt->p13, kk, jj, ii, zz, yy, xx, &c13, NULL, FUNC);
        eno3_apply(rt->p14, kk, jj, ii, zz, yy, xx, &c14, NULL, FUNC);
        eno3_apply(rt->p15, kk, jj, ii, zz, yy, xx, &c15, NULL, FUNC);
        eno3_apply(rt->p16, kk, jj, ii, zz, yy, xx, &c16, NULL, FUNC);
        eno3_apply(rt->p22, kk, jj, ii, zz, yy, xx, &c22, NULL, FUNC);
        eno3_apply(rt->p23, kk, jj, ii, zz, yy, xx, &c23, NULL, FUNC);
        eno3_apply(rt->p24, kk, jj, ii, zz, yy, xx, &c24, NULL, FUNC);
        eno3_apply(rt->p25, kk, jj, ii, zz, yy, xx, &c25, NULL, FUNC);
        eno3_apply(rt->p26, kk, jj, ii, zz, yy, xx, &c26, NULL, FUNC);
        eno3_apply(rt->p33, kk, jj, ii, zz, yy, xx, &c33, NULL, FUNC);
        eno3_apply(rt->p34, kk, jj, ii, zz, yy, xx, &c34, NULL, FUNC);
        eno3_apply(rt->p35, kk, jj, ii, zz, yy, xx, &c35, NULL, FUNC);
        eno3_apply(rt->p36, kk, jj, ii, zz, yy, xx, &c36, NULL, FUNC);
        eno3_apply(rt->p44, kk, jj, ii, zz, yy, xx, &c44, NULL, FUNC);
        eno3_apply(rt->p45, kk, jj, ii, zz, yy, xx, &c45, NULL, FUNC);
        eno3_apply(rt->p46, kk, jj, ii, zz, yy, xx, &c46, NULL, FUNC);
        eno3_apply(rt->p55, kk, jj, ii, zz, yy, xx, &c55, NULL, FUNC);
        eno3_apply(rt->p56, kk, jj, ii, zz, yy, xx, &c56, NULL, FUNC);
        eno3_apply(rt->p66, kk, jj, ii, zz, yy, xx, &c66, NULL, FUNC);
        eno3_apply(rt->pthe, kk, jj, ii, zz, yy, xx, &the, NULL, FUNC);
        eno3_apply(rt->pphi, kk, jj, ii, zz, yy, xx, &phi, NULL, FUNC);
        
	    double  Chr[9];  // Lapack SVD array 
	    
        float one = sqrt(traj[time][3]*traj[time][3] + traj[time][4]*traj[time][4] + traj[time][5]*traj[time][5]);
        n1 = traj[time][5]/one; n2 = traj[time][4]/one; n3 = traj[time][3]/one;
        
        /*Bond transformation*/
        bond(&phi,&the,&c11,&c12,&c13,&c14,&c15,&c16,&c22,&c23,&c24,&c25,&c26,&c33,&c34,&c35,&c36,&c44,&c45,&c46,&c55,&c56,&c66);
        
        /* Christoffel matrix */
	    Chr[0] = c11*n1*n1 + c66*n2*n2 + c55*n3*n3 + 2*c16*n1*n2 + 2*c15*n1*n3 + 2*c56*n2*n3;
	    Chr[4] = c66*n1*n1 + c22*n2*n2 + c44*n3*n3 + 2*c26*n1*n2 + 2*c46*n1*n3 + 2*c24*n2*n3;
	    Chr[8] = c55*n1*n1 + c44*n2*n2 + c33*n3*n3 + 2*c45*n1*n2 + 2*c35*n1*n3 + 2*c34*n2*n3;
	    Chr[1] = Chr[3] = c16*n1*n1 + c26*n2*n2 + c45*n3*n3 + (c12+c66)*n1*n2 + (c14+c56)*n1*n3 + (c25+c46)*n2*n3; 
	    Chr[2] = Chr[6] = c15*n1*n1 + c46*n2*n2 + c35*n3*n3 + (c14+c56)*n1*n2 + (c13+c55)*n1*n3 + (c36+c45)*n2*n3; 
	    Chr[5] = Chr[7] = c56*n1*n1 + c24*n2*n2 + c34*n3*n3 + (c25+c46)*n1*n2 + (c36+c45)*n1*n3 + (c23+c44)*n2*n3; 

        double A = -(Chr[0]+Chr[4]+Chr[8]);
	    double B = Chr[0]*Chr[4]+Chr[0]*Chr[8]+Chr[4]*Chr[8]-Chr[1]*Chr[1]-Chr[2]*Chr[2]-Chr[5]*Chr[5];
	    double C = Chr[0]*Chr[5]*Chr[5]+Chr[4]*Chr[2]*Chr[2]+Chr[8]*Chr[1]*Chr[1]-Chr[0]*Chr[4]*Chr[8]-2*Chr[1]*Chr[2]*Chr[5];
	    double D = -A*A/3 + B;
	    double Q = 2*pow(A/3,3)-A*B/3+C;
	    double nu = acos(-Q/(2*sqrt(pow(-D/3,3))));
	    double coss = cos(nu/3);
	    sine[time] = sqrt(1-coss*coss);
	    
        if (time < nt) { 
        	polar[time][0] = rt->oldg[time][2]; polar[time][1] = rt->oldg[time][1]; polar[time][2] = rt->oldg[time][0];
        } else {
        	polar[time][0] = rt->oldg[time-1][2]; polar[time][1] = rt->oldg[time-1][1]; polar[time][2] = rt->oldg[time-1][0]; // Duplicate the last two ray points
        }
    }
    }
    
    if (it > 0 && x[0] > rt->z0) {
	return (-it); /* exit through the side or bottom */
    } else {
	return it;
    }
}


/* 	$Id$	 */
