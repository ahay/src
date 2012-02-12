/* simple ray-tracing */
/*
  Copyright (C) 2011 KAUST

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

/* Reference : Theory of Seismic Imaging by Scales (Equation 11.1.4) */

#include <rsf.h>

#define NU 4

int nz; float z0,dz,*z;
int nx; float x0,dx,*x;
int nt; float dt;

float *v,*vz,*vx;

bool arc;

/* bilinear interpolation */
float int2(const float *f, /* [nx][nz] */
		   const float *c  /* {z,x} target */)
{
	int   iz,ix;
	float zo,xo,wx,wz;
	
	zo = c[0]; iz = (int)((zo - z0) / dz);
	xo = c[1]; ix = (int)((xo - x0) / dx);

	iz = (iz < 0   ) ? 0    : iz;
	iz = (iz > nz-2) ? nz-2 : iz;
	ix = (ix < 0   ) ? 0    : ix;
	ix = (ix > nx-2) ? nx-2 : ix;

	wx = (xo - x[ix]) / dx;
	wz = (zo - z[iz]) / dz;
	
	return												\
		(1.-wx)*(1.-wz) * f[ ix   *nz+ iz   ] +			\
		wx     *(1.-wz) * f[(ix+1)*nz+ iz   ] +			\
		(1.-wx)*    wz  * f[ ix   *nz+(iz+1)] +			\
		wx     *    wz  * f[(ix+1)*nz+(iz+1)];

}

/* rhs */
void evalf(float *f,      
		   const float *u)
{
	float vex,vez,ve,c[2];

	c[0] = u[2];
	c[1] = u[0];
	ve  = int2(v ,c);
	vex = int2(vx,c);
	vez = int2(vz,c);

	f[0] = u[1] * ve;
	f[2] = u[3] * ve;
	f[1] = (u[1]*u[1] - 1.) * vex + u[1]*u[3] * vez;
	f[3] = (u[3]*u[3] - 1.) * vez + u[1]*u[3] * vex;
	
	if (arc)
		for (int i=0; i < NU; i++) f[i] /= ve;
}

/* ode solver */
/* u[0] : x
   u[1] : dx/dt 
   u[2] : z
   u[3] : dz/dt 
*/
void rk4(float *u,       /* [NU*nt] */
		 const float *u0 /* [NU] */)
{
	int i,it;
	float uu[NU],*uo,w1,w2,w3,w4,k1[NU],k2[NU],k3[NU],k4[NU];

	w1 = w4 = dt / 6.;
	w2 = w3 = dt / 3.;

	for (i=0; i < NU; i++) u[i] = u0[i];

	for (uo=u, it=1; it < nt; it++, uo += NU) {
		evalf(k1,uo);
		for (i=0; i < NU; i++)
			uu[i] = uo[i] + .5*dt*k1[i];
		evalf(k2,uu);
		for (i=0; i < NU; i++)
			uu[i] = uo[i] + .5*dt*k2[i];
		evalf(k3,uu);
		for (i=0; i < NU; i++)
			uu[i] = uo[i] + dt*k3[i];
		evalf(k4,uu);

		for (i=0; i < NU; i++)
			u[it*NU + i] = uo[i] + w1*k1[i] + w2*k2[i] + w3*k3[i] + w4*k4[i];
	}
}

int main(int argc,char *argv[])
{
	int i,j,nxz,ns;
	float theta,s0,ds,zs,u0[NU],*u,*s;
	sf_file F,Fx,Fz,Fv,Fvx,Fvz;

	sf_init(argc,argv);

	Fv = sf_input("in");
	if (!sf_histint  (Fv,"n1",&nz) ||
		!sf_histint  (Fv,"n2",&nx) ||
		!sf_histfloat(Fv,"d1",&dz) ||
		!sf_histfloat(Fv,"d2",&dx))
		sf_error("Need n1= n2= in input");
	if (!sf_histfloat(Fv,"o1",&z0)) z0 = 0.;
	if (!sf_histfloat(Fv,"o2",&x0)) x0 = 0.;

	Fvx = sf_input("vx");
	if (!sf_histint(Fvx,"n1",&i) || i != nz ||
		!sf_histint(Fvx,"n2",&i) || i != nx)
		sf_error("Need n1=%d n2=%d in vx",nz,nx);

	Fvz = sf_input("vz");
	if (!sf_histint(Fvz,"n1",&i) || i != nz ||
		!sf_histint(Fvz,"n2",&i) || i != nx)
		sf_error("Need n1=%d n2=%d in vz",nz,nx);

	x = sf_floatalloc(nx);
	z = sf_floatalloc(nz);
	for (x[0]=x0, i=1; i < nx; i++) x[i]=x[i-1] + dx;
	for (z[0]=z0, i=1; i < nz; i++) z[i]=z[i-1] + dz;

	if (!sf_getfloat("ds",&ds)) ds = 1;       /* shot x */
	if (!sf_getint  ("ns",&ns)) ns = 1;       /* shot x */
	if (!sf_getfloat("s0",&s0)) s0 = x[nx/2]; /* shot x */
	if (!sf_getfloat("zs",&zs)) zs = z[nz/2]; /* shot z */
	if (!sf_getint  ("nt",&nt)) nt = 1;
	if (!sf_getfloat("dt",&dt)) dt = 1;
	if (!sf_getfloat("theta",&theta)) theta = 0; /* takeoff angle (degree) */
	if (!sf_getbool ("arc",&arc)) arc = false;
    /* if true, arclength; if false, traveltime */

	s = sf_floatalloc(ns);
	for (s[0]=s0, i=1; i < ns; i++) s[i]=s[i-1] + ds;

	Fx = sf_output("x");
	Fz = sf_output("z");
	for (i=0; i < 2; i++) {
		F = i ? Fx : Fz;
		sf_putint   (F,"n1",nt);
		sf_putint   (F,"n2",ns);
		sf_putfloat (F,"d1",dt);
		sf_putfloat (F,"d2",ds);
		sf_putfloat (F,"o1",0.);
		sf_putfloat (F,"o2",s0);
		sf_putstring(F,"label1",arc ? "arc" : "time");
		sf_putstring(F,"unit1" ,arc ? "m"   : "sec");
		sf_putstring(F,"label2","shot");
		sf_putstring(F,"unit2","");
	}

	nxz = nx * nz;
	v  = sf_floatalloc(nxz);
	vz = sf_floatalloc(nxz);
	vx = sf_floatalloc(nxz);
	sf_floatread(v ,nxz,Fv );
	sf_floatread(vz,nxz,Fvz);
	sf_floatread(vx,nxz,Fvx);

	u = sf_floatalloc(NU*nt);
	u0[2] = zs;
	u0[1] = sinf(theta * M_PI/180.);
	u0[3] = cosf(theta * M_PI/180.);
	for (i=0; i < ns; i++) {
		sf_warning("ray %d;",i+1);
		u0[0] = s[i];
		rk4(u,u0);

		for (j=0; j < nt; j++) sf_floatwrite(&u[j*NU]  ,1,Fx);
		for (j=0; j < nt; j++) sf_floatwrite(&u[j*NU+2],1,Fz);
	}
	sf_warning("\n");

	sf_fileclose(Fv);	sf_fileclose(Fvx);
	sf_fileclose(Fvz); 	sf_fileclose(Fx);
	sf_fileclose(Fz);
	free(u); free(x); free(z);
	free(v); free(vx); free(vz);
	return 0;
}
