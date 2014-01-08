/* 3D acoustic time-domain FD modeling
4th order in space, 2nd order in time. 
Sponge absorbing boundary condition.
*/
/*
  Copyright (C) 2013 Xi'an Jiaotong University, UT Austin (Pengliang Yang)
  Copyright (C) 2008 Colorado School of Mines
  
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
#ifdef _OPENMP
#include <omp.h>
#endif

#include "fd3dutil.h"
#include <time.h>


int main(int argc, char* argv[])
{
	int iz, ix, iy;
	sf_file Fv, Fw;

    	sf_init(argc,argv);
#ifdef _OPENMP
    	omp_init();
#endif
	Fv = sf_input("in");
	Fw = sf_output("out");

    	fdm3d fdm;
    	fdm = (fdm3d) sf_alloc(1,sizeof(*fdm));

    	if (!sf_histint(Fv,"n1",&fdm->nz)) sf_error("No n1= in input");
    	if (!sf_histint(Fv,"n2",&fdm->nx)) sf_error("No n2= in input");
    	if (!sf_histint(Fv,"n3",&fdm->ny)) sf_error("No n3= in input");
    	if (!sf_histfloat(Fv,"d1",&fdm->dz)) sf_error("No d1= in input");
    	if (!sf_histfloat(Fv,"d2",&fdm->dx)) sf_error("No d2= in input");
    	if (!sf_histfloat(Fv,"d3",&fdm->dy)) sf_error("No d3= in input");
	sf_putint(Fw,"n1",fdm->nz);
	sf_putint(Fw,"n2",fdm->nx);
	sf_putint(Fw,"n3",fdm->ny);


	int nt=300;
	float dt=0.001;
	float fm=20;
	int ns=1;

    	fdm->free=false;
    	fdm->verb=false;

    	fdm->nb=30;
    	fdm->oz=0.0;
    	fdm->ox=0.0;
    	fdm->oy=0.0;

    	fdm->nzpad=fdm->nz+2*fdm->nb;
    	fdm->nxpad=fdm->nx+2*fdm->nb;
    	fdm->nypad=fdm->ny+2*fdm->nb;
	
    	fdm->ozpad=fdm->oz-fdm->nb*fdm->dz;
    	fdm->oxpad=fdm->ox-fdm->nb*fdm->dx;
    	fdm->oypad=fdm->oy-fdm->nb*fdm->dy;

    	fdm->ompchunk=1;


	float *wlt=(float*)malloc(nt*sizeof(float));
	for(int it=0;it<nt;it++){
		float a=SF_PI*fm*(it*dt-1.0/fm);a*=a;
		wlt[it]=(1.0-2.0*a)*expf(-a);
	}

	sponge spo = sponge_make(fdm->nb);


	int **Szxy=sf_intalloc2(3,ns);// source position
	for(int is=0; is<ns; is++)
	{
		Szxy[is][0]=0;//iz, boundary excluded
		Szxy[is][1]=25;//ix, boundary excluded
		Szxy[is][2]=25;//iy, boundary excluded
	}

	float ***v0,***vv, ***u0, ***u1,***ptr=NULL;
    	v0=sf_floatalloc3(fdm->nz,fdm->nx,fdm->ny);
	vv=sf_floatalloc3(fdm->nzpad, fdm->nxpad, fdm->nypad);
	u0=sf_floatalloc3(fdm->nzpad, fdm->nxpad, fdm->nypad);
	u1=sf_floatalloc3(fdm->nzpad, fdm->nxpad, fdm->nypad);
	sf_floatread(v0[0][0],fdm->nz*fdm->nx*fdm->ny,Fv);

	expand3d(v0, vv, fdm);
	memset(u0[0][0],0,fdm->nzpad*fdm->nxpad*fdm->nypad*sizeof(float));
	memset(u1[0][0],0,fdm->nzpad*fdm->nxpad*fdm->nypad*sizeof(float));

    	/* precompute vv^2 * dt^2 */
    	for(iy=0; iy<fdm->nypad; iy++) 
	for(ix=0; ix<fdm->nxpad; ix++) 
	for(iz=0; iz<fdm->nzpad; iz++) 
	{
		float a= vv[iy][ix][iz] * dt;
		vv[iy][ix][iz] = a*a;
    	}

	fd3_init(fdm);
	for(int it=0; it<nt;it++)
	{
		add_source(u1, &wlt[it], Szxy, 1, fdm, true);
		step_forward(u0, u1, vv, fdm);
		ptr=u0;u0=u1;u1=ptr;

		sponge3d_apply(u0,spo,fdm);
		sponge3d_apply(u1,spo,fdm);
		free_surface(u0, u1, fdm);

		if (it==120)
		{
			window3d(v0,u0,fdm);
			sf_floatwrite(v0[0][0],fdm->nz*fdm->nx*fdm->ny,Fw);	
		}
	}

    
	free(wlt);
	free(spo);
	free(*Szxy);free(Szxy);
	free(**v0);free(*v0);free(v0);
	free(**vv);free(*vv);free(vv);
	free(**u0);free(*u0);free(u0);
	free(**u1);free(*u1);free(u1);

    	exit (0);
}
