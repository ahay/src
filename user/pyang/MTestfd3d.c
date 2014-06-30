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

bool frsf, verb;
int nz, nx, ny, nb, nzpad, nxpad, nypad, nt, ns;
float dz, dx, dy, fm, dt,c0, c11, c12, c21, c22, c31, c32;


void expand3d(float ***a, float ***b)
/*< expand domain a to b >*/
{
    int iz,ix,i3;

    for         (i3=0;i3<ny;i3++) {
	for     (ix=0;ix<nx;ix++) {
	    for (iz=0;iz<nz;iz++) {
		b[nb+i3][nb+ix][nb+iz] = a[i3][ix][iz];
	    }
	}
    }

    for         (i3=0; i3<nypad; i3++) {
	for     (ix=0; ix<nxpad; ix++) {
	    for (iz=0; iz<nb;    iz++) {
		b[i3][ix][  iz      ] = b[i3][ix][      nb  ];
		b[i3][ix][nzpad-iz-1] = b[i3][ix][nzpad-nb-1];
	    }
	}
    }


    for         (i3=0; i3<nypad; i3++) {
	for     (ix=0; ix<nb;    ix++) {
	    for (iz=0; iz<nzpad; iz++) {
		b[i3][      ix  ][iz] = b[i3][ nb  	][iz];
		b[i3][nxpad-ix-1][iz] = b[i3][nxpad-nb-1][iz];
	    }
	}
    }

    for         (i3=0; i3<nb;    i3++) {
	for     (ix=0; ix<nxpad; ix++) {
	    for (iz=0; iz<nzpad; iz++) {
		b[      i3  ][ix][iz] = b[      nb  ][ix][iz];
		b[nypad-i3-1][ix][iz] = b[nypad-nb-1][ix][iz];
	    }
	}
    }
}

void window3d(float ***a, float ***b)
/*< window b domain to a >*/
{
    int i3, ix, iz;

    for         (i3=0;i3<ny;i3++) {
	for     (ix=0;ix<nx;ix++) {
	    for (iz=0;iz<nz;iz++) {
		a[i3][ix][iz]=b[nb+i3][nb+ix][nb+iz];
	    }
	}
    }
}


void add_source(float ***u1, float *wlt, int **Szxy, int ns, bool add)
/*< add source >*/
{
	int is, sz, sx, sy;

	for (is=0; is<ns; is++)
	{
		sz=Szxy[is][0]+nb;
		sx=Szxy[is][1]+nb;
		sy=Szxy[is][2]+nb;

		if (add) u1[sy][sx][sz] += wlt[is];
		else 	u1[sy][sx][sz] -= wlt[is];
	}
}


void step_forward(float ***u0, float ***u1, float ***vv)
/*< step forward >*/
{
	int iz, ix, iy;
	float ua;

	for(iy=2; iy<nypad-2; iy++) 
	for(ix=2; ix<nxpad-2; ix++)
	for(iz=2; iz<nzpad-2; iz++) 
	{			    
		/* 4th order Laplacian operator */
		ua=c0 * u1[iy  ][ix  ][iz  ] + 
			c11*(u1[iy  ][ix  ][iz-1] + u1[iy  ][ix  ][iz+1]) +
			c12*(u1[iy  ][ix  ][iz-2] + u1[iy  ][ix  ][iz+2]) +
			c21*(u1[iy  ][ix-1][iz  ] + u1[iy  ][ix+1][iz  ]) +
			c22*(u1[iy  ][ix-2][iz  ] + u1[iy  ][ix+2][iz  ]) +
			c31*(u1[iy-1][ix  ][iz  ] + u1[iy+1][ix  ][iz  ]) +
			c32*(u1[iy-2][ix  ][iz  ] + u1[iy+2][ix  ][iz  ]) ;
		u0[iy][ix][iz] = 2*u1[iy][ix][iz]-u0[iy][ix][iz]+ua * vv[iy][ix][iz];
	}
}

void sponge3d_apply(float  ***uu, float *spo)
/*< apply boundary sponge >*/
{
    int iz,ix,iy,ib,ibz,ibx,iby;
    float w;

    for(ib=0; ib<nb; ib++) {
	w = spo[ib];

	ibz = nzpad-ib-1;
	for    (iy=0; iy<nypad; iy++) {
	    for(ix=0; ix<nxpad; ix++) {
		uu[iy][ix][ib ] *= w; /* z min */
		uu[iy][ix][ibz] *= w; /* z max */
	    }
	}

	ibx = nxpad-ib-1;
	for    (iy=0; iy<nypad; iy++) {
	    for(iz=0; iz<nzpad; iz++) {
		uu[iy][ib ][iz] *= w; /* x min */
		uu[iy][ibx][iz] *= w; /* x max */
	    }
	}
	
	iby = nypad-ib-1;
	for    (ix=0; ix<nxpad; ix++) {
	    for(iz=0; iz<nzpad; iz++) {
		uu[ib ][ix][iz] *= w; /* y min */
		uu[iby][ix][iz] *= w; /* y max */
	    }
	}

    }
}


void free_surface(float ***u0, float ***u1)
/*< handle free surface at the top >*/
{
	int iz,ix,iy;
	if(frsf) /* free surface */
	{ 
		for(iy=0; iy<nypad; iy++) 
		for(ix=0; ix<nxpad; ix++)
		for(iz=0; iz<nb; iz++) 
		{
			u0[iy][ix][iz]=0;
			u1[iy][ix][iz]=0;
		}
	}
}

int main(int argc, char* argv[])
{
	int iz, ix, iy, is, it, ib, kt;
	float t;
	float *wlt, *bndr;
	int **Szxy;
	float ***v0,***vv, ***u0, ***u1,***ptr;
	sf_file Fv, Fw;

    	sf_init(argc,argv);
	Fv = sf_input("in");
	Fw = sf_output("out");

    	if(!sf_getbool("verb",&verb)) verb=false;    /* verbosity */
    	if (!sf_histint(Fv,"n1",&nz)) sf_error("No n1= in input");
    	if (!sf_histint(Fv,"n2",&nx)) sf_error("No n2= in input");
    	if (!sf_histint(Fv,"n3",&ny)) sf_error("No n3= in input");
    	if (!sf_histfloat(Fv,"d1",&dz)) sf_error("No d1= in input");
    	if (!sf_histfloat(Fv,"d2",&dx)) sf_error("No d2= in input");
    	if (!sf_histfloat(Fv,"d3",&dy)) sf_error("No d3= in input");
    	if(!sf_getbool("verb",&verb)) verb=false;/* verbosity */
    	if(!sf_getbool("frsf",&frsf)) frsf=false;/* free surface or not */
   	if (!sf_getint("nt",&nt))  sf_error("nt required");
    	if (!sf_getint("kt",&kt)) sf_error("kt required");/* record wavefield at time kt */
	if (kt>nt) sf_error("make sure kt<=nt");
   	if (!sf_getint("ns",&ns))  ns=1;
   	if (!sf_getint("nb",&nb))  nb=30;
   	if (!sf_getfloat("dt",&dt))  sf_error("dt required");
   	if (!sf_getfloat("fm",&fm))  fm=20;

	sf_putint(Fw,"n1",nz);
	sf_putint(Fw,"n2",nx);
	sf_putint(Fw,"n3",ny);

	nzpad=nz+2*nb;
	nxpad=nx+2*nb;
	nypad=ny+2*nb;
	t = 1.0/(dz*dz);
	c11 = 4.0*t/3.0;
	c12=  -t/12.0;
	t = 1.0/(dx*dx);
	c21 = 4.0*t/3.0;
	c22=  -t/12.0;
	t = 1.0/(dy*dy);
	c31 = 4.0*t/3.0;
	c32=  -t/12.0;
	c0  = -2.0 * (c11 + c12 + c21 + c22 +c31 + c32);

	wlt=sf_floatalloc(nt);
	bndr=sf_floatalloc(nb);
	Szxy=sf_intalloc2(3,ns);// source position

    	v0=sf_floatalloc3(nz,nx,ny);
	vv=sf_floatalloc3(nzpad, nxpad, nypad);
	u0=sf_floatalloc3(nzpad, nxpad, nypad);
	u1=sf_floatalloc3(nzpad, nxpad, nypad);
	for(it=0;it<nt;it++){
		t=SF_PI*fm*(it*dt-1.0/fm);t=t*t;
		wlt[it]=(1.0-2.0*t)*expf(-t);
	}
	for(ib=0;ib<nb;ib++){
		t=0.015*(nb-ib);
		bndr[ib]=expf(-t*t);
	}
	for(is=0; is<ns; is++)
	{
		Szxy[is][0]=nz/2;//iz, boundary excluded
		Szxy[is][1]=nx/2;//ix, boundary excluded
		Szxy[is][2]=ny/2;//iy, boundary excluded
	}
	sf_floatread(v0[0][0],nz*nx*ny,Fv);
	expand3d(v0, vv);
	memset(u0[0][0],0,nzpad*nxpad*nypad*sizeof(float));
	memset(u1[0][0],0,nzpad*nxpad*nypad*sizeof(float));

    	/* precompute vv^2 * dt^2 */
    	for(iy=0; iy<nypad; iy++) 
	for(ix=0; ix<nxpad; ix++) 
	for(iz=0; iz<nzpad; iz++) 
	{
		t= vv[iy][ix][iz]*dt;
		vv[iy][ix][iz] = t*t;
    	}

	for(it=0; it<nt; it++)
	{
		add_source(u1, &wlt[it], Szxy, 1, true);
		step_forward(u0, u1, vv);
		ptr=u0;u0=u1;u1=ptr;

		sponge3d_apply(u0,bndr);
		sponge3d_apply(u1,bndr);
		free_surface(u0, u1);

		if (it==kt)
		{
			window3d(v0,u0);
			sf_floatwrite(v0[0][0],nz*nx*ny,Fw);	
		}
		if (verb) sf_warning("%d of %d;", it, nt);
	}

    
	free(wlt);
	free(bndr);
	free(*Szxy);free(Szxy);
	free(**v0);free(*v0);free(v0);
	free(**vv);free(*vv);free(vv);
	free(**u0);free(*u0);free(u0);
	free(**u1);free(*u1);free(u1);

    	exit (0);
}
