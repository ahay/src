/* 2-D Fourth-order Finite-difference wave extrapolation with ABC */
/*
  Copyright (C) 2013 University of Texas at Austin
  
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
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
int main(int argc, char* argv[]) 
{
    bool verb;                               /* verbose flag */
    sf_axis at,az,ax;                        /* cube axes */
    int nt, nz, nx, it, iz, ix;              /* dimension and index variables */
    float dt, dz, dx, idz, idx;              /* intervals */
    float **next,**curr,**lapl,**dcur,**dprv;/* tmp arrays */
    float *ww,**vv,**rr;                     /* I/O arrays*/
    float *wb, c;                            /* absorbing boundary parameters */
    int nb, nzb, nxb;                        /* absorbing boundary length */
    sf_file Fw=NULL,Fv=NULL,Fr=NULL,Fo=NULL; /* I/O files */
    /* Laplacian coefficients */
    float c0=-30./12.,c1=+16./12.,c2=- 1./12.;

    sf_init(argc,argv);

    if(! sf_getbool("verb",&verb)) verb=0;

    /* setup I/O files */
    Fw = sf_input("in");   /* source wavlet*/
    Fv = sf_input("vel");  /* velocity */
    Fr = sf_input("ref");  /* source location */
    Fo = sf_output("out"); /* output wavefield */

    /* Read/Write axes */
    at = sf_iaxa(Fw,1); nt = sf_n(at); dt = sf_d(at);
    az = sf_iaxa(Fv,1); nz = sf_n(az); dz = sf_d(az);
    ax = sf_iaxa(Fv,2); nx = sf_n(ax); dx = sf_d(ax);
    sf_oaxa(Fo,az,1); 
    sf_oaxa(Fo,ax,2); 
    sf_oaxa(Fo,at,3);

    idz = 1/(dz*dz);
    idx = 1/(dx*dx);
 
    if (!sf_getint("nb",&nb)) nb=30;    /* boundary length */
    if (!nb) sf_warning("No ABC applied! \n");
    if (!sf_getfloat("c",&c)) c = 0.01; /* decaying parameter */
    nzb = nz + 2*nb + 4;
    nxb = nx + 2*nb + 4;

    /* allocate memory for wavelet, velocity & source location */
    ww = sf_floatalloc(nt);     sf_floatread(ww,nt,Fw);
    vv = sf_floatalloc2(nzb,nxb);
    rr = sf_floatalloc2(nzb,nxb);

    /*input & extend velocity model*/
    for (ix=nb+2; ix<nx+nb+2; ix++){
        sf_floatread(vv[ix]+nb+2,nz,Fv);
	sf_floatread(rr[ix]+nb+2,nz,Fr);
         for (iz=0; iz<nb+2; iz++){
             vv[ix][iz] = vv[ix][nb+2];
             vv[ix][nz+nb+2+iz] = vv[ix][nz+nb+1]; //nb+1=nb+2-1
             rr[ix][iz] = rr[ix][nb+2];
             rr[ix][nz+nb+2+iz] = rr[ix][nz+nb+1];
         }
    }

    for (ix=0; ix<nb+2; ix++){
        for (iz=0; iz<nzb; iz++){
            vv[ix][iz] = vv[nb+2][iz];
            vv[nx+nb+2+ix][iz] = vv[nx+nb+1][iz];
            rr[ix][iz] = rr[nb+2][iz];
            rr[nx+nb+2+ix][iz] = rr[nx+nb+1][iz];
	}
    }

    /* allocate temporary arrays */
    next=sf_floatalloc2(nzb,nxb);
    curr=sf_floatalloc2(nzb,nxb);
    dprv=sf_floatalloc2(nzb,nxb);
    dcur=sf_floatalloc2(nzb,nxb);
    lapl=sf_floatalloc2(nzb,nxb);
 
    for (iz=0; iz<nzb; iz++) {
	for (ix=0; ix<nxb; ix++) {
	    next[ix][iz]=0;
	    curr[ix][iz]=0;
	    dprv[ix][iz]=0;
	    dcur[ix][iz]=0;
	    lapl[ix][iz]=0;
	}
    }
 
    wb =  nb? sf_floatalloc(nb): NULL;
    if (nb) { /* find absorbing coefficients */
	int ib;
	for(ib=0; ib<nb; ib++){
	    wb[ib]=exp(-c*c*(nb-1-ib)*(nb-1-ib));
	}  
    }

#ifdef _OPENMP
	sf_warning("Max num of threads=%d \n",omp_get_max_threads());
	//omp_set_num_threads(2);
	sf_warning("Num of cores available=%d \n",omp_get_num_procs());
#pragma omp parallel
	{
	    sf_warning("Number of threads in parallel region = %d \n",omp_get_num_threads());
	}
#endif

    /* Main loop: propagation in time */
    if(verb) fprintf(stderr,"\n");
    for (it=0; it < nt; it++) {
	if(verb) fprintf(stderr,"\b\b\b\b\b%d",it);

#ifdef _OPENMP
#pragma omp parallel for private(iz,ix) shared(curr,lapl) //num_threads(4)
#endif
	/* 4th order laplacian */
	for (iz=2; iz<nzb-2; iz++) {
	    for (ix=2; ix<nxb-2; ix++) {
		lapl[ix][iz] = 
		    c0* curr[ix  ][iz  ] * (idx+idz) + 
		    c1*(curr[ix-1][iz  ] + curr[ix+1][iz  ])*idx +
		    c2*(curr[ix-2][iz  ] + curr[ix+2][iz  ])*idx +
		    c1*(curr[ix  ][iz-1] + curr[ix  ][iz+1])*idz +
		    c2*(curr[ix  ][iz-2] + curr[ix  ][iz+2])*idz;
	  
	    }
	}

	/* time step */
#ifdef _OPENMP
#pragma omp parallel for private(iz,ix) shared(dcur,dprv,lapl,next,curr)
#endif
	for (iz=0; iz<nzb; iz++) {
	    for (ix=0; ix<nxb; ix++) {
		dcur[ix][iz] = dprv[ix][iz] + (lapl[ix][iz]+ww[it]*rr[ix][iz])*vv[ix][iz]*vv[ix][iz]*dt;
		next[ix][iz] = curr[ix][iz] + dcur[ix][iz]*dt;
//		next[ix][iz] = 2*curr[ix][iz] - prev[ix][iz] + (lapl[ix][iz]+ww[it]*rr[ix][iz])*vv[ix][iz]*vv[ix][iz]*dt2; 
//		if (!nb) {}
		    dprv[ix][iz] = dcur[ix][iz]; 
		    curr[ix][iz] = next[ix][iz];
 	    }
	}
	
	if (nb) {
	    /* absorbing boundary */
	    for (ix=0; ix < nb; ix++) {  
		for (iz=2; iz < nzb-2; iz++) {
		    curr[ix+2][iz] *= wb[ix];
		    curr[ix+nb+nx+2][iz] *= wb[nb-1-ix];
		    dprv[ix+2][iz] *= wb[ix];
		    dprv[ix+nb+nx+2][iz] *= wb[nb-1-ix];
		}
	    }
	    for (ix=2; ix < nxb-2; ix++) {  
		for (iz=0; iz < nb; iz++) {
		    curr[ix][iz+2] *= wb[iz];
		    curr[ix][iz+nz+nb+2] *= wb[nb-1-iz];
		    dprv[ix][iz+2] *= wb[iz];
		    dprv[ix][iz+nz+nb+2] *= wb[nb-1-iz];
		}
	    }
/*	    
	    for (iz=0; iz < nzb; iz++) {  
		for(ix=0; ix < nxb; ix++) {
		    dprv[ix][iz] = dcur[ix][iz]; 
		    curr[ix][iz] = next[ix][iz]; 
		}
	    }
*/
	}
	
	for (ix=nb+2; ix<nx+nb+2; ix++){
	    sf_floatwrite(curr[ix]+nb+2,nz,Fo);
	}  
    }
    
    if(verb) fprintf(stderr,"\n");    
    sf_close();
    exit(0); 
}           
           
