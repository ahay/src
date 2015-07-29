/*2D elastic time-domain FD modeling with GPU*/

/*
  Authors: Robin M. Weiss and Jeffrey Shragge

  Use of this code is freely available. In publications, please reference the paper: 
  Weiss and Shragge, "Solving 3D Anisotropic Elastic Wave Equations on Parallel 
  GPU Devices", GEOPHYSICS. http://software.seg.org/2012/0063

  This code is a GPU-enabled version of the ewefd2d module from the Madagascar
  software package (see: http://www.reproducibility.org).  It implements a 2D
  Finite-Difference Time Domain solver for the elastice wave equation with 
  2nd- and 8th- order temporal and spatial accuracy, respectively.  For more 
  information, see (Weiss and Shragge, "Solving 3D Anisotropic Elastic Wave 
  Equations on Parallel GPU Devices", GEOPHYSICS. http://software.seg.org/2012/0063)
*/

/*
  Copyright (C) 2012 University of Western Australia
  
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
#include <cuda.h>
#include <cuda_runtime_api.h>

extern "C" {
#include <rsf.h>
}

#include "fdutil.c"
#include "ewefd2d_kernels.cu"

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define NOP 4 /* derivative operator half-size */


// checks the current GPU device for an error flag and prints to stderr
static void sf_check_gpu_error (const char *msg) {
    cudaError_t err = cudaGetLastError ();
    if (cudaSuccess != err)
        sf_error ("Cuda error: %s: %s", msg, cudaGetErrorString (err));
}

// entry point
int main(int argc, char* argv[]) {
	
    bool verb,fsrf,snap,ssou,dabc;
    int  jsnap,ntsnap,jdata;

    /* I/O files */
    sf_file Fwav=NULL; /* wavelet   */
    sf_file Fsou=NULL; /* sources   */
    sf_file Frec=NULL; /* receivers */
    sf_file Fccc=NULL; /* velocity  */
    sf_file Fden=NULL; /* density   */
    sf_file Fdat=NULL; /* data      */
    sf_file Fwfl=NULL; /* wavefield */

    /* cube axes */
    sf_axis at,ax,az;
    sf_axis as,ar,ac;

    int     nt,nz,nx,ns,nr,nc,nb;
    int     it,iz,ix;
    float   dt,dz,dx,idz,idx;

    /* FDM structure */
    fdm2d    fdm=NULL;
    abcone2d abcs=NULL;
    sponge   spo=NULL;

    /* I/O arrays */
    float***ww=NULL;           /* wavelet   */
    pt2d   *ss=NULL;           /* sources   */
    pt2d   *rr=NULL;           /* receivers */
    float **dd=NULL;           /* data      */

    /*------------------------------------------------------------*/
    /* orthorombic footprint - 4 coefficients */
    /* c11 c13 
       .   c33 
               c55 */
	float *h_c11, *h_c33, *h_c55, *h_c13;
	float *d_c11, *d_c33, *d_c55, *d_c13;
    float **vs=NULL;

	// density
	float *h_ro, *d_ro;

    /*------------------------------------------------------------*/
    /* displacement: um = U @ t-1; uo = U @ t; up = U @ t+1 */

	float *d_umz, *d_uoz, *d_upz, *d_uaz, *d_utz;
	float *d_umx, *d_uox, *d_upx, *d_uax, *d_utx;
	
	// used for writing wavefield to output file
	float *h_uoz, *h_uox; 
	float **uoz, **uox;
	
    /* stress/strain tensor */ 
	float *d_tzz, *d_tzx, *d_txx;

    
    /*------------------------------------------------------------*/
    /* linear interpolation weights/indices */
    lint2d cs,cr;

    /* Gaussian bell */
    int nbell;

    /* wavefield cut params */
    sf_axis   acz=NULL,acx=NULL;
    int       nqz,nqx;
    float     oqz,oqx;
    float     dqz,dqx;
    float     **uc=NULL;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);


    /*------------------------------------------------------------*/
    /* execution flags */
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("snap",&snap)) snap=false; /* wavefield snapshots flag */
    if(! sf_getbool("free",&fsrf)) fsrf=false; /* free surface flag */
    if(! sf_getbool("ssou",&ssou)) ssou=false; /* stress source */
    if(! sf_getbool("dabc",&dabc)) dabc=false; /* absorbing BC */
    /*------------------------------------------------------------*/


    /*------------------------------------------------------------*/
    /* I/O files */
	Fwav = sf_input ("in" ); /* wavelet   */
    Fccc = sf_input ("ccc"); /* stiffness */
    Fden = sf_input ("den"); /* density   */
    Fsou = sf_input ("sou"); /* sources   */
    Frec = sf_input ("rec"); /* receivers */
    Fdat = sf_output("out"); /* data      */
    if(snap)
	Fwfl = sf_output("wfl"); /* wavefield */
    
    /*------------------------------------------------------------*/


	/*------------------------------------------------------------*/
    /* init GPU */
	int gpu;
	if (! sf_getint("gpu", &gpu)) gpu = 0;	/* ID of the GPU to be used */
	sf_warning("using GPU #%d", gpu);
	cudaSetDevice(gpu);
	cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
	cudaStream_t stream[6];
	for (int i = 0; i < 6; ++i) {
		cudaStreamCreate(&stream[i]);
	}
	

    /*------------------------------------------------------------*/
    /* axes */
    at = sf_iaxa(Fwav,3); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time */
    az = sf_iaxa(Fccc,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* depth */
    ax = sf_iaxa(Fccc,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* space */

    as = sf_iaxa(Fsou,2); sf_setlabel(as,"s"); if(verb) sf_raxa(as); /* sources */
    ar = sf_iaxa(Frec,2); sf_setlabel(ar,"r"); if(verb) sf_raxa(ar); /* receivers */

	nt = sf_n(at); dt = sf_d(at);
    nz = sf_n(az); dz = sf_d(az);
    nx = sf_n(ax); dx = sf_d(ax);

    ns = sf_n(as);
    nr = sf_n(ar);
    /*------------------------------------------------------------*/


	/*------------------------------------------------------------*/
    /* other execution parameters */
    if(! sf_getint("nbell",&nbell)) nbell=5;  /* bell size */
    if(verb) sf_warning("nbell=%d",nbell);
    if(! sf_getint("jdata",&jdata)) jdata=1;	/* extract receiver data every jdata time steps */
    if(snap) {  
		if(! sf_getint("jsnap",&jsnap)) jsnap=nt;	/* save wavefield every jsnap time steps */
    }
    /*------------------------------------------------------------*/


    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;

    fdm=fdutil_init(verb,fsrf,az,ax,nb,1);

	if (nbell * 2 + 1 > 32){
		sf_error("nbell must be <= 15\n"); 
	}

	float *h_bell;
	float *d_bell;
	h_bell = (float*)malloc((2*nbell+1)*(2*nbell+1)*sizeof(float));
	float s = 0.5*nbell;
	for (ix=-nbell;ix<=nbell;ix++) {
		for (iz=-nbell;iz<=nbell;iz++) {
	    	h_bell[(iz + nbell) * (2*nbell+1) + (ix + nbell)] = exp(-(iz*iz+ix*ix)/s);
		}
	}
	cudaMalloc((void**)&d_bell, (2*nbell+1)*(2*nbell+1)*sizeof(float));
	cudaMemcpy(d_bell, h_bell, (2*nbell+1)*(2*nbell+1)*sizeof(float), cudaMemcpyHostToDevice);

    sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if(verb) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if(verb) sf_raxa(ax);
    /*------------------------------------------------------------*/


	/*------------------------------------------------------------*/
    /* 2D vector components */
    nc=2;
    ac=sf_maxa(nc,0,1);
    /*------------------------------------------------------------*/


    /*------------------------------------------------------------*/
    /* setup output data file and arrays*/
    sf_oaxa(Fdat,ar,1);
    sf_oaxa(Fdat,ac,2);

    sf_setn(at,nt/jdata);
    sf_setd(at,dt*jdata);
    sf_oaxa(Fdat,at,3);

    /* setup output wavefield header and arrays*/
    if(snap) {
	
		uoz=sf_floatalloc2(fdm->nzpad,fdm->nxpad);
		uox=sf_floatalloc2(fdm->nzpad,fdm->nxpad);

		h_uoz = (float*)malloc(fdm->nzpad * fdm->nxpad * sizeof(float));
		h_uox = (float*)malloc(fdm->nzpad * fdm->nxpad * sizeof(float));
	
		nqz=sf_n(az);
		nqx=sf_n(ax);
	    
		oqz=sf_o(az);
		oqx=sf_o(ax);
	
		dqz=sf_d(az);
		dqx=sf_d(ax);
	
		acz = sf_maxa(nqz,oqz,dqz); if(verb)sf_raxa(acz);
		acx = sf_maxa(nqx,oqx,dqx); if(verb)sf_raxa(acx);
	
		uc=sf_floatalloc2(sf_n(acz),sf_n(acx));
	
		ntsnap=0;
		for(it=0; it<nt; it++) {
		    if(it%jsnap==0) ntsnap++;
		}
		sf_setn(at,  ntsnap);
		sf_setd(at,dt*jsnap);
		if(verb) sf_raxa(at);
	
		sf_oaxa(Fwfl,acz,1);
		sf_oaxa(Fwfl,acx,2);
		sf_oaxa(Fwfl,ac, 3);
		sf_oaxa(Fwfl,at, 4);
    }


    /*------------------------------------------------------------*/
    /* read source wavelet(s) and copy to GPU (into d_ww) */
    ww=sf_floatalloc3(ns,nc,nt); 
    sf_floatread(ww[0][0],nt*nc*ns,Fwav);

	float *h_ww;
	h_ww = (float*)malloc(ns*nc*nt*sizeof(float));
	for (int t = 0; t < nt; t++){
		for (int c = 0; c < nc; c++){
			for (int s = 0; s < ns; s++){
				h_ww[t * nc * ns + c * ns + s]=ww[t][c][s];
			}
		}
	}
	float *d_ww;
	cudaMalloc((void**)&d_ww, ns*nc*nt*sizeof(float));
	cudaMemcpy(d_ww, h_ww, ns*nc*nt*sizeof(float), cudaMemcpyHostToDevice);
    /*------------------------------------------------------------*/


	/*------------------------------------------------------------*/
	/* data array */
    dd=sf_floatalloc2(nr,nc);
	float *d_dd;
	float *h_dd;
	h_dd = (float*)malloc(nr * nc * sizeof(float));
	cudaMalloc((void**)&d_dd, nr*nc*sizeof(float));
	/*------------------------------------------------------------*/


    /*------------------------------------------------------------*/
    /* setup source/receiver coordinates */
    ss = (pt2d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt2d*) sf_alloc(nr,sizeof(*rr)); 

    pt2dread1(Fsou,ss,ns,2); /* read (x,z) coordinates */
    pt2dread1(Frec,rr,nr,2); /* read (x,z) coordinates */

	/* calculate 2d linear interpolation coefficients for source locations */
    cs = lint2d_make(ns,ss,fdm);
	float *d_Sw00, *d_Sw01, *d_Sw10, *d_Sw11;
	cudaMalloc((void**)&d_Sw00, ns * sizeof(float));
	cudaMalloc((void**)&d_Sw01, ns * sizeof(float));
	cudaMalloc((void**)&d_Sw10, ns * sizeof(float));
	cudaMalloc((void**)&d_Sw11, ns * sizeof(float));
	cudaMemcpy(d_Sw00, cs->w00, ns * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Sw01, cs->w01, ns * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Sw10, cs->w10, ns * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Sw11, cs->w11, ns * sizeof(float), cudaMemcpyHostToDevice);
	
	// z and x coordinates of each source
	int *d_Sjz, *d_Sjx;
	cudaMalloc((void**)&d_Sjz, ns * sizeof(int));
	cudaMalloc((void**)&d_Sjx, ns * sizeof(int));
	cudaMemcpy(d_Sjz, cs->jz, ns * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Sjx, cs->jx, ns * sizeof(int), cudaMemcpyHostToDevice);
	
	
	/* calculate 2d linear interpolation coefficients for receiver locations */
    cr = lint2d_make(nr,rr,fdm);
	float *d_Rw00, *d_Rw01, *d_Rw10, *d_Rw11;
	cudaMalloc((void**)&d_Rw00, nr * sizeof(float));
	cudaMalloc((void**)&d_Rw01, nr * sizeof(float));
	cudaMalloc((void**)&d_Rw10, nr * sizeof(float));
	cudaMalloc((void**)&d_Rw11, nr * sizeof(float));
	cudaMemcpy(d_Rw00, cr->w00, nr * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Rw01, cr->w01, nr * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Rw10, cr->w10, nr * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Rw11, cr->w11, nr * sizeof(float), cudaMemcpyHostToDevice);
	
	// z and x coordinates of each receiver
	int *d_Rjz, *d_Rjx;
	cudaMalloc((void**)&d_Rjz, nr * sizeof(int));
	cudaMalloc((void**)&d_Rjx, nr * sizeof(int));
	cudaMemcpy(d_Rjz, cr->jz, nr * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Rjx, cr->jx, nr * sizeof(int), cudaMemcpyHostToDevice);

    /*------------------------------------------------------------*/
    /* setup FD coefficients */
    idz = 1/dz;
    idx = 1/dx;
    /*------------------------------------------------------------*/ 
	
	
	/*------------------------------------------------------------*/ 
	/* Read density and stiffness model data and transfer to GPU */
	
	float *tt1 = (float*)malloc(nz * nx * sizeof(float)); 
	h_ro=(float*)malloc(fdm->nzpad * fdm->nxpad * sizeof(float));
	h_c11=(float*)malloc(fdm->nzpad * fdm->nxpad * sizeof(float));
	h_c33=(float*)malloc(fdm->nzpad * fdm->nxpad * sizeof(float));
	h_c55=(float*)malloc(fdm->nzpad * fdm->nxpad * sizeof(float));
	h_c13=(float*)malloc(fdm->nzpad * fdm->nxpad * sizeof(float));

    /* input density */
    sf_floatread(tt1,nz*nx,Fden);     expand_cpu(tt1,h_ro , fdm->nb, nx, fdm->nxpad, nz, fdm->nzpad);

    /* input stiffness */
    sf_floatread(tt1,nz*nx,Fccc );    expand_cpu(tt1,h_c11, fdm->nb, nx, fdm->nxpad, nz, fdm->nzpad);
    sf_floatread(tt1,nz*nx,Fccc );    expand_cpu(tt1,h_c33, fdm->nb, nx, fdm->nxpad, nz, fdm->nzpad);
    sf_floatread(tt1,nz*nx,Fccc );    expand_cpu(tt1,h_c55, fdm->nb, nx, fdm->nxpad, nz, fdm->nzpad);
    sf_floatread(tt1,nz*nx,Fccc );    expand_cpu(tt1,h_c13, fdm->nb, nx, fdm->nxpad, nz, fdm->nzpad);
	
	cudaMalloc((void **)&d_ro, fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_c11, fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_c33, fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_c55, fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_c13, fdm->nzpad * fdm->nxpad * sizeof(float));
	
	cudaMemcpy(d_ro, h_ro, fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_c11, h_c11, fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_c33, h_c33, fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_c55, h_c55, fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_c13, h_c13, fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);


    /*------------------------------------------------------------*/
	/* boundary condition setup */
   	float *d_bzl_s, *d_bzh_s;
	float *d_bxl_s, *d_bxh_s;

	float *d_spo;
    if(dabc) {
		/* one-way abc setup   */
		vs = sf_floatalloc2(fdm->nzpad,fdm->nxpad); 
		for (ix=0; ix<fdm->nxpad; ix++) {
		    for(iz=0; iz<fdm->nzpad; iz++) {
				vs[ix][iz] = sqrt(h_c55[iz * fdm->nxpad + ix]/h_ro[iz * fdm->nxpad + ix] );
		    }
		}
		abcs = abcone2d_make(NOP,dt,vs,fsrf,fdm);
		free(*vs); free(vs);

		cudaMalloc((void**)&d_bzl_s, fdm->nxpad * sizeof(float));
		cudaMalloc((void**)&d_bzh_s, fdm->nxpad * sizeof(float));
		cudaMalloc((void**)&d_bxl_s, fdm->nzpad * sizeof(float));
		cudaMalloc((void**)&d_bxh_s, fdm->nzpad * sizeof(float));
		
		cudaMemcpy(d_bzl_s, abcs->bzl, fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_bzh_s, abcs->bzh, fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_bxl_s, abcs->bxl, fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_bxh_s, abcs->bxh, fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
		
		/* sponge abc setup */
		spo = sponge_make(fdm->nb);
		
		// d_spo contains all of the sponge coefficients
		cudaMalloc((void**)&d_spo, fdm->nb * sizeof(float));
		cudaMemcpy(d_spo, spo->w, fdm->nb * sizeof(float), cudaMemcpyHostToDevice);
    }


    /*------------------------------------------------------------*/
    /* allocate wavefield arrays */
	cudaMalloc((void **)&d_umz, fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_uoz, fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_upz, fdm->nzpad * fdm->nxpad * sizeof(float));	
	cudaMalloc((void **)&d_uaz, fdm->nzpad * fdm->nxpad * sizeof(float));

	cudaMalloc((void **)&d_umx, fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_uox, fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_upx, fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_uax, fdm->nzpad * fdm->nxpad * sizeof(float));
	
	cudaMalloc((void **)&d_tzz, fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_tzx, fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_txx, fdm->nzpad * fdm->nxpad * sizeof(float));
	
	sf_check_gpu_error("allocate grid arrays");
	
	cudaMemset(d_umz, 0, fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_uoz, 0, fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_upz, 0, fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_uaz, 0, fdm->nzpad * fdm->nxpad * sizeof(float));
	
	cudaMemset(d_umx, 0, fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_uox, 0, fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_upx, 0, fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_uax, 0, fdm->nzpad * fdm->nxpad * sizeof(float));

	cudaMemset(d_tzz, 0, fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_tzx, 0, fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_txx, 0, fdm->nzpad * fdm->nxpad * sizeof(float));

	sf_check_gpu_error("initialize grid arrays");
	

	/*------------------------------------------------------------*/
    /* precompute 1/ro * dt^2									  */
	/*------------------------------------------------------------*/
	dim3 dimGrid5(ceil(fdm->nxpad/16.0f),ceil(fdm->nzpad/16.0f));
	dim3 dimBlock5(16,16);
	computeRo<<<dimGrid5, dimBlock5>>>(d_ro, dt, fdm->nxpad, fdm->nzpad, NOP);
	sf_check_gpu_error("computeRo Kernel");
	

    /*------------------------------------------------------------*/
    /* 
     *  MAIN LOOP
     */
    /*------------------------------------------------------------*/
    if(verb) fprintf(stderr,"\n");
    for (it=0; it<nt; it++) {
		if(verb) fprintf(stderr,"\b\b\b\b\b\b%d",it);
		/*------------------------------------------------------------*/
		/* from displacement to strain AND strain to stress           */
		/*		- Compute strains from displacements as in equation 1 */
		/*			- Step #1	(Steps denoted are as in Figure 2)	  */
		/*		- Compute stress from strain as in equation 2		  */
		/*			- Step #2										  */
		/*------------------------------------------------------------*/
			dim3 dimGrid9(ceil(fdm->nxpad/16.0f), ceil(fdm->nzpad/16.0f));
			dim3 dimBlock9(16,16);
			dispToStrain_strainToStress<<<dimGrid9, dimBlock9>>>(d_txx, d_tzz, d_tzx, d_uox, d_uoz, d_c11, d_c33, d_c55, d_c13, idx, idz, fdm->nxpad, fdm->nzpad, NOP);
			sf_check_gpu_error("dispToStrainToStress Kernel");


		/*------------------------------------------------------------*/
		/* free surface boundary condition							  */
		/*		- sets the z-component of stress tensor along the	  */
		/*			free surface boundary to 0						  */
		/*			- Step #3   									  */
		/*------------------------------------------------------------*/
			if(fsrf) {
				dim3 dimGrid3(ceil(fdm->nxpad/16.0f),ceil(fdm->nb/16.0f));
				dim3 dimBlock3(16,16);
				freeSurf<<<dimGrid3,dimBlock3>>>(d_tzz, d_tzx, fdm->nxpad, fdm->nb);
				sf_check_gpu_error("freeSurf Kernel");
			}


		/*------------------------------------------------------------*/
		/* inject stress source                                       */
		/*		- Step #4   										  */
		/*------------------------------------------------------------*/
			if(ssou) {
				dim3 dimGrid7(ns, 1, 1);
				dim3 dimBlock7(2 * nbell + 1, 2 * nbell + 1, 1);
				lint2d_bell_gpu<<<dimGrid7, dimBlock7>>>(d_tzz, d_ww, d_Sw00, d_Sw01, d_Sw10, d_Sw11, d_bell, d_Sjz, d_Sjx, it, nc, ns, 0, nbell, fdm->nxpad);
				lint2d_bell_gpu<<<dimGrid7, dimBlock7>>>(d_txx, d_ww, d_Sw00, d_Sw01, d_Sw10, d_Sw11, d_bell, d_Sjz, d_Sjx, it, nc, ns, 1, nbell, fdm->nxpad);			
				sf_check_gpu_error("lint2d_bell_gpu Kernel");
			}


		/*------------------------------------------------------------*/
		/* from stress to acceleration (first term in RHS of eq. 3)	  */
		/*		- Step #5											  */
		/*------------------------------------------------------------*/
			dim3 dimGrid4(ceil((fdm->nxpad-(2*NOP))/16.0f),ceil((fdm->nzpad-(2*NOP))/16.0f));
			dim3 dimBlock4(16,16);
			stressToAcceleration<<<dimGrid4, dimBlock4>>>(d_uax, d_uaz, d_txx, d_tzz, d_tzx, idx, idz, fdm->nxpad, fdm->nzpad);
			sf_check_gpu_error("stressToAcceleration Kernel");
		

		/*------------------------------------------------------------*/
		/* inject acceleration source  (second term in RHS of eq. 3)  */
		/*		- Step #6											  */
		/*------------------------------------------------------------*/
			if(!ssou) {
			    dim3 dimGrid8(ns, 1, 1);
				dim3 dimBlock8(2 * nbell + 1, 2 * nbell + 1, 1);
				lint2d_bell_gpu<<<dimGrid8, dimBlock8>>>(d_uaz, d_ww, d_Sw00, d_Sw01, d_Sw10, d_Sw11, d_bell, d_Sjz, d_Sjx, it, nc, ns, 0, nbell, fdm->nxpad);
				lint2d_bell_gpu<<<dimGrid8, dimBlock8>>>(d_uax, d_ww, d_Sw00, d_Sw01, d_Sw10, d_Sw11, d_bell, d_Sjz, d_Sjx, it, nc, ns, 1, nbell, fdm->nxpad);
				sf_check_gpu_error("lint2d_bell_gpu Kernel");
			}


		/*------------------------------------------------------------*/
		/* step forward in time                                       */
		/*		- Compute forward time step based on acceleration	  */
		/*			- Step #7
		/*------------------------------------------------------------*/
			dim3 dimGrid6(ceil(fdm->nxpad/16.0f),ceil(fdm->nzpad/12.0f));
			dim3 dimBlock6(16,12);
			stepTime<<<dimGrid6, dimBlock6>>>(d_upz, d_uoz, d_umz, d_uaz, d_upx, d_uox, d_umx, d_uax, d_ro, fdm->nxpad, fdm->nzpad);
			sf_check_gpu_error("stepTime Kernel");
		
		
		/* circulate wavefield arrays */
		d_utz=d_umz; d_utx=d_umx;
		d_umz=d_uoz; d_umx=d_uox;
		d_uoz=d_upz; d_uox=d_upx;
		d_upz=d_utz; d_upx=d_utx;
	
		/*------------------------------------------------------------*/
		/* apply boundary conditions                                  */
		/*		- Step #8
		/*------------------------------------------------------------*/
			if(dabc) {
				
				/*---------------------------------------------------------------*/
				/* apply One-way Absorbing BC as in (Clayton and Enquist, 1977)  */
				/*---------------------------------------------------------------*/
				/* One-way Absorbing BC */
				dim3 dimGrid_TB(ceil(fdm->nxpad/192.0f), 2, 1);
				dim3 dimBlock_TB(MIN(192, fdm->nxpad), 1, 1);

				dim3 dimGrid_LR(2, ceil(fdm->nzpad/192.0f), 1);
				dim3 dimBlock_LR(1, MIN(192, fdm->nzpad), 1);

				abcone2d_apply_TB_gpu<<<dimGrid_TB, dimBlock_TB, 0, stream[0]>>>(d_uoz, d_umz, d_bzl_s, d_bzh_s, fdm->nxpad, fdm->nzpad, fsrf);
				abcone2d_apply_LR_gpu<<<dimGrid_LR, dimBlock_LR, 0, stream[0]>>>(d_uoz, d_umz, d_bxl_s, d_bxh_s, fdm->nxpad, fdm->nzpad);

				abcone2d_apply_TB_gpu<<<dimGrid_TB, dimBlock_TB, 0, stream[1]>>>(d_uox, d_umx, d_bzl_s, d_bzh_s, fdm->nxpad, fdm->nzpad, fsrf);
				abcone2d_apply_LR_gpu<<<dimGrid_LR, dimBlock_LR, 0, stream[1]>>>(d_uox, d_umx, d_bxl_s, d_bxh_s, fdm->nxpad, fdm->nzpad);
	

				/*---------------------------------------------------------------*/
				/* apply Sponge BC as in (Cerjan, et al., 1985)                  */
				/*---------------------------------------------------------------*/
				dim3 dimGrid_TB2(ceil(fdm->nxpad/256.0f), (fdm->nb * 1), 1);
				dim3 dimBlock_TB2(256,1,1);

				dim3 dimGrid_LR2(ceil(fdm->nb/256.0f), fdm->nzpad, 1);
				dim3 dimBlock_LR2(MIN(256, fdm->nb),1,1);

				sponge2d_apply_LR_gpu<<<dimGrid_LR2, dimBlock_LR2, 0, stream[2]>>>(d_upz, d_spo, fdm->nxpad, fdm->nb, nx);
				sponge2d_apply_TB_gpu<<<dimGrid_TB2, dimBlock_TB2, 0, stream[2]>>>(d_upz, d_spo, fdm->nxpad, fdm->nb, nz);

				sponge2d_apply_LR_gpu<<<dimGrid_LR2, dimBlock_LR2, 0, stream[3]>>>(d_upx, d_spo, fdm->nxpad, fdm->nb, nx);
				sponge2d_apply_TB_gpu<<<dimGrid_TB2, dimBlock_TB2, 0, stream[3]>>>(d_upx, d_spo, fdm->nxpad, fdm->nb, nz);
	
				cudaStreamSynchronize(stream[0]);
				sponge2d_apply_LR_gpu<<<dimGrid_LR2, dimBlock_LR2, 0, stream[0]>>>(d_umz, d_spo, fdm->nxpad, fdm->nb, nx);			
				sponge2d_apply_TB_gpu<<<dimGrid_TB2, dimBlock_TB2, 0, stream[0]>>>(d_umz, d_spo, fdm->nxpad, fdm->nb, nz);

				sponge2d_apply_LR_gpu<<<dimGrid_LR2, dimBlock_LR2, 0, stream[4]>>>(d_uoz, d_spo, fdm->nxpad, fdm->nb, nx);
				sponge2d_apply_TB_gpu<<<dimGrid_TB2, dimBlock_TB2, 0, stream[4]>>>(d_uoz, d_spo, fdm->nxpad, fdm->nb, nz);
	
				cudaStreamSynchronize(stream[1]);
				sponge2d_apply_LR_gpu<<<dimGrid_LR2, dimBlock_LR2, 0, stream[1]>>>(d_umx, d_spo, fdm->nxpad, fdm->nb, nx);			
				sponge2d_apply_TB_gpu<<<dimGrid_TB2, dimBlock_TB2, 0, stream[1]>>>(d_umx, d_spo, fdm->nxpad, fdm->nb, nz);
	
				sponge2d_apply_LR_gpu<<<dimGrid_LR2, dimBlock_LR2, 0, stream[5]>>>(d_uox, d_spo, fdm->nxpad, fdm->nb, nx);			
				sponge2d_apply_TB_gpu<<<dimGrid_TB2, dimBlock_TB2, 0, stream[5]>>>(d_uox, d_spo, fdm->nxpad, fdm->nb, nz);
		
				cudaDeviceSynchronize();

				sf_check_gpu_error("boundary condition kernels");
			}	    

		/*------------------------------------------------------------*/
		/* cut wavefield and save 									  */
		/*		- Step #9
		/*------------------------------------------------------------*/
		    if(snap && it%jsnap==0) {
				cudaMemcpy(h_uox, d_uox, fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDeviceToHost);
				cudaMemcpy(h_uoz, d_uoz, fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDeviceToHost);

				for (int x = 0; x < fdm->nxpad; x++){
					for (int z = 0; z < fdm->nzpad; z++){
						uox[x][z] = h_uox[z * fdm->nxpad + x];
						uoz[x][z] = h_uoz[z * fdm->nxpad + x];
					}
				}

				cut2d(uoz,uc,fdm,acz,acx);
				sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);

				cut2d(uox,uc,fdm,acz,acx);
				sf_floatwrite(uc[0],sf_n(acz)*sf_n(acx),Fwfl);
		    }

		/*------------------------------------------------------------*/
		/* extract receiver data									  */
		/*------------------------------------------------------------*/
		    if(it%jdata==0) {
				dim3 dimGrid_extract(MIN(nr,ceil(nr/1024.0f)), 1, 1);
				dim3 dimBlock_extract(MIN(nr, 1024), 1, 1);
				lint2d_extract_gpu<<<dimGrid_extract, dimBlock_extract>>>(d_dd, nr, fdm->nxpad, d_uoz, d_uox, d_Rjz, d_Rjx, d_Rw00, d_Rw01, d_Rw10, d_Rw11);
				sf_check_gpu_error("lint2d_extract kernel");

				cudaMemcpy(h_dd, d_dd, nr * nc * sizeof(float), cudaMemcpyDeviceToHost);
				sf_floatwrite(h_dd, nr*nc, Fdat);
				
			}
	}
	
    if(verb) fprintf(stderr,"\n");
    
    /*------------------------------------------------------------*/
    /* deallocate host arrays */

    free(**ww); free(*ww); free(ww);
    free(ss);
    free(rr);
    free(*dd);  free(dd);
	
	if (snap){
		free(*uoz); free(uoz);
	    free(*uox); free(uox);
	    free(h_uoz);  
		free(h_uox);
		free(*uc);  free(uc);    
	}
	
	free(h_c11); free(h_c33); free(h_c55); free(h_c13); free(h_ro);
	free(h_bell);
	free(h_dd);
	

	/*------------------------------------------------------------*/
    /* deallocate GPU arrays */

	cudaFree(d_ro);
	cudaFree(d_c11); 	cudaFree(d_c33); 	cudaFree(d_c55); 	cudaFree(d_c13);
	cudaFree(d_umz); 	cudaFree(d_uoz); 	cudaFree(d_upz); 	cudaFree(d_uaz); 	cudaFree(d_utz);
	cudaFree(d_umx); 	cudaFree(d_uox); 	cudaFree(d_upx); 	cudaFree(d_uax); 	cudaFree(d_utx);
	cudaFree(d_tzz); 	cudaFree(d_tzx); 	cudaFree(d_txx);
	cudaFree(d_Sw00); 	cudaFree(d_Sw01); 	cudaFree(d_Sw10); 	cudaFree(d_Sw11);
	cudaFree(d_Sjz); 	cudaFree(d_Sjx);
	cudaFree(d_Rjz);	cudaFree(d_Rjx);
	cudaFree(d_bell);
	cudaFree(d_ww);
	cudaFree(d_dd);

	if (dabc){
		cudaFree(d_bzl_s); 	cudaFree(d_bzh_s);
		cudaFree(d_bxl_s); 	cudaFree(d_bxh_s);
		cudaFree(d_spo);
	}

    sf_close();
	exit(0);
}

