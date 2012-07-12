/* 3D elastic time-domain FD modeling with multiple GPUs coordinated via MPI*/

/*
  Authors: Robin M. Weiss and Jeffrey Shragge

  Use of this code is freely avaiable. In publications, please reference the paper: 
  Weiss and Shragge, "Solving 3D Anisotropic Elastic Wave Equations on Parallel 
  GPU Devices", GEOPHYSICS. http://software.seg.org/2012/0063

  This code is a GPU-enabled version of the ewefd3d module from the Madagascar
  software package (see: http://www.reproducibility.org).  It implements a 3D
  Finite-Difference Time Domain solver for the elastice wave equation with 
  2nd- and 8th- order temporal and spatial accuracy, respectively.  Computation
  is distributed across an arbitrary number of GPU devices and coordinted by MPI.  
  For more information, see (Weiss and Shragge, "Solving 3D Anisotropic Elastic Wave 
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
#include <stdio.h>
#include <mpi.h>
#include <cuda.h>
#include <cuda_runtime_api.h>

extern "C" {
#include <rsf.h>
}

#include "fdutil.c"
#include "ewefd3d_kernels.cu"

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define NOP 4 /* derivative operator half-size */


// checks the current GPU device for an error flag and prints to stderr
static void sf_check_gpu_error (int gpu, const char *msg) {
    cudaError_t err = cudaGetLastError ();
     if (cudaSuccess != err)
        sf_error ("Cuda error on GPU %d: %s: %s", gpu, msg, cudaGetErrorString (err));
}

// entry point
int main (int argc, char* argv[]) {
	
	// Initialize MPI
	int rank, size;
	MPI_Init (&argc, &argv);				/* start MPI */
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &size);	/* get number of processes */
	MPI_Status status;
	MPI_Request request;
	
	// Initialize RSF
	sf_init(argc, argv);
	
	bool verb,fsrf,snap,ssou,dabc,interp;
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
    sf_axis at,ax,ay,az;
    sf_axis as,ar,ac;

    int     nt,nz,nx,ny,ns,nr,nc,nb;
    int     it,iz,ix,iy;
    float   dt,dz,dx,dy,idz,idx,idy;

    /* FDM structure */
    fdm3d    fdm=NULL;

    /* I/O arrays */
    float***ww=NULL;           /* wavelet   */
    pt3d   *ss=NULL;           /* sources   */
    pt3d   *rr=NULL;           /* receivers */

    /*------------------------------------------------------------*/
    /* orthorombic footprint - 9 coefficients */
    /* c11 c12 c13 
       .   c22 c23 
       .   .   c33 
                  c44
                     c55
                        c66 */
	float *d_c11, *d_c22, *d_c33, *d_c44, *d_c55, *d_c66, *d_c12, *d_c13, *d_c23;
	float *h_c11, *h_c22, *h_c33, *h_c44, *h_c55, *h_c66, *h_c12, *h_c13, *h_c23;
	
	// density
	float *h_ro, *d_ro;	
	
    /*------------------------------------------------------------*/
    /* displacement: um = U @ t-1; uo = U @ t; up = U @ t+1 */
	float *d_umz,*d_uoz,*d_upz,*d_uaz,*d_utz; 
    float *d_umx,*d_uox,*d_upx,*d_uax,*d_utx;
    float *d_umy,*d_uoy,*d_upy,*d_uay,*d_uty;

	// used for writing wavefield to file, only needed if snap=y
	float ***uox, ***uoy, ***uoz;
	float *h_uox, *h_uoy, *h_uoz;

    /* stress/strain tensor */ 
	float *d_tzz,*d_txx,*d_tyy,*d_txy,*d_tyz,*d_tzx;

    /*------------------------------------------------------------*/
    /* linear interpolation weights/indices */
    lint3d cs,cr;

    /* Gaussian bell */
    int nbell;
    
    /* wavefield cut params */
    sf_axis   acz=NULL,acx=NULL,acy=NULL;
    int       nqz,nqx,nqy;
    float     oqz,oqx,oqy;
    float     dqz,dqx,dqy;
    float     ***uc=NULL;

	/*------------------------------------------------------------*/
    /* execution flags */
    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("snap",&snap)) snap=false; /* wavefield snapshots flag */
    if(! sf_getbool("free",&fsrf)) fsrf=false; /* free surface flag */
    if(! sf_getbool("ssou",&ssou)) ssou=false; /* stress source */
    if(! sf_getbool("dabc",&dabc)) dabc=false; /* absorbing BC */
	if(! sf_getbool("interp",&interp)) interp=true; /* perform linear interpolation on receiver locations */
    /*------------------------------------------------------------*/


    /*------------------------------------------------------------*/
    /* I/O files */
    Fwav = sf_input ("wav"); /* wavelet   */
    Fccc = sf_input ("ccc"); /* stiffness */
    Fden = sf_input ("den"); /* density   */
    Fsou = sf_input ("sou"); /* sources   */
    Frec = sf_input ("rec"); /* receivers */
    Fwfl = sf_output("wfl"); /* wavefield */
    Fdat = sf_output("out"); /* data      */
	/*------------------------------------------------------------*/


	/*------------------------------------------------------------*/
    /* initialize GPU */
	int ngpu;
	if(!sf_getint("ngpu",&ngpu)) ngpu=1; 	/* Number of GPUs in each node, must be set to lowest common number of GPUs*/
	
	cudaSetDevice(rank % ngpu);	// attach this rank to a GPU
	cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
	sf_warning("rank %d using GPU %d\n", rank, rank % ngpu);
	
	/*------------------------------------------------------------*/
	
	
	/*------------------------------------------------------------*/
    /* axes */
    at = sf_iaxa(Fwav,3); sf_setlabel(at,"t"); if(verb && rank == 0) sf_raxa(at); /* time */
    az = sf_iaxa(Fccc,1); sf_setlabel(az,"z"); if(verb && rank == 0) sf_raxa(az); /* depth */
    ax = sf_iaxa(Fccc,2); sf_setlabel(ax,"x"); if(verb && rank == 0) sf_raxa(ax); /* space x */
    ay = sf_iaxa(Fccc,3); sf_setlabel(ay,"y"); if(verb && rank == 0) sf_raxa(ay); /* space y */

    as = sf_iaxa(Fsou,2); sf_setlabel(as,"s"); if(verb && rank == 0) sf_raxa(as); /* sources */
    ar = sf_iaxa(Frec,2); sf_setlabel(ar,"r"); if(verb && rank == 0) sf_raxa(ar); /* receivers */

    nt = sf_n(at); dt = sf_d(at);
    nz = sf_n(az); dz = sf_d(az);
    nx = sf_n(ax); dx = sf_d(ax);
    ny = sf_n(ay); dy = sf_d(ay);

    ns = sf_n(as);
    nr = sf_n(ar);
    /*------------------------------------------------------------*/


	/*------------------------------------------------------------*/
    /* other execution parameters */
    if(! sf_getint("nbell",&nbell)) nbell=5;  /* bell size */
    if(verb && rank == 0) sf_warning("nbell=%d",nbell);
    if(! sf_getint("jdata",&jdata)) jdata=1;	/* extract receiver data every jdata time steps */
    if(snap) {  
	if(! sf_getint("jsnap",&jsnap)) jsnap=nt; /* save wavefield every jsnap time steps */
    }
    /*------------------------------------------------------------*/


    /*------------------------------------------------------------*/
    /* expand domain for FD operators and ABC */
    if( !sf_getint("nb",&nb) || nb<NOP) nb=NOP;

	fdm=fdutil3d_init(verb,fsrf,az,ax,ay,nb,1);
	
	sf_setn(az,fdm->nzpad); sf_seto(az,fdm->ozpad); if(verb && rank == 0) sf_raxa(az);
    sf_setn(ax,fdm->nxpad); sf_seto(ax,fdm->oxpad); if(verb && rank == 0) sf_raxa(ax);
    sf_setn(ay,fdm->nypad); sf_seto(ay,fdm->oypad); if(verb && rank == 0) sf_raxa(ay);
	/*------------------------------------------------------------*/
	
	
	/*------------------------------------------------------------*/
    /* compute sub-domain dimmensions (domain decomposition) */
	
	int nyinterior = (fdm->nypad / size); // size of sub-domains in y-dimension EXCLUDING any ghost cells from adjacent GPUs
	
	int nylocal = fdm->nypad / size; // size of sub-domains in y-dimension INCLUDING any ghost cells from adjacent GPUs
	
	if (rank == 0 || rank == size-1){	
		nylocal += 4;	// exterior nodes require 4 additional ghost slices from neighbor
	}
	else {								
		nylocal += 8;	// interior nodes require 8 additional ghost slices from neighbors
	}
	
	if (size == 1){
		nylocal = fdm->nypad;	// Running in 1-GPU mode, this node will contain all of the grid in the y-direction
	}
	
	
	// check that all dimmeionsons are ok for FD kernels
	if ((fdm->nzpad - 8) % 24 != 0){
		sf_error("nz + 2*nb - 8 is not a multiple of 24");
	}
	if ((fdm->nxpad - 8) % 24 != 0){
		sf_error("nx + 2*nb - 8 is not a multiple of 24");
	}
	if ((fdm->nypad % size) != 0){
		sf_error("You are using %d GPUs.\n(ny + 2*nb) must me a multiple of %d\nChange model dimensions or select a different number of GPUs", size, size);
	}
	/*------------------------------------------------------------*/
	

	/*------------------------------------------------------------*/
    /* setup bell for source injection smoothing */
	if (nbell * 2 + 1 > 32){
		sf_error("nbell must be <= 15\n");
	}
	
	float *h_bell;
	h_bell = (float*)malloc((2*nbell+1)*(2*nbell+1)*(2*nbell+1)*sizeof(float));
	
	float *d_bell;
	cudaMalloc((void**)&d_bell, (2*nbell+1)*(2*nbell+1)*(2*nbell+1)*sizeof(float));
	
	float s = 0.5*nbell;
    for (iy=-nbell;iy<=nbell;iy++) {
		for (ix=-nbell;ix<=nbell;ix++) {
	    	for(iz=-nbell;iz<=nbell;iz++) {
				h_bell[(iy + nbell) * (2*nbell+1) * (2*nbell+1) + (iz + nbell) * (2*nbell+1) + (ix + nbell)] = exp(-(iz*iz+ix*ix+iy*iy)/s);
	    	}
		}    
    }

	cudaMemcpy(d_bell, h_bell, (2*nbell+1)*(2*nbell+1)*(2*nbell+1)*sizeof(float), cudaMemcpyHostToDevice);
	/*------------------------------------------------------------*/
	
	
	/*------------------------------------------------------------*/
	/* 3D vector components */
    nc=3;
	ac=sf_maxa(nc  ,0,1);
	/*------------------------------------------------------------*/

	
	/*------------------------------------------------------------*/
    /* setup output data files and arrays */
	if (snap){
		h_uoz = (float*)malloc(nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
		h_uox = (float*)malloc(nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
		h_uoy = (float*)malloc(nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
	}
	
	if (rank == 0){
	    sf_oaxa(Fdat,ar,1);
	    sf_oaxa(Fdat,ac,2);

	    sf_setn(at,nt/jdata);
	    sf_setd(at,dt*jdata);
	    sf_oaxa(Fdat,at,3);

	    /* setup output wavefield header */
	    if(snap) {
			
			// Used to accumulate wavefield data from other GPUs
			uoz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
			uox=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
			uoy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
			
			nqz=sf_n(az);
			nqx=sf_n(ax);
			nqy=sf_n(ay);
            
			oqz=sf_o(az);
			oqx=sf_o(ax);
			oqy=sf_o(ay);

			dqz=sf_d(az);
			dqx=sf_d(ax);
			dqy=sf_d(ay);

			acz = sf_maxa(nqz,oqz,dqz); sf_raxa(acz);
			acx = sf_maxa(nqx,oqx,dqx); sf_raxa(acx);
			acy = sf_maxa(nqy,oqy,dqy); sf_raxa(acy);

			uc=sf_floatalloc3(sf_n(acz),sf_n(acx),sf_n(acy));

			ntsnap=0;
			for(it=0; it<nt; it++) {
			    if(it%jsnap==0) ntsnap++;
			}
			sf_setn(at,  ntsnap);
			sf_setd(at,dt*jsnap);
			if(verb && rank == 0) sf_raxa(at);

			sf_oaxa(Fwfl,acz,1);
			sf_oaxa(Fwfl,acx,2);
			sf_oaxa(Fwfl,acy,3);
			sf_oaxa(Fwfl,ac, 4);
			sf_oaxa(Fwfl,at, 5);
	    }
	}
	/*------------------------------------------------------------*/
	
	
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
	float *d_dd;
	float *h_dd;
	float *h_dd_combined;
	h_dd = (float*)malloc(nr * nc * sizeof(float));
	cudaMalloc((void**)&d_dd, nr*nc*sizeof(float));
	
	if (rank == 0){		// rank 0 accumulates data from other ranks
		h_dd_combined = (float*)malloc(nr * nc * sizeof(float));
	}
	/*------------------------------------------------------------*/
	
	
	/*------------------------------------------------------------*/
    /* get source/receiver coordinates */
    ss = (pt3d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt3d*) sf_alloc(nr,sizeof(*rr)); 

    pt3dread1(Fsou,ss,ns,3); /* read (x,y,z) coordinates */
    pt3dread1(Frec,rr,nr,3); /* read (x,y,z) coordinates */
    
	/* calculate 3d linear interpolation coefficients for source locations */
    cs = lint3d_make(ns,ss,fdm);
	float *d_Sw000, *d_Sw001, *d_Sw010, *d_Sw011, *d_Sw100, *d_Sw101, *d_Sw110, *d_Sw111;
	cudaMalloc((void**)&d_Sw000, ns * sizeof(float));
	cudaMalloc((void**)&d_Sw001, ns * sizeof(float));
	cudaMalloc((void**)&d_Sw010, ns * sizeof(float));
	cudaMalloc((void**)&d_Sw011, ns * sizeof(float));
	cudaMalloc((void**)&d_Sw100, ns * sizeof(float));
	cudaMalloc((void**)&d_Sw101, ns * sizeof(float));
	cudaMalloc((void**)&d_Sw110, ns * sizeof(float));
	cudaMalloc((void**)&d_Sw111, ns * sizeof(float));
	cudaMemcpy(d_Sw000, cs->w000, ns * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Sw001, cs->w001, ns * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Sw010, cs->w010, ns * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Sw011, cs->w011, ns * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Sw100, cs->w100, ns * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Sw101, cs->w101, ns * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Sw110, cs->w110, ns * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Sw111, cs->w111, ns * sizeof(float), cudaMemcpyHostToDevice);
	
	// z, x, and y coordinates of each source
	int *d_Sjz, *d_Sjx, *d_Sjy;	
	cudaMalloc((void**)&d_Sjz, ns * sizeof(int));
	cudaMalloc((void**)&d_Sjx, ns * sizeof(int));
	cudaMalloc((void**)&d_Sjy, ns * sizeof(int));
	cudaMemcpy(d_Sjz, cs->jz, ns * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Sjx, cs->jx, ns * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Sjy, cs->jy, ns * sizeof(int), cudaMemcpyHostToDevice);


	/* calculate 3d linear interpolation coefficients for receiver locations */
    cr = lint3d_make(nr,rr,fdm);
	float *d_Rw000, *d_Rw001, *d_Rw010, *d_Rw011, *d_Rw100, *d_Rw101, *d_Rw110, *d_Rw111;
	if (interp){
		cudaMalloc((void**)&d_Rw000, nr * sizeof(float));
		cudaMalloc((void**)&d_Rw001, nr * sizeof(float));
		cudaMalloc((void**)&d_Rw010, nr * sizeof(float));
		cudaMalloc((void**)&d_Rw011, nr * sizeof(float));
		cudaMalloc((void**)&d_Rw100, nr * sizeof(float));
		cudaMalloc((void**)&d_Rw101, nr * sizeof(float));
		cudaMalloc((void**)&d_Rw110, nr * sizeof(float));
		cudaMalloc((void**)&d_Rw111, nr * sizeof(float));
		cudaMemcpy(d_Rw000, cr->w000, nr * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_Rw001, cr->w001, nr * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_Rw010, cr->w010, nr * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_Rw011, cr->w011, nr * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_Rw100, cr->w100, nr * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_Rw101, cr->w101, nr * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_Rw110, cr->w110, nr * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_Rw111, cr->w111, nr * sizeof(float), cudaMemcpyHostToDevice);
	}
	
	// z, x, and y coordinates of each receiver
	int *d_Rjz, *d_Rjx, *d_Rjy;
	cudaMalloc((void**)&d_Rjz, nr * sizeof(int));
	cudaMalloc((void**)&d_Rjx, nr * sizeof(int));
	cudaMalloc((void**)&d_Rjy, nr * sizeof(int));
	cudaMemcpy(d_Rjz, cr->jz, nr * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Rjx, cr->jx, nr * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Rjy, cr->jy, nr * sizeof(int), cudaMemcpyHostToDevice);
	/*------------------------------------------------------------*/
    

    /*------------------------------------------------------------*/
    /* setup FD coefficients */
    idz = 1/dz;
    idx = 1/dx;
    idy = 1/dy;
	/*------------------------------------------------------------*/
    
	
	/*------------------------------------------------------------*/ 
	/* Read density and stiffness model data */

	// allocate space on GPU for this GPU's portion of the domain
	cudaMalloc((void **)&d_c11, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_c22, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_c33, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_c44, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_c55, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_c66, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_c12, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_c13, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_c23, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_ro, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
	
	if (rank == 0){	// read in density and stiffness arrays and expand to padded dimensions
		float *tt1 = (float*)malloc(nz*nx*ny*sizeof(float));
		
	    /* density */
	    h_ro = (float*)malloc(fdm->nzpad * fdm->nxpad * fdm->nypad * sizeof(float));
	    sf_floatread(tt1,nz*nx*ny,Fden);     expand_cpu(tt1, h_ro, fdm->nb, nx, fdm->nxpad, ny, fdm->nypad, nz, fdm->nzpad);
	
	    /* stiffness */
	    h_c11 = (float*)malloc(fdm->nzpad * fdm->nxpad * fdm->nypad * sizeof(float));
	    h_c22 = (float*)malloc(fdm->nzpad * fdm->nxpad * fdm->nypad * sizeof(float)); 
	    h_c33 = (float*)malloc(fdm->nzpad * fdm->nxpad * fdm->nypad * sizeof(float)); 
	    h_c44 = (float*)malloc(fdm->nzpad * fdm->nxpad * fdm->nypad * sizeof(float)); 
	    h_c55 = (float*)malloc(fdm->nzpad * fdm->nxpad * fdm->nypad * sizeof(float)); 
	    h_c66 = (float*)malloc(fdm->nzpad * fdm->nxpad * fdm->nypad * sizeof(float)); 
	    h_c12 = (float*)malloc(fdm->nzpad * fdm->nxpad * fdm->nypad * sizeof(float)); 
	    h_c13 = (float*)malloc(fdm->nzpad * fdm->nxpad * fdm->nypad * sizeof(float)); 
	    h_c23 = (float*)malloc(fdm->nzpad * fdm->nxpad * fdm->nypad * sizeof(float));
	    sf_floatread(tt1,nz*nx*ny,Fccc);    expand_cpu(tt1,h_c11,fdm->nb, nx, fdm->nxpad, ny, fdm->nypad, nz, fdm->nzpad);
	    sf_floatread(tt1,nz*nx*ny,Fccc);    expand_cpu(tt1,h_c22,fdm->nb, nx, fdm->nxpad, ny, fdm->nypad, nz, fdm->nzpad);    
	    sf_floatread(tt1,nz*nx*ny,Fccc);    expand_cpu(tt1,h_c33,fdm->nb, nx, fdm->nxpad, ny, fdm->nypad, nz, fdm->nzpad);    
	    sf_floatread(tt1,nz*nx*ny,Fccc);    expand_cpu(tt1,h_c44,fdm->nb, nx, fdm->nxpad, ny, fdm->nypad, nz, fdm->nzpad);
	    sf_floatread(tt1,nz*nx*ny,Fccc);    expand_cpu(tt1,h_c55,fdm->nb, nx, fdm->nxpad, ny, fdm->nypad, nz, fdm->nzpad);    
	    sf_floatread(tt1,nz*nx*ny,Fccc);    expand_cpu(tt1,h_c66,fdm->nb, nx, fdm->nxpad, ny, fdm->nypad, nz, fdm->nzpad);
	    sf_floatread(tt1,nz*nx*ny,Fccc);    expand_cpu(tt1,h_c12,fdm->nb, nx, fdm->nxpad, ny, fdm->nypad, nz, fdm->nzpad);
	    sf_floatread(tt1,nz*nx*ny,Fccc);    expand_cpu(tt1,h_c13,fdm->nb, nx, fdm->nxpad, ny, fdm->nypad, nz, fdm->nzpad);
	    sf_floatread(tt1,nz*nx*ny,Fccc);    expand_cpu(tt1,h_c23,fdm->nb, nx, fdm->nxpad, ny, fdm->nypad, nz, fdm->nzpad);
	    free(tt1);
	}
	
		
	if (rank == 0){
		// put GPU 0's sub-domain into GPU memory
		cudaMemcpy(d_ro, h_ro, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_c11, h_c11, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_c22, h_c22, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_c33, h_c33, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_c44, h_c44, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_c55, h_c55, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_c66, h_c66, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_c12, h_c12, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_c13, h_c13, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_c23, h_c23, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
				
		// send remaining sub-domains to other GPUs
		for (int dst = 1; dst < size; dst++){
			MPI_Isend(h_ro  + (dst * nyinterior * fdm->nzpad * fdm->nxpad), (nyinterior * fdm->nzpad * fdm->nxpad), MPI_FLOAT, dst, 0, MPI_COMM_WORLD, &request);
			MPI_Isend(h_c11 + (dst * nyinterior * fdm->nzpad * fdm->nxpad), (nyinterior * fdm->nzpad * fdm->nxpad), MPI_FLOAT, dst, 1, MPI_COMM_WORLD, &request);
			MPI_Isend(h_c22 + (dst * nyinterior * fdm->nzpad * fdm->nxpad), (nyinterior * fdm->nzpad * fdm->nxpad), MPI_FLOAT, dst, 2, MPI_COMM_WORLD, &request);
			MPI_Isend(h_c33 + (dst * nyinterior * fdm->nzpad * fdm->nxpad), (nyinterior * fdm->nzpad * fdm->nxpad), MPI_FLOAT, dst, 3, MPI_COMM_WORLD, &request);
			MPI_Isend(h_c44 + (dst * nyinterior * fdm->nzpad * fdm->nxpad), (nyinterior * fdm->nzpad * fdm->nxpad), MPI_FLOAT, dst, 4, MPI_COMM_WORLD, &request);
			MPI_Isend(h_c55 + (dst * nyinterior * fdm->nzpad * fdm->nxpad), (nyinterior * fdm->nzpad * fdm->nxpad), MPI_FLOAT, dst, 5, MPI_COMM_WORLD, &request);
			MPI_Isend(h_c66 + (dst * nyinterior * fdm->nzpad * fdm->nxpad), (nyinterior * fdm->nzpad * fdm->nxpad), MPI_FLOAT, dst, 6, MPI_COMM_WORLD, &request);
			MPI_Isend(h_c12 + (dst * nyinterior * fdm->nzpad * fdm->nxpad), (nyinterior * fdm->nzpad * fdm->nxpad), MPI_FLOAT, dst, 7, MPI_COMM_WORLD, &request);
			MPI_Isend(h_c13 + (dst * nyinterior * fdm->nzpad * fdm->nxpad), (nyinterior * fdm->nzpad * fdm->nxpad), MPI_FLOAT, dst, 8, MPI_COMM_WORLD, &request);
			MPI_Isend(h_c23 + (dst * nyinterior * fdm->nzpad * fdm->nxpad), (nyinterior * fdm->nzpad * fdm->nxpad), MPI_FLOAT, dst, 9, MPI_COMM_WORLD, &request);	
		}
	}
	else {
		// allocate space to receive sub-domain from rank0
		h_ro=(float*)malloc(nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
		h_c11=(float*)malloc(nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
		h_c22=(float*)malloc(nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
		h_c33=(float*)malloc(nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
		h_c44=(float*)malloc(nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
		h_c55=(float*)malloc(nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
		h_c66=(float*)malloc(nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
		h_c12=(float*)malloc(nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
		h_c13=(float*)malloc(nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
		h_c23=(float*)malloc(nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
				
		// receive stiffness and density sub-domain from rank0
		MPI_Recv(h_ro,  (nyinterior * fdm->nzpad * fdm->nxpad), MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
		MPI_Recv(h_c11, (nyinterior * fdm->nzpad * fdm->nxpad), MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);
		MPI_Recv(h_c22, (nyinterior * fdm->nzpad * fdm->nxpad), MPI_FLOAT, 0, 2, MPI_COMM_WORLD, &status);
		MPI_Recv(h_c33, (nyinterior * fdm->nzpad * fdm->nxpad), MPI_FLOAT, 0, 3, MPI_COMM_WORLD, &status);
		MPI_Recv(h_c44, (nyinterior * fdm->nzpad * fdm->nxpad), MPI_FLOAT, 0, 4, MPI_COMM_WORLD, &status);
		MPI_Recv(h_c55, (nyinterior * fdm->nzpad * fdm->nxpad), MPI_FLOAT, 0, 5, MPI_COMM_WORLD, &status);
		MPI_Recv(h_c66, (nyinterior * fdm->nzpad * fdm->nxpad), MPI_FLOAT, 0, 6, MPI_COMM_WORLD, &status);
		MPI_Recv(h_c12, (nyinterior * fdm->nzpad * fdm->nxpad), MPI_FLOAT, 0, 7, MPI_COMM_WORLD, &status);
		MPI_Recv(h_c13, (nyinterior * fdm->nzpad * fdm->nxpad), MPI_FLOAT, 0, 8, MPI_COMM_WORLD, &status);
		MPI_Recv(h_c23, (nyinterior * fdm->nzpad * fdm->nxpad), MPI_FLOAT, 0, 9, MPI_COMM_WORLD, &status);
		
		// copy stiffness and density sub-domain to GPU
		cudaMemcpy(d_ro, h_ro, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_c11, h_c11, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_c22, h_c22, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_c33, h_c33, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_c44, h_c44, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_c55, h_c55, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_c66, h_c66, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_c12, h_c12, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_c13, h_c13, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_c23, h_c23, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		
	}
	
	/*------------------------------------------------------------*/
	/* Boundary condition setup */
	float *d_bzl_s, *d_bzh_s;
	float *d_bxl_s, *d_bxh_s;
	float *d_byl_s, *d_byh_s;
	
	float *h_bzl_s, *h_bzh_s;
	float *h_bxl_s, *h_bxh_s;
	float *h_byl_s, *h_byh_s;
	
	float spo;
	if(dabc) {
		/* one-way abc setup   */
		float *vs1 = (float*)malloc(fdm->nzpad * fdm->nxpad * nyinterior * sizeof(float));
		for (iy = 0; iy < nyinterior; iy++) {
		    for (ix = 0; ix < fdm->nxpad; ix++) {
				for (iz = 0; iz < fdm->nzpad; iz++) {
					vs1[iy * fdm->nzpad * fdm->nxpad + iz * fdm->nxpad + ix] = sqrt(h_c55[iy * fdm->nxpad * fdm->nzpad + iz * fdm->nxpad + ix] / h_ro[iy * fdm->nxpad * fdm->nzpad + iz * fdm->nxpad + ix]);
				}
		    }
		}
		
		float d;
		
		h_bzl_s = (float*)malloc(fdm->nxpad * nyinterior * sizeof(float));
		h_bzh_s = (float*)malloc(fdm->nxpad * nyinterior * sizeof(float));
		for (int ix = 0; ix < fdm->nxpad; ix++){
			for (int iy = 0; iy < nyinterior; iy++){
				d = vs1[iy * fdm->nzpad * fdm->nxpad + NOP * fdm->nxpad + ix] *dt/fdm->dz;
				h_bzl_s[iy * fdm->nxpad + ix] = (1-d)/(1+d);
			    d = vs1[iy * fdm->nzpad * fdm->nxpad + (fdm->nzpad-NOP-1) * fdm->nxpad + ix] *dt/fdm->dz;
			 	h_bzh_s[iy * fdm->nxpad + ix] = (1-d)/(1+d);
			}
		}
		cudaMalloc((void**)&d_bzl_s, fdm->nxpad * nyinterior * sizeof(float));
		cudaMalloc((void**)&d_bzh_s, fdm->nxpad * nyinterior * sizeof(float));
		cudaMemcpy(d_bzl_s, h_bzl_s, fdm->nxpad * nyinterior * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_bzh_s, h_bzh_s, fdm->nxpad * nyinterior * sizeof(float), cudaMemcpyHostToDevice);
		
		
		h_bxl_s = (float*)malloc(fdm->nzpad * nyinterior * sizeof(float));
		h_bxh_s = (float*)malloc(fdm->nzpad * nyinterior * sizeof(float));
		for (int iz = 0; iz < fdm->nzpad; iz++){
			for (int iy = 0; iy < nyinterior; iy++){
				d = vs1[iy * fdm->nzpad * fdm->nxpad + iz * fdm->nxpad + NOP] *dt/fdm->dx;
				h_bxl_s[iy * fdm->nzpad + iz] = (1-d)/(1+d);
			    d = vs1[iy * fdm->nzpad * fdm->nxpad + iz * fdm->nxpad + (fdm->nxpad-NOP-1)] *dt/fdm->dx;
			 	h_bxh_s[iy * fdm->nzpad + iz] = (1-d)/(1+d);
			}
		}
		cudaMalloc((void**)&d_bxl_s, fdm->nzpad * nyinterior * sizeof(float));
		cudaMalloc((void**)&d_bxh_s, fdm->nzpad * nyinterior * sizeof(float));
		cudaMemcpy(d_bxl_s, h_bxl_s, fdm->nzpad * nyinterior * sizeof(float), cudaMemcpyHostToDevice);
		cudaMemcpy(d_bxh_s, h_bxh_s, fdm->nzpad * nyinterior * sizeof(float), cudaMemcpyHostToDevice);
		
		
		if (rank == 0){
			h_byl_s = (float*)malloc(fdm->nzpad * fdm->nxpad * sizeof(float));
			for (int ix = 0; ix < fdm->nxpad; ix++){
				for (int iz = 0; iz < fdm->nzpad; iz++){
					d = vs1[NOP * fdm->nzpad * fdm->nxpad + iz * fdm->nxpad + ix] *dt/fdm->dy;
					h_byl_s[ix * fdm->nzpad + iz] = (1-d)/(1+d);
				}
			}
			cudaMalloc((void**)&d_byl_s, fdm->nzpad * fdm->nxpad * sizeof(float));
			cudaMemcpy(d_byl_s, h_byl_s, fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		}
		
		
		if (rank == size-1){
			h_byh_s = (float*)malloc(fdm->nzpad * fdm->nxpad * sizeof(float));
			for (int ix = 0; ix < fdm->nxpad; ix++){
				for (int iz = 0; iz < fdm->nzpad; iz++){
					d = vs1[(nyinterior-NOP-1) * fdm->nzpad * fdm->nxpad + iz * fdm->nxpad + ix] *dt/fdm->dy;
					h_byh_s[ix * fdm->nzpad + iz] = (1-d)/(1+d);
				}
			}
			cudaMalloc((void**)&d_byh_s, fdm->nzpad * fdm->nxpad * sizeof(float));
			cudaMemcpy(d_byh_s, h_byh_s, fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyHostToDevice);
		}
		
		/* sponge set up */
		// sponge coefficients are calculated inside the sponge kernel on GPU based on spo
		spo = (sqrt(2.0) * 4.0f * nb);
    }
	/*------------------------------------------------------------*/
	
	
	/*------------------------------------------------------------*/	
	/* allocate wavefield arrays */
    
	cudaMalloc((void **)&d_umz, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_uoz, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_upz, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));	
	cudaMalloc((void **)&d_uaz, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	
	cudaMalloc((void **)&d_umx, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_uox, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_upx, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));	
	cudaMalloc((void **)&d_uax, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	
	cudaMalloc((void **)&d_umy, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_uoy, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_upy, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));	
	cudaMalloc((void **)&d_uay, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	
	cudaMalloc((void **)&d_tzz, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_tyy, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_txx, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_txy, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_tyz, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMalloc((void **)&d_tzx, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	
	sf_check_gpu_error(rank % ngpu, "allocate grid arrays");
	
	
	cudaMemset(d_umz, 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_uoz, 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_upz, 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_uaz, 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));

	cudaMemset(d_umx, 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_uox, 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_upx, 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_uax, 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	
	cudaMemset(d_umy, 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_uoy, 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_upy, 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_uay, 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	
	cudaMemset(d_tzz, 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_tyy, 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_txx, 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_txy, 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_tyz, 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	cudaMemset(d_tzx, 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
	
	sf_check_gpu_error(rank % ngpu, "initialize grid arrays");
	
	
	// Used for exchanging halo regions between neighboring GPUs
	float *h_tzz_l_send, *h_tzz_l_recv, *h_tzz_h_send, *h_tzz_h_recv;
	float *h_tyy_l_send, *h_tyy_l_recv, *h_tyy_h_send, *h_tyy_h_recv;
	float *h_txx_l_send, *h_txx_l_recv, *h_txx_h_send, *h_txx_h_recv;
	float *h_txy_l_send, *h_txy_l_recv, *h_txy_h_send, *h_txy_h_recv;
	float *h_tyz_l_send, *h_tyz_l_recv, *h_tyz_h_send, *h_tyz_h_recv;
	float *h_tzx_l_send, *h_tzx_l_recv, *h_tzx_h_send, *h_tzx_h_recv;
	
	float *h_uoz_l_send, *h_uoz_l_recv, *h_uoz_h_send, *h_uoz_h_recv;
	float *h_uoy_l_send, *h_uoy_l_recv, *h_uoy_h_send, *h_uoy_h_recv;
	float *h_uox_l_send, *h_uox_l_recv, *h_uox_h_send, *h_uox_h_recv;
	
	if (size > 1){
		h_tzz_l_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_tzz_l_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_tzz_h_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_tzz_h_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));

		h_tyy_l_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_tyy_l_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_tyy_h_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_tyy_h_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));

		h_txx_l_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_txx_l_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_txx_h_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_txx_h_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));

		h_txy_l_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_txy_l_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_txy_h_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_txy_h_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));

		h_tyz_l_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_tyz_l_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_tyz_h_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_tyz_h_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));

		h_tzx_l_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_tzx_l_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_tzx_h_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_tzx_h_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));

		h_uoz_l_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_uoz_l_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_uoz_h_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_uoz_h_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));

		h_uoy_l_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_uoy_l_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_uoy_h_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_uoy_h_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));

		h_uox_l_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_uox_l_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_uox_h_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_uox_h_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
	}
	

	/*------------------------------------------------------------*/
    /* precompute 1/ro * dt^2 */
	/*------------------------------------------------------------*/
	dim3 dimGrid1(ceil(fdm->nxpad/8.0f),ceil(fdm->nzpad/8.0f),ceil(nyinterior/8.0f));
	dim3 dimBlock1(8,8,8);
	computeRo<<<dimGrid1, dimBlock1>>>(d_ro, dt, fdm->nxpad, fdm->nzpad, nyinterior);
	sf_check_gpu_error(rank % ngpu, "computeRo Kernel");
		
		
	/*------------------------------------------------------------*/
    /* 
     *  MAIN LOOP
     */
    /*------------------------------------------------------------*/
	if(verb && rank == 0) fprintf(stderr,"\n");
	for (it=0; it<nt; it++) {
		if(verb && rank == 0) fprintf(stderr,"\b\b\b\b\b%d",it);
	
		/*------------------------------------------------------------*/
		/* from displacement to strain                                */
		/*		- Compute strains from displacements as in equation 1 */
		/*			- Step #1	(Steps denoted are as in Figure 2)	  */
 		/*------------------------------------------------------------*/
			dim3 dimGrid2(ceil((fdm->nxpad-2*NOP)/24.0f), ceil((fdm->nzpad-2*NOP)/24.0f));
			dim3 dimBlock2(24,24,1);
			dispToStrain<<<dimGrid2, dimBlock2, 32*32*3*sizeof(float)>>>(fdm->nxpad, nylocal, fdm->nzpad, d_uox, d_uoy, d_uoz, d_txx, d_tyy, d_tzz, d_txy, d_tyz, d_tzx, idx, idy, idz);
			sf_check_gpu_error(rank % ngpu, "dispToStrain Kernel");
		
		
		/*------------------------------------------------------------*/
		/* from strain to stress                                      */
		/*		- Compute stress from strain as in equation 2		  */
		/*			- Step #2										  */
		/*------------------------------------------------------------*/
			dim3 dimGrid3(ceil(fdm->nxpad/192.0f), ceil(fdm->nzpad/1.0f), ceil(nyinterior/1.0f));
			dim3 dimBlock3(192,1,1);
			strainToStress<<<dimGrid3, dimBlock3>>>(rank, fdm->nxpad, fdm->nzpad, nyinterior, d_c11, d_c12, d_c13, d_c22, d_c23, d_c33, d_c44, d_c55, d_c66, d_txx, d_tyy, d_tzz, d_txy, d_tyz, d_tzx);
			sf_check_gpu_error(rank % ngpu, "strainToStress Kernel");


		/*------------------------------------------------------------*/
		/* free surface                                               */
		/*		- sets the z-component of stress tensor along the	  */
		/*			free surface boundary to 0						  */
		/*			- Step #3										  */
		/*------------------------------------------------------------*/
			if(fsrf) {
				dim3 dimGrid4(ceil(fdm->nxpad/8.0f), ceil(fdm->nb/8.0f), ceil(nyinterior/8.0f));
				dim3 dimBlock4(8,8,8);
				freeSurf<<<dimGrid4, dimBlock4>>>(rank, fdm->nxpad, nyinterior, fdm->nzpad, fdm->nb, d_tzz, d_tyz, d_tzx);
				sf_check_gpu_error(rank % ngpu, "freeSurf Kernel");
			}
		
		
		/*------------------------------------------------------------*/
		/* inject stress source                                       */
		/*		- Step #4											  */
		/*------------------------------------------------------------*/
			if(ssou) {
				dim3 dimGrid5(ns, 1, 1);
				dim3 dimBlock5(2 * nbell + 1, 2 * nbell + 1, 1);
				lint3d_bell_gpu<<<dimGrid5, dimBlock5>>>(rank, it, nc, ns, 0, nbell, fdm->nxpad, nyinterior, fdm->nzpad, d_tzz, d_bell, d_Sjx, d_Sjz, d_Sjy, d_ww, d_Sw000, d_Sw001, d_Sw010, d_Sw011, d_Sw100, d_Sw101, d_Sw110, d_Sw111);
				lint3d_bell_gpu<<<dimGrid5, dimBlock5>>>(rank, it, nc, ns, 1, nbell, fdm->nxpad, nyinterior, fdm->nzpad, d_txx, d_bell, d_Sjx, d_Sjz, d_Sjy, d_ww, d_Sw000, d_Sw001, d_Sw010, d_Sw011, d_Sw100, d_Sw101, d_Sw110, d_Sw111);
				lint3d_bell_gpu<<<dimGrid5, dimBlock5>>>(rank, it, nc, ns, 2, nbell, fdm->nxpad, nyinterior, fdm->nzpad, d_tyy, d_bell, d_Sjx, d_Sjz, d_Sjy, d_ww, d_Sw000, d_Sw001, d_Sw010, d_Sw011, d_Sw100, d_Sw101, d_Sw110, d_Sw111);
				sf_check_gpu_error(rank % ngpu, "lint3d_bell_gpu Kernel");
			}
	
		/*------------------------------------------------------------*/
		/* exchange halo regions of d_t arrays between GPUs           */
		/*------------------------------------------------------------*/
		if (size > 1){	// using multiple GPUs, must exchange halo regions between neighboring GPUs
			if (rank == 0){
				// get high halo region from d_t arrays on GPU and send to rank+1
				cudaMemcpy(h_tzz_h_send, d_tzz + (fdm->nxpad * fdm->nzpad * (nyinterior - 4)), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_tzz_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tyy_h_send, d_tyy + (fdm->nxpad * fdm->nzpad * (nyinterior - 4)), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_tyy_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_txx_h_send, d_txx + (fdm->nxpad * fdm->nzpad * (nyinterior - 4)), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_txx_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_txy_h_send, d_txy + (fdm->nxpad * fdm->nzpad * (nyinterior - 4)), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_txy_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 3, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tyz_h_send, d_tyz + (fdm->nxpad * fdm->nzpad * (nyinterior - 4)), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_tyz_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 4, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tzx_h_send, d_tzx + (fdm->nxpad * fdm->nzpad * (nyinterior - 4)), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_tzx_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 5, MPI_COMM_WORLD, &request);

				// receive low halo region of d_t arrays from from rank+1 and copy to GPU
				MPI_Recv(h_tzz_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_tzz + (fdm->nxpad * fdm->nzpad * nyinterior), h_tzz_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				MPI_Recv(h_tyy_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_tyy + (fdm->nxpad * fdm->nzpad * nyinterior), h_tyy_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				MPI_Recv(h_txx_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_txx + (fdm->nxpad * fdm->nzpad * nyinterior), h_txx_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				MPI_Recv(h_txy_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 3, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_txy + (fdm->nxpad * fdm->nzpad * nyinterior), h_txy_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				MPI_Recv(h_tyz_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 4, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_tyz + (fdm->nxpad * fdm->nzpad * nyinterior), h_tyz_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				MPI_Recv(h_tzx_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 5, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_tzx + (fdm->nxpad * fdm->nzpad * nyinterior), h_tzx_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
			}
			else if (rank == size-1){
				// get low halo region from d_t arrays on GPU and send to rank-1
				cudaMemcpy(h_tzz_l_send, d_tzz + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_tzz_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tyy_l_send, d_tyy + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_tyy_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_txx_l_send, d_txx + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_txx_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_txy_l_send, d_txy + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_txy_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 3, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tyz_l_send, d_tyz + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_tyz_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 4, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tzx_l_send, d_tzx + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_tzx_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 5, MPI_COMM_WORLD, &request);

				// receive high halo region of d_t arrays from from rank-1 and copy to GPU
				MPI_Recv(h_tzz_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_tzz, h_tzz_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				MPI_Recv(h_tyy_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_tyy, h_tyy_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				MPI_Recv(h_txx_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_txx, h_txx_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				MPI_Recv(h_txy_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 3, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_txy, h_txy_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				MPI_Recv(h_tyz_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 4, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_tyz, h_tyz_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				MPI_Recv(h_tzx_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 5, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_tzx, h_tzx_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
			}
			else {
				// get low halo region from d_t arrays on GPU and send to rank-1
				cudaMemcpy(h_tzz_l_send, d_tzz + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_tzz_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tyy_l_send, d_tyy + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_tyy_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_txx_l_send, d_txx + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_txx_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_txy_l_send, d_txy + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_txy_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 3, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tyz_l_send, d_tyz + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_tyz_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 4, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tzx_l_send, d_tzx + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_tzx_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 5, MPI_COMM_WORLD, &request);

				// get high halo region from d_t arrays on GPU and send to rank+1
				cudaMemcpy(h_tzz_h_send, d_tzz + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_tzz_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tyy_h_send, d_tyy + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_tyy_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_txx_h_send, d_txx + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_txx_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_txy_h_send, d_txy + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_txy_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 3, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tyz_h_send, d_tyz + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_tyz_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 4, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tzx_h_send, d_tzx + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
				MPI_Isend(h_tzx_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 5, MPI_COMM_WORLD, &request);

				// receive high halo region of d_t arrays from from rank-1 and copy to GPU
				MPI_Recv(h_tzz_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_tzz, h_tzz_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				MPI_Recv(h_tyy_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_tyy, h_tyy_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				MPI_Recv(h_txx_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_txx, h_txx_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				MPI_Recv(h_txy_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 3, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_txy, h_txy_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				MPI_Recv(h_tyz_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 4, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_tyz, h_tyz_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				MPI_Recv(h_tzx_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 5, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_tzx, h_tzx_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);

				// receive low halo region of d_t arrays from from rank+1 and copy to GPU
				MPI_Recv(h_tzz_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_tzz + (fdm->nxpad * fdm->nzpad * 4 + nyinterior * fdm->nxpad * fdm->nzpad), h_tzz_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				MPI_Recv(h_tyy_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_tyy + (fdm->nxpad * fdm->nzpad * 4 + nyinterior * fdm->nxpad * fdm->nzpad), h_tyy_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				MPI_Recv(h_txx_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_txx + (fdm->nxpad * fdm->nzpad * 4 + nyinterior * fdm->nxpad * fdm->nzpad), h_txx_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				MPI_Recv(h_txy_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 3, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_txy + (fdm->nxpad * fdm->nzpad * 4 + nyinterior * fdm->nxpad * fdm->nzpad), h_txy_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				MPI_Recv(h_tyz_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 4, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_tyz + (fdm->nxpad * fdm->nzpad * 4 + nyinterior * fdm->nxpad * fdm->nzpad), h_tyz_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				MPI_Recv(h_tzx_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 5, MPI_COMM_WORLD, &status);
				cudaMemcpy(d_tzx + (fdm->nxpad * fdm->nzpad * 4 + nyinterior * fdm->nxpad * fdm->nzpad), h_tzx_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
			}
		}
	
	
		/*------------------------------------------------------------*/
		/* from stress to acceleration  (first term in RHS of eq. 3)  */
		/*		- Step #5
		/*------------------------------------------------------------*/
			dim3 dimGrid6((fdm->nxpad-2*NOP)/24.0f, (fdm->nzpad-2*NOP)/24.0f);
			dim3 dimBlock6(24,24,1);
			stressToAccel<<<dimGrid6, dimBlock6, 32*32*5*sizeof(float)>>>(fdm->nxpad, fdm->nzpad, nylocal, idx, idy, idz, d_txx, d_tyy, d_tzz, d_txy, d_tzx, d_tyz, d_uax, d_uay, d_uaz);
			sf_check_gpu_error(rank % ngpu, "stressToAccel Kernel");

	
		/*------------------------------------------------------------*/
		/* inject acceleration source  (second term in RHS of eq. 3)  */
		/*		- Step #6											  */
		/*------------------------------------------------------------*/
			if(!ssou) {
				dim3 dimGrid7(ns, 1, 1);
				dim3 dimBlock7(2 * nbell + 1, 2 * nbell + 1, 1);
				lint3d_bell_gpu<<<dimGrid7, dimBlock7>>>(rank, it, nc, ns, 0, nbell, fdm->nxpad, nyinterior, fdm->nzpad, d_uaz, d_bell, d_Sjx, d_Sjz, d_Sjy, d_ww, d_Sw000, d_Sw001, d_Sw010, d_Sw011, d_Sw100, d_Sw101, d_Sw110, d_Sw111);
				lint3d_bell_gpu<<<dimGrid7, dimBlock7>>>(rank, it, nc, ns, 1, nbell, fdm->nxpad, nyinterior, fdm->nzpad, d_uax, d_bell, d_Sjx, d_Sjz, d_Sjy, d_ww, d_Sw000, d_Sw001, d_Sw010, d_Sw011, d_Sw100, d_Sw101, d_Sw110, d_Sw111);
				lint3d_bell_gpu<<<dimGrid7, dimBlock7>>>(rank, it, nc, ns, 2, nbell, fdm->nxpad, nyinterior, fdm->nzpad, d_uay, d_bell, d_Sjx, d_Sjz, d_Sjy, d_ww, d_Sw000, d_Sw001, d_Sw010, d_Sw011, d_Sw100, d_Sw101, d_Sw110, d_Sw111);
				sf_check_gpu_error(rank % ngpu, "lint3d_bell_gpu Kernel");
			
			}
		

		/*------------------------------------------------------------*/
		/* step forward in time                                       */
		/*		- Compute forward time step based on acceleration	  */
		/*			- Step #7
		/*------------------------------------------------------------*/
			dim3 dimGrid8(ceil(fdm->nxpad/192.0f), ceil(fdm->nzpad/1.0f), ceil(nyinterior/1.0f));
			dim3 dimBlock8(192,1,1);
			stepTime<<<dimGrid8, dimBlock8>>>(rank, fdm->nxpad, nyinterior, fdm->nzpad, d_ro, d_uox, d_umx, d_uax, d_upx, d_uoy, d_umy, d_uay, d_upy, d_uoz, d_umz, d_uaz, d_upz);
			sf_check_gpu_error(rank % ngpu, "stepTime Kernel");
		
		
		/* circulate wavefield arrays */
		d_utz=d_umz; d_uty=d_umy; d_utx=d_umx;
		d_umz=d_uoz; d_umy=d_uoy; d_umx=d_uox;
		d_uoz=d_upz; d_uoy=d_upy; d_uox=d_upx;
		d_upz=d_utz; d_upy=d_uty; d_upx=d_utx;
	
	
		/*------------------------------------------------------------*/
		/* apply boundary conditions                                  */
		/*		- Step #8											  */
		/*------------------------------------------------------------*/
		if (dabc){
			
			/*---------------------------------------------------------------*/
			/* apply One-way Absorbing BC as in (Clayton and Enquist, 1977)  */
			/*---------------------------------------------------------------*/
			dim3 dimGrid_abc_XY(ceil(fdm->nxpad/32.0f),ceil(nyinterior/32.0f),2);
			dim3 dimBlock_abc_XY(32,32,1);
			abcone3d_apply_XY<<<dimGrid_abc_XY,dimBlock_abc_XY>>>(rank, fdm->nxpad, nyinterior, fdm->nzpad, d_uox, d_umx, d_bzl_s, d_bzh_s);
			abcone3d_apply_XY<<<dimGrid_abc_XY,dimBlock_abc_XY>>>(rank, fdm->nxpad, nyinterior, fdm->nzpad, d_uoy, d_umy, d_bzl_s, d_bzh_s);
			abcone3d_apply_XY<<<dimGrid_abc_XY,dimBlock_abc_XY>>>(rank, fdm->nxpad, nyinterior, fdm->nzpad, d_uoz, d_umz, d_bzl_s, d_bzh_s);

			dim3 dimGrid_abc_ZY(2, ceil(nyinterior/32.0f), ceil(fdm->nzpad/32.0f));
			dim3 dimBlock_abc_ZY(1,32,32);                                                                  
			abcone3d_apply_ZY<<<dimGrid_abc_ZY,dimBlock_abc_ZY>>>(rank, fdm->nxpad, nyinterior, fdm->nzpad, d_uox, d_umx, d_bxl_s, d_bxh_s);
			abcone3d_apply_ZY<<<dimGrid_abc_ZY,dimBlock_abc_ZY>>>(rank, fdm->nxpad, nyinterior, fdm->nzpad, d_uoy, d_umy, d_bxl_s, d_bxh_s);
			abcone3d_apply_ZY<<<dimGrid_abc_ZY,dimBlock_abc_ZY>>>(rank, fdm->nxpad, nyinterior, fdm->nzpad, d_uoz, d_umz, d_bxl_s, d_bxh_s);
		
			dim3 dimGrid_abc_XZ(ceil(fdm->nxpad/32.0f),1,ceil(fdm->nzpad/32.0f));
			dim3 dimBlock_abc_XZ(32,1,32);
			if (rank == 0){
				abcone3d_apply_XZ_low<<<dimGrid_abc_XZ,dimBlock_abc_XZ>>>(fdm->nxpad, fdm->nzpad, d_uox, d_umx, d_byl_s);
				abcone3d_apply_XZ_low<<<dimGrid_abc_XZ,dimBlock_abc_XZ>>>(fdm->nxpad, fdm->nzpad, d_uoy, d_umy, d_byl_s);
				abcone3d_apply_XZ_low<<<dimGrid_abc_XZ,dimBlock_abc_XZ>>>(fdm->nxpad, fdm->nzpad, d_uoz, d_umz, d_byl_s);
			}
			if (rank == size-1){
				abcone3d_apply_XZ_high<<<dimGrid_abc_XZ,dimBlock_abc_XZ>>>(fdm->nxpad, nylocal, fdm->nzpad, d_uox, d_umx, d_byh_s);
				abcone3d_apply_XZ_high<<<dimGrid_abc_XZ,dimBlock_abc_XZ>>>(fdm->nxpad, nylocal, fdm->nzpad, d_uoy, d_umy, d_byh_s);
				abcone3d_apply_XZ_high<<<dimGrid_abc_XZ,dimBlock_abc_XZ>>>(fdm->nxpad, nylocal, fdm->nzpad, d_uoz, d_umz, d_byh_s);
			}
		
		
			/*---------------------------------------------------------------*/
			/* apply Sponge BC as in (Cerjan, et al., 1985)                  */
			/*---------------------------------------------------------------*/
			dim3 dimGrid_spng_XY(ceil(fdm->nxpad/192.0f),nyinterior,1);
			dim3 dimBlock_spng_XY(192,1,1);                                            
			sponge3d_apply_XY<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(rank, d_umz, fdm->nxpad, nyinterior, fdm->nzpad, nb);
			sponge3d_apply_XY<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(rank, d_uoz, fdm->nxpad, nyinterior, fdm->nzpad, nb);
			sponge3d_apply_XY<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(rank, d_upz, fdm->nxpad, nyinterior, fdm->nzpad, nb);
		                                                                 
			sponge3d_apply_XY<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(rank, d_umx, fdm->nxpad, nyinterior, fdm->nzpad, nb);
			sponge3d_apply_XY<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(rank, d_uox, fdm->nxpad, nyinterior, fdm->nzpad, nb);
			sponge3d_apply_XY<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(rank, d_upx, fdm->nxpad, nyinterior, fdm->nzpad, nb);
		                                                                 
			sponge3d_apply_XY<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(rank, d_umy, fdm->nxpad, nyinterior, fdm->nzpad, nb);
			sponge3d_apply_XY<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(rank, d_uoy, fdm->nxpad, nyinterior, fdm->nzpad, nb);
			sponge3d_apply_XY<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(rank, d_upy, fdm->nxpad, nyinterior, fdm->nzpad, nb);
		
		
			dim3 dimGrid_spng_ZY(ceil(nb/8.0f),ceil(fdm->nzpad/8.0f),ceil(nyinterior/8.0f));
			dim3 dimBlock_spng_ZY(8,8,8);
			sponge3d_apply_ZY<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(rank, d_umz, fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
			sponge3d_apply_ZY<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(rank, d_uoz, fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
			sponge3d_apply_ZY<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(rank, d_upz, fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
		
			sponge3d_apply_ZY<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(rank, d_umx, fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
			sponge3d_apply_ZY<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(rank, d_uox, fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
			sponge3d_apply_ZY<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(rank, d_upx, fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
		
			sponge3d_apply_ZY<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(rank, d_umy, fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
			sponge3d_apply_ZY<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(rank, d_uoy, fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
			sponge3d_apply_ZY<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(rank, d_upy, fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
		
		
			dim3 dimGrid_spng_XZ(ceil(fdm->nxpad/192.0f),1,fdm->nzpad);
			dim3 dimBlock_spng_XZ(192,1,1);
			if (rank == 0){
				sponge3d_apply_XZ_low<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_umz, fdm->nxpad, fdm->nzpad, nb);
				sponge3d_apply_XZ_low<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_uoz, fdm->nxpad, fdm->nzpad, nb);
				sponge3d_apply_XZ_low<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_upz, fdm->nxpad, fdm->nzpad, nb);
                                                                                           
				sponge3d_apply_XZ_low<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_umx, fdm->nxpad, fdm->nzpad, nb);
				sponge3d_apply_XZ_low<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_uox, fdm->nxpad, fdm->nzpad, nb);
				sponge3d_apply_XZ_low<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_upx, fdm->nxpad, fdm->nzpad, nb);
                                                                                           
				sponge3d_apply_XZ_low<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_umy, fdm->nxpad, fdm->nzpad, nb);
				sponge3d_apply_XZ_low<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_uoy, fdm->nxpad, fdm->nzpad, nb);
				sponge3d_apply_XZ_low<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_upy, fdm->nxpad, fdm->nzpad, nb);
			}
			if (rank == size-1){
				sponge3d_apply_XZ_high<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_umz, fdm->nxpad, nylocal, fdm->nzpad, nb);
				sponge3d_apply_XZ_high<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_uoz, fdm->nxpad, nylocal, fdm->nzpad, nb);
				sponge3d_apply_XZ_high<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_upz, fdm->nxpad, nylocal, fdm->nzpad, nb);
                                                                                            
				sponge3d_apply_XZ_high<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_umx, fdm->nxpad, nylocal, fdm->nzpad, nb);
				sponge3d_apply_XZ_high<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_uox, fdm->nxpad, nylocal, fdm->nzpad, nb);
				sponge3d_apply_XZ_high<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_upx, fdm->nxpad, nylocal, fdm->nzpad, nb);
                                                                                            
				sponge3d_apply_XZ_high<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_umy, fdm->nxpad, nylocal, fdm->nzpad, nb);
				sponge3d_apply_XZ_high<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_uoy, fdm->nxpad, nylocal, fdm->nzpad, nb);
				sponge3d_apply_XZ_high<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_upy, fdm->nxpad, nylocal, fdm->nzpad, nb);
			}
			
			sf_check_gpu_error(rank % ngpu, "Boundary Condition Kernels");
			
		}
	
	
		/*------------------------------------------------------------*/
		/* exchange halo regions of d_uo arrays between GPUs          */
		/*------------------------------------------------------------*/
		if (size > 1){
				if (rank == 0){
					// get high halo region from d_uo arrays on GPU and send to rank+1
					cudaMemcpy(h_uoz_h_send, d_uoz + (fdm->nxpad * fdm->nzpad * (nyinterior - 4)), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
					MPI_Isend(h_uoz_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &request);
					cudaMemcpy(h_uoy_h_send, d_uoy + (fdm->nxpad * fdm->nzpad * (nyinterior - 4)), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
					MPI_Isend(h_uoy_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &request);
					cudaMemcpy(h_uox_h_send, d_uox + (fdm->nxpad * fdm->nzpad * (nyinterior - 4)), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
					MPI_Isend(h_uox_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &request);
					sf_check_gpu_error(rank % ngpu, "memcpy high to host");
				
					// receive low halo region of d_uo arrays from from rank+1 and copy to GPU
					MPI_Recv(h_uoz_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &status);
					cudaMemcpy(d_uoz + (fdm->nxpad * fdm->nzpad * nyinterior), h_uoz_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
					MPI_Recv(h_uoy_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &status);
					cudaMemcpy(d_uoy + (fdm->nxpad * fdm->nzpad * nyinterior), h_uoy_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
					MPI_Recv(h_uox_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &status);
					cudaMemcpy(d_uox + (fdm->nxpad * fdm->nzpad * nyinterior), h_uox_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				
				}
				else if (rank == size-1){
					// get low halo region from d_uo arrays on GPU and send to rank-1
					cudaMemcpy(h_uoz_l_send, d_uoz + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
					MPI_Isend(h_uoz_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &request);
					cudaMemcpy(h_uoy_l_send, d_uoy + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
					MPI_Isend(h_uoy_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &request);
					cudaMemcpy(h_uox_l_send, d_uox + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
					MPI_Isend(h_uox_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, &request);
					sf_check_gpu_error(rank % ngpu, "memcpy low to host");
				
					// receive high halo region of d_uo arrays from from rank-1 and copy to GPU
					MPI_Recv(h_uoz_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &status);
					cudaMemcpy(d_uoz, h_uoz_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
					MPI_Recv(h_uoy_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &status);
					cudaMemcpy(d_uoy, h_uoy_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
					MPI_Recv(h_uox_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, &status);
					cudaMemcpy(d_uox, h_uox_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				}
				else {
					// get low halo region from d_uo arrays on GPU and send to rank-1
					cudaMemcpy(h_uoz_l_send, d_uoz + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
					MPI_Isend(h_uoz_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &request);
					cudaMemcpy(h_uoy_l_send, d_uoy + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
					MPI_Isend(h_uoy_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &request);
					cudaMemcpy(h_uox_l_send, d_uox + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
					MPI_Isend(h_uox_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, &request);
				
					// get high halo region from d_uo arrays on GPU and send to rank+1
					cudaMemcpy(h_uoz_h_send, d_uoz + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
					MPI_Isend(h_uoz_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &request);
					cudaMemcpy(h_uoy_h_send, d_uoy + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
					MPI_Isend(h_uoy_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &request);
					cudaMemcpy(h_uox_h_send, d_uox + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDeviceToHost);
					MPI_Isend(h_uox_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &request);
				
					// receive high halo region of d_uo arrays from from rank-1 and copy to GPU
					MPI_Recv(h_uoz_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &status);
					cudaMemcpy(d_uoz, h_uoz_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
					MPI_Recv(h_uoy_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &status);
					cudaMemcpy(d_uoy, h_uoy_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
					MPI_Recv(h_uox_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, &status);
					cudaMemcpy(d_uox, h_uox_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				
					// receive low halo region of d_uo arrays from from rank+1 and copy to GPU
					MPI_Recv(h_uoz_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &status);
					cudaMemcpy(d_uoz + (fdm->nxpad * fdm->nzpad * 4 + nyinterior * fdm->nxpad * fdm->nzpad), h_uoz_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
					MPI_Recv(h_uoy_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &status);
					cudaMemcpy(d_uoy + (fdm->nxpad * fdm->nzpad * 4 + nyinterior * fdm->nxpad * fdm->nzpad), h_uoy_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
					MPI_Recv(h_uox_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &status);
					cudaMemcpy(d_uox + (fdm->nxpad * fdm->nzpad * 4 + nyinterior * fdm->nxpad * fdm->nzpad), h_uox_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyHostToDevice);
				}
		}
	
	
		/*------------------------------------------------------------*/
		/* cut wavefield and save 									  */
		/*		- Step #9											  */
		/*------------------------------------------------------------*/
		if(snap && it%jsnap==0) {
		
			if (rank == 0){	// accumulate wavefield data and write to file
			
				// Copy displacements back to CPU excluding the ghost cells from neighbor GPU
				cudaMemcpy(h_uox, d_uox, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDeviceToHost);
				cudaMemcpy(h_uoy, d_uoy, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDeviceToHost);
				cudaMemcpy(h_uoz, d_uoz, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDeviceToHost);
				sf_check_gpu_error(rank % ngpu, "Memcpy wavefield");
			
				// write rank0's wavefield into the output arrays
				for (int y = 0; y < nyinterior; y++){
					for (int z = 0; z < fdm->nzpad; z++){
						for (int x = 0; x < fdm->nxpad; x++){
							uox[y][x][z] = h_uox[y * fdm->nzpad * fdm->nxpad + z * fdm->nxpad + x];
							uoy[y][x][z] = h_uoy[y * fdm->nzpad * fdm->nxpad + z * fdm->nxpad + x];
							uoz[y][x][z] = h_uoz[y * fdm->nzpad * fdm->nxpad + z * fdm->nxpad + x];
						}
					}
				}
			
				// receive wavefield from other GPUs and write to output arrays
				for (int r = 1; r < size; r++){
					MPI_Recv(h_uox, nyinterior * fdm->nzpad * fdm->nxpad, MPI_FLOAT, r, 0, MPI_COMM_WORLD, &status);
					MPI_Recv(h_uoy, nyinterior * fdm->nzpad * fdm->nxpad, MPI_FLOAT, r, 1, MPI_COMM_WORLD, &status);
					MPI_Recv(h_uoz, nyinterior * fdm->nzpad * fdm->nxpad, MPI_FLOAT, r, 2, MPI_COMM_WORLD, &status);
				
					for (int y = 0; y < nyinterior; y++){
						for (int z = 0; z < fdm->nzpad; z++){
							for (int x = 0; x < fdm->nxpad; x++){
								uox[r * nyinterior + y][x][z] = h_uox[y * fdm->nzpad * fdm->nxpad + z * fdm->nxpad + x];
								uoy[r * nyinterior + y][x][z] = h_uoy[y * fdm->nzpad * fdm->nxpad + z * fdm->nxpad + x];
								uoz[r * nyinterior + y][x][z] = h_uoz[y * fdm->nzpad * fdm->nxpad + z * fdm->nxpad + x];
							}
						}
					}
				}

				// Write complete wavefield to output file
				cut3d(uoz,uc,fdm,acz,acx,acy);
				sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);

				cut3d(uox,uc,fdm,acz,acx,acy);
				sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);

				cut3d(uoy,uc,fdm,acz,acx,acy);
				sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);
			
			}
			else {	// send wavefield data to rank0
			
				// Copy displacements back to CPU excluding the ghost cells from neighbor GPU
				cudaMemcpy(h_uox, d_uox + 4 * fdm->nxpad * fdm->nzpad, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDeviceToHost);
				cudaMemcpy(h_uoy, d_uoy + 4 * fdm->nxpad * fdm->nzpad, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDeviceToHost);
				cudaMemcpy(h_uoz, d_uoz + 4 * fdm->nxpad * fdm->nzpad, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDeviceToHost);
			
				MPI_Send(h_uox, nyinterior * fdm->nzpad * fdm->nxpad, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
				MPI_Send(h_uoy, nyinterior * fdm->nzpad * fdm->nxpad, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
				MPI_Send(h_uoz, nyinterior * fdm->nzpad * fdm->nxpad, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
			}
		}
		
		/*------------------------------------------------------------*/
		/* extract receiver data									  */
		/*------------------------------------------------------------*/
		if(it%jdata==0) {
			if (interp){	// if using linear interpolation
				cudaMemset(d_dd, 0, nr*nc*sizeof(float));
				dim3 dimGrid_extract(MIN(nr,ceil(nr/1024.0f)), 1, 1);
				dim3 dimBlock_extract(MIN(nr, 1024), 1, 1);
				lint3d_extract_gpu<<<dimGrid_extract, dimBlock_extract>>>(rank, d_dd, nr, fdm->nxpad, nyinterior, fdm->nzpad, d_uoz, d_uox, d_uoy, d_Rjz, d_Rjx, d_Rjy, d_Rw000, d_Rw001, d_Rw010, d_Rw011, d_Rw100, d_Rw101, d_Rw110, d_Rw111);
				sf_check_gpu_error(rank % ngpu, "lint3d_extract kernel");
			
				cudaMemcpy(h_dd, d_dd, nr * nc * sizeof(float), cudaMemcpyDeviceToHost);
			}
			else {
				cudaMemset(d_dd, 0, nr*nc*sizeof(float));
				dim3 dimGrid_extract(MIN(nr,ceil(nr/1024.0f)), 1, 1);
				dim3 dimBlock_extract(MIN(nr, 1024), 1, 1);
				extract_gpu<<<dimGrid_extract, dimBlock_extract>>>(rank, d_dd, nr, fdm->nxpad, nyinterior, fdm->nzpad, d_uoz, d_uox, d_uoy, d_Rjz, d_Rjx, d_Rjy);
				sf_check_gpu_error(rank % ngpu, "extract_gpu kernel");
			
				cudaMemcpy(h_dd, d_dd, nr * nc * sizeof(float), cudaMemcpyDeviceToHost);
			}
		
			// MPI Reduce on the h_dd arrays from each GPU to accumulate reciever data
			MPI_Reduce(h_dd, h_dd_combined, nr * nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
		
			if (rank == 0){
				// write receiver data to ourput file
				sf_floatwrite(h_dd_combined, nr*nc, Fdat);
			}
		
		}

	} // end of MAIN LOOP
	
	
	/*------------------------------------------------------------*/
    /* deallocate host arrays */

	free(**ww); free(*ww); free(ww);
    free(ss);
    free(rr);

	if (snap && rank == 0){
		free(**uc);  free(*uc);  free(uc);
		free(**uoz); free(*uoz); free(uoz);
	    free(**uox); free(*uox); free(uox);
	    free(**uoy); free(*uoy); free(uoy);
	}
	if (snap){
		free(h_uox); free(h_uoy); free(h_uoz);
	}
	
	
	/*------------------------------------------------------------*/
    /* deallocate GPU arrays */

	cudaFree(d_bell); cudaFree(d_ww); cudaFree(d_dd);
	
	cudaFree(d_Sw000); cudaFree(d_Sw001); cudaFree(d_Sw010); cudaFree(d_Sw011); cudaFree(d_Sw100); cudaFree(d_Sw101); cudaFree(d_Sw110); cudaFree(d_Sw111);
	
	cudaFree(d_Sjz); cudaFree(d_Sjx); cudaFree(d_Sjy);
	
	if (interp){
		cudaFree(d_Rw000); cudaFree(d_Rw001); cudaFree(d_Rw010); cudaFree(d_Rw011); cudaFree(d_Rw100); cudaFree(d_Rw101); cudaFree(d_Rw110); cudaFree(d_Rw111);
	}
	
	cudaFree(d_Rjz); cudaFree(d_Rjx); cudaFree(d_Rjy);
	
	cudaFree(d_c11); cudaFree(d_c22); cudaFree(d_c33); cudaFree(d_c44); cudaFree(d_c55); cudaFree(d_c66); cudaFree(d_c12); cudaFree(d_c13); cudaFree(d_c23);
	
	if (dabc){
		cudaFree(d_bzl_s); cudaFree(d_bzh_s);
		cudaFree(d_bxl_s); cudaFree(d_bxh_s);
		cudaFree(d_byl_s); cudaFree(d_byh_s);
	}
	
	cudaFree(d_ro);
	
	cudaFree(d_umz); cudaFree(d_uoz); cudaFree(d_upz); cudaFree(d_uaz);
	cudaFree(d_umx); cudaFree(d_uox); cudaFree(d_upx); cudaFree(d_uax);
	cudaFree(d_umy); cudaFree(d_uoy); cudaFree(d_upy); cudaFree(d_uay);
	
	cudaFree(d_tzz); cudaFree(d_tyy); cudaFree(d_txx); cudaFree(d_txy); cudaFree(d_tyz); cudaFree(d_tzx);
	
	
	MPI_Finalize();
	sf_close();
	exit(0);
	
}