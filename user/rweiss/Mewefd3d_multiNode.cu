/* 3D elastic time-domain FD modeling with multiple GPUs coordinated via MPI and p2p*/

/*
  Authors: Robin M. Weiss and Jeffrey Shragge

  This code is a GPU-enabled version of the ewefd3d module from the Madagascar
  software package (see: http://www.reproducibility.org).  It implements a 3D
  Finite-Difference Time Domain solver for the elastice wave equation with 
  2nd- and 8th- order temporal and spatial accuracy, respectively.  Computation
  is distributed across an arbitrary number of GPU devices and coordinted by MPI.
  Inter-GPU communication that occurs within a node is carried out with direct p2p
  communication over the PCIe bus via cudaMemCpy(DeviceToDevice).  Communication
  between GPUs located in seperate nodes occurs via MPI send/receive calls. 
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

/*

SOMETHING YOU SHOULD PROBABLY KNOW ABOUT THIS CODE:

All arrays used internally in this program order axes Y, Z, X (slowest to fastest)
However, madagascar wants axes to be ordered Y, X, Z and this is the standard in RSF files (i think)

Therefore, whenever reading arrays from a file, or writing data to a file, the X and Z axes must be transposed.

I appologize for the inconvenience and encourage someone to take the time to fix this and make the internal representation the same as the external.

sorry for being difficult,
---Robin

*/



#include <stdio.h>
#include <mpi.h>
#include <cuda.h>
#include <cuda_runtime_api.h>

extern "C" {
#include <rsf.h>
}

#include "fdutil.c"
#include "ewefd3d_multiNode_kernels.cu"

#define MIN(x, y) (((x) < (y)) ? (x) : (y))
#define NOP 4 /* derivative operator half-size */


// checks the current GPU device for an error flag and prints to stderr
static void sf_check_gpu_error (int rank, const char *msg) {
    cudaError_t err = cudaGetLastError ();
     if (cudaSuccess != err)
        sf_error ("Cuda error on node %d: %s: %s", rank, msg, cudaGetErrorString (err));
}




// entry point
int main (int argc, char* argv[]) {
	
	// Initialize MPI
	int rank, nNode;
	MPI_Init (&argc, &argv);				/* start MPI */
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);	/* get current process id */
	MPI_Comm_size (MPI_COMM_WORLD, &nNode);	/* get number of nodes */
	MPI_Status status;
	MPI_Request request;
	
	MPI_Request *recvRequests = (MPI_Request*)malloc(12 * sizeof(MPI_Request));
	MPI_Status *recvStatuses = (MPI_Status*)malloc(12 * sizeof(MPI_Status));
	
	// Initialize RSF
	sf_init(argc, argv);

	bool verb,fsrf,snap,ssou,dabc,interp,wavSrc;
    int  jsnap,ntsnap,jdata;

    /* I/O files */
    sf_file Fwav=NULL; /* wavelet   */
    sf_file Fsou=NULL; /* sources   */
    sf_file Frec=NULL; /* receivers */
    sf_file Fccc=NULL; /* velocity  */
    sf_file Fden=NULL; /* density   */
    sf_file Fdat=NULL; /* data      */
    sf_file Fwfl=NULL; /* wavefield */

	sf_file Fum=NULL;
	sf_file Fuo=NULL;

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
    /* displacement: um = U @ t-1; uo = U @ t; up = U @ t+1 */

	// used for writing wavefield to file, only needed if snap=y
	float ***uox, ***uoy, ***uoz;
	float *h_uox, *h_uoy, *h_uoz;

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
	if(! sf_getbool("wavSrc",&wavSrc)) wavSrc=true; /* if yes, look for a source wavelet.  if no, look for initial displacement fields (uo and um) */
    /*------------------------------------------------------------*/


    /*------------------------------------------------------------*/
    /* I/O files */
    Fwav = sf_input ("wav"); /* source wavelet.  must always be present to establish time but contents is ignored if wavSrc=n   */
    Fccc = sf_input ("ccc"); /* stiffness */
    Fden = sf_input ("den"); /* density   */
	Fsou = sf_input ("sou"); /* source locations   */
    Frec = sf_input ("rec"); /* receiver locations */
    Fwfl = sf_output("wfl"); /* wavefield */
    Fdat = sf_output("out"); /* data      */
    
	if (!wavSrc){
		Fum  = sf_input ("um"); /* if wavSrc=n, looks for um file containing the previous displacement timestep. axes: z, x, then y */
		Fuo  = sf_input ("uo"); /* if wavSrc=n, looks for uo file containing the current displacement timestep. axes: z, x, then y */
	}
	/*------------------------------------------------------------*/
	
	/*------------------------------------------------------------*/
    /* initialize GPU */

	int availGPU;
	cudaGetDeviceCount(&availGPU);
	
	int ngpu;
	if(!sf_getint("ngpu",&ngpu)) ngpu=1; 	/* Number of GPUs in each node, must be set to lowest common number of GPUs*/

	if (availGPU < ngpu){
		sf_error("ngpu parameter (%d) is greater than %d GPUs in node %d", ngpu, availGPU, rank);
	}
	
	for (int g = 0; g < ngpu; g++){
		cudaSetDevice(g);
		cudaDeviceSetCacheConfig(cudaFuncCachePreferL1);
	}
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
    /* check dimensions of all other files */
	if (sf_n(sf_iaxa(Fden,1)) != nz) sf_error("Dimension missmatch on z-axis of density model");
	if (sf_n(sf_iaxa(Fden,2)) != nx) sf_error("Dimension missmatch on x-axis of density model");
	if (sf_n(sf_iaxa(Fden,3)) != ny) sf_error("Dimension missmatch on y-axis of density model");
	
	if (!wavSrc){
		if (sf_n(sf_iaxa(Fuo,1)) != nz+2*fdm->nb) sf_error("Dimension missmatch on z-axis of initial displacement field (uo)");
		if (sf_n(sf_iaxa(Fuo,2)) != nx+2*fdm->nb) sf_error("Dimension missmatch on x-axis of initial displacement field (uo)");
		if (sf_n(sf_iaxa(Fuo,3)) != ny+2*fdm->nb) sf_error("Dimension missmatch on y-axis of initial displacement field (uo)");
		
		if (sf_n(sf_iaxa(Fum,1)) != nz+2*fdm->nb) sf_error("Dimension missmatch on z-axis of initial displacement field (um)");
		if (sf_n(sf_iaxa(Fum,2)) != nx+2*fdm->nb) sf_error("Dimension missmatch on x-axis of initial displacement field (um)");
		if (sf_n(sf_iaxa(Fum,3)) != ny+2*fdm->nb) sf_error("Dimension missmatch on y-axis of initial displacement field (um)");
	}
	
	
	/*------------------------------------------------------------*/
    /* compute sub-domain dimmensions (domain decomposition) */


	// check that dimmeionsons are ok for FD kernels
	if ((fdm->nzpad - 8) % 24 != 0){
		sf_error("nz + 2*nb - 8 is not a multiple of 24");
	}
	if ((fdm->nxpad - 8) % 24 != 0){
		sf_error("nx + 2*nb - 8 is not a multiple of 24");
	}
	if ((fdm->nypad % (ngpu * nNode)) != 0){
		sf_error("You are using %d GPUs.\n(ny + 2*nb) must me a multiple of %d\nChange model dimensions or select a different number of GPUs", ngpu, ngpu);
	}
	
	
	int nyinterior = (fdm->nypad / (nNode * ngpu) ); // size of sub-domains in y-dimension EXCLUDING any ghost cells from adjacent GPUs
	
	// nodes require 8 additional ghost slices
	int nylocal = nyinterior + 8;
	
	// if using 1 GPU, this GPU holds the entire domain
	if (ngpu == 1 && nNode == 1){
		nylocal = fdm->nypad;
	}
	
	/*------------------------------------------------------------*/


	/*------------------------------------------------------------*/
    /* setup bell for source injection smoothing */
	if (nbell * 2 + 1 > 32){
		sf_error("nbell must be <= 15\n");
	}
	
	float *h_bell;
	h_bell = (float*)malloc((2*nbell+1)*(2*nbell+1)*(2*nbell+1)*sizeof(float));
	
	float s = 0.5*nbell;
    for (iy=-nbell;iy<=nbell;iy++) {
		for (ix=-nbell;ix<=nbell;ix++) {
	    	for(iz=-nbell;iz<=nbell;iz++) {
				h_bell[(iy + nbell) * (2*nbell+1) * (2*nbell+1) + (iz + nbell) * (2*nbell+1) + (ix + nbell)] = exp(-(iz*iz+ix*ix+iy*iy)/s);
	    	}
		}    
    }

	// copy bell coeficients to the GPUs
	float **d_bell = (float**)malloc(ngpu*sizeof(float*));
	for (int g = 0; g < ngpu; g++){
		cudaSetDevice(g);
		cudaMalloc(&d_bell[g], (2*nbell+1)*(2*nbell+1)*(2*nbell+1)*sizeof(float));
		cudaMemcpy(d_bell[g], h_bell, (2*nbell+1)*(2*nbell+1)*(2*nbell+1)*sizeof(float), cudaMemcpyDefault);
	}
	/*------------------------------------------------------------*/
	
	
	/*------------------------------------------------------------*/
	/* 3D vector components */
    nc=3;
	ac=sf_maxa(nc  ,0,1);
	/*------------------------------------------------------------*/
	
	
	/*------------------------------------------------------------*/
     /* setup output data files and arrays */
    sf_oaxa(Fdat,ar,1);
    sf_oaxa(Fdat,ac,2);

    sf_setn(at,nt/jdata);
    sf_setd(at,dt*jdata);
    sf_oaxa(Fdat,at,3);

    if(snap) {

		// Used to accumulate wavefield data from other GPUs
		uoz=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
		uox=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
		uoy=sf_floatalloc3(fdm->nzpad,fdm->nxpad,fdm->nypad);
		h_uoz = (float*)malloc(fdm->nypad * fdm->nzpad * fdm->nxpad * sizeof(float));
		h_uox = (float*)malloc(fdm->nypad * fdm->nzpad * fdm->nxpad * sizeof(float));
		h_uoy = (float*)malloc(fdm->nypad * fdm->nzpad * fdm->nxpad * sizeof(float));

		nqz=sf_n(az);
		nqx=sf_n(ax);
		nqy=sf_n(ay);

		oqz=sf_o(az);
		oqx=sf_o(ax);
		oqy=sf_o(ay);

		dqz=sf_d(az);
		dqx=sf_d(ax);
		dqy=sf_d(ay);

		acz = sf_maxa(nqz,oqz,dqz); if (rank==0) sf_raxa(acz);
		acx = sf_maxa(nqx,oqx,dqx); if (rank==0) sf_raxa(acx);
		acy = sf_maxa(nqy,oqy,dqy); if (rank==0) sf_raxa(acy);

		uc=sf_floatalloc3(sf_n(acz),sf_n(acx),sf_n(acy));

		ntsnap=0;
		for(it=0; it<nt; it++) {
		    if(it%jsnap==0) ntsnap++;
		}
		sf_setn(at,  ntsnap);
		sf_setd(at,dt*jsnap);
		if(verb && rank==0) sf_raxa(at);

		sf_oaxa(Fwfl,acz,1);
		sf_oaxa(Fwfl,acx,2);
		sf_oaxa(Fwfl,acy,3);
		sf_oaxa(Fwfl,ac, 4);
		sf_oaxa(Fwfl,at, 5);
    }
    /*------------------------------------------------------------*/


	/*------------------------------------------------------------*/
    /* read source wavelet(s) and copy to each GPU (into d_ww) */
    
	float *h_ww;
	float **d_ww = (float**)malloc(ngpu*sizeof(float*));

	if (wavSrc){
		ww=sf_floatalloc3(ns,nc,nt); 
	    sf_floatread(ww[0][0],nt*nc*ns,Fwav);

		h_ww = (float*)malloc(ns*nc*nt*sizeof(float));
		for (int t = 0; t < nt; t++){
			for (int c = 0; c < nc; c++){
				for (int s = 0; s < ns; s++){
					h_ww[t * nc * ns + c * ns + s]=ww[t][c][s];
				}
			}
		}

		for (int g = 0; g < ngpu; g++){
			cudaSetDevice(g);
			cudaMalloc(&d_ww[g], ns*nc*nt*sizeof(float));
			cudaMemcpy(d_ww[g], h_ww, ns*nc*nt*sizeof(float), cudaMemcpyDefault);
		}
	}
    /*------------------------------------------------------------*/


	/*------------------------------------------------------------*/
	/* data array */
	float *h_dd = (float*)malloc(nr * nc * sizeof(float));
	
	float *h_dd_combined = (float*)malloc(nr * nc * sizeof(float));		
	float *h_dd_locals = (float*)malloc(nr * nc * sizeof(float));
	float *h_dd_allNodes = (float*)malloc(nr * nc * sizeof(float));
	
	float **d_dd = (float**)malloc(ngpu*sizeof(float*));
	for (int g = 0; g < ngpu; g++){
		cudaSetDevice(g);
		cudaMalloc(&d_dd[g], nr*nc*sizeof(float));
	}
    /*------------------------------------------------------------*/


	/*------------------------------------------------------------*/
    /* setup source/receiver coordinates */
    ss = (pt3d*) sf_alloc(ns,sizeof(*ss)); 
    rr = (pt3d*) sf_alloc(nr,sizeof(*rr)); 

    pt3dread1(Fsou,ss,ns,3); /* read (x,y,z) coordinates */
    pt3dread1(Frec,rr,nr,3); /* read (x,y,z) coordinates */

	/* calculate 3d linear interpolation coefficients for source locations and copy to each GPU*/
    cs = lint3d_make(ns,ss,fdm);
	float **d_Sw000 = (float**)malloc(ngpu*sizeof(float*));
	float **d_Sw001 = (float**)malloc(ngpu*sizeof(float*));
	float **d_Sw010 = (float**)malloc(ngpu*sizeof(float*));
	float **d_Sw011 = (float**)malloc(ngpu*sizeof(float*));
	float **d_Sw100 = (float**)malloc(ngpu*sizeof(float*));
	float **d_Sw101 = (float**)malloc(ngpu*sizeof(float*));
	float **d_Sw110 = (float**)malloc(ngpu*sizeof(float*));
	float **d_Sw111 = (float**)malloc(ngpu*sizeof(float*));

	for (int g = 0; g < ngpu; g++){
		cudaSetDevice(g);
		cudaMalloc(&d_Sw000[g], ns * sizeof(float));
		cudaMalloc(&d_Sw001[g], ns * sizeof(float));
		cudaMalloc(&d_Sw010[g], ns * sizeof(float));
		cudaMalloc(&d_Sw011[g], ns * sizeof(float));
		cudaMalloc(&d_Sw100[g], ns * sizeof(float));
		cudaMalloc(&d_Sw101[g], ns * sizeof(float));
		cudaMalloc(&d_Sw110[g], ns * sizeof(float));
		cudaMalloc(&d_Sw111[g], ns * sizeof(float));
		cudaMemcpy(d_Sw000[g], cs->w000, ns * sizeof(float), cudaMemcpyDefault);
		cudaMemcpy(d_Sw001[g], cs->w001, ns * sizeof(float), cudaMemcpyDefault);
		cudaMemcpy(d_Sw010[g], cs->w010, ns * sizeof(float), cudaMemcpyDefault);
		cudaMemcpy(d_Sw011[g], cs->w011, ns * sizeof(float), cudaMemcpyDefault);
		cudaMemcpy(d_Sw100[g], cs->w100, ns * sizeof(float), cudaMemcpyDefault);
		cudaMemcpy(d_Sw101[g], cs->w101, ns * sizeof(float), cudaMemcpyDefault);
		cudaMemcpy(d_Sw110[g], cs->w110, ns * sizeof(float), cudaMemcpyDefault);
		cudaMemcpy(d_Sw111[g], cs->w111, ns * sizeof(float), cudaMemcpyDefault);		
	}

	// z, x, and y coordinates of each source
	int **d_Sjz = (int**)malloc(ngpu*sizeof(int*));
	int **d_Sjx = (int**)malloc(ngpu*sizeof(int*));
	int **d_Sjy = (int**)malloc(ngpu*sizeof(int*));
	for (int g = 0; g < ngpu; g++){
		cudaSetDevice(g);
		cudaMalloc(&d_Sjz[g], ns * sizeof(int));
		cudaMalloc(&d_Sjx[g], ns * sizeof(int));
		cudaMalloc(&d_Sjy[g], ns * sizeof(int));
		cudaMemcpy(d_Sjz[g], cs->jz, ns * sizeof(int), cudaMemcpyDefault);
		cudaMemcpy(d_Sjx[g], cs->jx, ns * sizeof(int), cudaMemcpyDefault);
		cudaMemcpy(d_Sjy[g], cs->jy, ns * sizeof(int), cudaMemcpyDefault);	
	}


	/* calculate 3d linear interpolation coefficients for receiver locations and copy to each GPU*/
	cr = lint3d_make(nr,rr,fdm);
	float **d_Rw000 = (float**)malloc(ngpu*sizeof(float*));
	float **d_Rw001 = (float**)malloc(ngpu*sizeof(float*));
	float **d_Rw010 = (float**)malloc(ngpu*sizeof(float*));
	float **d_Rw011 = (float**)malloc(ngpu*sizeof(float*));
	float **d_Rw100 = (float**)malloc(ngpu*sizeof(float*));
	float **d_Rw101 = (float**)malloc(ngpu*sizeof(float*));
	float **d_Rw110 = (float**)malloc(ngpu*sizeof(float*));
	float **d_Rw111 = (float**)malloc(ngpu*sizeof(float*));
	if (interp){
		for (int g = 0; g < ngpu; g++){
			cudaSetDevice(g);
			cudaMalloc(&d_Rw000[g], nr * sizeof(float));
			cudaMalloc(&d_Rw001[g], nr * sizeof(float));
			cudaMalloc(&d_Rw010[g], nr * sizeof(float));
			cudaMalloc(&d_Rw011[g], nr * sizeof(float));
			cudaMalloc(&d_Rw100[g], nr * sizeof(float));
			cudaMalloc(&d_Rw101[g], nr * sizeof(float));
			cudaMalloc(&d_Rw110[g], nr * sizeof(float));
			cudaMalloc(&d_Rw111[g], nr * sizeof(float));
			cudaMemcpy(d_Rw000[g], cr->w000, nr * sizeof(float), cudaMemcpyDefault);
			cudaMemcpy(d_Rw001[g], cr->w001, nr * sizeof(float), cudaMemcpyDefault);
			cudaMemcpy(d_Rw010[g], cr->w010, nr * sizeof(float), cudaMemcpyDefault);
			cudaMemcpy(d_Rw011[g], cr->w011, nr * sizeof(float), cudaMemcpyDefault);
			cudaMemcpy(d_Rw100[g], cr->w100, nr * sizeof(float), cudaMemcpyDefault);
			cudaMemcpy(d_Rw101[g], cr->w101, nr * sizeof(float), cudaMemcpyDefault);
			cudaMemcpy(d_Rw110[g], cr->w110, nr * sizeof(float), cudaMemcpyDefault);
			cudaMemcpy(d_Rw111[g], cr->w111, nr * sizeof(float), cudaMemcpyDefault);
		}
	}

	// z, x, and y coordinates of each receiver
	int **d_Rjz = (int**)malloc(ngpu*sizeof(int*));
	int **d_Rjx = (int**)malloc(ngpu*sizeof(int*));
	int **d_Rjy = (int**)malloc(ngpu*sizeof(int*));
	for (int g = 0; g < ngpu; g++){
		cudaSetDevice(g);
		cudaMalloc(&d_Rjz[g], nr * sizeof(int));
		cudaMalloc(&d_Rjx[g], nr * sizeof(int));
		cudaMalloc(&d_Rjy[g], nr * sizeof(int));
		cudaMemcpy(d_Rjz[g], cr->jz, nr * sizeof(int), cudaMemcpyDefault);
		cudaMemcpy(d_Rjx[g], cr->jx, nr * sizeof(int), cudaMemcpyDefault);
		cudaMemcpy(d_Rjy[g], cr->jy, nr * sizeof(int), cudaMemcpyDefault);	
	}
	/*------------------------------------------------------------*/
	
	
	/*------------------------------------------------------------*/	
    /* setup FD coefficients */
    idz = 1/dz;
    idx = 1/dx;
    idy = 1/dy;
	/*------------------------------------------------------------*/
	
	
	/*------------------------------------------------------------*/ 
	/* read in model density and stiffness arrays */
	
	float *h_ro, *h_c11, *h_c22, *h_c33, *h_c44, *h_c55, *h_c66, *h_c12, *h_c13, *h_c23;
	
	if (nNode == 1){
		float *tt = (float*)malloc(nz*nx*ny*sizeof(float));

	    /* input density */
	    h_ro = (float*)malloc(fdm->nzpad * fdm->nxpad * fdm->nypad * sizeof(float));
	    sf_floatread(tt,nz*nx*ny,Fden);     expand_all(tt, h_ro, fdm->nb, nx, fdm->nxpad, ny, fdm->nypad, nz, fdm->nzpad);

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
	    sf_floatread(tt,nz*nx*ny,Fccc);    expand_all(tt,h_c11,fdm->nb, nx, fdm->nxpad, ny, fdm->nypad, nz, fdm->nzpad);
	    sf_floatread(tt,nz*nx*ny,Fccc);    expand_all(tt,h_c22,fdm->nb, nx, fdm->nxpad, ny, fdm->nypad, nz, fdm->nzpad);    
	    sf_floatread(tt,nz*nx*ny,Fccc);    expand_all(tt,h_c33,fdm->nb, nx, fdm->nxpad, ny, fdm->nypad, nz, fdm->nzpad);    
	    sf_floatread(tt,nz*nx*ny,Fccc);    expand_all(tt,h_c44,fdm->nb, nx, fdm->nxpad, ny, fdm->nypad, nz, fdm->nzpad);
	    sf_floatread(tt,nz*nx*ny,Fccc);    expand_all(tt,h_c55,fdm->nb, nx, fdm->nxpad, ny, fdm->nypad, nz, fdm->nzpad);    
	    sf_floatread(tt,nz*nx*ny,Fccc);    expand_all(tt,h_c66,fdm->nb, nx, fdm->nxpad, ny, fdm->nypad, nz, fdm->nzpad);
	    sf_floatread(tt,nz*nx*ny,Fccc);    expand_all(tt,h_c12,fdm->nb, nx, fdm->nxpad, ny, fdm->nypad, nz, fdm->nzpad);
	    sf_floatread(tt,nz*nx*ny,Fccc);    expand_all(tt,h_c13,fdm->nb, nx, fdm->nxpad, ny, fdm->nypad, nz, fdm->nzpad);
	    sf_floatread(tt,nz*nx*ny,Fccc);    expand_all(tt,h_c23,fdm->nb, nx, fdm->nxpad, ny, fdm->nypad, nz, fdm->nzpad);
	    free(tt);
	}
	else {

		if (rank == 0){
			float *tt = (float*)malloc(((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz*sizeof(float));

		    /* input density */
		    h_ro = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float));
		    sf_floatread(tt,((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz,Fden);	expand_low(tt, h_ro, fdm->nb, nx, fdm->nxpad, (nyinterior - nb) + (nyinterior * (ngpu-1)), (nyinterior * ngpu), nz, fdm->nzpad);
		
		    /* stiffness */
		    h_c11 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float));
		    h_c22 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float)); 
		    h_c33 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float)); 
		    h_c44 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float)); 
		    h_c55 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float)); 
		    h_c66 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float)); 
		    h_c12 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float)); 
		    h_c13 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float)); 
		    h_c23 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float));
		
		    sf_floatread(tt,((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz,Fccc);    expand_low(tt, h_c11, fdm->nb, nx, fdm->nxpad, (nyinterior - nb) + (nyinterior * (ngpu-1)), (nyinterior * ngpu), nz, fdm->nzpad);
			sf_seek(Fccc, 1*ny*nx*nz*sizeof(float), 0);
			sf_floatread(tt,((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz,Fccc);    expand_low(tt, h_c22, fdm->nb, nx, fdm->nxpad, (nyinterior - nb) + (nyinterior * (ngpu-1)), (nyinterior * ngpu), nz, fdm->nzpad);    
			sf_seek(Fccc, 2*ny*nx*nz*sizeof(float), 0);
			sf_floatread(tt,((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz,Fccc);    expand_low(tt, h_c33, fdm->nb, nx, fdm->nxpad, (nyinterior - nb) + (nyinterior * (ngpu-1)), (nyinterior * ngpu), nz, fdm->nzpad);    
			sf_seek(Fccc, 3*ny*nx*nz*sizeof(float), 0);
			sf_floatread(tt,((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz,Fccc);    expand_low(tt, h_c44, fdm->nb, nx, fdm->nxpad, (nyinterior - nb) + (nyinterior * (ngpu-1)), (nyinterior * ngpu), nz, fdm->nzpad);
			sf_seek(Fccc, 4*ny*nx*nz*sizeof(float), 0);
			sf_floatread(tt,((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz,Fccc);    expand_low(tt, h_c55, fdm->nb, nx, fdm->nxpad, (nyinterior - nb) + (nyinterior * (ngpu-1)), (nyinterior * ngpu), nz, fdm->nzpad);    
		    sf_seek(Fccc, 5*ny*nx*nz*sizeof(float), 0);
			sf_floatread(tt,((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz,Fccc);    expand_low(tt, h_c66, fdm->nb, nx, fdm->nxpad, (nyinterior - nb) + (nyinterior * (ngpu-1)), (nyinterior * ngpu), nz, fdm->nzpad);
		    sf_seek(Fccc, 6*ny*nx*nz*sizeof(float), 0);
			sf_floatread(tt,((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz,Fccc);    expand_low(tt, h_c12, fdm->nb, nx, fdm->nxpad, (nyinterior - nb) + (nyinterior * (ngpu-1)), (nyinterior * ngpu), nz, fdm->nzpad);
		    sf_seek(Fccc, 7*ny*nx*nz*sizeof(float), 0);
			sf_floatread(tt,((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz,Fccc);    expand_low(tt, h_c13, fdm->nb, nx, fdm->nxpad, (nyinterior - nb) + (nyinterior * (ngpu-1)), (nyinterior * ngpu), nz, fdm->nzpad);
		    sf_seek(Fccc, 8*ny*nx*nz*sizeof(float), 0);
			sf_floatread(tt,((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz,Fccc);    expand_low(tt, h_c23, fdm->nb, nx, fdm->nxpad, (nyinterior - nb) + (nyinterior * (ngpu-1)), (nyinterior * ngpu), nz, fdm->nzpad);
			free(tt);
		}
		else if (rank = nNode-1){

			float *tt = (float*)malloc(((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz*sizeof(float));

		    /* input density */
		    h_ro = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float));

			sf_seek(Fden, ((nyinterior * ngpu * (rank-1) + (nyinterior - nb) + (nyinterior * (ngpu-1))) * nx * nz)*sizeof(float), 0); 
		    sf_floatread(tt,((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz,Fden);	expand_high(tt, h_ro, fdm->nb, nx, fdm->nxpad, (nyinterior - nb) + (nyinterior * (ngpu-1)), (nyinterior * ngpu), nz, fdm->nzpad);

			/* stiffness */
		    h_c11 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float));
		    h_c22 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float)); 
		    h_c33 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float)); 
		    h_c44 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float)); 
		    h_c55 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float)); 
		    h_c66 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float)); 
		    h_c12 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float)); 
		    h_c13 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float)); 
		    h_c23 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float));
		
		    sf_floatread(tt,((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz,Fccc);     expand_high(tt, h_c11, fdm->nb, nx, fdm->nxpad, (nyinterior - nb) + (nyinterior * (ngpu-1)), (nyinterior * ngpu), nz, fdm->nzpad);
		    sf_seek(Fccc, (1*ny*nx*nz + (nyinterior * ngpu * (rank-1) + (nyinterior - nb) + (nyinterior * (ngpu-1))) * nx * nz)*sizeof(float), 0);
			sf_floatread(tt,((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz,Fccc);     expand_high(tt, h_c22, fdm->nb, nx, fdm->nxpad, (nyinterior - nb) + (nyinterior * (ngpu-1)), (nyinterior * ngpu), nz, fdm->nzpad);    
		    sf_seek(Fccc, (2*ny*nx*nz + (nyinterior * ngpu * (rank-1) + (nyinterior - nb) + (nyinterior * (ngpu-1))) * nx * nz)*sizeof(float), 0);
			sf_floatread(tt,((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz,Fccc);     expand_high(tt, h_c33, fdm->nb, nx, fdm->nxpad, (nyinterior - nb) + (nyinterior * (ngpu-1)), (nyinterior * ngpu), nz, fdm->nzpad);    
		    sf_seek(Fccc, (3*ny*nx*nz + (nyinterior * ngpu * (rank-1) + (nyinterior - nb) + (nyinterior * (ngpu-1))) * nx * nz)*sizeof(float), 0);
			sf_floatread(tt,((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz,Fccc);     expand_high(tt, h_c44, fdm->nb, nx, fdm->nxpad, (nyinterior - nb) + (nyinterior * (ngpu-1)), (nyinterior * ngpu), nz, fdm->nzpad);
		    sf_seek(Fccc, (4*ny*nx*nz + (nyinterior * ngpu * (rank-1) + (nyinterior - nb) + (nyinterior * (ngpu-1))) * nx * nz)*sizeof(float), 0);
			sf_floatread(tt,((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz,Fccc);     expand_high(tt, h_c55, fdm->nb, nx, fdm->nxpad, (nyinterior - nb) + (nyinterior * (ngpu-1)), (nyinterior * ngpu), nz, fdm->nzpad);    
		    sf_seek(Fccc, (5*ny*nx*nz + (nyinterior * ngpu * (rank-1) + (nyinterior - nb) + (nyinterior * (ngpu-1))) * nx * nz)*sizeof(float), 0);
			sf_floatread(tt,((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz,Fccc);     expand_high(tt, h_c66, fdm->nb, nx, fdm->nxpad, (nyinterior - nb) + (nyinterior * (ngpu-1)), (nyinterior * ngpu), nz, fdm->nzpad);
		    sf_seek(Fccc, (6*ny*nx*nz + (nyinterior * ngpu * (rank-1) + (nyinterior - nb) + (nyinterior * (ngpu-1))) * nx * nz)*sizeof(float), 0);
			sf_floatread(tt,((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz,Fccc);     expand_high(tt, h_c12, fdm->nb, nx, fdm->nxpad, (nyinterior - nb) + (nyinterior * (ngpu-1)), (nyinterior * ngpu), nz, fdm->nzpad);
		    sf_seek(Fccc, (7*ny*nx*nz + (nyinterior * ngpu * (rank-1) + (nyinterior - nb) + (nyinterior * (ngpu-1))) * nx * nz)*sizeof(float), 0);
			sf_floatread(tt,((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz,Fccc);     expand_high(tt, h_c13, fdm->nb, nx, fdm->nxpad, (nyinterior - nb) + (nyinterior * (ngpu-1)), (nyinterior * ngpu), nz, fdm->nzpad);
		    sf_seek(Fccc, (8*ny*nx*nz + (nyinterior * ngpu * (rank-1) + (nyinterior - nb) + (nyinterior * (ngpu-1))) * nx * nz)*sizeof(float), 0);
			sf_floatread(tt,((nyinterior - nb) + (nyinterior * (ngpu-1)))*nx*nz,Fccc);     expand_high(tt, h_c23, fdm->nb, nx, fdm->nxpad, (nyinterior - nb) + (nyinterior * (ngpu-1)), (nyinterior * ngpu), nz, fdm->nzpad);
			free(tt);
		}
		else {
			
			float *tt = (float*)malloc((nyinterior * ngpu)*nx*nz*sizeof(float));

		    /* input density */
		    h_ro = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float));

			sf_seek(Fden, ((nyinterior * ngpu * (rank-1) + (nyinterior - nb) + (nyinterior * (ngpu-1))) * nx * nz)*sizeof(float), 0); 
		    sf_floatread(tt,(nyinterior * ngpu)*nx*nz,Fden);	expand_interior(tt, h_ro, fdm->nb, nx, fdm->nxpad, (nyinterior * ngpu), (nyinterior * ngpu), nz, fdm->nzpad);

			/* stiffness */
		    h_c11 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float));
		    h_c22 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float)); 
		    h_c33 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float)); 
		    h_c44 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float)); 
		    h_c55 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float)); 
		    h_c66 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float)); 
		    h_c12 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float)); 
		    h_c13 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float)); 
		    h_c23 = (float*)malloc(fdm->nzpad * fdm->nxpad * (nyinterior * ngpu) * sizeof(float));
			
			sf_seek(Fccc, (0*ny*nx*nz + (nyinterior * ngpu * (rank-1) + (nyinterior - nb) + (nyinterior * (ngpu-1))) * nx * nz)*sizeof(float), 0);
		    sf_floatread(tt,nz*nx*(nyinterior * ngpu),Fccc);    expand_interior(tt,h_c11,fdm->nb, nx, fdm->nxpad, (nyinterior * ngpu), (nyinterior * ngpu), nz, fdm->nzpad);
		    sf_seek(Fccc, (1*ny*nx*nz + (nyinterior * ngpu * (rank-1) + (nyinterior - nb) + (nyinterior * (ngpu-1))) * nx * nz)*sizeof(float), 0);
			sf_floatread(tt,nz*nx*(nyinterior * ngpu),Fccc);    expand_interior(tt,h_c22,fdm->nb, nx, fdm->nxpad, (nyinterior * ngpu), (nyinterior * ngpu), nz, fdm->nzpad);    
		    sf_seek(Fccc, (2*ny*nx*nz + (nyinterior * ngpu * (rank-1) + (nyinterior - nb) + (nyinterior * (ngpu-1))) * nx * nz)*sizeof(float), 0);
			sf_floatread(tt,nz*nx*(nyinterior * ngpu),Fccc);    expand_interior(tt,h_c33,fdm->nb, nx, fdm->nxpad, (nyinterior * ngpu), (nyinterior * ngpu), nz, fdm->nzpad);    
		    sf_seek(Fccc, (3*ny*nx*nz + (nyinterior * ngpu * (rank-1) + (nyinterior - nb) + (nyinterior * (ngpu-1))) * nx * nz)*sizeof(float), 0);
			sf_floatread(tt,nz*nx*(nyinterior * ngpu),Fccc);    expand_interior(tt,h_c44,fdm->nb, nx, fdm->nxpad, (nyinterior * ngpu), (nyinterior * ngpu), nz, fdm->nzpad);
		    sf_seek(Fccc, (4*ny*nx*nz + (nyinterior * ngpu * (rank-1) + (nyinterior - nb) + (nyinterior * (ngpu-1))) * nx * nz)*sizeof(float), 0);
			sf_floatread(tt,nz*nx*(nyinterior * ngpu),Fccc);    expand_interior(tt,h_c55,fdm->nb, nx, fdm->nxpad, (nyinterior * ngpu), (nyinterior * ngpu), nz, fdm->nzpad);    
		    sf_seek(Fccc, (5*ny*nx*nz + (nyinterior * ngpu * (rank-1) + (nyinterior - nb) + (nyinterior * (ngpu-1))) * nx * nz)*sizeof(float), 0);
			sf_floatread(tt,nz*nx*(nyinterior * ngpu),Fccc);    expand_interior(tt,h_c66,fdm->nb, nx, fdm->nxpad, (nyinterior * ngpu), (nyinterior * ngpu), nz, fdm->nzpad);
		    sf_seek(Fccc, (6*ny*nx*nz + (nyinterior * ngpu * (rank-1) + (nyinterior - nb) + (nyinterior * (ngpu-1))) * nx * nz)*sizeof(float), 0);
			sf_floatread(tt,nz*nx*(nyinterior * ngpu),Fccc);    expand_interior(tt,h_c12,fdm->nb, nx, fdm->nxpad, (nyinterior * ngpu), (nyinterior * ngpu), nz, fdm->nzpad);
		    sf_seek(Fccc, (7*ny*nx*nz + (nyinterior * ngpu * (rank-1) + (nyinterior - nb) + (nyinterior * (ngpu-1))) * nx * nz)*sizeof(float), 0);
			sf_floatread(tt,nz*nx*(nyinterior * ngpu),Fccc);    expand_interior(tt,h_c13,fdm->nb, nx, fdm->nxpad, (nyinterior * ngpu), (nyinterior * ngpu), nz, fdm->nzpad);
		    sf_seek(Fccc, (8*ny*nx*nz + (nyinterior * ngpu * (rank-1) + (nyinterior - nb) + (nyinterior * (ngpu-1))) * nx * nz)*sizeof(float), 0);
			sf_floatread(tt,nz*nx*(nyinterior * ngpu),Fccc);    expand_interior(tt,h_c23,fdm->nb, nx, fdm->nxpad, (nyinterior * ngpu), (nyinterior * ngpu), nz, fdm->nzpad);
		    free(tt);
		}
	}
	
	// device arrays for local GPUs
	float **d_ro = (float**)malloc(ngpu*sizeof(float*));
	float **d_c11 = (float**)malloc(ngpu*sizeof(float*));
	float **d_c22 = (float**)malloc(ngpu*sizeof(float*));
	float **d_c33 = (float**)malloc(ngpu*sizeof(float*));
	float **d_c44 = (float**)malloc(ngpu*sizeof(float*));
	float **d_c55 = (float**)malloc(ngpu*sizeof(float*));
	float **d_c66 = (float**)malloc(ngpu*sizeof(float*));
	float **d_c12 = (float**)malloc(ngpu*sizeof(float*));
	float **d_c13 = (float**)malloc(ngpu*sizeof(float*));
	float **d_c23 = (float**)malloc(ngpu*sizeof(float*));

	// allocate density and stiffness sub-domain arrays on each GPU and copy the data
	for (int g = 0; g < ngpu; g++){
		cudaSetDevice(g);
		cudaMalloc(&d_ro[g] , nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMalloc(&d_c11[g], nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMalloc(&d_c22[g], nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMalloc(&d_c33[g], nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMalloc(&d_c44[g], nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMalloc(&d_c55[g], nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMalloc(&d_c66[g], nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMalloc(&d_c12[g], nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMalloc(&d_c13[g], nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMalloc(&d_c23[g], nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float));
		
		cudaMemcpy(d_ro[g] , h_ro  + g * nyinterior * fdm->nzpad * fdm->nxpad, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDefault);
		cudaMemcpy(d_c11[g], h_c11 + g * nyinterior * fdm->nzpad * fdm->nxpad, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDefault);
		cudaMemcpy(d_c22[g], h_c22 + g * nyinterior * fdm->nzpad * fdm->nxpad, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDefault);
		cudaMemcpy(d_c33[g], h_c33 + g * nyinterior * fdm->nzpad * fdm->nxpad, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDefault);
		cudaMemcpy(d_c44[g], h_c44 + g * nyinterior * fdm->nzpad * fdm->nxpad, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDefault);
		cudaMemcpy(d_c55[g], h_c55 + g * nyinterior * fdm->nzpad * fdm->nxpad, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDefault);
		cudaMemcpy(d_c66[g], h_c66 + g * nyinterior * fdm->nzpad * fdm->nxpad, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDefault);
		cudaMemcpy(d_c12[g], h_c12 + g * nyinterior * fdm->nzpad * fdm->nxpad, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDefault);
		cudaMemcpy(d_c13[g], h_c13 + g * nyinterior * fdm->nzpad * fdm->nxpad, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDefault);
		cudaMemcpy(d_c23[g], h_c23 + g * nyinterior * fdm->nzpad * fdm->nxpad, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDefault);	
	}
	
	
	
	/*------------------------------------------------------------*/
	/* Boundary condition setup */
	
	float spo;
	
	float **d_bzl_s = (float**)malloc(ngpu*sizeof(float*));
	float **d_bzh_s = (float**)malloc(ngpu*sizeof(float*));
	float **d_bxl_s = (float**)malloc(ngpu*sizeof(float*));
	float **d_bxh_s = (float**)malloc(ngpu*sizeof(float*));
	float **d_byl_s = (float**)malloc(ngpu*sizeof(float*));
	float **d_byh_s = (float**)malloc(ngpu*sizeof(float*));
	
	if (dabc){
		
		/* ABC set up */
		for (int g = 0; g < ngpu; g++){
			cudaSetDevice(g);
			cudaMalloc(&d_bzl_s[g], fdm->nxpad * nyinterior * sizeof(float));
			cudaMalloc(&d_bzh_s[g], fdm->nxpad * nyinterior * sizeof(float));
			cudaMalloc(&d_bxl_s[g], fdm->nzpad * nyinterior * sizeof(float));
			cudaMalloc(&d_bxh_s[g], fdm->nzpad * nyinterior * sizeof(float));
		}
		
		if (rank==0){
			cudaSetDevice(0);
			cudaMalloc(&d_byl_s[0], fdm->nzpad * fdm->nxpad * sizeof(float));
		}
		if (rank == nNode-1){
			cudaSetDevice(ngpu-1);
			cudaMalloc(&d_byh_s[ngpu-1], fdm->nzpad * fdm->nxpad * sizeof(float));
		}
		
		
		for (int g = 0; g < ngpu; g++){
			cudaSetDevice(g);
			dim3 dimGrid_init_XY(ceil(fdm->nxpad/32.0f), ceil(nyinterior/32.0f));
			dim3 dimBlock_init_XY(32,32);
			initABC_XY<<<dimGrid_init_XY, dimBlock_init_XY>>>(fdm->nxpad, fdm->nzpad, nyinterior, d_bzl_s[g], d_bzh_s[g], d_c55[g], d_ro[g], dt, fdm->dz);
			
			dim3 dimGrid_init_ZY(ceil(fdm->nzpad/32.0f), ceil(nyinterior/32.0f));
			dim3 dimBlock_init_ZY(32,32);
			initABC_ZY<<<dimGrid_init_ZY, dimBlock_init_ZY>>>(fdm->nxpad, fdm->nzpad, nyinterior, d_bxl_s[g], d_bxh_s[g], d_c55[g], d_ro[g], dt, fdm->dx);
		}
		
		if (rank==0){
			cudaSetDevice(0);
			dim3 dimGrid_init_XZ(ceil(fdm->nxpad/32.0f), ceil(fdm->nzpad/32.0f));
			dim3 dimBlock_init_XZ(32,32);
			initABC_XZ_low<<<dimGrid_init_XZ, dimBlock_init_XZ>>>(fdm->nxpad, fdm->nzpad, nyinterior, d_byl_s[0], d_c55[0], d_ro[0], dt, fdm->dy);
		}
		if (rank == nNode-1){
			cudaSetDevice(ngpu-1);
			dim3 dimGrid_init_XZ(ceil(fdm->nxpad/32.0f), ceil(fdm->nzpad/32.0f));
			dim3 dimBlock_init_XZ(32,32);
			initABC_XZ_high<<<dimGrid_init_XZ, dimBlock_init_XZ>>>(fdm->nxpad, fdm->nzpad, nyinterior, d_byh_s[ngpu-1], d_c55[ngpu-1], d_ro[ngpu-1], dt, fdm->dy);
		}
		
		
		/* sponge set up */
		// sponge coefficients are calculated inside the sponge kernel on GPU based on spo
		spo = (sqrt(2.0) * 4.0f * nb);
	}
	
	/*------------------------------------------------------------*/
	
	
	
	/*------------------------------------------------------------*/
    /* displacement: um = U @ t-1; uo = U @ t; up = U @ t+1 */
	float **d_umx = (float **)malloc(ngpu*sizeof(float*));
	float **d_uox = (float **)malloc(ngpu*sizeof(float*));
	float **d_upx = (float **)malloc(ngpu*sizeof(float*));
	float **d_uax = (float **)malloc(ngpu*sizeof(float*));
	float **d_utx = (float **)malloc(ngpu*sizeof(float*));

	float **d_umy = (float **)malloc(ngpu*sizeof(float*));
	float **d_uoy = (float **)malloc(ngpu*sizeof(float*));
	float **d_upy = (float **)malloc(ngpu*sizeof(float*));
	float **d_uay = (float **)malloc(ngpu*sizeof(float*));
	float **d_uty = (float **)malloc(ngpu*sizeof(float*));

	float **d_umz = (float **)malloc(ngpu*sizeof(float*));
	float **d_uoz = (float **)malloc(ngpu*sizeof(float*));
	float **d_upz = (float **)malloc(ngpu*sizeof(float*));
	float **d_uaz = (float **)malloc(ngpu*sizeof(float*));
	float **d_utz = (float **)malloc(ngpu*sizeof(float*));
	
	float **d_tzz = (float **)malloc(ngpu*sizeof(float*));
	float **d_txx = (float **)malloc(ngpu*sizeof(float*));
	float **d_tyy = (float **)malloc(ngpu*sizeof(float*));
	float **d_txy = (float **)malloc(ngpu*sizeof(float*));
	float **d_tyz = (float **)malloc(ngpu*sizeof(float*));
	float **d_tzx = (float **)malloc(ngpu*sizeof(float*));
	
	float *h_uTemp = (float *)malloc(nyinterior * fdm->nxpad * fdm->nzpad * sizeof(float));
	

	// allocate and initialize displacement, accel, and stress/strain arrasys to 0 on each GPU
	for (int g = 0; g < ngpu; g++){
		cudaSetDevice(g);
		
		cudaMalloc(&d_umx[g], nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMalloc(&d_uox[g], nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMalloc(&d_upx[g], nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));	
		cudaMalloc(&d_uax[g], nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		
		cudaMalloc(&d_umy[g], nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMalloc(&d_uoy[g], nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMalloc(&d_upy[g], nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));	
		cudaMalloc(&d_uay[g], nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		
		cudaMalloc(&d_umz[g], nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMalloc(&d_uoz[g], nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMalloc(&d_upz[g], nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));	
		cudaMalloc(&d_uaz[g], nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		
		cudaMalloc(&d_tzz[g], nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMalloc(&d_tyy[g], nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMalloc(&d_txx[g], nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMalloc(&d_txy[g], nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMalloc(&d_tyz[g], nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMalloc(&d_tzx[g], nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		
		sf_check_gpu_error(rank, "allocate grid arrays");
		
		cudaMemset(d_upz[g], 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMemset(d_uaz[g], 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		
		cudaMemset(d_upx[g], 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMemset(d_uax[g], 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		
		cudaMemset(d_upy[g], 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMemset(d_uay[g], 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		
		if (wavSrc){
			cudaMemset(d_umz[g], 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
			cudaMemset(d_uoz[g], 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
			
			cudaMemset(d_umx[g], 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
			cudaMemset(d_uox[g], 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));

			cudaMemset(d_umy[g], 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
			cudaMemset(d_uoy[g], 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		}
		else {
			// get umz interior chunk
			sf_seek(Fum, (rank * nyinterior * ngpu + g * nyinterior)*fdm->nxpad*fdm->nzpad*sizeof(float), 0);
			sf_floatreadtransp(h_uTemp, Fum, fdm->nxpad, nyinterior, fdm->nzpad);
			cudaMemcpy(d_umz[g] + 4 * fdm->nxpad * fdm->nzpad, h_uTemp, nyinterior * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			
			// get uoz interior chunk
			sf_seek(Fuo, (rank * nyinterior * ngpu + g * nyinterior)*fdm->nxpad*fdm->nzpad*sizeof(float), 0);
			sf_floatreadtransp(h_uTemp, Fuo, fdm->nxpad, nyinterior, fdm->nzpad);
			cudaMemcpy(d_uoz[g] + 4 * fdm->nxpad * fdm->nzpad, h_uTemp, nyinterior * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
						
			// get umx interior chunk
			sf_seek(Fum, ((fdm->nypad * fdm->nzpad * fdm->nxpad) + (rank * nyinterior * ngpu + g * nyinterior)*fdm->nxpad*fdm->nzpad)*sizeof(float), 0);
			sf_floatreadtransp(h_uTemp, Fum, fdm->nxpad, nyinterior, fdm->nzpad);
			cudaMemcpy(d_umx[g] + 4 * fdm->nxpad * fdm->nzpad, h_uTemp, nyinterior * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			
			// get uox interior chunk
			sf_seek(Fuo, ((fdm->nypad * fdm->nzpad * fdm->nxpad) + (rank * nyinterior * ngpu + g * nyinterior)*fdm->nxpad*fdm->nzpad)*sizeof(float), 0);
			sf_floatreadtransp(h_uTemp, Fuo, fdm->nxpad, nyinterior, fdm->nzpad);
			cudaMemcpy(d_uox[g] + 4 * fdm->nxpad * fdm->nzpad, h_uTemp, nyinterior * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
						
			// get umy interior chunk
			sf_seek(Fum, (2*(fdm->nypad * fdm->nzpad * fdm->nxpad) + (rank * nyinterior * ngpu + g * nyinterior)*fdm->nxpad*fdm->nzpad)*sizeof(float), 0);
			sf_floatreadtransp(h_uTemp, Fum, fdm->nxpad, nyinterior, fdm->nzpad);
			cudaMemcpy(d_umy[g] + 4 * fdm->nxpad * fdm->nzpad, h_uTemp, nyinterior * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			
			// get uoy interior chunk
			sf_seek(Fuo, (2*(fdm->nypad * fdm->nzpad * fdm->nxpad) + (rank * nyinterior * ngpu + g * nyinterior)*fdm->nxpad*fdm->nzpad)*sizeof(float), 0);
			sf_floatreadtransp(h_uTemp, Fuo, fdm->nxpad, nyinterior, fdm->nzpad);
			cudaMemcpy(d_uoy[g] + 4 * fdm->nxpad * fdm->nzpad, h_uTemp, nyinterior * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
						
			// get high halos
			if (rank != nNode-1 || g != ngpu-1){
				
				// umz
				sf_seek(Fum, (rank * nyinterior * ngpu + g * nyinterior + nyinterior)*fdm->nxpad*fdm->nzpad*sizeof(float), 0);
				sf_floatreadtransp(h_uTemp, Fum, fdm->nxpad, 4, fdm->nzpad);
				cudaMemcpy(d_umz[g] + (4 + nyinterior) * fdm->nxpad * fdm->nzpad, h_uTemp, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				
				// uoz
				sf_seek(Fuo, (rank * nyinterior * ngpu + g * nyinterior + nyinterior)*fdm->nxpad*fdm->nzpad*sizeof(float), 0);
				sf_floatreadtransp(h_uTemp, Fuo, fdm->nxpad, 4, fdm->nzpad);
				cudaMemcpy(d_uoz[g] + (4 + nyinterior) * fdm->nxpad * fdm->nzpad, h_uTemp, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				
				
				// umx 
				sf_seek(Fum, ((fdm->nypad * fdm->nzpad * fdm->nxpad) + (rank * nyinterior * ngpu + g * nyinterior + nyinterior)*fdm->nxpad*fdm->nzpad)*sizeof(float), 0);
				sf_floatreadtransp(h_uTemp, Fum, fdm->nxpad, 4, fdm->nzpad);
				cudaMemcpy(d_umx[g] + (4 + nyinterior) * fdm->nxpad * fdm->nzpad, h_uTemp, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
		
				// uox
				sf_seek(Fuo, ((fdm->nypad * fdm->nzpad * fdm->nxpad) + (rank * nyinterior * ngpu + g * nyinterior + nyinterior)*fdm->nxpad*fdm->nzpad)*sizeof(float), 0);
				sf_floatreadtransp(h_uTemp, Fuo, fdm->nxpad, 4, fdm->nzpad);
				cudaMemcpy(d_uox[g] + (4 + nyinterior) * fdm->nxpad * fdm->nzpad, h_uTemp, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				
				
				// umy
				sf_seek(Fum, (2*(fdm->nypad * fdm->nzpad * fdm->nxpad) + (rank * nyinterior * ngpu + g * nyinterior + nyinterior)*fdm->nxpad*fdm->nzpad)*sizeof(float), 0);
				sf_floatreadtransp(h_uTemp, Fum, fdm->nxpad, 4, fdm->nzpad);
				cudaMemcpy(d_umy[g] + (4 + nyinterior) * fdm->nxpad * fdm->nzpad, h_uTemp, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
		
				// uoy
				sf_seek(Fuo, (2*(fdm->nypad * fdm->nzpad * fdm->nxpad) + (rank * nyinterior * ngpu + g * nyinterior + nyinterior)*fdm->nxpad*fdm->nzpad)*sizeof(float), 0);
				sf_floatreadtransp(h_uTemp, Fuo, fdm->nxpad, 4, fdm->nzpad);
				cudaMemcpy(d_uoy[g] + (4 + nyinterior) * fdm->nxpad * fdm->nzpad, h_uTemp, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				
			}
			
			// get low halos
			if (rank != 0 || g != 0){
				
				// umz
				sf_seek(Fum, (rank * nyinterior * ngpu + g * nyinterior - 4)*fdm->nxpad*fdm->nzpad*sizeof(float), 0);
				sf_floatreadtransp(h_uTemp, Fum, fdm->nxpad, 4, fdm->nzpad);
				cudaMemcpy(d_umz[g], h_uTemp, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				
				// uoz
				sf_seek(Fuo, (rank * nyinterior * ngpu + g * nyinterior - 4)*fdm->nxpad*fdm->nzpad*sizeof(float), 0);
				sf_floatreadtransp(h_uTemp, Fuo, fdm->nxpad, 4, fdm->nzpad);
				cudaMemcpy(d_uoz[g], h_uTemp, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				
				
				// umx 
				sf_seek(Fum, ((fdm->nypad * fdm->nzpad * fdm->nxpad) + (rank * nyinterior * ngpu + g * nyinterior - 4)*fdm->nxpad*fdm->nzpad)*sizeof(float), 0);
				sf_floatreadtransp(h_uTemp, Fum, fdm->nxpad, 4, fdm->nzpad);
				cudaMemcpy(d_umx[g], h_uTemp, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
					
				// uox
				sf_seek(Fuo, ((fdm->nypad * fdm->nzpad * fdm->nxpad) + (rank * nyinterior * ngpu + g * nyinterior - 4)*fdm->nxpad*fdm->nzpad)*sizeof(float), 0);
				sf_floatreadtransp(h_uTemp, Fuo, fdm->nxpad, 4, fdm->nzpad);
				cudaMemcpy(d_uox[g], h_uTemp, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				
				
				// umy
				sf_seek(Fum, (2*(fdm->nypad * fdm->nzpad * fdm->nxpad) + (rank * nyinterior * ngpu + g * nyinterior - 4)*fdm->nxpad*fdm->nzpad)*sizeof(float), 0);
				sf_floatreadtransp(h_uTemp, Fum, fdm->nxpad, 4, fdm->nzpad);
				cudaMemcpy(d_umy[g], h_uTemp, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
					
				// uoy
				sf_seek(Fuo, (2*(fdm->nypad * fdm->nzpad * fdm->nxpad) + (rank * nyinterior * ngpu + g * nyinterior - 4)*fdm->nxpad*fdm->nzpad)*sizeof(float), 0);
				sf_floatreadtransp(h_uTemp, Fuo, fdm->nxpad, 4, fdm->nzpad);
				cudaMemcpy(d_uoy[g], h_uTemp, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				
			}
		}

		cudaMemset(d_tzz[g], 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMemset(d_tyy[g], 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMemset(d_txx[g], 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMemset(d_txy[g], 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMemset(d_tyz[g], 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		cudaMemset(d_tzx[g], 0, nylocal * fdm->nzpad * fdm->nxpad * sizeof(float));
		
		sf_check_gpu_error(rank, "initialize grid arrays");
	}
	
	// Used for exchanging halo regions between remote GPUs
	// TODO: all of these probably aren't really needed, could probably re-use some, not sure tho...
	float *h_tzz_l_send, *h_tzz_l_recv, *h_tzz_h_send, *h_tzz_h_recv;
	float *h_tyy_l_send, *h_tyy_l_recv, *h_tyy_h_send, *h_tyy_h_recv;
	float *h_txx_l_send, *h_txx_l_recv, *h_txx_h_send, *h_txx_h_recv;
	float *h_txy_l_send, *h_txy_l_recv, *h_txy_h_send, *h_txy_h_recv;
	float *h_tyz_l_send, *h_tyz_l_recv, *h_tyz_h_send, *h_tyz_h_recv;
	float *h_tzx_l_send, *h_tzx_l_recv, *h_tzx_h_send, *h_tzx_h_recv;
	
	float *h_upz_l_send, *h_upz_l_recv, *h_upz_h_send, *h_upz_h_recv;
	float *h_upy_l_send, *h_upy_l_recv, *h_upy_h_send, *h_upy_h_recv;
	float *h_upx_l_send, *h_upx_l_recv, *h_upx_h_send, *h_upx_h_recv;
	
	if (nNode > 1){
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

		h_upz_l_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_upz_l_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_upz_h_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_upz_h_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));

		h_upy_l_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_upy_l_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_upy_h_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_upy_h_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));

		h_upx_l_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_upx_l_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_upx_h_send = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
		h_upx_h_recv = (float*)malloc(4 * fdm->nxpad * fdm->nzpad * sizeof(float));
	}
	
	
	/*------------------------------------------------------------*/
    /* precompute 1/ro * dt^2 									  */	
	/*------------------------------------------------------------*/
	for (int g = 0; g < ngpu; g++){
		cudaSetDevice(g);
		dim3 dimGrid1(ceil(fdm->nxpad/8.0f),ceil(fdm->nzpad/8.0f),ceil(nyinterior/8.0f));
		dim3 dimBlock1(8,8,8);
		computeRo<<<dimGrid1, dimBlock1>>>(d_ro[g], dt, fdm->nxpad, fdm->nzpad, nyinterior);
	}
	sf_check_gpu_error(rank, "computeRo Kernel");
		
	int first;
	int last;
	
	/*------------------------------------------------------------*/
    /* 
     *  MAIN LOOP
     */
    /*------------------------------------------------------------*/
    if(verb && rank==0) fprintf(stderr,"\n");
    for (it=0; it<nt; it++) {
		if(verb && rank==0) fprintf(stderr,"\b\b\b\b\b%d",it);
		
		/*------------------------------------------------------------*/
		/* from displacement to strain     OUTSIDE                    */
		/*------------------------------------------------------------*/
		for (int g = 0; g < ngpu; g++){
			cudaSetDevice(g);
			dim3 dimGrid2((fdm->nxpad-2*NOP)/24.0f, (fdm->nzpad-2*NOP)/24.0f);
			dim3 dimBlock2(24,24,1);
			if (rank == 0 && g == 0){
				first = 1;
				last = 0;
			}
			else if (rank == nNode-1 && g == ngpu-1){
				first = 0;
				last = 1;
			}
			else {
				first = 0;
				last = 0;
			}
			
			// NOTE: These all have to be seperate kernels because nvcc does something REALLY stupid with register usage when you add if-statements......
			if (!first){	// do the low region
				dispToStrainOutside<<<dimGrid2, dimBlock2, 32*32*3*sizeof(float)>>>(0, fdm->nxpad, nylocal, fdm->nzpad, d_uox[g], d_uoy[g], d_uoz[g], d_txx[g], d_tyy[g], d_tzz[g], d_txy[g], d_tyz[g], d_tzx[g], idx, idy, idz);
			}
			if (!last){		// do the high region
				dispToStrainOutside<<<dimGrid2, dimBlock2, 32*32*3*sizeof(float)>>>(1, fdm->nxpad, nylocal, fdm->nzpad, d_uox[g], d_uoy[g], d_uoz[g], d_txx[g], d_tyy[g], d_tzz[g], d_txy[g], d_tyz[g], d_tzx[g], idx, idy, idz);			
			}
		}
		
		sf_check_gpu_error(rank, "dispToStrain Kernel");
		
		
		/*------------------------------------------------------------*/
		/* from strain to stress        OUTSIDE                       */
		/*------------------------------------------------------------*/
		for (int g = 0; g < ngpu; g++){
			cudaSetDevice(g);
			dim3 dimGrid3_outside(ceil(fdm->nxpad/192.0f), fdm->nzpad, 8);
			dim3 dimBlock3_outside(192,1,1);
			strainToStressOutside<<<dimGrid3_outside, dimBlock3_outside>>>(fdm->nxpad, fdm->nzpad, nyinterior, d_c11[g], d_c12[g], d_c13[g], d_c22[g], d_c23[g], d_c33[g], d_c44[g], d_c55[g], d_c66[g], d_txx[g], d_tyy[g], d_tzz[g], d_txy[g], d_tyz[g], d_tzx[g]);
		}
		sf_check_gpu_error(rank, "strainToStress Kernel");
		
			
		/*------------------------------------------------------------*/
		/* free surface       OUTSIDE                                 */
		/*------------------------------------------------------------*/
		if(fsrf) {
			for (int g = 0; g < ngpu; g++){
				cudaSetDevice(g);
				dim3 dimGrid4(ceil(fdm->nxpad/8.0f), ceil(fdm->nb/8.0f), 1);
				dim3 dimBlock4(8,8,8);
				freeSurfOutside<<<dimGrid4, dimBlock4>>>(fdm->nxpad, nyinterior, fdm->nzpad, fdm->nb, d_tzz[g], d_tyz[g], d_tzx[g]);
			}
			sf_check_gpu_error(rank, "freeSurf Kernel");
		}
		
		
		/*------------------------------------------------------------*/
		/* inject stress source           OUTSIDE                     */
		/*------------------------------------------------------------*/
		if(ssou && wavSrc) {
			for (int g = 0; g < ngpu; g++){
				cudaSetDevice(g);
				dim3 dimGrid5(ns, 1, 1);
				dim3 dimBlock5(2 * nbell + 1, 2 * nbell + 1, 1);
				lint3d_bell_gpuOutside<<<dimGrid5, dimBlock5>>>(rank, ngpu, g, it, nc, ns, 0, nbell, fdm->nxpad, nyinterior, fdm->nzpad, d_tzz[g], d_bell[g], d_Sjx[g], d_Sjz[g], d_Sjy[g], d_ww[g], d_Sw000[g], d_Sw001[g], d_Sw010[g], d_Sw011[g], d_Sw100[g], d_Sw101[g], d_Sw110[g], d_Sw111[g]);
				lint3d_bell_gpuOutside<<<dimGrid5, dimBlock5>>>(rank, ngpu, g, it, nc, ns, 1, nbell, fdm->nxpad, nyinterior, fdm->nzpad, d_txx[g], d_bell[g], d_Sjx[g], d_Sjz[g], d_Sjy[g], d_ww[g], d_Sw000[g], d_Sw001[g], d_Sw010[g], d_Sw011[g], d_Sw100[g], d_Sw101[g], d_Sw110[g], d_Sw111[g]);
				lint3d_bell_gpuOutside<<<dimGrid5, dimBlock5>>>(rank, ngpu, g, it, nc, ns, 2, nbell, fdm->nxpad, nyinterior, fdm->nzpad, d_tyy[g], d_bell[g], d_Sjx[g], d_Sjz[g], d_Sjy[g], d_ww[g], d_Sw000[g], d_Sw001[g], d_Sw010[g], d_Sw011[g], d_Sw100[g], d_Sw101[g], d_Sw110[g], d_Sw111[g]);					
			}
			sf_check_gpu_error(rank, "lint3d_bell_gpu Kernel");
		}
		
		/*------------------------------------------------------------*/
		/* send halo regions of d_t arrays to remote GPU	          */
		/*------------------------------------------------------------*/
		if (nNode > 1){	// using multiple nodes, must exchange halo regions between neighboring exterior GPUs
			if (rank == 0){
				// get high halo region from d_t array on GPU ngpu-1 and send to rank+1
				cudaMemcpy(h_tzz_h_send, d_tzz[ngpu-1] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_tzz_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tyy_h_send, d_tyy[ngpu-1] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_tyy_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_txx_h_send, d_txx[ngpu-1] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_txx_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_txy_h_send, d_txy[ngpu-1] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_txy_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 3, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tyz_h_send, d_tyz[ngpu-1] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_tyz_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 4, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tzx_h_send, d_tzx[ngpu-1] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_tzx_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 5, MPI_COMM_WORLD, &request);

				sf_check_gpu_error(rank, "send get high");

			}
			else if (rank == nNode-1){
				// get low halo region from d_t arrays on GPU 0 and send to rank-1
				cudaMemcpy(h_tzz_l_send, d_tzz[0] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_tzz_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tyy_l_send, d_tyy[0] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_tyy_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_txx_l_send, d_txx[0] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_txx_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_txy_l_send, d_txy[0] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_txy_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 3, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tyz_l_send, d_tyz[0] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_tyz_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 4, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tzx_l_send, d_tzx[0] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_tzx_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 5, MPI_COMM_WORLD, &request);

				sf_check_gpu_error(rank, "send get low");
				
			}
			else {
				// get low halo region from d_t arrays on GPU 0 and send to rank-1
				cudaMemcpy(h_tzz_l_send, d_tzz[0] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_tzz_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tyy_l_send, d_tyy[0] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_tyy_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_txx_l_send, d_txx[0] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_txx_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_txy_l_send, d_txy[0] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_txy_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 3, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tyz_l_send, d_tyz[0] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_tyz_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 4, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tzx_l_send, d_tzx[0] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_tzx_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 5, MPI_COMM_WORLD, &request);

				// get high halo region from d_t arrays on GPU ngpu-1 and send to rank+1
				cudaMemcpy(h_tzz_h_send, d_tzz[ngpu-1] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_tzz_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tyy_h_send, d_tyy[ngpu-1] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_tyy_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_txx_h_send, d_txx[ngpu-1] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_txx_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_txy_h_send, d_txy[ngpu-1] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_txy_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 3, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tyz_h_send, d_tyz[ngpu-1] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_tyz_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 4, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_tzx_h_send, d_tzx[ngpu-1] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_tzx_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 5, MPI_COMM_WORLD, &request);

			}

			sf_check_gpu_error(rank, "send remote halos");
		}
		
		/*------------------------------------------------------------*/
		/* start recv'ing halo regions of d_t arrays from remote GPU  */
		/*------------------------------------------------------------*/
		if (nNode > 1){
			if (rank == 0){
				
				// receive low halo region of d_t arrays from rank+1 and copy into GPU ngpu-1 high halo
				MPI_Irecv(h_tzz_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &recvRequests[0]);
				MPI_Irecv(h_tyy_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &recvRequests[1]);
				MPI_Irecv(h_txx_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &recvRequests[2]);
				MPI_Irecv(h_txy_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 3, MPI_COMM_WORLD, &recvRequests[3]);
				MPI_Irecv(h_tyz_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 4, MPI_COMM_WORLD, &recvRequests[4]);
				MPI_Irecv(h_tzx_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 5, MPI_COMM_WORLD, &recvRequests[5]);
				
			}
			else if (rank == nNode - 1){
				
				// receive high halo region of d_t arrays from rank-1 and copy into GPU 0 low halo
				MPI_Irecv(h_tzz_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &recvRequests[0]);
				MPI_Irecv(h_tyy_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &recvRequests[1]);
				MPI_Irecv(h_txx_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, &recvRequests[2]);
				MPI_Irecv(h_txy_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 3, MPI_COMM_WORLD, &recvRequests[3]);
				MPI_Irecv(h_tyz_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 4, MPI_COMM_WORLD, &recvRequests[4]);
				MPI_Irecv(h_tzx_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 5, MPI_COMM_WORLD, &recvRequests[5]);
				
			}
			
			else {
				
				// receive high halo region of d_t arrays from from rank-1 and copy to GPU 0
				MPI_Irecv(h_tzz_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &recvRequests[0]);
				MPI_Irecv(h_tyy_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &recvRequests[1]);
				MPI_Irecv(h_txx_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, &recvRequests[2]);
				MPI_Irecv(h_txy_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 3, MPI_COMM_WORLD, &recvRequests[3]);
				MPI_Irecv(h_tyz_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 4, MPI_COMM_WORLD, &recvRequests[4]);
				MPI_Irecv(h_tzx_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 5, MPI_COMM_WORLD, &recvRequests[5]);

				// receive low halo region of d_t arrays from from rank+1 and copy to GPU ngpu-1
				MPI_Irecv(h_tzz_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &recvRequests[6]);
				MPI_Irecv(h_tyy_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &recvRequests[7]);
				MPI_Irecv(h_txx_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &recvRequests[8]);
				MPI_Irecv(h_txy_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 3, MPI_COMM_WORLD, &recvRequests[9]);
				MPI_Irecv(h_tyz_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 4, MPI_COMM_WORLD, &recvRequests[10]);
				MPI_Irecv(h_tzx_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 5, MPI_COMM_WORLD, &recvRequests[11]);
								
			}
			
		}
		
		
		/*------------------------------------------------------------*/
		/* exchange halo regions of d_t arrays between local GPUs     */
		/*------------------------------------------------------------*/		
		if (ngpu > 1){ // using multiple GPUs, must exchange halo regions between neighboring GPUs
			// high halo region of d_t arrays on GPU 0 to GPU 1
			cudaMemcpy(d_tzz[1], d_tzz[0] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			cudaMemcpy(d_tyy[1], d_tyy[0] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			cudaMemcpy(d_txx[1], d_txx[0] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			cudaMemcpy(d_txy[1], d_txy[0] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			cudaMemcpy(d_tyz[1], d_tyz[0] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			cudaMemcpy(d_tzx[1], d_tzx[0] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			
			// exchange halo regions of d_t arrays between all internal GPUs
			for (int g = 1; g < ngpu-1; g++){
				// high halo region of GPU g to low halo region of GPU g+1
				cudaMemcpy(d_tzz[g+1], d_tzz[g] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_tyy[g+1], d_tyy[g] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_txx[g+1], d_txx[g] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_txy[g+1], d_txy[g] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_tyz[g+1], d_tyz[g] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_tzx[g+1], d_tzx[g] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				
				// low halo region of GPU g to high halo region of GPU g-1
				cudaMemcpy(d_tzz[g-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), d_tzz[g] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_tyy[g-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), d_tyy[g] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_txx[g-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), d_txx[g] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_txy[g-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), d_txy[g] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_tyz[g-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), d_tyz[g] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_tzx[g-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), d_tzx[g] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			}
			
			// low halo region of d_t arrays on GPU (ngpu-1) to GPU (ngpu-2)
			cudaMemcpy(d_tzz[ngpu-2] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), d_tzz[ngpu-1] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			cudaMemcpy(d_tyy[ngpu-2] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), d_tyy[ngpu-1] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			cudaMemcpy(d_txx[ngpu-2] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), d_txx[ngpu-1] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			cudaMemcpy(d_txy[ngpu-2] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), d_txy[ngpu-1] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			cudaMemcpy(d_tyz[ngpu-2] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), d_tyz[ngpu-1] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			cudaMemcpy(d_tzx[ngpu-2] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), d_tzx[ngpu-1] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
		
			sf_check_gpu_error(rank, "exchange local halos");
		}
		
		
		/*------------------------------------------------------------*/
		/* from displacement to strain       INSIDE                   */
		/*------------------------------------------------------------*/
		for (int g = 0; g < ngpu; g++){
			cudaSetDevice(g);
			dim3 dimGrid2((fdm->nxpad-2*NOP)/24.0f, (fdm->nzpad-2*NOP)/24.0f);
			dim3 dimBlock2(24,24,1);
			dispToStrainInside<<<dimGrid2, dimBlock2, 32*32*3*sizeof(float)>>>(fdm->nxpad, nylocal, fdm->nzpad, d_uox[g], d_uoy[g], d_uoz[g], d_txx[g], d_tyy[g], d_tzz[g], d_txy[g], d_tyz[g], d_tzx[g], idx, idy, idz);						
		}
		sf_check_gpu_error(rank, "dispToStrain Kernel");
		
		
		/*------------------------------------------------------------*/
		/* from strain to stress     INSIDE                           */
		/*------------------------------------------------------------*/
		for (int g = 0; g < ngpu; g++){
			cudaSetDevice(g);
			dim3 dimGrid3_inside(ceil(fdm->nxpad/192.0f), fdm->nzpad, nyinterior-8);
			dim3 dimBlock3_inside(192,1,1);
			strainToStressInside<<<dimGrid3_inside, dimBlock3_inside>>>(fdm->nxpad, fdm->nzpad, nyinterior, d_c11[g], d_c12[g], d_c13[g], d_c22[g], d_c23[g], d_c33[g], d_c44[g], d_c55[g], d_c66[g], d_txx[g], d_tyy[g], d_tzz[g], d_txy[g], d_tyz[g], d_tzx[g]);
		}
		sf_check_gpu_error(rank, "strainToStress Kernel");
		
			
		/*------------------------------------------------------------*/
		/* free surface         INSIDE                                */
		/*------------------------------------------------------------*/
		if(fsrf) {
			for (int g = 0; g < ngpu; g++){
				cudaSetDevice(g);
				dim3 dimGrid4(ceil(fdm->nxpad/8.0f), ceil(fdm->nb/8.0f), ceil(nyinterior/8.0f));
				dim3 dimBlock4(8,8,8);
				freeSurfInside<<<dimGrid4, dimBlock4>>>(fdm->nxpad, nyinterior, fdm->nzpad, fdm->nb, d_tzz[g], d_tyz[g], d_tzx[g]);
				
			}
			sf_check_gpu_error(rank, "freeSurf Kernel");
		}
		
		
		/*------------------------------------------------------------*/
		/* inject stress source     INSIDE                            */
		/*------------------------------------------------------------*/
		if(ssou && wavSrc) {
			for (int g = 0; g < ngpu; g++){
				cudaSetDevice(g);
				dim3 dimGrid5(ns, 1, 1);
				dim3 dimBlock5(2 * nbell + 1, 2 * nbell + 1, 1);
				lint3d_bell_gpuInside<<<dimGrid5, dimBlock5>>>(rank, ngpu, g, it, nc, ns, 0, nbell, fdm->nxpad, nyinterior, fdm->nzpad, d_tzz[g], d_bell[g], d_Sjx[g], d_Sjz[g], d_Sjy[g], d_ww[g], d_Sw000[g], d_Sw001[g], d_Sw010[g], d_Sw011[g], d_Sw100[g], d_Sw101[g], d_Sw110[g], d_Sw111[g]);
				lint3d_bell_gpuInside<<<dimGrid5, dimBlock5>>>(rank, ngpu, g, it, nc, ns, 1, nbell, fdm->nxpad, nyinterior, fdm->nzpad, d_txx[g], d_bell[g], d_Sjx[g], d_Sjz[g], d_Sjy[g], d_ww[g], d_Sw000[g], d_Sw001[g], d_Sw010[g], d_Sw011[g], d_Sw100[g], d_Sw101[g], d_Sw110[g], d_Sw111[g]);
				lint3d_bell_gpuInside<<<dimGrid5, dimBlock5>>>(rank, ngpu, g, it, nc, ns, 2, nbell, fdm->nxpad, nyinterior, fdm->nzpad, d_tyy[g], d_bell[g], d_Sjx[g], d_Sjz[g], d_Sjy[g], d_ww[g], d_Sw000[g], d_Sw001[g], d_Sw010[g], d_Sw011[g], d_Sw100[g], d_Sw101[g], d_Sw110[g], d_Sw111[g]);	
			}
			sf_check_gpu_error(rank, "lint3d_bell_gpu Kernel");
		}
		
		
		/*------------------------------------------------------------*/
		/* copy halos from d_t arrays from remote GPU into local GPU  */
		/*------------------------------------------------------------*/
		// Wait for all MPI transfers to complete
		if (nNode > 1){
			if (rank == 0 || rank == nNode-1){
				MPI_Waitall(6, recvRequests, recvStatuses);
			}
			else {
				MPI_Waitall(12, recvRequests, recvStatuses);
			}
		}
		if (nNode > 1){
			if (rank == 0){
				
				cudaMemcpy(d_tzz[ngpu-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), h_tzz_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_tyy[ngpu-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), h_tyy_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_txx[ngpu-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), h_txx_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_txy[ngpu-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), h_txy_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_tyz[ngpu-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), h_tyz_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_tzx[ngpu-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), h_tzx_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);

				sf_check_gpu_error(rank, "exchange insert low");
			}
			else if (rank == nNode - 1){
				
				cudaMemcpy(d_tzz[0], h_tzz_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_tyy[0], h_tyy_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_txx[0], h_txx_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_txy[0], h_txy_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_tyz[0], h_tyz_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_tzx[0], h_tzx_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				
				sf_check_gpu_error(rank, "exchange insert high");
			}
			else {
				
				cudaMemcpy(d_tzz[0], h_tzz_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_tyy[0], h_tyy_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_txx[0], h_txx_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_txy[0], h_txy_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_tyz[0], h_tyz_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_tzx[0], h_tzx_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				
				cudaMemcpy(d_tzz[ngpu-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), h_tzz_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_tyy[ngpu-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), h_tyy_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_txx[ngpu-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), h_txx_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_txy[ngpu-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), h_txy_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_tyz[ngpu-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), h_tyz_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_tzx[ngpu-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), h_tzx_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				
			}
		}
		

		
		/*------------------------------------------------------------*/
		/* from stress to acceleration               OUTSIDE          */
		/*------------------------------------------------------------*/
		for (int g = 0; g < ngpu; g++){
			cudaSetDevice(g);
			dim3 dimGrid6((fdm->nxpad-2*NOP)/24.0f, (fdm->nzpad-2*NOP)/24.0f);
			dim3 dimBlock6(24,24,1);
			if (rank == 0 && g == 0){
				first = 1;
				last = 0;
			}
			else if (rank == nNode-1 && g == ngpu-1){
				first = 0;
				last = 1;
			}
			else {
				first = 0;
				last = 0;
			}
			if (!first){	// do the low region
				stressToAccelOutside<<<dimGrid6, dimBlock6, 32*32*5*sizeof(float)>>>(0, fdm->nxpad, fdm->nzpad, nylocal, idx, idy, idz, d_txx[g], d_tyy[g], d_tzz[g], d_txy[g], d_tzx[g], d_tyz[g], d_uax[g], d_uay[g], d_uaz[g]);
			}
			if (!last){
				stressToAccelOutside<<<dimGrid6, dimBlock6, 32*32*5*sizeof(float)>>>(1, fdm->nxpad, fdm->nzpad, nylocal, idx, idy, idz, d_txx[g], d_tyy[g], d_tzz[g], d_txy[g], d_tzx[g], d_tyz[g], d_uax[g], d_uay[g], d_uaz[g]);
			}
		}
		sf_check_gpu_error(rank, "stressToAccel Kernel");
		
		
		/*------------------------------------------------------------*/
		/* inject acceleration source      OUTSIDE                    */
		/*------------------------------------------------------------*/
		if(!ssou && wavSrc) {
			for (int g = 0; g < ngpu; g++){
				cudaSetDevice(g);
				dim3 dimGrid7(ns, 1, 1);
				dim3 dimBlock7(2 * nbell + 1, 2 * nbell + 1, 1);
				lint3d_bell_gpuOutside<<<dimGrid7, dimBlock7>>>(rank, ngpu, g, it, nc, ns, 0, nbell, fdm->nxpad, nyinterior, fdm->nzpad, d_uaz[g], d_bell[g], d_Sjx[g], d_Sjz[g], d_Sjy[g], d_ww[g], d_Sw000[g], d_Sw001[g], d_Sw010[g], d_Sw011[g], d_Sw100[g], d_Sw101[g], d_Sw110[g], d_Sw111[g]);
				lint3d_bell_gpuOutside<<<dimGrid7, dimBlock7>>>(rank, ngpu, g, it, nc, ns, 1, nbell, fdm->nxpad, nyinterior, fdm->nzpad, d_uax[g], d_bell[g], d_Sjx[g], d_Sjz[g], d_Sjy[g], d_ww[g], d_Sw000[g], d_Sw001[g], d_Sw010[g], d_Sw011[g], d_Sw100[g], d_Sw101[g], d_Sw110[g], d_Sw111[g]);
				lint3d_bell_gpuOutside<<<dimGrid7, dimBlock7>>>(rank, ngpu, g, it, nc, ns, 2, nbell, fdm->nxpad, nyinterior, fdm->nzpad, d_uay[g], d_bell[g], d_Sjx[g], d_Sjz[g], d_Sjy[g], d_ww[g], d_Sw000[g], d_Sw001[g], d_Sw010[g], d_Sw011[g], d_Sw100[g], d_Sw101[g], d_Sw110[g], d_Sw111[g]);	
			}
		
			sf_check_gpu_error(rank, "lint3d_bell_gpu Kernel");
		}
		
		
		
		/*------------------------------------------------------------*/
		/* step forward in time        OUTSIDE                        */
		/*------------------------------------------------------------*/
		for (int g = 0; g < ngpu; g++){
			cudaSetDevice(g);
			dim3 dimGrid8_outside(ceil(fdm->nxpad/192.0f), fdm->nzpad, 8);
			dim3 dimBlock8_outside(192,1,1);
			stepTimeOutside<<<dimGrid8_outside, dimBlock8_outside>>>(fdm->nxpad, nyinterior, fdm->nzpad, d_ro[g], d_uox[g], d_umx[g], d_uax[g], d_upx[g], d_uoy[g], d_umy[g], d_uay[g], d_upy[g], d_uoz[g], d_umz[g], d_uaz[g], d_upz[g]);
		}
		sf_check_gpu_error(rank, "stepTime Kernel");
		

		/*------------------------------------------------------------*/
		/* apply boundary conditions          OUTSIDE                 */
		/*------------------------------------------------------------*/
		if(dabc){
			
			/* One-way Absorbing BC */
			for (int g = 0; g < ngpu; g++){
				cudaSetDevice(g);
				
				dim3 dimGrid_abc_XY(ceil(fdm->nxpad/32.0f),2,2);
				dim3 dimBlock_abc_XY(32,4,1);
				abcone3d_apply_XY_outside<<<dimGrid_abc_XY,dimBlock_abc_XY>>>(fdm->nxpad, nyinterior, fdm->nzpad, d_upx[g], d_uox[g], d_bzl_s[g], d_bzh_s[g]);
				abcone3d_apply_XY_outside<<<dimGrid_abc_XY,dimBlock_abc_XY>>>(fdm->nxpad, nyinterior, fdm->nzpad, d_upy[g], d_uoy[g], d_bzl_s[g], d_bzh_s[g]);
				abcone3d_apply_XY_outside<<<dimGrid_abc_XY,dimBlock_abc_XY>>>(fdm->nxpad, nyinterior, fdm->nzpad, d_upz[g], d_uoz[g], d_bzl_s[g], d_bzh_s[g]);
				
				dim3 dimGrid_abc_ZY(2, 2, ceil(fdm->nzpad/32.0f));
				dim3 dimBlock_abc_ZY(1,4,32);
				abcone3d_apply_ZY_outside<<<dimGrid_abc_ZY,dimBlock_abc_ZY>>>(fdm->nxpad, nyinterior, fdm->nzpad, d_upx[g], d_uox[g], d_bxl_s[g], d_bxh_s[g]);
				abcone3d_apply_ZY_outside<<<dimGrid_abc_ZY,dimBlock_abc_ZY>>>(fdm->nxpad, nyinterior, fdm->nzpad, d_upy[g], d_uoy[g], d_bxl_s[g], d_bxh_s[g]);
				abcone3d_apply_ZY_outside<<<dimGrid_abc_ZY,dimBlock_abc_ZY>>>(fdm->nxpad, nyinterior, fdm->nzpad, d_upz[g], d_uoz[g], d_bxl_s[g], d_bxh_s[g]);
			}
			
			
			/* sponge BC */
			for (int g = 0; g < ngpu; g++){
				cudaSetDevice(g);
				dim3 dimGrid_spng_XY(ceil(fdm->nxpad/192.0f),4,2);
				dim3 dimBlock_spng_XY(192,1,1);                                            
				sponge3d_apply_XY_outside<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(d_umz[g], fdm->nxpad, nyinterior, fdm->nzpad, nb);
				sponge3d_apply_XY_outside<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(d_uoz[g], fdm->nxpad, nyinterior, fdm->nzpad, nb);
				sponge3d_apply_XY_outside<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(d_upz[g], fdm->nxpad, nyinterior, fdm->nzpad, nb);
			                                                                    
				sponge3d_apply_XY_outside<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(d_umx[g], fdm->nxpad, nyinterior, fdm->nzpad, nb);
				sponge3d_apply_XY_outside<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(d_uox[g], fdm->nxpad, nyinterior, fdm->nzpad, nb);
				sponge3d_apply_XY_outside<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(d_upx[g], fdm->nxpad, nyinterior, fdm->nzpad, nb);
			                                                                    
				sponge3d_apply_XY_outside<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(d_umy[g], fdm->nxpad, nyinterior, fdm->nzpad, nb);
				sponge3d_apply_XY_outside<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(d_uoy[g], fdm->nxpad, nyinterior, fdm->nzpad, nb);
				sponge3d_apply_XY_outside<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(d_upy[g], fdm->nxpad, nyinterior, fdm->nzpad, nb);
			
			
				dim3 dimGrid_spng_ZY(ceil(nb/8.0f),ceil(fdm->nzpad/8.0f),2);
				dim3 dimBlock_spng_ZY(8,8,4);
				sponge3d_apply_ZY_outside<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(d_umz[g], fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
				sponge3d_apply_ZY_outside<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(d_uoz[g], fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
				sponge3d_apply_ZY_outside<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(d_upz[g], fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
							                                                                    
				sponge3d_apply_ZY_outside<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(d_umx[g], fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
				sponge3d_apply_ZY_outside<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(d_uox[g], fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
				sponge3d_apply_ZY_outside<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(d_upx[g], fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
							                                                                    
				sponge3d_apply_ZY_outside<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(d_umy[g], fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
				sponge3d_apply_ZY_outside<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(d_uoy[g], fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
				sponge3d_apply_ZY_outside<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(d_upy[g], fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
			}
			
			sf_check_gpu_error(rank, "Boundary Condition Kernels");
			
		}

		

		/*------------------------------------------------------------*/
		/* start sending halo regions of d_up arrays to remote GPUs   */
		/*------------------------------------------------------------*/
		if (nNode > 1){	// using multiple nodes, must exchange halo regions between neighboring nodes GPUs
			if (rank == 0){
				// get high halo region from d_t array on GPU ngpu-1 and send to rank+1
				cudaMemcpy(h_upx_h_send, d_upx[ngpu-1] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_upx_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_upy_h_send, d_upy[ngpu-1] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_upy_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_upz_h_send, d_upz[ngpu-1] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_upz_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &request);
				
			}
			else if (rank == nNode-1){
				// get low halo region from d_t arrays on GPU 0 and send to rank-1
				cudaMemcpy(h_upx_l_send, d_upx[0] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_upx_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_upy_l_send, d_upy[0] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_upy_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_upz_l_send, d_upz[0] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_upz_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, &request);
				
			}
			else {
				// get low halo region from d_t arrays on GPU 0 and send to rank-1
				cudaMemcpy(h_upx_l_send, d_upx[0] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_upx_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_upy_l_send, d_upy[0] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_upy_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_upz_l_send, d_upz[0] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_upz_l_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, &request);
				
				// get high halo region from d_t arrays on GPU ngpu-1 and send to rank+1
				cudaMemcpy(h_upx_h_send, d_upx[ngpu-1] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_upx_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_upy_h_send, d_upy[ngpu-1] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_upy_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &request);
				cudaMemcpy(h_upz_h_send, d_upz[ngpu-1] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				MPI_Isend(h_upz_h_send, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &request);
			}
		}
		
		
		/*------------------------------------------------------------*/
		/* start recv'ing halo regions of d_up arrays from remote GPU */
		/*------------------------------------------------------------*/
		if (nNode > 1){	// using multiple nodes, must exchange halo regions between neighboring nodes GPUs
			if (rank == 0){
				// receive low halo region of d_t arrays from from rank+1 and copy to GPU ngpu-1
				MPI_Irecv(h_upx_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &recvRequests[0]);
				MPI_Irecv(h_upy_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &recvRequests[1]);
				MPI_Irecv(h_upz_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &recvRequests[2]);
			}
			else if (rank == nNode-1){				
				// receive high halo region of d_t arrays from from rank-1 and copy to GPU 0
				MPI_Irecv(h_upx_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &recvRequests[0]);
				MPI_Irecv(h_upy_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &recvRequests[1]);
				MPI_Irecv(h_upz_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, &recvRequests[2]);
			}
			else {				
				// receive high halo region of d_t arrays from from rank-1 and copy to GPU 0
				MPI_Irecv(h_upx_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &recvRequests[0]);
				MPI_Irecv(h_upy_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 1, MPI_COMM_WORLD, &recvRequests[1]);
				MPI_Irecv(h_upz_h_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank-1, 2, MPI_COMM_WORLD, &recvRequests[2]);
		
				// receive low halo region of d_t arrays from from rank+1 and copy to GPU ngpu-1
				MPI_Irecv(h_upx_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &recvRequests[3]);
				MPI_Irecv(h_upy_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 1, MPI_COMM_WORLD, &recvRequests[4]);
				MPI_Irecv(h_upz_l_recv, 4 * fdm->nxpad * fdm->nzpad, MPI_FLOAT, rank+1, 2, MPI_COMM_WORLD, &recvRequests[5]);
			}
		}
		
		
		/*------------------------------------------------------------*/
		/* exchange halo regions of d_up arrays between local GPUs    */
		/*------------------------------------------------------------*/
		if (ngpu > 1){ // using multiple local GPUs, must exchange halo regions between neighboring local GPUs
			// high halo region of d_up arrays on GPU 0 to GPU 1
			cudaMemcpy(d_upx[1], d_upx[0] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			cudaMemcpy(d_upy[1], d_upy[0] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			cudaMemcpy(d_upz[1], d_upz[0] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			
			// exchange halo regions of d_up arrays between all internal GPUs
			for (int g = 1; g < ngpu-1; g++){
				// high halo region of GPU g to low halo region of GPU g+1
				cudaMemcpy(d_upx[g+1], d_upx[g] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_upy[g+1], d_upy[g] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_upz[g+1], d_upz[g] + (fdm->nxpad * fdm->nzpad * nyinterior), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				
				// low halo region of GPU g to high halo region of GPU g-1
				cudaMemcpy(d_upx[g-1] + (fdm->nxpad * fdm->nzpad * (nylocal - 4)), d_upx[g] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_upy[g-1] + (fdm->nxpad * fdm->nzpad * (nylocal - 4)), d_upy[g] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_upz[g-1] + (fdm->nxpad * fdm->nzpad * (nylocal - 4)), d_upz[g] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			}
			
			// low halo region of d_up arrays on GPU (ngpu-1) to GPU (ngpu-2)
			cudaMemcpy(d_upx[ngpu-2] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), d_upx[ngpu-1] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			cudaMemcpy(d_upy[ngpu-2] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), d_upy[ngpu-1] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			cudaMemcpy(d_upz[ngpu-2] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), d_upz[ngpu-1] + (4 * fdm->nxpad * fdm->nzpad), 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
		}
		
		
		/*------------------------------------------------------------*/
		/* from stress to acceleration         INSIDE                 */
		/*------------------------------------------------------------*/
		for (int g = 0; g < ngpu; g++){
			cudaSetDevice(g);
			dim3 dimGrid6((fdm->nxpad-2*NOP)/24.0f, (fdm->nzpad-2*NOP)/24.0f);
			dim3 dimBlock6(24,24,1);
			stressToAccelInside<<<dimGrid6, dimBlock6, 32*32*5*sizeof(float)>>>(fdm->nxpad, fdm->nzpad, nylocal, idx, idy, idz, d_txx[g], d_tyy[g], d_tzz[g], d_txy[g], d_tzx[g], d_tyz[g], d_uax[g], d_uay[g], d_uaz[g]);
		}
		sf_check_gpu_error(rank, "stressToAccel Kernel");
		
		
		/*------------------------------------------------------------*/
		/* inject acceleration source      INSIDE                     */
		/*------------------------------------------------------------*/
		if(!ssou && wavSrc) {
			for (int g = 0; g < ngpu; g++){
				cudaSetDevice(g);
				dim3 dimGrid7(ns, 1, 1);
				dim3 dimBlock7(2 * nbell + 1, 2 * nbell + 1, 1);
				lint3d_bell_gpuInside<<<dimGrid7, dimBlock7>>>(rank, ngpu, g, it, nc, ns, 0, nbell, fdm->nxpad, nyinterior, fdm->nzpad, d_uaz[g], d_bell[g], d_Sjx[g], d_Sjz[g], d_Sjy[g], d_ww[g], d_Sw000[g], d_Sw001[g], d_Sw010[g], d_Sw011[g], d_Sw100[g], d_Sw101[g], d_Sw110[g], d_Sw111[g]);
				lint3d_bell_gpuInside<<<dimGrid7, dimBlock7>>>(rank, ngpu, g, it, nc, ns, 1, nbell, fdm->nxpad, nyinterior, fdm->nzpad, d_uax[g], d_bell[g], d_Sjx[g], d_Sjz[g], d_Sjy[g], d_ww[g], d_Sw000[g], d_Sw001[g], d_Sw010[g], d_Sw011[g], d_Sw100[g], d_Sw101[g], d_Sw110[g], d_Sw111[g]);
				lint3d_bell_gpuInside<<<dimGrid7, dimBlock7>>>(rank, ngpu, g, it, nc, ns, 2, nbell, fdm->nxpad, nyinterior, fdm->nzpad, d_uay[g], d_bell[g], d_Sjx[g], d_Sjz[g], d_Sjy[g], d_ww[g], d_Sw000[g], d_Sw001[g], d_Sw010[g], d_Sw011[g], d_Sw100[g], d_Sw101[g], d_Sw110[g], d_Sw111[g]);	
			}
			sf_check_gpu_error(rank, "lint3d_bell_gpu Kernel");
		}
		
		
		/*------------------------------------------------------------*/
		/* step forward in time        INSIDE                         */
		/*------------------------------------------------------------*/
		for (int g = 0; g < ngpu; g++){
			cudaSetDevice(g);
			dim3 dimGrid8_inside(ceil(fdm->nxpad/192.0f), fdm->nzpad, nyinterior-8);
			dim3 dimBlock8_inside(192,1,1);
			stepTimeInside<<<dimGrid8_inside, dimBlock8_inside>>>(fdm->nxpad, nyinterior, fdm->nzpad, d_ro[g], d_uox[g], d_umx[g], d_uax[g], d_upx[g], d_uoy[g], d_umy[g], d_uay[g], d_upy[g], d_uoz[g], d_umz[g], d_uaz[g], d_upz[g]);
		}
		sf_check_gpu_error(rank, "stepTime Kernel");
		
		
		/*------------------------------------------------------------*/
		/* apply boundary conditions          INSIDE                  */
		/*------------------------------------------------------------*/
		if(dabc){
			
			/* One-way Absorbing BC */
			for (int g = 0; g < ngpu; g++){
				cudaSetDevice(g);
				
				dim3 dimGrid_abc_XY(ceil(fdm->nxpad/32.0f),ceil((nyinterior - 8)/32.0f),2);
				dim3 dimBlock_abc_XY(32,32,1);
				abcone3d_apply_XY_inside<<<dimGrid_abc_XY,dimBlock_abc_XY>>>(fdm->nxpad, nyinterior, fdm->nzpad, d_upx[g], d_uox[g], d_bzl_s[g], d_bzh_s[g]);
				abcone3d_apply_XY_inside<<<dimGrid_abc_XY,dimBlock_abc_XY>>>(fdm->nxpad, nyinterior, fdm->nzpad, d_upy[g], d_uoy[g], d_bzl_s[g], d_bzh_s[g]);
				abcone3d_apply_XY_inside<<<dimGrid_abc_XY,dimBlock_abc_XY>>>(fdm->nxpad, nyinterior, fdm->nzpad, d_upz[g], d_uoz[g], d_bzl_s[g], d_bzh_s[g]);
				
				dim3 dimGrid_abc_ZY(2, ceil((nyinterior-8)/32.0f), ceil(fdm->nzpad/32.0f));
				dim3 dimBlock_abc_ZY(1,32,32);
				abcone3d_apply_ZY_inside<<<dimGrid_abc_ZY,dimBlock_abc_ZY>>>(fdm->nxpad, nyinterior, fdm->nzpad, d_upx[g], d_uox[g], d_bxl_s[g], d_bxh_s[g]);
				abcone3d_apply_ZY_inside<<<dimGrid_abc_ZY,dimBlock_abc_ZY>>>(fdm->nxpad, nyinterior, fdm->nzpad, d_upy[g], d_uoy[g], d_bxl_s[g], d_bxh_s[g]);
				abcone3d_apply_ZY_inside<<<dimGrid_abc_ZY,dimBlock_abc_ZY>>>(fdm->nxpad, nyinterior, fdm->nzpad, d_upz[g], d_uoz[g], d_bxl_s[g], d_bxh_s[g]);
			}
			
			if (rank == 0){
				cudaSetDevice(0);
				dim3 dimGrid_abc_XZ(ceil(fdm->nxpad/32.0f),1,ceil(fdm->nzpad/32.0f));
				dim3 dimBlock_abc_XZ(32,1,32);
				abcone3d_apply_XZ_low<<<dimGrid_abc_XZ,dimBlock_abc_XZ>>>(fdm->nxpad, fdm->nzpad, d_upx[0], d_uox[0], d_byl_s[0]);
				abcone3d_apply_XZ_low<<<dimGrid_abc_XZ,dimBlock_abc_XZ>>>(fdm->nxpad, fdm->nzpad, d_upy[0], d_uoy[0], d_byl_s[0]);
				abcone3d_apply_XZ_low<<<dimGrid_abc_XZ,dimBlock_abc_XZ>>>(fdm->nxpad, fdm->nzpad, d_upz[0], d_uoz[0], d_byl_s[0]);
				
			}
			
			if (rank == nNode-1){
				cudaSetDevice(ngpu-1);
				dim3 dimGrid_abc_XZ(ceil(fdm->nxpad/32.0f),1,ceil(fdm->nzpad/32.0f));
				dim3 dimBlock_abc_XZ(32,1,32);
				abcone3d_apply_XZ_high<<<dimGrid_abc_XZ,dimBlock_abc_XZ>>>(fdm->nxpad, nylocal, fdm->nzpad, d_upx[ngpu-1], d_uox[ngpu-1], d_byh_s[ngpu-1]);
				abcone3d_apply_XZ_high<<<dimGrid_abc_XZ,dimBlock_abc_XZ>>>(fdm->nxpad, nylocal, fdm->nzpad, d_upy[ngpu-1], d_uoy[ngpu-1], d_byh_s[ngpu-1]);
				abcone3d_apply_XZ_high<<<dimGrid_abc_XZ,dimBlock_abc_XZ>>>(fdm->nxpad, nylocal, fdm->nzpad, d_upz[ngpu-1], d_uoz[ngpu-1], d_byh_s[ngpu-1]);
			}
			
			/* sponge BC */
			for (int g = 0; g < ngpu; g++){
				cudaSetDevice(g);
				dim3 dimGrid_spng_XY(ceil(fdm->nxpad/192.0f),(nyinterior - 8),1);
				dim3 dimBlock_spng_XY(192,1,1);                                            
				sponge3d_apply_XY_inside<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(d_umz[g], fdm->nxpad, nyinterior, fdm->nzpad, nb);
				sponge3d_apply_XY_inside<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(d_uoz[g], fdm->nxpad, nyinterior, fdm->nzpad, nb);
				sponge3d_apply_XY_inside<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(d_upz[g], fdm->nxpad, nyinterior, fdm->nzpad, nb);
			                                                                   
				sponge3d_apply_XY_inside<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(d_umx[g], fdm->nxpad, nyinterior, fdm->nzpad, nb);
				sponge3d_apply_XY_inside<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(d_uox[g], fdm->nxpad, nyinterior, fdm->nzpad, nb);
				sponge3d_apply_XY_inside<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(d_upx[g], fdm->nxpad, nyinterior, fdm->nzpad, nb);
			                                                                   
				sponge3d_apply_XY_inside<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(d_umy[g], fdm->nxpad, nyinterior, fdm->nzpad, nb);
				sponge3d_apply_XY_inside<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(d_uoy[g], fdm->nxpad, nyinterior, fdm->nzpad, nb);
				sponge3d_apply_XY_inside<<<dimGrid_spng_XY,dimBlock_spng_XY>>>(d_upy[g], fdm->nxpad, nyinterior, fdm->nzpad, nb);
			
			
				dim3 dimGrid_spng_ZY(ceil(nb/8.0f),ceil(fdm->nzpad/8.0f),ceil((nyinterior - 8)/8.0f));
				dim3 dimBlock_spng_ZY(8,8,8);
				sponge3d_apply_ZY_inside<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(d_umz[g], fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
				sponge3d_apply_ZY_inside<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(d_uoz[g], fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
				sponge3d_apply_ZY_inside<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(d_upz[g], fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
							                                                                   
				sponge3d_apply_ZY_inside<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(d_umx[g], fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
				sponge3d_apply_ZY_inside<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(d_uox[g], fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
				sponge3d_apply_ZY_inside<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(d_upx[g], fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
						                                                                       
				sponge3d_apply_ZY_inside<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(d_umy[g], fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
				sponge3d_apply_ZY_inside<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(d_uoy[g], fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
				sponge3d_apply_ZY_inside<<<dimGrid_spng_ZY,dimBlock_spng_ZY>>>(d_upy[g], fdm->nxpad, nyinterior, fdm->nzpad, nb, nx, spo);
			}
			
			if (rank == 0){
				cudaSetDevice(0);
				dim3 dimGrid_spng_XZ(ceil(fdm->nxpad/192.0f),1,fdm->nzpad);
				dim3 dimBlock_spng_XZ(192,1,1);
				sponge3d_apply_XZ_low<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_umz[0], fdm->nxpad, fdm->nzpad, nb);
				sponge3d_apply_XZ_low<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_uoz[0], fdm->nxpad, fdm->nzpad, nb);
				sponge3d_apply_XZ_low<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_upz[0], fdm->nxpad, fdm->nzpad, nb);
					                                                                          
				sponge3d_apply_XZ_low<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_umx[0], fdm->nxpad, fdm->nzpad, nb);
				sponge3d_apply_XZ_low<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_uox[0], fdm->nxpad, fdm->nzpad, nb);
				sponge3d_apply_XZ_low<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_upx[0], fdm->nxpad, fdm->nzpad, nb);
					                                                                          
				sponge3d_apply_XZ_low<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_umy[0], fdm->nxpad, fdm->nzpad, nb);
				sponge3d_apply_XZ_low<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_uoy[0], fdm->nxpad, fdm->nzpad, nb);
				sponge3d_apply_XZ_low<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_upy[0], fdm->nxpad, fdm->nzpad, nb);
			}
			
			if (rank == nNode-1){
				cudaSetDevice(ngpu-1);
				dim3 dimGrid_spng_XZ(ceil(fdm->nxpad/192.0f),1,fdm->nzpad);
				dim3 dimBlock_spng_XZ(192,1,1);                                                                      
				sponge3d_apply_XZ_high<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_umz[ngpu-1], fdm->nxpad, nylocal, fdm->nzpad, nb);
				sponge3d_apply_XZ_high<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_uoz[ngpu-1], fdm->nxpad, nylocal, fdm->nzpad, nb);
				sponge3d_apply_XZ_high<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_upz[ngpu-1], fdm->nxpad, nylocal, fdm->nzpad, nb);
					                                                                                                               
				sponge3d_apply_XZ_high<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_umx[ngpu-1], fdm->nxpad, nylocal, fdm->nzpad, nb);
				sponge3d_apply_XZ_high<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_uox[ngpu-1], fdm->nxpad, nylocal, fdm->nzpad, nb);
				sponge3d_apply_XZ_high<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_upx[ngpu-1], fdm->nxpad, nylocal, fdm->nzpad, nb);
					                                                                                                               
				sponge3d_apply_XZ_high<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_umy[ngpu-1], fdm->nxpad, nylocal, fdm->nzpad, nb);
				sponge3d_apply_XZ_high<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_uoy[ngpu-1], fdm->nxpad, nylocal, fdm->nzpad, nb);
				sponge3d_apply_XZ_high<<<dimGrid_spng_XZ,dimBlock_spng_XZ>>>(d_upy[ngpu-1], fdm->nxpad, nylocal, fdm->nzpad, nb);
			}
			
			sf_check_gpu_error(rank, "Boundary Condition Kernels");
			
		}
		
		
		/*------------------------------------------------------------*/
		/* copy halos from d_up arrays from remote GPU into local GPU */
		/*------------------------------------------------------------*/
		// Wait for all MPI transfers to complete
		if (nNode > 1){
			if (rank == 0 || rank == nNode-1){
				MPI_Waitall(3, recvRequests, recvStatuses);
			}
			else {
				MPI_Waitall(6, recvRequests, recvStatuses);
			}
		}
		if (nNode > 1){	// using multiple nodes, must exchange halo regions between neighboring nodes GPUs
			if (rank == 0){
				cudaMemcpy(d_upx[ngpu-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), h_upx_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_upy[ngpu-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), h_upy_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_upz[ngpu-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), h_upz_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			}
			else if (rank == nNode-1){				
				cudaMemcpy(d_upx[0], h_upx_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_upy[0], h_upy_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_upz[0], h_upz_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
		
			}
			else {				
				cudaMemcpy(d_upx[0], h_upx_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_upy[0], h_upy_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_upz[0], h_upz_h_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				
				cudaMemcpy(d_upx[ngpu-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), h_upx_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_upy[ngpu-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), h_upy_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
				cudaMemcpy(d_upz[ngpu-1] + (fdm->nxpad * fdm->nzpad * (nyinterior + 4)), h_upz_l_recv, 4 * fdm->nxpad * fdm->nzpad * sizeof(float), cudaMemcpyDefault);
			}
		}
		 
					
		/* circulate wavefield arrays */
		for (int g = 0; g < ngpu; g++){
			d_utz[g]=d_umz[g]; d_uty[g]=d_umy[g]; d_utx[g]=d_umx[g];
			d_umz[g]=d_uoz[g]; d_umy[g]=d_uoy[g]; d_umx[g]=d_uox[g];
			d_uoz[g]=d_upz[g]; d_uoy[g]=d_upy[g]; d_uox[g]=d_upx[g];
			d_upz[g]=d_utz[g]; d_upy[g]=d_uty[g]; d_upx[g]=d_utx[g];
		}
		
		
		/*------------------------------------------------------------*/
		/* cut wavefield and save */
		/*------------------------------------------------------------*/
		if(snap && it%jsnap==0) {
			
			if (rank == 0){	// accumulate wavefield data and write to file
				
				// write local GPU's portions of wavefield data into output arrays
				if (nNode == 1 && ngpu == 1){
					cudaMemcpy(h_uox, d_uox[0], fdm->nypad * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDefault);
					cudaMemcpy(h_uoy, d_uoy[0], fdm->nypad * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDefault);
					cudaMemcpy(h_uoz, d_uoz[0], fdm->nypad * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDefault);
				}
				else {
					for (int g = 0; g < ngpu; g++){
						cudaMemcpy(h_uox + g * (nyinterior * fdm->nzpad * fdm->nxpad), d_uox[g] + 4 * fdm->nzpad * fdm->nxpad, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDefault);
						cudaMemcpy(h_uoy + g * (nyinterior * fdm->nzpad * fdm->nxpad), d_uoy[g] + 4 * fdm->nzpad * fdm->nxpad, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDefault);
						cudaMemcpy(h_uoz + g * (nyinterior * fdm->nzpad * fdm->nxpad), d_uoz[g] + 4 * fdm->nzpad * fdm->nxpad, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDefault);
					}
				}
			
				// receive wavefield from other GPUs and write to output arrays
				for (int r = 1; r < nNode; r++){
					MPI_Recv(h_uox + r * ngpu * nyinterior * fdm->nzpad * fdm->nxpad, ngpu * nyinterior * fdm->nzpad * fdm->nxpad, MPI_FLOAT, r, 0, MPI_COMM_WORLD, &status);
					MPI_Recv(h_uoy + r * ngpu * nyinterior * fdm->nzpad * fdm->nxpad, ngpu * nyinterior * fdm->nzpad * fdm->nxpad, MPI_FLOAT, r, 1, MPI_COMM_WORLD, &status);
					MPI_Recv(h_uoz + r * ngpu * nyinterior * fdm->nzpad * fdm->nxpad, ngpu * nyinterior * fdm->nzpad * fdm->nxpad, MPI_FLOAT, r, 2, MPI_COMM_WORLD, &status);
				}
			
			
				for (int y = 0; y < fdm->nypad; y++){
					for (int z = 0; z < fdm->nzpad; z++){
						for (int x = 0; x < fdm->nxpad; x++){
							uox[y][x][z] = h_uox[y * fdm->nzpad * fdm->nxpad + z * fdm->nxpad + x];
							uoy[y][x][z] = h_uoy[y * fdm->nzpad * fdm->nxpad + z * fdm->nxpad + x];
							uoz[y][x][z] = h_uoz[y * fdm->nzpad * fdm->nxpad + z * fdm->nxpad + x];
						}
					}
				}
				
				// Write complete wavefield to output file
				// sf_floatwrite(h_uoz,fdm->nypad * fdm->nzpad * fdm->nxpad,Fwfl);
				// sf_floatwrite(h_uox,fdm->nypad * fdm->nzpad * fdm->nxpad,Fwfl);
				// sf_floatwrite(h_uoy,fdm->nypad * fdm->nzpad * fdm->nxpad,Fwfl);

				// Write wavefield arrays to output file
				cut3d(uoz,uc,fdm,acz,acx,acy);	// flip axes back
				sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);
				
				cut3d(uox,uc,fdm,acz,acx,acy);	// flip axes back
				sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);
				
				cut3d(uoy,uc,fdm,acz,acx,acy);	// flip axes back
				sf_floatwrite(uc[0][0],sf_n(acx)*sf_n(acy)*sf_n(acz),Fwfl);
			
			}
			else {	// send wavefield data to rank0
			
				// write local GPU's portions of wavefield data into output arrays
				for (int g = 0; g < ngpu; g++){
					cudaMemcpy(h_uox + g * nyinterior * fdm->nzpad * fdm->nxpad, d_uox[g] + 4 * fdm->nzpad * fdm->nxpad, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDefault);
					cudaMemcpy(h_uoy + g * nyinterior * fdm->nzpad * fdm->nxpad, d_uoy[g] + 4 * fdm->nzpad * fdm->nxpad, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDefault);
					cudaMemcpy(h_uoz + g * nyinterior * fdm->nzpad * fdm->nxpad, d_uoz[g] + 4 * fdm->nzpad * fdm->nxpad, nyinterior * fdm->nzpad * fdm->nxpad * sizeof(float), cudaMemcpyDefault);
				}
			
				MPI_Send(h_uox, ngpu * nyinterior * fdm->nzpad * fdm->nxpad, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
				MPI_Send(h_uoy, ngpu * nyinterior * fdm->nzpad * fdm->nxpad, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
				MPI_Send(h_uoz, ngpu * nyinterior * fdm->nzpad * fdm->nxpad, MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
			}

		}
		
		/*------------------------------------------------------------*/
		/* extract receiver data									  */
		/*------------------------------------------------------------*/
		if(it%jdata==0) {
			for (int g = 0; g < ngpu; g++){
				cudaSetDevice(g);
				cudaMemset(d_dd[g], 0, nr*nc*sizeof(float));
				dim3 dimGrid_extract(MIN(nr,ceil(nr/1024.0f)), 1, 1);
				dim3 dimBlock_extract(MIN(nr, 1024), 1, 1);
				if (interp){
					lint3d_extract_gpu<<<dimGrid_extract, dimBlock_extract>>>(rank, ngpu, g, d_dd[g], nr, fdm->nxpad, nyinterior, fdm->nzpad, d_uoz[g], d_uox[g], d_uoy[g], d_Rjz[g], d_Rjx[g], d_Rjy[g], d_Rw000[g], d_Rw001[g], d_Rw010[g], d_Rw011[g], d_Rw100[g], d_Rw101[g], d_Rw110[g], d_Rw111[g]);
					sf_check_gpu_error(rank, "lint3d_extract kernel");
				}
				else {
					extract_gpu<<<dimGrid_extract, dimBlock_extract>>>(rank, ngpu, g, d_dd[g], nr, fdm->nxpad, nyinterior, fdm->nzpad, d_uoz[g], d_uox[g], d_uoy[g], d_Rjz[g], d_Rjx[g], d_Rjy[g]);
					sf_check_gpu_error(rank, "extract_gpu kernel");
				}
			}
			
			// copy local GPU 0's receiver data into h_dd_locals
			cudaMemcpy(h_dd_locals, d_dd[0], nr * nc * sizeof(float), cudaMemcpyDefault);
			
			// add all other GPU's recever data to h_dd_combined
			for (int g = 1; g < ngpu; g++){
				cudaMemcpy(h_dd, d_dd[g], nr * nc * sizeof(float), cudaMemcpyDefault);
				for (int i = 0; i < nr * nc; i++){
					h_dd_locals[i] += h_dd[i];
				}
			}
			
			// MPI Reduce with SUM on the h_dd_locals arrays from each GPU to accumulate reciever data into h_dd_allNodes
			MPI_Reduce(h_dd_locals, h_dd_allNodes, nr * nc, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			
			if (rank == 0){
				// write receiver data to output file
				sf_floatwrite(h_dd_allNodes, nr*nc, Fdat);
			}
			
		}
		
		
	} // END MAIN LOOP
	
	
	/*------------------------------------------------------------*/
    /* deallocate host arrays */

	if (wavSrc){
		free(**ww); free(*ww); free(ww); free(h_ww);
	}
	free(h_dd); free(h_dd_combined);
	free(ss); free(rr);
	free(h_bell);
	free(h_ro);
	free(h_c11); free(h_c22); free(h_c33); free(h_c44); free(h_c55); free(h_c66); free(h_c12); free(h_c13); free(h_c23);
	
	if (snap){
		free(h_uoz); free(h_uox); free(h_uoy);
		free(**uc);  free(*uc);  free(uc);
		free(**uoz); free(*uoz); free(uoz);
	    free(**uox); free(*uox); free(uox);
	    free(**uoy); free(*uoy); free(uoy);
	}
	
	
	/*------------------------------------------------------------*/
    /* deallocate GPU arrays */

	for (int g = 0; g < ngpu; g++){
	
		if (wavSrc){
			cudaFree(&d_ww[g]);
		}
		cudaFree(&d_dd[g]);
		cudaFree(&d_bell[g]);
	
		cudaFree(&d_ro[g]);
		cudaFree(&d_c11[g]);
		cudaFree(&d_c22[g]);
		cudaFree(&d_c33[g]);
		cudaFree(&d_c44[g]);
		cudaFree(&d_c55[g]);
		cudaFree(&d_c66[g]);
		cudaFree(&d_c12[g]);
		cudaFree(&d_c13[g]);
		cudaFree(&d_c23[g]);
	
		if (dabc){
			cudaFree(&d_bzl_s[g]);
			cudaFree(&d_bzh_s[g]);
			cudaFree(&d_bxl_s[g]);
			cudaFree(&d_bxh_s[g]);
			cudaFree(&d_byl_s[0]);
			cudaFree(&d_byh_s[ngpu-1]);
		}
	
		cudaFree(&d_umx[g]); cudaFree(&d_umy[g]); cudaFree(&d_umz[g]);
		cudaFree(&d_uox[g]); cudaFree(&d_uoy[g]); cudaFree(&d_uoz[g]);
		cudaFree(&d_upx[g]); cudaFree(&d_upy[g]); cudaFree(&d_upz[g]);
		cudaFree(&d_uax[g]); cudaFree(&d_uay[g]); cudaFree(&d_uaz[g]);
	
		cudaFree(&d_tzz[g]); cudaFree(&d_tyy[g]); cudaFree(&d_txx[g]); 
		cudaFree(&d_txy[g]); cudaFree(&d_tyz[g]); cudaFree(&d_tzx[g]);
	
		cudaFree(&d_Sjz[g]);
		cudaFree(&d_Sjx[g]);
		cudaFree(&d_Sjy[g]);
		cudaFree(&d_Sw000[g]);
		cudaFree(&d_Sw001[g]);
		cudaFree(&d_Sw010[g]);
		cudaFree(&d_Sw011[g]);
		cudaFree(&d_Sw100[g]);
		cudaFree(&d_Sw101[g]);
		cudaFree(&d_Sw110[g]);
		cudaFree(&d_Sw111[g]);
	
		cudaFree(&d_Rjz[g]);
		cudaFree(&d_Rjx[g]);
		cudaFree(&d_Rjy[g]);
		if (interp){
			cudaFree(&d_Rw000[g]);
			cudaFree(&d_Rw001[g]);
			cudaFree(&d_Rw010[g]);
			cudaFree(&d_Rw011[g]);
			cudaFree(&d_Rw100[g]);
			cudaFree(&d_Rw101[g]);
			cudaFree(&d_Rw110[g]);
			cudaFree(&d_Rw111[g]);
		}
	}
	
	
	MPI_Finalize();
	sf_close();
	exit(0);
	
}








