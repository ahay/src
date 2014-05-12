/*
  Copyright (C) 2013  Xi'an Jiaotong University (Pengliang Yang)
    Email: ypl.2100@gmail.com	
    Acknowledgement: This code is written with the help of Baoli Wang.

References:
    [1] Clayton, Robert, and Björn Engquist. "Absorbing boundary 
	conditions for acoustic and elastic wave equations." Bulletin 
	of the Seismological Society of America 67.6 (1977): 1529-1540.
    [2] Tarantola, Albert. "Inversion of seismic reflection data in the 
	acoustic approximation." Geophysics 49.8 (1984): 1259-1266.
    [3] Pica, A., J. P. Diet, and A. Tarantola. "Nonlinear inversion 
	of seismic reflection data in a laterally invariant medium." 
	Geophysics 55.3 (1990): 284-292.
    [4] Hager, William W., and Hongchao Zhang. "A survey of nonlinear
	conjugate gradient methods." Pacific journal of Optimization 
	2.1 (2006): 35-58.
*/
__global__ void cuda_set_sg(int *sxz, int sxbeg, int szbeg, int jsx, int jsz, int ns, int nz)
/*< set the positions of sources/geophones >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<ns) sxz[id]=(szbeg+id*jsz)+nz*(sxbeg+id*jsx);
}

__global__ void cuda_ricker_wavelet(float *wlt, float amp, float fm, float dt, int nt)
/*< generate ricker wavelet with time deley >*/
{
	int it=threadIdx.x+blockDim.x*blockIdx.x;
    	float tmp = PI*fm*(it*dt-1.0/fm);
    	tmp *=tmp;
    	if (it<nt) wlt[it]= amp*(1.0-2.0*tmp)*expf(-tmp);
}

__global__ void cuda_add_source(float *p, float *source, int *sxz, int ns, bool add)
/*< add==true, add (inject) the source; add==false, subtract the source >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<ns)
	{
		if (add)	p[sxz[id]]+=source[id];
		else 		p[sxz[id]]-=source[id];
	}	
}


__global__ void cuda_step_forward(float *p0, float *p1, float *vv, float dtz, float dtx, int nz, int nx)
{
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int id=i1+i2*nz;

	__shared__ float s_p0[Block_Size1+2][Block_Size2+2];
	__shared__ float s_p1[Block_Size1+2][Block_Size2+2];
	s_p0[threadIdx.x+1][threadIdx.y+1]=p0[id];
	s_p1[threadIdx.x+1][threadIdx.y+1]=p1[id];
	if(threadIdx.x==0)
	{
		if(blockIdx.x>0)		{ s_p0[threadIdx.x][threadIdx.y+1]=p0[id-1];	s_p1[threadIdx.x][threadIdx.y+1]=p1[id-1];}
		else				{ s_p0[threadIdx.x][threadIdx.y+1]=0.0f;	s_p1[threadIdx.x][threadIdx.y+1]=0.0f;}
	}
	if(threadIdx.x==blockDim.x-1)
	{
		if(blockIdx.x<gridDim.x-1)	{ s_p0[threadIdx.x+2][threadIdx.y+1]=p0[id+1];	s_p1[threadIdx.x+2][threadIdx.y+1]=p1[id+1];}
		else				{ s_p0[threadIdx.x+2][threadIdx.y+1]=0.0f;	s_p1[threadIdx.x+2][threadIdx.y+1]=0.0f;}
	}
	if(threadIdx.y==0)
	{
		if(blockIdx.y>0)		{ s_p0[threadIdx.x+1][threadIdx.y]=p1[id-nz]; 	s_p1[threadIdx.x+1][threadIdx.y]=p1[id-nz];}
		else				{ s_p0[threadIdx.x+1][threadIdx.y]=0.0f;	s_p1[threadIdx.x+1][threadIdx.y]=0.0f;}
	}
	if(threadIdx.y==blockDim.y-1)
	{
		if(blockIdx.y<gridDim.y-1)	{ s_p0[threadIdx.x+1][threadIdx.y+2]=p1[id+nz];	s_p1[threadIdx.x+1][threadIdx.y+2]=p1[id+nz];}
		else				{ s_p0[threadIdx.x+1][threadIdx.y+2]=0.0f;	s_p1[threadIdx.x+1][threadIdx.y+2]=0.0f;}
	}
	__syncthreads();

	float v1=vv[id]*dtz;
	float v2=vv[id]*dtx; 
	float c1=v1*v1*(s_p1[threadIdx.x+2][threadIdx.y+1]-2.0*s_p1[threadIdx.x+1][threadIdx.y+1]+s_p1[threadIdx.x][threadIdx.y+1]);
	float c2=v2*v2*(s_p1[threadIdx.x+1][threadIdx.y+2]-2.0*s_p1[threadIdx.x+1][threadIdx.y+1]+s_p1[threadIdx.x+1][threadIdx.y]);
/*
	if(i1==0)// top boundary: free surface, no absorbing
	{
		c1=v1*(-s_p1[threadIdx.x+1][threadIdx.y+1]+s_p1[threadIdx.x+2][threadIdx.y+1]
					+s_p0[threadIdx.x+1][threadIdx.y+1]-s_p0[threadIdx.x+2][threadIdx.y+1]);
		if(i2>0 && i2<nx-1) c2=0.5*c2;
	}
*/
	if(i1==nz-1) /* bottom boundary*/
	{
		c1=v1*(s_p1[threadIdx.x][threadIdx.y+1]-s_p1[threadIdx.x+1][threadIdx.y+1]
					-s_p0[threadIdx.x][threadIdx.y+1]+s_p0[threadIdx.x+1][threadIdx.y+1]);
		if(i2>0 && i2<nx-1) c2=0.5*c2;
	}
	if(i2==0)/* left boundary*/
	{
		if(i1>0 && i1<nz-1) c1=0.5*c1;
		c2=v2*(-s_p1[threadIdx.x+1][threadIdx.y+1]+s_p1[threadIdx.x+1][threadIdx.y+2]
					+s_p0[threadIdx.x+1][threadIdx.y+1]-s_p0[threadIdx.x+1][threadIdx.y+2]);

	}
	if(i2==nx-1) /* right boundary*/
	{
		if(i1>0 && i1<nz-1) c1=0.5*c1;
		c2=v2*(s_p1[threadIdx.x+1][threadIdx.y]-s_p1[threadIdx.x+1][threadIdx.y+1]
					-s_p0[threadIdx.x+1][threadIdx.y]+s_p0[threadIdx.x+1][threadIdx.y+1]);
	}
	
	
	if (i1<nz && i2<nx) p0[id]=2*s_p1[threadIdx.x+1][threadIdx.y+1]-s_p0[threadIdx.x+1][threadIdx.y+1]+c1+c2;
}


__global__ void cuda_step_backward(float *p0, float *p1, float *vv, float *lap, float *sillum, float dtz, float dtx, int nz, int nx)
{
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int id=i1+i2*nz;

	__shared__ float s_p1[Block_Size1+2][Block_Size2+2];
	s_p1[threadIdx.x+1][threadIdx.y+1]=p1[id];
	if(threadIdx.x==0)
	{
		if(blockIdx.x>0)		{ s_p1[threadIdx.x][threadIdx.y+1]=p1[id-1];}
		else				{ s_p1[threadIdx.x][threadIdx.y+1]=0.0f;}
	}
	if(threadIdx.x==blockDim.x-1)
	{
		if(blockIdx.x<gridDim.x-1)	{s_p1[threadIdx.x+2][threadIdx.y+1]=p1[id+1];}
		else				{s_p1[threadIdx.x+2][threadIdx.y+1]=0.0f;}
	}
	if(threadIdx.y==0)
	{
		if(blockIdx.y>0)		{s_p1[threadIdx.x+1][threadIdx.y]=p1[id-nz];}
		else				{s_p1[threadIdx.x+1][threadIdx.y]=0.0f;}
	}
	if(threadIdx.y==blockDim.y-1)
	{
		if(blockIdx.y<gridDim.y-1)	{s_p1[threadIdx.x+1][threadIdx.y+2]=p1[id+nz];}
		else				{s_p1[threadIdx.x+1][threadIdx.y+2]=0.0f;}
	}
	__syncthreads();

	float v1=vv[id]*dtz;
	float v2=vv[id]*dtx; 
	float c1=v1*v1*(s_p1[threadIdx.x+2][threadIdx.y+1]-2.0*s_p1[threadIdx.x+1][threadIdx.y+1]+s_p1[threadIdx.x][threadIdx.y+1]);
	float c2=v2*v2*(s_p1[threadIdx.x+1][threadIdx.y+2]-2.0*s_p1[threadIdx.x+1][threadIdx.y+1]+s_p1[threadIdx.x+1][threadIdx.y]);

	if (i1<nz && i2<nx) 
	{
		p0[id]=2*s_p1[threadIdx.x+1][threadIdx.y+1]-p0[id]+c1+c2;
		lap[id]=c1+c2;/* laplacian(p) */
		sillum[id]+=s_p1[threadIdx.x+1][threadIdx.y+1]*s_p1[threadIdx.x+1][threadIdx.y+1];
	}
}

__global__ void cuda_rw_bndr(float *bndr, float *p1, int nz, int nx, bool read)
/*< read/write boundaries >*/
{
	int id=threadIdx.x+blockIdx.x*blockDim.x;
	if(read){
		if(id<nz) p1[id]=bndr[id];//left boundary
		else if (id<2*nz) p1[(id-nz)+nz*(nx-1)]=bndr[id];//right boundary
		else if (id<2*nz+nx) p1[nz-1+nz*(id-2*nz)]=bndr[id];//bottom boundary
	}else{
		if(id<nz) bndr[id]=p1[id];//left boundary
		else if (id<2*nz) bndr[id]=p1[(id-nz)+nz*(nx-1)];//right boundary
		else if (id<2*nz+nx) bndr[id]=p1[nz-1+nz*(id-2*nz)];//bottom boundary
	}
}

__global__ void cuda_record(float*p, float *seis, int *gxz, int ng)
/*< record the seismogram at time it >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<ng) seis[id]=p[gxz[id]];
}

__global__ void cuda_cal_residual(float *dobs, float *derr, float *p, int *gxz, int ng)
/*< calculate residual wavefield at the receiver positions
   dcal: p[gxz[id]]
   dobs: d_{obs}
   derr: d_{err}=d_{cal}-d_{obs} >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	if(id<ng) derr[id]=p[gxz[id]]-dobs[id];
}

__global__ void cuda_cal_objective(float *obj, float *err, int ng)
/*< calculate the value of objective/misfit function: obj >*/
{
  	__shared__ float  sdata[Block_Size];
    	int tid=threadIdx.x;
	int s,id;
	float a;
    	sdata[tid]=0.0f;
	for(s=0; s<(ng+Block_Size-1)/Block_Size; s++)
	{
		id=s*blockDim.x+threadIdx.x;
		a=(id<ng)?err[id]:0.0f;
		sdata[tid] += a*a;	
	} 
    	__syncthreads();

    	/* do reduction in shared mem */
    	for(s=blockDim.x/2; s>32; s>>=1) 
    	{/* synchronization is required between different blocks */	
		if (threadIdx.x < s) sdata[tid] += sdata[tid + s]; 
		__syncthreads();
    	}
   	if (tid < 32)
   	{/* no synchronization required in the same block */
		if (blockDim.x >=  64) { sdata[tid] += sdata[tid + 32]; }
		if (blockDim.x >=  32) { sdata[tid] += sdata[tid + 16]; }
		if (blockDim.x >=  16) { sdata[tid] += sdata[tid +  8]; }
		if (blockDim.x >=   8) { sdata[tid] += sdata[tid +  4]; }
		if (blockDim.x >=   4) { sdata[tid] += sdata[tid +  2]; }
		if (blockDim.x >=   2) { sdata[tid] += sdata[tid +  1]; }
    	}
     
    	if (tid == 0) { *obj=sdata[0]; }
}


__global__ void cuda_cal_gradient(float *g1, float *lap, float *gp, int nz, int nx)
/*< calculate gradient >*/
{
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int id=i1+nz*i2;	

	if (i1<nz && i2<nx) 
	{
		g1[id]+=lap[id]*gp[id];
	}
}

__global__ void cuda_scale_gradient(float *g1, float *vv, float *illum, int nz, int nx)
/*< scale gradient >*/
{
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int id=i1+nz*i2;
	if (i1>=1 && i1<nz-1 && i2>=1 && i2<nx-1)	g1[id]*=2.0f/(vv[id]);//g1[id]*=2.0f/(vv[id]*illum[id]);
	__syncthreads();
	// handling the outliers at the boundary
	if (i1==0) 	g1[id]=g1[1+nz*i2];
	else if (i1==nz-1)	g1[id]=g1[nz-2+nz*i2];
	__syncthreads();
	if (i2==0)	g1[id]=g1[i1+nz];
	else if (i2==nx-1)	g1[id]=g1[i1+nz*(nx-2)];
	__syncthreads();
}
__global__ void cuda_init_bell(float *bell)
/*< initialize 1D bell function:<<<1, 2*nbell+1>>> >*/
{
    	int tid=threadIdx.x;
	bell[tid]=expf(-(tid-nbell)*(tid-nbell)/(0.5*nbell));
}


__global__ void cuda_bell_smoothz(float *g, float *smg, float *bell, int nz, int nx)
/*< smoothing with gaussian function >*/
{
	int i;
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*blockDim.x;
	int id=i1+i2*nz;
	float s=0;
	for(i=-nbell; i<=nbell; i++) if(i1+i>=0 && i1+i<nz) s+=bell[i+nbell]*g[id+i];
	if(i1<nz && i2<nx) smg[id]=s;
}
__global__ void cuda_bell_smoothx(float *g, float *smg, float *bell, int nz, int nx)
/*< smoothing with gaussian function >*/
{
	int i;
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*blockDim.x;
	int id=i1+i2*nz;
	float s=0;
	for(i=-nbell; i<=nbell; i++) if(i2+i>=0 && i2+i<nx) s+=bell[i+nbell]*g[id+nz*i];
	if(i1<nz && i2<nx) smg[id]=s;
}

__global__ void cuda_cal_beta(float *beta, float *g0, float *g1, float *cg, int N)
/*< calculate beta for nonlinear conjugate gradient algorithm 
configuration requirement: <<<1,Block_Size>>> >*/
{
    	__shared__ float sdata[Block_Size];
	__shared__ float tdata[Block_Size];
	__shared__ float rdata[Block_Size];
    	int tid = threadIdx.x;
	int s, id; 
	float a,b,c, beta_HS, beta_DY;
    	sdata[tid] = 0.0f;
	tdata[tid] = 0.0f;
	rdata[tid] = 0.0f;
	for(s=0; s<(N+Block_Size-1)/Block_Size; s++)
	{
		id=s*blockDim.x+threadIdx.x;
		a=(id<N)?g0[id]:0.0f;
		b=(id<N)?g1[id]:0.0f;
		c=(id<N)?cg[id]:0.0f;

		// HS: Hestenses-Stiefel NLCG algorithm 
		sdata[tid] += b*(b-a);	// numerator of HS
		tdata[tid] += c*(b-a);	// denominator of HS,DY
		rdata[tid] += b*b;	// numerator of DY
		
/*   	
		// PRP: Polark-Ribiere-Polyar NLCG algorithm 
		sdata[tid] += b*(b-a);	// numerator
		tdata[tid] += a*a;	// denominator
		// HS: Hestenses-Stiefel NLCG algorithm 
		sdata[tid] += b*(b-a);	// numerator
		tdata[tid] += c*(b-a);	// denominator
		// FR: Fletcher-Reeves NLCG algorithm 
		sdata[tid] += b*b;	// numerator
		tdata[tid] += a*a;	// denominator
		// PRP: Polark-Ribiere-Polyar NLCG algorithm 
		sdata[tid] += b*(b-a);	// numerator
		tdata[tid] += a*a;	// denominator
		// CD: Fletcher NLCG algorithm  
		sdata[tid] += b*b;	// numerator
		tdata[tid] -= c*a;	// denominator
		// DY: Dai-Yuan NLCG algorithm 
		sdata[tid] += b*b;	// numerator
		tdata[tid] += c*(b-a);	// denominator
*/
	} 
    	__syncthreads();

    	// do reduction in shared mem
    	for(s=blockDim.x/2; s>32; s>>=1) 
    	{
		if (tid < s)	{ sdata[tid]+=sdata[tid+s]; tdata[tid]+=tdata[tid+s]; rdata[tid]+=rdata[tid+s];}
		__syncthreads();
    	}     
   	if (tid < 32)
   	{
		if (blockDim.x >=64) { sdata[tid]+=sdata[tid+32]; tdata[tid]+=tdata[tid+32]; rdata[tid]+=rdata[tid+32];}
		if (blockDim.x >=32) { sdata[tid]+=sdata[tid+16]; tdata[tid]+=tdata[tid+16]; rdata[tid]+=rdata[tid+16];}
		if (blockDim.x >=16) { sdata[tid]+=sdata[tid+ 8]; tdata[tid]+=tdata[tid+ 8]; rdata[tid]+=rdata[tid+ 8];}
		if (blockDim.x >= 8) { sdata[tid]+=sdata[tid+ 4]; tdata[tid]+=tdata[tid+ 4]; rdata[tid]+=rdata[tid+ 4];}
		if (blockDim.x >= 4) { sdata[tid]+=sdata[tid+ 2]; tdata[tid]+=tdata[tid+ 2]; rdata[tid]+=rdata[tid+ 2];}
		if (blockDim.x >= 2) { sdata[tid]+=sdata[tid+ 1]; tdata[tid]+=tdata[tid+ 1]; rdata[tid]+=rdata[tid+ 1];}
    	}
     
	if (tid == 0) 
	{ 
		beta_HS=beta_DY=0;
		if(fabsf(tdata[0])>EPS) 
		{
			beta_HS=sdata[0]/tdata[0]; 
			beta_DY=rdata[0]/tdata[0];
		} 
		*beta=MAX(0, MIN(beta_HS, beta_DY));// Hybrid HS-DY method combined with iteration restart
	}	
}

__global__ void cuda_cal_conjgrad(float *g1, float *cg, float beta, int nz, int nx)
/*< calculate nonlinear conjugate gradient >*/
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+i2*nz;

	if (i1<nz && i2<nx) cg[id]=-g1[id]+beta*cg[id];
}

__global__ void cuda_cal_vtmp(float *vtmp, float *vv, float *cg, float *epsil, int nz, int nx)
/*< calculate temporary velocity >*/ 
{
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*blockDim.x;
	int id=i1+i2*nz;
	if (i1<nz && i2<nx)	vtmp[id]=vv[id]+(*epsil)*cg[id];
}

__global__ void cuda_cal_epsilon(float *vv, float *cg, float *epsil, int N)
/*< calculate estimated stepsize (epsil) according to Taratola's method
configuration requirement: <<<1, Block_Size>>> >*/ 
{
	int s,id;
	float a,b;
    	__shared__ float sdata[Block_Size];/* find max(|vv(:)|) */
	__shared__ float tdata[Block_Size];/* find max(|cg(:)|) */
    	int tid = threadIdx.x;
    	sdata[tid] = 0.0f;
    	tdata[tid] = 0.0f;
	for(s=0; s<(N+Block_Size-1)/Block_Size; s++)
	{
		id=s*blockDim.x+threadIdx.x;
		a=(id<N)?fabsf(vv[id]):0.0f;
		b=(id<N)?fabsf(cg[id]):0.0f;
		sdata[tid]=MAX(sdata[tid], a);
		tdata[tid]=MAX(tdata[tid], b);
	} 
    	__syncthreads();

    	/* do reduction in shared mem */
    	for(s=blockDim.x/2; s>32; s>>=1) 
    	{
		if (threadIdx.x < s)	{sdata[tid]=MAX(sdata[tid], sdata[tid+s]);tdata[tid]=MAX(tdata[tid], tdata[tid+s]);} 
		__syncthreads();
    	}  
   	if (tid < 32)
   	{
		if (blockDim.x >=  64) { sdata[tid] =MAX(sdata[tid],sdata[tid + 32]);tdata[tid]=MAX(tdata[tid], tdata[tid+32]);}
		if (blockDim.x >=  32) { sdata[tid] =MAX(sdata[tid],sdata[tid + 16]);tdata[tid]=MAX(tdata[tid], tdata[tid+16]);}
		if (blockDim.x >=  16) { sdata[tid] =MAX(sdata[tid],sdata[tid + 8]);tdata[tid]=MAX(tdata[tid], tdata[tid+8]);}
		if (blockDim.x >=   8) { sdata[tid] =MAX(sdata[tid],sdata[tid + 4]);tdata[tid]=MAX(tdata[tid], tdata[tid+4]);}
		if (blockDim.x >=   4) { sdata[tid] =MAX(sdata[tid],sdata[tid + 2]);tdata[tid]=MAX(tdata[tid], tdata[tid+2]);}
		if (blockDim.x >=   2) { sdata[tid] =MAX(sdata[tid],sdata[tid + 1]);tdata[tid]=MAX(tdata[tid], tdata[tid+1]);}
    	}

    	if (tid == 0) { if(fabsf(tdata[0])>EPS) *epsil=0.01*sdata[0]/tdata[0];	else *epsil=0.0f; }
}


__global__ void cuda_sum_alpha12(float *alpha1, float *alpha2, float *p, float *dobs, float *derr, int *gxz, int ng)
/*< calculate the numerator and denominator of alpha
	alpha1: numerator; length=ng
	alpha2: denominator; length=ng >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	float a,b;
	if(id<ng)
	{
		a=dobs[id]+derr[id];/* since f(mk)-dobs[id]=derr[id], thus f(mk)=b+c; */ 
		b=p[gxz[id]]-a;/* f(mk+epsil*cg)-f(mk) */ 
		alpha1[id]-=b*derr[id]; 
		alpha2[id]+=b*b; 
	}
}


__global__ void cuda_cal_alpha(float *alpha, float *alpha1, float *alpha2, float epsil, int ng)
/*< calculate searched stepsize (alpha) according to Taratola's method
configuration requirement: <<<1, Block_Size>>> >*/ 
{
	int s,id;
	float a,b;
  	__shared__ float sdata[Block_Size];
	__shared__ float tdata[Block_Size];
    	int tid=threadIdx.x;
    	sdata[tid]=0.0f;
	tdata[tid]=0.0f;
	for(int s=0; s<(ng+Block_Size-1)/Block_Size; s++)
	{
		id=s*blockDim.x+threadIdx.x;
		a=(id<ng)?alpha1[id]:0.0f;
		b=(id<ng)?alpha2[id]:0.0f;
		sdata[tid] +=a;	
		tdata[tid] +=b;	
	} 
    	__syncthreads();

    	/* do reduction in shared mem */
    	for(s=blockDim.x/2; s>32; s>>=1) 
    	{
		if (threadIdx.x < s) { sdata[tid] += sdata[tid + s];tdata[tid] += tdata[tid + s]; } __syncthreads();
    	}
   	if (tid < 32)
   	{
		if (blockDim.x >=  64) { sdata[tid] += sdata[tid + 32]; tdata[tid] += tdata[tid + 32];}
		if (blockDim.x >=  32) { sdata[tid] += sdata[tid + 16]; tdata[tid] += tdata[tid + 16];}
		if (blockDim.x >=  16) { sdata[tid] += sdata[tid +  8]; tdata[tid] += tdata[tid +  8];}
		if (blockDim.x >=   8) { sdata[tid] += sdata[tid +  4]; tdata[tid] += tdata[tid +  4];}
		if (blockDim.x >=   4) { sdata[tid] += sdata[tid +  2]; tdata[tid] += tdata[tid +  2];}
		if (blockDim.x >=   2) { sdata[tid] += sdata[tid +  1]; tdata[tid] += tdata[tid +  1];}
    	}
     
    	if (tid == 0) *alpha=epsil*sdata[0]/tdata[0];
}

__global__ void cuda_update_vel(float *vv, float *cg, float alpha, int nz, int nx)
/*< update velocity model with obtained stepsize (alpha) >*/
{
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*blockDim.x;
	int id=i1+i2*nz;
	if (i1>=1 && i1<nz-1 && i2>=1 && i2<nx-1) vv[id]=vv[id]+alpha*cg[id];
	__syncthreads();
	// handling the outliers at the boundary
	if (i1==0) 	vv[id]=vv[1+nz*i2];
	else if (i1==nz-1)	vv[id]=vv[nz-2+nz*i2];
	__syncthreads();
	if (i2==0)	vv[id]=vv[i1+nz];
	else if (i2==nx-1)	vv[id]=vv[i1+nz*(nx-2)];
	__syncthreads();
}


