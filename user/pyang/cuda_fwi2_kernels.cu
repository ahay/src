
// set the positions of sources and geophones
__global__ void cuda_set_sg(int *sxz, int sxbeg, int szbeg, int jsx, int jsz, int ns, int npml, int nnz)
{
    	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<ns) sxz[id]=nnz*(sxbeg+id*jsx+npml)+(szbeg+id*jsz+npml);
}

// generate ricker wavelet with time deley
__global__ void cuda_ricker_wavelet(float *wavelet, float fm, float dt, int nt)
{
	int it=threadIdx.x+blockDim.x*blockIdx.x;
    	float tmp = PI*fm*(it*dt-1.0/fm);	//delay the wavelet to exhibit all waveform
    	tmp *=tmp;
    	if (it<nt) wavelet[it]= (1.0-2.0*tmp)*exp(-tmp);	// ricker wavelet at time: t=nt*dt
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


// record the seismogram at time kt
__global__ void cuda_record(float*p, float *seis_kt, int *Gxz, int ng)
{
    	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<ng) seis_kt[id]=p[Gxz[id]];
}


// initialize the absorbing boundary condition (ABC) coefficients along x direction
__global__ void cuda_init_abcx(float *vel, float *bx1, float *bx2, float dx, float dz, float dt, int npml, int nnz, int nnx)
{
	// bx1: left and right PML ABC coefficients, decay p (px,pz) along x direction
	// bx2: left and right PML ABC coefficients, decay v (vx,vz) along x direction
	// only 2 blocks used horizontally, blockIdx.x=0, 1

	// id: position in top or bottom PML zone itself
	// blockIdx.x==0, left PML zone; blockIdx.x==1, right PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=blockIdx.y*npml+threadIdx.y;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+threadIdx.y);
	int ik=i1+nnz*i2;

	float Rc=1.0e-5f;
	float d=npml*MAX(dx,dz);
	float d0=3.0f*vel[id]*logf(Rc)/d/2.0f;
	float tmp1, tmp2;

	if (i2<npml) 	// left PML zone
	{
		tmp1=(float)(npml-i2);
		tmp2=tmp1-0.5f;
	
	}
	else		// right PML zone
	{
		tmp1=i2-npml+0.5f;
		tmp2=(tmp1+0.5f);	
	}
	tmp1=tmp1/npml;
	tmp2=tmp2/npml;
	tmp1=tmp1*tmp1;
	tmp2=tmp2*tmp2;
	bx1[ik]=expf(d0*tmp1*dt);
	bx2[ik]=expf(d0*tmp2*dt);
}


// initialize the absorbing boundary condition (ABC) coefficients along z direction
__global__ void cuda_init_abcz(float *vel, float *bz1, float *bz2, float dx, float dz, float dt, int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone

	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int id=nnz*i2+(blockIdx.x*(nnz-npml)+threadIdx.x);
	int ik=i1+2*npml*i2;

	float Rc=1.0e-5f;
	float d=npml*MAX(dx,dz);
	float d0=3.0f*vel[id]*logf(Rc)/d/2.0f;
	float tmp1, tmp2;

	if (i1<npml) 	// top PML zone
	{	
		tmp1=(float)(npml-i1);
		tmp2=tmp1-0.5f;	
	}
	else		// bottom PML zone
	{
		tmp1=i1-npml+0.5f;
		tmp2=(float)(tmp1+0.5f);	
	}
	tmp1=tmp1/npml;
	tmp2=tmp2/npml;
	tmp1=tmp1*tmp1;
	tmp2=tmp2*tmp2;
	bz1[ik]=expf(d0*tmp1*dt);
	bz2[ik]=expf(d0*tmp2*dt);
}


__global__ void cuda_cal_grad(float *grad, float *sp1, float *rp1, float *power, float _dz2, float _dx2, int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+nnz*i2;

	__shared__ float s_p1[Block_Size1+2][Block_Size2+2];
	s_p1[threadIdx.x+1][threadIdx.y+1]=sp1[id];
	if(threadIdx.x==0)
	{
		if(blockIdx.x>0)	{ s_p1[threadIdx.x][threadIdx.y+1]=sp1[id-1];}
		else			{ s_p1[threadIdx.x][threadIdx.y+1]=0.0f;}
	}
	if(threadIdx.x==blockDim.x-1)
	{
		if(blockIdx.x<gridDim.x-1)	{s_p1[threadIdx.x+2][threadIdx.y+1]=sp1[id+1];}
		else				{s_p1[threadIdx.x+2][threadIdx.y+1]=0.0f;}
	}
	if(threadIdx.y==0)
	{
		if(blockIdx.y>0)	{s_p1[threadIdx.x+1][threadIdx.y]=sp1[id-nnz];}
		else			{s_p1[threadIdx.x+1][threadIdx.y]=0.0f;}
	}
	if(threadIdx.y==blockDim.y-1)
	{
		if(blockIdx.y<gridDim.y-1)	{s_p1[threadIdx.x+1][threadIdx.y+2]=sp1[id+nnz];}
		else				{s_p1[threadIdx.x+1][threadIdx.y+2]=0.0f;}
	}
	__syncthreads();
	if (i1>=npml && i1<nnz-npml && i2>=npml && i2<nnx-npml) 	
	{
		float diff1=s_p1[threadIdx.x+2][threadIdx.y+1]-2.0*s_p1[threadIdx.x+1][threadIdx.y+1]+s_p1[threadIdx.x][threadIdx.y+1];	
		float diff2=s_p1[threadIdx.x+1][threadIdx.y+2]-2.0*s_p1[threadIdx.x+1][threadIdx.y+1]+s_p1[threadIdx.x+1][threadIdx.y];
		float tmp=_dz2*diff1+_dx2*diff2;
		grad[id]+=tmp*rp1[id];	
		power[id]+=tmp*tmp;	
	}
}

__global__ void cuda_scale_grad(float *grad, float *vel, float *power, int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+nnz*i2;
	if (i1>=npml && i1<nnz-npml && i2>=npml && i2<nnx-npml) grad[id]*=2.0/(vel[id]*power[id]);//
	else grad[id]=0;
}



// d_dcal: d_{cal}
// d_dobs: d_{obs}
// d_derr: d_{err}=d_{cal}-d_{obs}
__global__ void cuda_cal_residuals(float *d_dcal, float *d_dobs, float *d_derr, int ng)
{
	int id=blockIdx.x*blockDim.x+threadIdx.x;
	if (id<ng) d_derr[id]=d_dcal[id]-d_dobs[id];
}

// calculate the value of objective function
// d_obj: obj (objective function value)
__global__ void cuda_cal_objective(float *obj,  float *res, int ng)
{
  	__shared__ float  sdata[Block_Size];
    	int tid=threadIdx.x;
    	sdata[tid]=0.0f;
	float a;
	for(int s=0; s<(ng+Block_Size-1)/Block_Size; s++)
	{
		int id=s*blockDim.x+threadIdx.x;
		a=(id<ng)?res[id]:0.0f;
		sdata[tid] += a*a;	
	} 
    	__syncthreads();

    	// do reduction in shared mem
    	for(int s=blockDim.x/2; s>32; s>>=1) 
    	{
		if (threadIdx.x < s) sdata[tid] += sdata[tid + s]; __syncthreads();
    	}
   	if (tid < 32)
   	{
		if (blockDim.x >=  64) { sdata[tid] += sdata[tid + 32]; }
		if (blockDim.x >=  32) { sdata[tid] += sdata[tid + 16]; }
		if (blockDim.x >=  16) { sdata[tid] += sdata[tid +  8]; }
		if (blockDim.x >=   8) { sdata[tid] += sdata[tid +  4]; }
		if (blockDim.x >=   4) { sdata[tid] += sdata[tid +  2]; }
		if (blockDim.x >=   2) { sdata[tid] += sdata[tid +  1]; }
    	}
     
    	if (tid == 0) { *obj=sdata[0]; }
}



__global__ void cuda_cal_beta(float *beta, float *g0, float *g1, float *cg, int N)
/*< calculate beta for nonlinear conjugate gradient algorithm 
configuration requirement: <<<1,Block_Size>>> >*/
{
    	__shared__ float sdata[Block_Size];
	__shared__ float tdata[Block_Size];
	__shared__ float rdata[Block_Size];
    	int tid = threadIdx.x;
    	sdata[tid] = 0.0f;
	tdata[tid] = 0.0f;
	rdata[tid] = 0.0f;
	for(int s=0; s<(N+Block_Size-1)/Block_Size; s++)
	{
		int id=s*blockDim.x+threadIdx.x;
		float a=(id<N)?g0[id]:0.0f;
		float b=(id<N)?g1[id]:0.0f;
		float c=(id<N)?cg[id]:0.0f;

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
    	for(int s=blockDim.x/2; s>32; s>>=1) 
    	{
		if (threadIdx.x < s)	{ sdata[tid]+=sdata[tid+s]; tdata[tid]+=tdata[tid+s]; rdata[tid]+=rdata[tid+s];}
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
		float beta_HS=0;
		float beta_DY=0;
		if(fabsf(tdata[0])>EPS) 
		{
			beta_HS=sdata[0]/tdata[0]; 
			beta_DY=rdata[0]/tdata[0];
		} 
		*beta=MAX(0, MIN(beta_HS, beta_DY));// Hybrid HS-DY method incorprating iteration restart
	}	
}


__global__ void cuda_cal_conjgrad(float *g1, float *cg, float beta, int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+i2*nnz;

	if (i1>=npml && i1<nnz-npml && i2>=npml && i2<nnx-npml) cg[id]=-g1[id]+beta*cg[id];
	else cg[id]=0;
}


__global__ void cuda_cal_vtmp(float *vtmp, float *vv, float *cg, float epsil, int npml, int nnz, int nnx)
/*< calculate temporary velocity >*/ 
{
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*blockDim.x;
	int id=i1+i2*nnz;
	if (i1>=npml && i1<nnz-npml && i2>=npml && i2<nnx-npml) vtmp[id]=vv[id]+epsil*cg[id];
	__syncthreads();
	if (i1<npml) vtmp[id]=vtmp[npml+i2*nnz];
	else if (i1>=nnz-npml) vtmp[id]=vtmp[nnz-npml-1+i2*nnz];
	__syncthreads();
	if (i2<npml) vtmp[id]=vtmp[i1+nnz*npml];
	else if (i2>nnx-npml) vtmp[id]=vtmp[i1+nnz*(nnx-npml-1)];
}


__global__ void cuda_update_vel(float *vv, float *cg, float alpha, int npml, int nnz, int nnx)
/*< update velocity model with obtained stepsize (alpha) >*/
{
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*blockDim.x;
	int id=i1+i2*nnz;
	if (i1>=npml && i1<nnz-npml && i2>=npml && i2<nnx-npml) vv[id]=vv[id]+alpha*cg[id];
	__syncthreads();
	if (i1<npml) vv[id]=vv[npml+i2*nnz];
	else if (i1>=nnz-npml) vv[id]=vv[nnz-npml-1+i2*nnz];
	__syncthreads();
	if (i2<npml) vv[id]=vv[i1+nnz*npml];
	else if (i2>nnx-npml) vv[id]=vv[i1+nnz*(nnx-npml-1)];
}

__global__ void cuda_cal_epsilon(float *vv, float *cg, float *epsil, int N)
/*< calculate estimated stepsize (epsil) according to Taratola's method
configuration requirement: <<<1, Block_Size>>> >*/ 
{
    	__shared__ float sdata[Block_Size];// find max(|vv(:)|)
	__shared__ float tdata[Block_Size];// find max(|cg(:)|)
    	int tid = threadIdx.x;
    	sdata[tid] = 0.0f;
    	tdata[tid] = 0.0f;
	for(int s=0; s<(N+Block_Size-1)/Block_Size; s++)
	{
		int id=s*blockDim.x+threadIdx.x;
		float a=(id<N)?fabsf(vv[id]):0.0f;
		float b=(id<N)?fabsf(cg[id]):0.0f;
		sdata[tid]= MAX(sdata[tid], a);
		tdata[tid]= MAX(tdata[tid], b);
	} 
    	__syncthreads();

    	// do reduction in shared mem
    	for(int s=blockDim.x/2; s>32; s>>=1) 
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


// calculate the numerator and denominator of alpha
// alpha1: numerator; length=ng
// alpha2: denominator; length=ng
__global__ void cuda_sum_alpha12(float *alpha1, float *alpha2, float *dcaltmp, float *dobs, float *derr, int ng)
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	float a=(id<ng)?dcaltmp[id]:0.0f;
	float b=(id<ng)?dobs[id]:0.0f;
	float c=(id<ng)?derr[id]:0.0f;
	float d=b+c;// dcal[id]-dobs[id]=derr[id]; here d=dcal[id]=b+c;
	if(id<ng) { alpha1[id]-=(a-d)*c; alpha2[id]+=(a-d)*(a-d); }
}

//<<<1,Block_Size>>>
__global__ void cuda_cal_alpha(float *alpha, float *alpha1, float *alpha2, float epsilon, int ng)
{
  	__shared__ float sdata[Block_Size];
	__shared__ float tdata[Block_Size];
    	int tid=threadIdx.x;
    	sdata[tid]=0.0f;
	tdata[tid]=0.0f;
	for(int s=0; s<(ng+Block_Size-1)/Block_Size; s++)
	{
		int id=s*blockDim.x+threadIdx.x;
		float a=(id<ng)?alpha1[id]:0.0f;
		float b=(id<ng)?alpha2[id]:0.0f;
		sdata[tid] +=a;	
		tdata[tid] +=b;	
	} 
    	__syncthreads();

    	// do reduction in shared mem
    	for(int s=blockDim.x/2; s>32; s>>=1) 
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
     
    	if (tid == 0) { if(fabs(tdata[0])>1e-15f) *alpha=epsilon*sdata[0]/tdata[0];	else *alpha=0.0f;}
}


///////////////////////////////////////////// NJ=2 //////////////////////////////////////////////////
__global__ void cuda_forward_v_2(float *p, float *vx, float *vz, float _dx, float _dz, 
	int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+i2*nnz;

	__shared__ float s_p[Block_Size1+1][Block_Size2+1];
    	s_p[threadIdx.x][threadIdx.y]=p[id];
	if (threadIdx.x>blockDim.x-2)
	{
		if (blockIdx.x<gridDim.x-1)	s_p[threadIdx.x+1][threadIdx.y]=p[id+1];
		else				s_p[threadIdx.x+1][threadIdx.y]=0.0f;
	} 
	if (threadIdx.y>blockDim.y-2)
	{
		if (blockIdx.y<gridDim.y-1)	s_p[threadIdx.x][threadIdx.y+1]=p[id+nnz];
		else				s_p[threadIdx.x][threadIdx.y+1]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml))
	{
		float diff1=(s_p[threadIdx.x+1][threadIdx.y]-s_p[threadIdx.x][threadIdx.y]);// .x-->1st dim--> i1
		float diff2=(s_p[threadIdx.x][threadIdx.y+1]-s_p[threadIdx.x][threadIdx.y]);// .y-->2nd dim--> i2
		vz[id]=_dz*diff1;
		vx[id]=_dx*diff2;
	}
}

__global__ void cuda_PML_vz_2(float *p, float *convpz, float *bz, float *vz, float _dz, 
	int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0,1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int ik=i1+2*npml*i2;
	int id=(blockIdx.x*(nnz-npml)+threadIdx.x)+nnz*i2;

	__shared__ float s_p[33][Block_Size2];	
    	s_p[threadIdx.x][threadIdx.y]=p[id]; 
	if (threadIdx.x>30)
	{
		if (blockIdx.x<gridDim.x-1)	s_p[threadIdx.x+1][threadIdx.y]=p[id+1];
		else				s_p[threadIdx.x+1][threadIdx.y]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml)) 
	{
		float diff1=(s_p[threadIdx.x+1][threadIdx.y]-s_p[threadIdx.x][threadIdx.y]);
		convpz[ik]=bz[ik]*convpz[ik]+(bz[ik]-1.0f)*_dz*diff1;	
		vz[id]+=convpz[ik];
	}
}

__global__ void cuda_PML_vx_2(float *p,  float *convpx, float *bx, float *vx, float _dx, 
	int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*npml;
	int ik=i1+nnz*i2;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+	threadIdx.y);

	__shared__ float s_p[Block_Size1][33];
    	s_p[threadIdx.x][threadIdx.y]=p[id]; 
	if (threadIdx.y>30)
	{
		if (blockIdx.y<gridDim.y-1)	s_p[threadIdx.x][threadIdx.y+1]=p[id+nnz];
		else				s_p[threadIdx.x][threadIdx.y+1]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml))
	{
		float diff2=(s_p[threadIdx.x][threadIdx.y+1]-s_p[threadIdx.x][threadIdx.y]);
		convpx[ik]=bx[ik]*convpx[ik]+(bx[ik]-1.0f)*_dx*diff2;	
		vx[id]+=convpx[ik];
	}
}

__global__ void cuda_forward_p_2(float *vel, float *p0, float *p1, float *vx, float *vz, float dt, float _dx, float _dz, 
	int npml, int nnz, int nnx)
{
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int id=i1+i2*nnz;

	__shared__ float s_v1[Block_Size1+1][Block_Size2];
	__shared__ float s_v2[Block_Size1][Block_Size2+1];
    	s_v1[threadIdx.x+1][threadIdx.y]=vz[id]; 
    	s_v2[threadIdx.x][threadIdx.y+1]=vx[id]; 
	if (threadIdx.x<1)
	{
		if (blockIdx.x) 		s_v1[threadIdx.x][threadIdx.y]=vz[id-1];
		else				s_v1[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y<1)
	{
		if (blockIdx.y)			s_v2[threadIdx.x][threadIdx.y]=vx[id-nnz];
		else 				s_v2[threadIdx.x][threadIdx.y]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml))
	{
		float diff1=(s_v1[threadIdx.x+1][threadIdx.y]-s_v1[threadIdx.x][threadIdx.y]);
		float diff2=(s_v2[threadIdx.x][threadIdx.y+1]-s_v2[threadIdx.x][threadIdx.y]);
		p0[id]=2*p1[id]-p0[id]+dt*dt*vel[id]*vel[id]*(_dz*diff1+_dx*diff2);
	}
}



__global__ void cuda_PML_pz_2(float *vel, float *p0,  float *convvz, float *bz, float *vz, float dt, float _dz, 
	int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int ik=i1+2*npml*i2;
	int id=(blockIdx.x*(nnz-npml)+	threadIdx.x)+	nnz*	i2;

	__shared__ float s_v1[33][Block_Size2];
    	s_v1[threadIdx.x+1][threadIdx.y]=vz[id]; 
	if (threadIdx.x<1)
	{
		if (blockIdx.x)			s_v1[threadIdx.x][threadIdx.y]=vz[id-1];
		else 				s_v1[threadIdx.x][threadIdx.y]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml))
	{
		float diff1=(s_v1[threadIdx.x+1][threadIdx.y]-s_v1[threadIdx.x][threadIdx.y]);
		convvz[ik]=bz[ik]*convvz[ik]+(bz[ik]-1.0f)*_dz*diff1;	
		p0[id]+=dt*dt*vel[id]*vel[id]*convvz[ik];
	}
}

__global__ void cuda_PML_px_2(float *vel, float *p0,  float *convvx, float *bx, float *vx, float dt, float _dx, 
	int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*npml;
	int ik=i1+nnz*i2;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+threadIdx.y);

	__shared__ float s_v2[Block_Size1][33];
    	s_v2[threadIdx.x][threadIdx.y+1]=vx[id]; 
	if (threadIdx.y<1)
	{
		if (blockIdx.y) 		s_v2[threadIdx.x][threadIdx.y]=vx[id-nnz];
		else				s_v2[threadIdx.x][threadIdx.y]=0.0f;
	}	
	__syncthreads();

	if (!( i1<npml))
	{
		float diff2=(s_v2[threadIdx.x][threadIdx.y+1]-s_v2[threadIdx.x][threadIdx.y]);
		convvx[ik]=bx[ik]*convvx[ik]+(bx[ik]-1.0f)*_dx*diff2;	
		p0[id]+=dt*dt*vel[id]*vel[id]*convvx[ik];
	}
}


////////////////////////////////////////// NJ=4 ////////////////////////////////////////////
__global__ void cuda_forward_v_4(float *p, float *vx, float *vz, float _dx, float _dz, 
	int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+i2*nnz;

	__shared__ float s_p[Block_Size1+3][Block_Size2+3];
    	s_p[threadIdx.x+1][threadIdx.y+1]=p[id]; 
	if (threadIdx.x<1)
	{
		if (blockIdx.x)			s_p[threadIdx.x][threadIdx.y+1]=p[id-1];
		else 				s_p[threadIdx.x][threadIdx.y+1]=0.0f;
	}
	if (threadIdx.x>blockDim.x-3)
	{
		if (blockIdx.x<gridDim.x-1)	s_p[threadIdx.x+3][threadIdx.y+1]=p[id+2];
		else				s_p[threadIdx.x+3][threadIdx.y+1]=0.0f;
	}
	if (threadIdx.y<1)
	{
		if (blockIdx.y) 		s_p[threadIdx.x+1][threadIdx.y]=p[id-nnz];
		else				s_p[threadIdx.x+1][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>blockDim.y-3)
	{
		if (blockIdx.y<gridDim.y-1)	s_p[threadIdx.x+1][threadIdx.y+3]=p[id+2*nnz];
		else				s_p[threadIdx.x+1][threadIdx.y+3]=0.0f;
	}
	__syncthreads();   
  
	if (!( i1<npml))
	{
		float diff1=1.125f*(s_p[threadIdx.x+2][threadIdx.y+1]-s_p[threadIdx.x+1][threadIdx.y+1])
			-0.041666666666667f*(s_p[threadIdx.x+3][threadIdx.y+1]-s_p[threadIdx.x][threadIdx.y+1]);
		float diff2=1.125f*(s_p[threadIdx.x+1][threadIdx.y+2]-s_p[threadIdx.x+1][threadIdx.y+1])
			-0.041666666666667f*(s_p[threadIdx.x+1][threadIdx.y+3]-s_p[threadIdx.x+1][threadIdx.y]);

		vz[id]=_dz*diff1;
		vx[id]=_dx*diff2;
	}
}

__global__ void cuda_PML_vz_4(float *p,  float *convpz, float *bz, float *vz, float _dz, 
	int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int ik=i1+2*npml*i2;
	int id=blockIdx.x*(nnz-npml)+threadIdx.x+nnz*i2;

	__shared__ float s_p[35][Block_Size2];
    	s_p[threadIdx.x+1][threadIdx.y]=p[id]; 
	if (threadIdx.x<1)
	{
		if (blockIdx.x)			s_p[threadIdx.x][threadIdx.y]=p[id-1];
		else 				s_p[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>29)
	{
		if (blockIdx.x<gridDim.x-1)	s_p[threadIdx.x+3][threadIdx.y]=p[id+2];
		else				s_p[threadIdx.x+3][threadIdx.y]=0.0f;
	}
	__syncthreads();


	if (!( i1<npml)) 
	{
		float diff1=1.125f*(s_p[threadIdx.x+2][threadIdx.y]-s_p[threadIdx.x+1][threadIdx.y])
			-0.041666666666667f*(s_p[threadIdx.x+3][threadIdx.y]-s_p[threadIdx.x][threadIdx.y]);
		convpz[ik]=bz[ik]*convpz[ik]+(bz[ik]-1.0f)*_dz*diff1;	
		vz[id]+=convpz[ik];
	}
}
__global__ void cuda_PML_vx_4(float *p,  float *convpx, float *bx, float *vx, float _dx, 
	int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*npml;
	int ik=i1+i2*nnz;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+threadIdx.y);

	__shared__ float s_p[Block_Size1][35];// npml+3=35; Block_SizeX=32; Block_SizeY=8;
    	s_p[threadIdx.x][threadIdx.y+1]=p[id]; 
	if (threadIdx.y<1)
	{
		if (blockIdx.y) 		s_p[threadIdx.x][threadIdx.y]=p[id-nnz];
		else				s_p[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>29)
	{
		if (blockIdx.y<gridDim.y-1)	s_p[threadIdx.x][threadIdx.y+3]=p[id+2*nnz];
		else				s_p[threadIdx.x][threadIdx.y+3]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml))
	{
		float diff2=1.125f*(s_p[threadIdx.x][threadIdx.y+2]-s_p[threadIdx.x][threadIdx.y+1])
			-0.041666666666667f*(s_p[threadIdx.x][threadIdx.y+3]-s_p[threadIdx.x][threadIdx.y]);
		convpx[ik]=bx[ik]*convpx[ik]+(bx[ik]-1.0f)*_dx*diff2;	
		vx[id]+=convpx[ik];
	}
}

__global__ void cuda_forward_p_4(float *vel, float *p0, float *p1, float *vx, float *vz, float dt, float _dx, float _dz, 
	int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+i2*nnz;

	__shared__ float s_vx[Block_Size1][Block_Size2+3];
	__shared__ float s_vz[Block_Size1+3][Block_Size2];
    	s_vx[threadIdx.x][threadIdx.y+2]=vx[id]; 
    	s_vz[threadIdx.x+2][threadIdx.y]=vz[id]; 

	if (threadIdx.x<2)
	{
		if (blockIdx.x)			s_vz[threadIdx.x][threadIdx.y]=vz[id-2];
		else 				s_vz[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>blockDim.x-2)
	{
		if (blockIdx.x<gridDim.x-1)	s_vz[threadIdx.x+3][threadIdx.y]=vz[id+1];
		else				s_vz[threadIdx.x+3][threadIdx.y]=0.0f;
	}
	if (threadIdx.y<2)
	{
		if (blockIdx.y) 		s_vx[threadIdx.x][threadIdx.y]=vx[id-2*nnz];
		else				s_vx[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>blockDim.y-2)
	{
		if (blockIdx.y<gridDim.y-1)	s_vx[threadIdx.x][threadIdx.y+3]=vx[id+nnz];
		else				s_vx[threadIdx.x][threadIdx.y+3]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml))
	{
		float diff2=1.125f*(s_vx[threadIdx.x][threadIdx.y+2]-s_vx[threadIdx.x][threadIdx.y+1])
			 -0.041666666666667f*(s_vx[threadIdx.x][threadIdx.y+3]-s_vx[threadIdx.x][threadIdx.y]);
		float diff1=1.125f*(s_vz[threadIdx.x+2][threadIdx.y]-s_vz[threadIdx.x+1][threadIdx.y])+
			 -0.041666666666667f*(s_vz[threadIdx.x+3][threadIdx.y]-s_vz[threadIdx.x][threadIdx.y]);
		p0[id]=2.0f*p1[id]-p0[id]+dt*dt*vel[id]*vel[id]*(_dz*diff1+_dx*diff2);
	}
}

__global__ void cuda_PML_pz_4(float *vel, float *p,  float *convvz, float *bz, float *vz, float dt, float _dz, 
	int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int ik=i1+2*npml*i2;
	int id=blockIdx.x*(nnz-npml)+threadIdx.x+nnz*i2;

	__shared__ float s_vz[35][Block_Size2];
    	s_vz[threadIdx.x+2][threadIdx.y]=vz[id]; 
	if (threadIdx.x<2)
	{
		if (blockIdx.x)			s_vz[threadIdx.x][threadIdx.y]=vz[id-2];
		else 				s_vz[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>30)
	{
		if (blockIdx.x<gridDim.x-1)	s_vz[threadIdx.x+3][threadIdx.y]=vz[id+1];
		else				s_vz[threadIdx.x+3][threadIdx.y]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml))
	{
		float diff1=1.125f*(s_vz[threadIdx.x+2][threadIdx.y]-s_vz[threadIdx.x+1][threadIdx.y])+
			-0.041666666666667f*(s_vz[threadIdx.x+3][threadIdx.y]-s_vz[threadIdx.x][threadIdx.y]);
		convvz[ik]=bz[ik]*convvz[ik]+(bz[ik]-1.0f)*_dz*diff1;	
		p[id]+=dt*dt*vel[id]*vel[id]*convvz[ik];
	}
}
__global__ void cuda_PML_px_4(float *vel, float *p,  float *convvx, float *bx, float *vx, float dt, float _dx, 
	int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*npml;
	int ik=i1+i2*nnz;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+threadIdx.y);

	__shared__ float s_vx[Block_Size1][35];
    	s_vx[threadIdx.x][threadIdx.y+2]=vx[id]; 
	if (threadIdx.y<2)
	{
		if (blockIdx.y) 		s_vx[threadIdx.x][threadIdx.y]=vx[id-2*nnz];
		else				s_vx[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>30)
	{
		if (blockIdx.y<gridDim.y-1)	s_vx[threadIdx.x][threadIdx.y+3]=vx[id+nnz];
		else				s_vx[threadIdx.x][threadIdx.y+3]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml))
	{
		float diff2=1.125f*(s_vx[threadIdx.x][threadIdx.y+2]-s_vx[threadIdx.x][threadIdx.y+1])
			-0.041666666666667f*(s_vx[threadIdx.x][threadIdx.y+3]-s_vx[threadIdx.x][threadIdx.y]);
		convvx[ik]=bx[ik]*convvx[ik]+(bx[ik]-1.0f)*_dx*diff2;	
		p[id]+=dt*dt*vel[id]*vel[id]*convvx[ik];
	}
}


////////////////////////////////////////// NJ=6 ////////////////////////////////////////
__global__ void cuda_forward_v_6(float *p, float *vx, float *vz, float _dx, float _dz, 
	int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+i2*nnz;

	__shared__ float s_p[Block_Size1+5][Block_Size2+5];
    	s_p[threadIdx.x+2][threadIdx.y+2]=p[id]; 
	if (threadIdx.x<2)
	{
		if (blockIdx.x)			s_p[threadIdx.x][threadIdx.y+2]=p[id-2];
		else 				s_p[threadIdx.x][threadIdx.y+2]=0.0f;
	}
	if (threadIdx.x>blockDim.x-4)
	{
		if (blockIdx.x<gridDim.x-1)	s_p[threadIdx.x+5][threadIdx.y+2]=p[id+3];
		else				s_p[threadIdx.x+5][threadIdx.y+2]=0.0f;
	}
	if (threadIdx.y<2)
	{
		if (blockIdx.y) 		s_p[threadIdx.x+2][threadIdx.y]=p[id-2*nnz];
		else				s_p[threadIdx.x+2][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>blockDim.y-4)
	{
		if (blockIdx.y<gridDim.y-1)	s_p[threadIdx.x+2][threadIdx.y+5]=p[id+3*nnz];
		else				s_p[threadIdx.x+2][threadIdx.y+5]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml))
	{
		float diff1=1.171875f*(s_p[threadIdx.x+3][threadIdx.y+2]-s_p[threadIdx.x+2][threadIdx.y+2])
			-0.065104166666667f*(s_p[threadIdx.x+4][threadIdx.y+2]-s_p[threadIdx.x+1][threadIdx.y+2])
			+0.0046875f*(s_p[threadIdx.x+5][threadIdx.y+2]-s_p[threadIdx.x][threadIdx.y+2]);
		float diff2=1.171875f*(s_p[threadIdx.x+2][threadIdx.y+3]-s_p[threadIdx.x+2][threadIdx.y+2])
			-0.065104166666667f*(s_p[threadIdx.x+2][threadIdx.y+4]-s_p[threadIdx.x+2][threadIdx.y+1])
			+0.0046875f*(s_p[threadIdx.x+2][threadIdx.y+5]-s_p[threadIdx.x+2][threadIdx.y]);

		vz[id]=_dz*diff1;
		vx[id]=_dx*diff2;
	}
}

__global__ void cuda_PML_vz_6(float *p,  float *convpz, float *bz, float *vz, float _dz, 
	int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int ik=i1+2*npml*i2;
	int id=blockIdx.x*(nnz-npml)+threadIdx.x+nnz*i2;

	__shared__ float s_p[37][Block_Size2];
    	s_p[threadIdx.x+2][threadIdx.y]=p[id]; 
	if (threadIdx.x<2)
	{
		if (blockIdx.x)			s_p[threadIdx.x][threadIdx.y]=p[id-2];
		else 				s_p[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>28)
	{
		if (blockIdx.x<gridDim.x-1)	s_p[threadIdx.x+5][threadIdx.y]=p[id+3];
		else				s_p[threadIdx.x+5][threadIdx.y]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml)) 
	{
		float diff1=1.171875f*(s_p[threadIdx.x+3][threadIdx.y]-s_p[threadIdx.x+2][threadIdx.y])
			-0.065104166666667f*(s_p[threadIdx.x+4][threadIdx.y]-s_p[threadIdx.x+1][threadIdx.y])
			+0.0046875f*(s_p[threadIdx.x+5][threadIdx.y]-s_p[threadIdx.x][threadIdx.y]);
		convpz[ik]=bz[ik]*convpz[ik]+(bz[ik]-1.0f)*_dz*diff1;	
		vz[id]+=convpz[ik];
	}
}
__global__ void cuda_PML_vx_6(float *p,  float *convpx, float *bx, float *vx, float _dx,
	int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*npml;
	int ik=i1+i2*nnz;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+threadIdx.y);

	__shared__ float s_p[Block_Size1][37];
    	s_p[threadIdx.x][threadIdx.y+2]=p[id]; 
	if (threadIdx.y<2)
	{
		if (blockIdx.y) 		s_p[threadIdx.x][threadIdx.y]=p[id-2*nnz];
		else				s_p[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>28)
	{
		if (blockIdx.y<gridDim.y-1)	s_p[threadIdx.x][threadIdx.y+5]=p[id+3*nnz];
		else				s_p[threadIdx.x][threadIdx.y+5]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml))
	{
		float diff2=1.171875f*(s_p[threadIdx.x][threadIdx.y+3]-s_p[threadIdx.x][threadIdx.y+2])
			-0.065104166666667f*(s_p[threadIdx.x][threadIdx.y+4]-s_p[threadIdx.x][threadIdx.y+1])
			+0.0046875f*(s_p[threadIdx.x][threadIdx.y+5]-s_p[threadIdx.x][threadIdx.y]);
		convpx[ik]=bx[ik]*convpx[ik]+(bx[ik]-1.0f)*_dx*diff2;	
		vx[id]+=convpx[ik];
	}
}
__global__ void cuda_forward_p_6(float *vel, float *p0, float *p1, float *vx, float *vz, float dt, float _dx, float _dz, 
	int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+i2*nnz;

	__shared__ float s_vx[Block_Size1][Block_Size2+5];
	__shared__ float s_vz[Block_Size1+5][Block_Size2];
    	s_vx[threadIdx.x][threadIdx.y+3]=vx[id]; 
    	s_vz[threadIdx.x+3][threadIdx.y]=vz[id]; 

	if (threadIdx.x<3)
	{
		if (blockIdx.x)			s_vz[threadIdx.x][threadIdx.y]=vz[id-3];
		else 				s_vz[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>blockDim.x-3)
	{
		if (blockIdx.x<gridDim.x-1)	s_vz[threadIdx.x+5][threadIdx.y]=vz[id+2];
		else				s_vz[threadIdx.x+5][threadIdx.y]=0.0f;
	}
	if (threadIdx.y<3)
	{
		if (blockIdx.y) 		s_vx[threadIdx.x][threadIdx.y]=vx[id-3*nnz];
		else				s_vx[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>blockDim.y-3)
	{
		if (blockIdx.y<gridDim.y-1)	s_vx[threadIdx.x][threadIdx.y+5]=vx[id+2*nnz];
		else				s_vx[threadIdx.x][threadIdx.y+5]=0.0f;
	}
	__syncthreads();


	if (!( i1<npml))
	{
		float diff2=1.171875f*(s_vx[threadIdx.x][threadIdx.y+3]-s_vx[threadIdx.x][threadIdx.y+2])
			-0.065104166666667f*(s_vx[threadIdx.x][threadIdx.y+4]-s_vx[threadIdx.x][threadIdx.y+1])+
			0.0046875f*(s_vx[threadIdx.x][threadIdx.y+5]-s_vx[threadIdx.x][threadIdx.y]);
		float diff1=1.171875f*(s_vz[threadIdx.x+3][threadIdx.y]-s_vz[threadIdx.x+2][threadIdx.y])+
			-0.065104166666667f*(s_vz[threadIdx.x+4][threadIdx.y]-s_vz[threadIdx.x+1][threadIdx.y])+
			0.0046875f*(s_vz[threadIdx.x+5][threadIdx.y]-s_vz[threadIdx.x][threadIdx.y]);
		p0[id]=2.0f*p1[id]-p0[id]+dt*dt*vel[id]*vel[id]*(_dz*diff1+_dx*diff2);
	}
}

__global__ void cuda_PML_pz_6(float *vel, float *p0,  float *convvz, float *bz, float *vz, float dt, float _dz,
	int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int ik=i1+2*npml*i2;
	int id=blockIdx.x*(nnz-npml)+threadIdx.x+nnz*i2;

	__shared__ float s_vz[37][Block_Size2];// npml+5=37; Block_SizeX=32; Block_SizeY=8;
    	s_vz[threadIdx.x+3][threadIdx.y]=vz[id]; 
	if (threadIdx.x<3)
	{
		if (blockIdx.x)			s_vz[threadIdx.x][threadIdx.y]=vz[id-3];
		else 				s_vz[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>29)
	{
		if (blockIdx.x<gridDim.x-1)	s_vz[threadIdx.x+5][threadIdx.y]=vz[id+2];
		else				s_vz[threadIdx.x+5][threadIdx.y]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml))
	{
		float diff1=1.171875f*(s_vz[threadIdx.x+3][threadIdx.y]-s_vz[threadIdx.x+2][threadIdx.y])
			-0.065104166666667f*(s_vz[threadIdx.x+4][threadIdx.y]-s_vz[threadIdx.x+1][threadIdx.y])
			+0.0046875f*(s_vz[threadIdx.x+5][threadIdx.y]-s_vz[threadIdx.x][threadIdx.y]);
		convvz[ik]=bz[ik]*convvz[ik]+(bz[ik]-1.0f)*_dz*diff1;	
		p0[id]+=dt*dt*vel[id]*vel[id]*convvz[ik];
	}
}
__global__ void cuda_PML_px_6(float *vel, float *p,  float *convvx, float *bx, float *vx, float dt, float _dx, 
	int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*npml;
	int ik=i1+i2*nnz;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+threadIdx.y);

	__shared__ float s_vx[Block_Size1][37];
    	s_vx[threadIdx.x][threadIdx.y+3]=vx[id]; 
	if (threadIdx.y<3)
	{
		if (blockIdx.y) 		s_vx[threadIdx.x][threadIdx.y]=vx[id-3*nnz];
		else				s_vx[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>29)
	{
		if (blockIdx.y<gridDim.y-1)	s_vx[threadIdx.x][threadIdx.y+5]=vx[id+2*nnz];
		else				s_vx[threadIdx.x][threadIdx.y+5]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml))
	{
		float diff2=1.171875f*(s_vx[threadIdx.x][threadIdx.y+3]-s_vx[threadIdx.x][threadIdx.y+2])
			-0.065104166666667f*(s_vx[threadIdx.x][threadIdx.y+4]-s_vx[threadIdx.x][threadIdx.y+1])
			+0.0046875f*(s_vx[threadIdx.x][threadIdx.y+5]-s_vx[threadIdx.x][threadIdx.y]);
		convvx[ik]=bx[ik]*convvx[ik]+(bx[ik]-1.0f)*_dx*diff2;	
		p[id]+=dt*dt*vel[id]*vel[id]*convvx[ik];
	}
}


//////////////////////////////////////////// NJ=8 ////////////////////////////////////
__global__ void cuda_forward_v_8(float *p, float *vx, float *vz, float _dx, float _dz, 
	int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+nnz*i2;

	__shared__ float s_p[Block_Size1+7][Block_Size2+7];
    	s_p[threadIdx.x+3][threadIdx.y+3]=p[id]; 
	if (threadIdx.x<3)
	{
		if (blockIdx.x)			s_p[threadIdx.x][threadIdx.y+3]=p[id-3];
		else 				s_p[threadIdx.x][threadIdx.y+3]=0.0f;
	}
	if (threadIdx.x>blockDim.x-5)
	{
		if (blockIdx.x<gridDim.x-1)	s_p[threadIdx.x+7][threadIdx.y+3]=p[id+4];
		else				s_p[threadIdx.x+7][threadIdx.y+3]=0.0f;
	}
	if (threadIdx.y<3)
	{
		if (blockIdx.y) 		s_p[threadIdx.x+3][threadIdx.y]=p[id-3*nnz];
		else				s_p[threadIdx.x+3][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>blockDim.y-5)
	{
		if (blockIdx.y<gridDim.y-1)	s_p[threadIdx.x+3][threadIdx.y+7]=p[id+4*nnz];
		else				s_p[threadIdx.x+3][threadIdx.y+7]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml))
	{
		float diff1=1.1962890625000f*(s_p[threadIdx.x+4][threadIdx.y+3]-s_p[threadIdx.x+3][threadIdx.y+3])
			-0.0797526041667f*(s_p[threadIdx.x+5][threadIdx.y+3]-s_p[threadIdx.x+2][threadIdx.y+3])
			+0.0095703125000f*(s_p[threadIdx.x+6][threadIdx.y+3]-s_p[threadIdx.x+1][threadIdx.y+3])
			-0.0006975446429f*(s_p[threadIdx.x+7][threadIdx.y+3]-s_p[threadIdx.x][threadIdx.y+3]);
		float diff2=1.1962890625000f*(s_p[threadIdx.x+3][threadIdx.y+4]-s_p[threadIdx.x+3][threadIdx.y+3])
			-0.0797526041667f*(s_p[threadIdx.x+3][threadIdx.y+5]-s_p[threadIdx.x+3][threadIdx.y+2])
			+0.0095703125000f*(s_p[threadIdx.x+3][threadIdx.y+6]-s_p[threadIdx.x+3][threadIdx.y+1])
			-0.0006975446429f*(s_p[threadIdx.x+3][threadIdx.y+7]-s_p[threadIdx.x+3][threadIdx.y]);

		vz[id]=_dz*diff1;
		vx[id]=_dx*diff2;
	}
}


__global__ void cuda_PML_vz_8(float *p,  float *convpz, float *bz, float *vz, float _dz, 
	int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int ik=i1+2*npml*i2;
	int id=blockIdx.x*(nnz-npml)+threadIdx.x+nnz*i2;

	__shared__ float s_p[39][Block_Size2];
    	s_p[threadIdx.x+3][threadIdx.y]=p[id]; 
	if (threadIdx.x<3)
	{
		if (blockIdx.x)			s_p[threadIdx.x][threadIdx.y]=p[id-3];
		else 				s_p[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>27)
	{
		if (blockIdx.x<gridDim.x-1)	s_p[threadIdx.x+7][threadIdx.y]=p[id+4];
		else				s_p[threadIdx.x+7][threadIdx.y]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml)) 
	{
		float diff1=1.1962890625000f*(s_p[threadIdx.x+4][threadIdx.y]-s_p[threadIdx.x+3][threadIdx.y])
			-0.0797526041667f*(s_p[threadIdx.x+5][threadIdx.y]-s_p[threadIdx.x+2][threadIdx.y])
			+0.0095703125000f*(s_p[threadIdx.x+6][threadIdx.y]-s_p[threadIdx.x+1][threadIdx.y])
			-0.0006975446429f*(s_p[threadIdx.x+7][threadIdx.y]-s_p[threadIdx.x][threadIdx.y]);
		convpz[ik]=bz[ik]*convpz[ik]+(bz[ik]-1.0f)*_dz*diff1;	
		vz[id]+=convpz[ik];
	}
}
__global__ void cuda_PML_vx_8(float *p,  float *convpx, float *bx, float *vx, float _dx, 
	int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*npml;
	int ik=i1+i2*nnz;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+threadIdx.y);

	__shared__ float s_p[Block_Size1][39];// npml+7=39; Block_SizeX=32; Block_SizeY=8;
    	s_p[threadIdx.x][threadIdx.y+3]=p[id]; 
	if (threadIdx.y<3)
	{
		if (blockIdx.y) 		s_p[threadIdx.x][threadIdx.y]=p[id-3*nnz];
		else				s_p[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>27)//32-4
	{
		if (blockIdx.y<gridDim.y-1)	s_p[threadIdx.x][threadIdx.y+7]=p[id+4*nnz];
		else				s_p[threadIdx.x][threadIdx.y+7]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml))
	{
		float diff2=1.1962890625000f*(s_p[threadIdx.x][threadIdx.y+4]-s_p[threadIdx.x][threadIdx.y+3])
			-0.0797526041667f*(s_p[threadIdx.x][threadIdx.y+5]-s_p[threadIdx.x][threadIdx.y+2])+
			0.0095703125000f*(s_p[threadIdx.x][threadIdx.y+6]-s_p[threadIdx.x][threadIdx.y+1])+
			-0.0006975446429f*(s_p[threadIdx.x][threadIdx.y+7]-s_p[threadIdx.x][threadIdx.y]);
		convpx[ik]=bx[ik]*convpx[ik]+(bx[ik]-1.0f)*_dx*diff2;	
		vx[id]+=convpx[ik];
	}
}

__global__ void cuda_forward_p_8(float *vel, float *p0, float *p1, float *vx, float *vz, float dt, float _dx, float _dz, 
	int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+i2*nnz;

	__shared__ float s_vx[Block_Size1][Block_Size2+7];
	__shared__ float s_vz[Block_Size1+7][Block_Size2];
    	s_vx[threadIdx.x][threadIdx.y+4]=vx[id]; 
    	s_vz[threadIdx.x+4][threadIdx.y]=vz[id]; 

	if (threadIdx.x<4)
	{
		if (blockIdx.x)			s_vz[threadIdx.x][threadIdx.y]=vz[id-4];
		else 				s_vz[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>blockDim.x-4)
	{
		if (blockIdx.x<gridDim.x-1)	s_vz[threadIdx.x+7][threadIdx.y]=vz[id+3];
		else				s_vz[threadIdx.x+7][threadIdx.y]=0.0f;
	}
	if (threadIdx.y<4)
	{
		if (blockIdx.y) 		s_vx[threadIdx.x][threadIdx.y]=vx[id-4*nnz];
		else				s_vx[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>blockDim.y-4)
	{
		if (blockIdx.y<gridDim.y-1)	s_vx[threadIdx.x][threadIdx.y+7]=vx[id+3*nnz];
		else				s_vx[threadIdx.x][threadIdx.y+7]=0.0f;
	}
	__syncthreads();


	if (!( i1<npml))
	{
		float diff2=1.1962890625000f*(s_vx[threadIdx.x][threadIdx.y+4]-s_vx[threadIdx.x][threadIdx.y+3])
			-0.0797526041667f*(s_vx[threadIdx.x][threadIdx.y+5]-s_vx[threadIdx.x][threadIdx.y+2])+
			0.0095703125000f*(s_vx[threadIdx.x][threadIdx.y+6]-s_vx[threadIdx.x][threadIdx.y+1])+
			-0.0006975446429f*(s_vx[threadIdx.x][threadIdx.y+7]-s_vx[threadIdx.x][threadIdx.y]);
		float diff1=1.1962890625000f*(s_vz[threadIdx.x+4][threadIdx.y]-s_vz[threadIdx.x+3][threadIdx.y])+
			-0.0797526041667f*(s_vz[threadIdx.x+5][threadIdx.y]-s_vz[threadIdx.x+2][threadIdx.y])+
			0.0095703125000f*(s_vz[threadIdx.x+6][threadIdx.y]-s_vz[threadIdx.x+1][threadIdx.y])+
			-0.0006975446429f*(s_vz[threadIdx.x+7][threadIdx.y]-s_vz[threadIdx.x][threadIdx.y]);
		p0[id]=2.0f*p1[id]-p0[id]+dt*dt*vel[id]*vel[id]*(_dz*diff1+_dx*diff2);
	}
}
__global__ void cuda_PML_pz_8(float *vel, float *p0,  float *convvz, float *bz, float *vz, float dt, float _dz, 
	int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int ik=i1+2*npml*i2;
	int id=blockIdx.x*(nnz-npml)+threadIdx.x+nnz*i2;

	__shared__ float s_vz[39][Block_Size2];
    	s_vz[threadIdx.x+4][threadIdx.y]=vz[id]; 
	if (threadIdx.x<4)
	{
		if (blockIdx.x)			s_vz[threadIdx.x][threadIdx.y]=vz[id-4];
		else 				s_vz[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>28)
	{
		if (blockIdx.x<gridDim.x-1)	s_vz[threadIdx.x+7][threadIdx.y]=vz[id+3];
		else				s_vz[threadIdx.x+7][threadIdx.y]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml))
	{
		float diff1=1.1962890625000f*(s_vz[threadIdx.x+4][threadIdx.y]-s_vz[threadIdx.x+3][threadIdx.y])
			-0.0797526041667f*(s_vz[threadIdx.x+5][threadIdx.y]-s_vz[threadIdx.x+2][threadIdx.y])
			+0.0095703125000f*(s_vz[threadIdx.x+6][threadIdx.y]-s_vz[threadIdx.x+1][threadIdx.y])
			-0.0006975446429f*(s_vz[threadIdx.x+7][threadIdx.y]-s_vz[threadIdx.x][threadIdx.y]);
		convvz[ik]=bz[ik]*convvz[ik]+(bz[ik]-1.0f)*_dz*diff1;	
		p0[id]+=dt*dt*vel[id]*vel[id]*convvz[ik];
	}
}
__global__ void cuda_PML_px_8(float *vel, float *p0,  float *convvx, float *bx, float *vx, float dt, float _dx, 
	int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*npml;
	int ik=i1+i2*nnz;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+threadIdx.y);

	__shared__ float s_vx[Block_Size1][39];
    	s_vx[threadIdx.x][threadIdx.y+4]=vx[id]; 
	if (threadIdx.y<4)
	{
		if (blockIdx.y) 		s_vx[threadIdx.x][threadIdx.y]=vx[id-4*nnz];
		else				s_vx[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>28)
	{
		if (blockIdx.y<gridDim.y-1)	s_vx[threadIdx.x][threadIdx.y+7]=vx[id+3*nnz];
		else				s_vx[threadIdx.x][threadIdx.y+7]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml))
	{
		float diff2=1.1962890625000f*(s_vx[threadIdx.x][threadIdx.y+4]-s_vx[threadIdx.x][threadIdx.y+3])
			-0.0797526041667f*(s_vx[threadIdx.x][threadIdx.y+5]-s_vx[threadIdx.x][threadIdx.y+2])
			+0.0095703125000f*(s_vx[threadIdx.x][threadIdx.y+6]-s_vx[threadIdx.x][threadIdx.y+1])
			-0.0006975446429f*(s_vx[threadIdx.x][threadIdx.y+7]-s_vx[threadIdx.x][threadIdx.y]);
		convvx[ik]=bx[ik]*convvx[ik]+(bx[ik]-1.0f)*_dx*diff2;
		p0[id]+=dt*dt*vel[id]*vel[id]*convvx[ik];	
	}
}


/////////////////////////////////// NJ=10 ///////////////////////////////////////
__global__ void cuda_forward_v_10(float *p, float *vx, float *vz, float _dx, float _dz,
	int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+i2*nnz;

	__shared__ float s_p[Block_Size1+9][Block_Size2+9];
    	s_p[threadIdx.x+4][threadIdx.y+4]=p[id]; 
	if (threadIdx.x<4)
	{
		if (blockIdx.x)			s_p[threadIdx.x][threadIdx.y+4]=p[id-4];
		else 				s_p[threadIdx.x][threadIdx.y+4]=0.0f;
	}
	if (threadIdx.x>blockDim.x-6)
	{
		if (blockIdx.x<gridDim.x-1)	s_p[threadIdx.x+9][threadIdx.y+4]=p[id+5];
		else				s_p[threadIdx.x+9][threadIdx.y+4]=0.0f;
	}
	if (threadIdx.y<4)
	{
		if (blockIdx.y) 		s_p[threadIdx.x+4][threadIdx.y]=p[id-4*nnz];
		else				s_p[threadIdx.x+4][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>blockDim.y-6)
	{
		if (blockIdx.y<gridDim.y-1)	s_p[threadIdx.x+4][threadIdx.y+9]=p[id+5*nnz];
		else				s_p[threadIdx.x+4][threadIdx.y+9]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml))
	{
		float diff1=1.211242675781250f*(s_p[threadIdx.x+5][threadIdx.y+4]-s_p[threadIdx.x+4][threadIdx.y+4])
			-0.089721679687500f*(s_p[threadIdx.x+6][threadIdx.y+4]-s_p[threadIdx.x+3][threadIdx.y+4])
			+0.013842773437500f*(s_p[threadIdx.x+7][threadIdx.y+4]-s_p[threadIdx.x+2][threadIdx.y+4])
			-0.001765659877232f*(s_p[threadIdx.x+8][threadIdx.y+4]-s_p[threadIdx.x+1][threadIdx.y+4])
			+0.000118679470486f*(s_p[threadIdx.x+9][threadIdx.y+4]-s_p[threadIdx.x][threadIdx.y+4]);
		float diff2= 1.211242675781250f*(s_p[threadIdx.x+4][threadIdx.y+5]-s_p[threadIdx.x+4][threadIdx.y+4])
			-0.089721679687500f*(s_p[threadIdx.x+4][threadIdx.y+6]-s_p[threadIdx.x+4][threadIdx.y+3])
			+0.013842773437500f*(s_p[threadIdx.x+4][threadIdx.y+7]-s_p[threadIdx.x+4][threadIdx.y+2])
			-0.001765659877232f*(s_p[threadIdx.x+4][threadIdx.y+8]-s_p[threadIdx.x+4][threadIdx.y+1])
			+0.000118679470486f*(s_p[threadIdx.x+4][threadIdx.y+9]-s_p[threadIdx.x+4][threadIdx.y]);
		vz[id]=_dz*diff1;
		vx[id]=_dx*diff2;
	}
}



__global__ void cuda_PML_vz_10(float *p,  float *convpz, float *bz, float *vz, float _dz, 
	int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int ik=i1+2*npml*i2;
	int id=blockIdx.x*(nnz-npml)+threadIdx.x+nnz*i2;

	__shared__ float s_p[41][Block_Size2];	// npml+9=41; Block_SizeX=32; Block_SizeY=8;
    	s_p[threadIdx.x+4][threadIdx.y]=p[id]; 
	if (threadIdx.x<4)
	{
		if (blockIdx.x)			s_p[threadIdx.x][threadIdx.y]=p[id-4];
		else 				s_p[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>26)// (npml-1)-5
	{
		if (blockIdx.x<gridDim.x-1)	s_p[threadIdx.x+9][threadIdx.y]=p[id+5];
		else				s_p[threadIdx.x+9][threadIdx.y]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml)) 
	{
		float diff1=1.1962890625000f*(s_p[threadIdx.x+5][threadIdx.y]-s_p[threadIdx.x+4][threadIdx.y])
			-0.089721679687500f*(s_p[threadIdx.x+6][threadIdx.y]-s_p[threadIdx.x+3][threadIdx.y])
			+0.013842773437500f*(s_p[threadIdx.x+7][threadIdx.y]-s_p[threadIdx.x+2][threadIdx.y])
			-0.001765659877232f*(s_p[threadIdx.x+8][threadIdx.y]-s_p[threadIdx.x+1][threadIdx.y])
			+0.000118679470486f*(s_p[threadIdx.x+9][threadIdx.y]-s_p[threadIdx.x][threadIdx.y]);
		convpz[ik]=bz[ik]*convpz[ik]+(bz[ik]-1.0f)*_dz*diff1;	
		vz[id]+=convpz[ik];
	}
}
__global__ void cuda_PML_vx_10(float *p,  float *convpx, float *bx, float *vx, float _dx, 
	int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*npml;
	int ik=i1+i2*nnz;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+threadIdx.y);
	__shared__ float s_p[Block_Size1][41];
    	s_p[threadIdx.x][threadIdx.y+4]=p[id]; 
	if (threadIdx.y<4)
	{
		if (blockIdx.y) 		s_p[threadIdx.x][threadIdx.y]=p[id-4*nnz];
		else				s_p[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>26)
	{
		if (blockIdx.y<gridDim.y-1)	s_p[threadIdx.x][threadIdx.y+9]=p[id+5*nnz];
		else				s_p[threadIdx.x][threadIdx.y+9]=0.0f;
	}
	__syncthreads();

	if (!(i1<npml))
	{
		float diff2=1.1962890625000f*(s_p[threadIdx.x][threadIdx.y+5]-s_p[threadIdx.x][threadIdx.y+4])
			-0.089721679687500f*(s_p[threadIdx.x][threadIdx.y+6]-s_p[threadIdx.x][threadIdx.y+3])
			+0.013842773437500f*(s_p[threadIdx.x][threadIdx.y+7]-s_p[threadIdx.x][threadIdx.y+2])
			-0.001765659877232f*(s_p[threadIdx.x][threadIdx.y+8]-s_p[threadIdx.x][threadIdx.y+1])
			+0.000118679470486f*(s_p[threadIdx.x][threadIdx.y+9]-s_p[threadIdx.x][threadIdx.y]);
		convpx[ik]=bx[ik]*convpx[ik]+(bx[ik]-1.0f)*_dx*diff2;	
		vx[id]+=convpx[ik];
	}
}

__global__ void cuda_forward_p_10(float *vel, float *p0, float *p1, float *vx, float *vz, float dt, float _dx, float _dz, 
	int npml, int nnz, int nnx)
{
	int i1=blockIdx.x*blockDim.x+threadIdx.x;
	int i2=blockIdx.y*blockDim.y+threadIdx.y;
	int id=i1+i2*nnz;

	__shared__ float s_vx[Block_Size1][Block_Size2+9];
	__shared__ float s_vz[Block_Size1+9][Block_Size2];
    	s_vx[threadIdx.x][threadIdx.y+5]=vx[id]; 
    	s_vz[threadIdx.x+5][threadIdx.y]=vz[id]; 

	if (threadIdx.x<5)
	{
		if (blockIdx.x)			s_vz[threadIdx.x][threadIdx.y]=vz[id-5];
		else 				s_vz[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>blockDim.x-5)
	{
		if (blockIdx.x<gridDim.x-1)	s_vz[threadIdx.x+9][threadIdx.y]=vz[id+4];
		else				s_vz[threadIdx.x+9][threadIdx.y]=0.0f;
	}
	if (threadIdx.y<5)
	{
		if (blockIdx.y) 		s_vx[threadIdx.x][threadIdx.y]=vx[id-5*nnz];
		else				s_vx[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>blockDim.y-5)
	{
		if (blockIdx.y<gridDim.y-1)	s_vx[threadIdx.x][threadIdx.y+9]=vx[id+4*nnz];
		else				s_vx[threadIdx.x][threadIdx.y+9]=0.0f;
	}
	__syncthreads();	

	if (!(i1<npml))
	{
		float diff1=1.1962890625000f*(s_vz[threadIdx.x+5][threadIdx.y]-s_vz[threadIdx.x+4][threadIdx.y])
			-0.089721679687500f*(s_vz[threadIdx.x+6][threadIdx.y]-s_vz[threadIdx.x+3][threadIdx.y])
			+0.013842773437500f*(s_vz[threadIdx.x+7][threadIdx.y]-s_vz[threadIdx.x+2][threadIdx.y])
			-0.001765659877232f*(s_vz[threadIdx.x+8][threadIdx.y]-s_vz[threadIdx.x+1][threadIdx.y])
			+0.000118679470486f*(s_vz[threadIdx.x+9][threadIdx.y]-s_vz[threadIdx.x][threadIdx.y]);
		float diff2=1.1962890625000f*(s_vx[threadIdx.x][threadIdx.y+5]-s_vx[threadIdx.x][threadIdx.y+4])
			-0.089721679687500f*(s_vx[threadIdx.x][threadIdx.y+6]-s_vx[threadIdx.x][threadIdx.y+3])
			+0.013842773437500f*(s_vx[threadIdx.x][threadIdx.y+7]-s_vx[threadIdx.x][threadIdx.y+2])
			-0.001765659877232f*(s_vx[threadIdx.x][threadIdx.y+8]-s_vx[threadIdx.x][threadIdx.y+1])
			+0.000118679470486f*(s_vx[threadIdx.x][threadIdx.y+9]-s_vx[threadIdx.x][threadIdx.y]);
		p0[id]=2.0f*p1[id]-p0[id]+dt*dt*vel[id]*vel[id]*(_dz*diff1+_dx*diff2);
	}
}
__global__ void cuda_PML_pz_10(float *vel, float *p0,  float *convvz, float *bz, float *vz, float dt, float _dz, 
	int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*npml;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int ik=i1+2*npml*i2;
	int id=blockIdx.x*(nnz-npml)+threadIdx.x+nnz*i2;

	__shared__ float s_vz[41][Block_Size2];
    	s_vz[threadIdx.x+5][threadIdx.y]=vz[id]; 
	if (threadIdx.x<5)
	{
		if (blockIdx.x)			s_vz[threadIdx.x][threadIdx.y]=vz[id-5];
		else 				s_vz[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.x>27)
	{
		if (blockIdx.x<gridDim.x-1)	s_vz[threadIdx.x+9][threadIdx.y]=vz[id+4];
		else				s_vz[threadIdx.x+9][threadIdx.y]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml))
	{
		float diff1=1.1962890625000f*(s_vz[threadIdx.x+5][threadIdx.y]-s_vz[threadIdx.x+4][threadIdx.y])		 				-0.089721679687500f*(s_vz[threadIdx.x+6][threadIdx.y]-s_vz[threadIdx.x+3][threadIdx.y])
			+0.013842773437500f*(s_vz[threadIdx.x+7][threadIdx.y]-s_vz[threadIdx.x+2][threadIdx.y])
			-0.001765659877232f*(s_vz[threadIdx.x+8][threadIdx.y]-s_vz[threadIdx.x+1][threadIdx.y])
			+0.000118679470486f*(s_vz[threadIdx.x+9][threadIdx.y]-s_vz[threadIdx.x][threadIdx.y]);
		convvz[ik]=bz[ik]*convvz[ik]+(bz[ik]-1.0f)*_dz*diff1;	
		p0[id]+=dt*dt*vel[id]*vel[id]*convvz[ik];
	}
}
__global__ void cuda_PML_px_10(float *vel, float *p0,  float *convvx, float *bx, float *vx, float dt, float _dx,
	int npml, int nnz, int nnx)
{
	// bz1: top and bottom PML ABC coefficients, decay p (px,pz) along z direction
	// bz2: top and bottom PML ABC coefficients, decay v (vx,vz) along z direction
	// only 2 blocks used vertically, blockIdx.y=0, 1

	// id: position in whole zone(including PML)
	// ik: position in top or bottom PML zone itself
	// blockIdx.y==0, top PML zone; blockIdx.y==1, bottom PML zone
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*npml;
	int ik=i1+i2*nnz;
	int id=i1+nnz*(blockIdx.y*(nnx-npml)+threadIdx.y);

	__shared__ float s_vx[Block_Size1][41];
    	s_vx[threadIdx.x][threadIdx.y+5]=vx[id]; 
	if (threadIdx.y<5)
	{
		if (blockIdx.y) 		s_vx[threadIdx.x][threadIdx.y]=vx[id-5*nnz];
		else				s_vx[threadIdx.x][threadIdx.y]=0.0f;
	}
	if (threadIdx.y>27)
	{
		if (blockIdx.y<gridDim.y-1)	s_vx[threadIdx.x][threadIdx.y+9]=vx[id+4*nnz];
		else				s_vx[threadIdx.x][threadIdx.y+9]=0.0f;
	}
	__syncthreads();

	if (!( i1<npml))
	{
		float diff2=1.1962890625000f*(s_vx[threadIdx.x][threadIdx.y+5]-s_vx[threadIdx.x][threadIdx.y+4])
			 -0.089721679687500f*(s_vx[threadIdx.x][threadIdx.y+6]-s_vx[threadIdx.x][threadIdx.y+3])+
			0.013842773437500f*(s_vx[threadIdx.x][threadIdx.y+7]-s_vx[threadIdx.x][threadIdx.y+2])
			 -0.001765659877232f*(s_vx[threadIdx.x][threadIdx.y+8]-s_vx[threadIdx.x][threadIdx.y+1])+
			0.000118679470486f*(s_vx[threadIdx.x][threadIdx.y+9]-s_vx[threadIdx.x][threadIdx.y]);
		convvx[ik]=bx[ik]*convvx[ik]+(bx[ik]-1.0f)*_dx*diff2;	
		p0[id]+=dt*dt*vel[id]*vel[id]*convvx[ik];
	}
}
/////////////////////////////////// read and save the boundary //////////////////////////////////
// read and write the inner computation zone boundary coefficients from and into RAM along z direction
// adj==flase, write and save boundary; adj==true, read the boundary
__global__ void cuda_rw_boundary_z(float *boundary_z, float *p, int npml, int nnz, int nnx, int NJ, bool adj)
{
	int nx=nnx-2*npml;
	int nz=nnz-2*npml;
    	int i1=threadIdx.x+blockDim.x*blockIdx.x;
	int i2=threadIdx.y+blockDim.y*blockIdx.y;
	int id=i1+2*(NJ-1)*i2;
	int i1p=i1+npml;
	int i2p=i2+npml;
	int idp=i1p+nnz*i2p;

	if (i1<NJ-1 && i2<nx)
	{
		if(adj)
		{
			p[idp]=boundary_z[id];
			p[idp+(nz-NJ+1)]=boundary_z[id+(NJ-1)];
		}
		else	
		{
			boundary_z[id]=p[idp];	
			boundary_z[id+(NJ-1)]=p[idp+(nz-NJ+1)];
		}
	}
}

// read and write the inner computation zone boundary coefficients from and into RAM along x direction
// adj==flase, write and save boundary; adj==true, read the boundary
__global__ void cuda_rw_boundary_x(float *boundary_x, float *p, int npml, int nnz, int nnx, int NJ, bool adj)
{
	int nx=nnx-2*npml;
	int nz=nnz-2*npml;
    	int i1=threadIdx.x+blockDim.x*blockIdx.x;
	int i2=threadIdx.y+blockDim.y*blockIdx.y;
	int id=i1+nz*i2;
	int i1p=i1+npml;
	int i2p=i2+npml;
	int idp=i1p+nnz*i2p;

	if (i1<nz && i2<NJ-1)
	{
		if (adj)
		{
			p[idp]=boundary_x[id];
			p[idp+nnz*(nx-NJ+1)]=boundary_x[id+nz*(NJ-1)];
		}
		else
		{
			boundary_x[id]=p[idp];
			boundary_x[id+nz*(NJ-1)]=p[idp+nnz*(nx-NJ+1)];
		}
	}
}

