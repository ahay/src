/* 3D 8-th order elastic wave propagation with sponge ABC and moment tensor source
By Yangkang Chen, 2020
Revised in July, 2021
Currently still in a draft version

DEMO:
https://github.com/chenyk1990/tutorials/blob/main/demo/efd3dmt/SConstruct
*/

#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static bool iffree;
static int nb, nz, nx, ny, nt, nzpad, nxpad, nypad;
static float dz, dx, dy, _dz, _dx, _dy, dt, fm;

void expand3d(float*** b, float*** a)
/*< expand domain of 'a' to 'b': source(a)-->destination(b) >*/
{
    int iz,ix,iy;

#ifdef _OPENMP
#pragma omp parallel for default(none)	\
	private(iy,ix,iz)			\
	shared(b,a,nb,nz,nx,ny)
#endif
	for     (iy=0;iy<ny;iy++) {
    for     (ix=0;ix<nx;ix++) {
	for (iz=0;iz<nz;iz++) {
	    b[nb+iy][nb+ix][nb+iz] = a[iy][ix][iz];
	}
    }
	}
	
	for     (iy=0; iy<nypad; iy++) {
    for     (ix=0; ix<nxpad; ix++) {
	for (iz=0; iz<nb;    iz++) {
	    b[iy][ix][      iz  ] = b[iy][ix][nb  ];
	    b[iy][ix][nzpad-iz-1] = b[iy][ix][nzpad-nb-1];
	}
    }
    }

	for     (iy=0; iy<nypad;    iy++) {
    for     (ix=0; ix<nb;    ix++) {
	for (iz=0; iz<nzpad; iz++) {
	    b[iy][ix 	 ][iz] = b[iy][nb  		][iz];
	    b[iy][nxpad-ix-1 ][iz] = b[iy][nxpad-nb-1	][iz];
	}
    }
    }
    
	for     (iy=0; iy<nb;    iy++) {
    for     (ix=0; ix<nxpad;    ix++) {
	for (iz=0; iz<nzpad; iz++) {
	    b[iy		][ix][iz] = b[nb			][ix][iz];
	    b[nypad-iy-1][ix][iz] = b[nypad-nb-1][ix][iz];
	}
    }
    }
}


void window3d(float ***a, float ***b)
/*< window 'b' to 'a': source(b)-->destination(a) >*/
{
    int iz,ix,iy;

#ifdef _OPENMP
#pragma omp parallel for default(none)	\
	private(iy,ix,iz)			\
	shared(b,a,nb,nz,nx,ny)
#endif
	for     (iy=0;iy<ny;iy++) {
    for     (ix=0;ix<nx;ix++) {
	for (iz=0;iz<nz;iz++) {
	    a[iy][ix][iz]=b[nb+iy][nb+ix][nb+iz] ;
	}
    }
    }
}


void forward_uvy_uvx_uvz(float ***uvx, float ***uvy, float ***uvz, float ***txx, float ***tyy, float ***tzz, float ***txz, float ***txy, float ***tyz, float ***rho)
/*< forward step: update uvx, uvz >*/
{
	int i1, i2, i3;
	float diff1, diff2, diff3, diff4, diff5, diff6, diff7, diff8, diff9;

// #ifdef _OPENMP
// #pragma omp parallel for default(none) 	\
// 	private(i1,i2,i3,diff1,diff2,diff3,diff4,diff5,diff6,diff7,diff8,diff9)		\
// 	shared(uvx,uvy,uvz,txx,tyy,tzz,txz,txy,tyz,rho,nxpad,nypad,nzpad,dt,_dx,_dy,_dz)
// #endif
	int nnz=0, nnx=0, nny=0, nsum=0;
	for(i3=4; i3<nypad-3; i3++)
	for(i2=4; i2<nxpad-3; i2++)
	for(i1=4; i1<nzpad-3; i1++)
	{
		diff1 = 1.1962890625000f*(tzz[i3][i2][i1]-tzz[i3][i2][i1-1])
			-0.0797526041667f*(tzz[i3][i2][i1+1]-tzz[i3][i2][i1-2])
			+0.0095703125000f*(tzz[i3][i2][i1+2]-tzz[i3][i2][i1-3])
			-0.0006975446429f*(tzz[i3][i2][i1+3]-tzz[i3][i2][i1-4]);
		diff2 = 1.1962890625000f*(txz[i3][i2][i1]-txz[i3][i2-1][i1])
			-0.0797526041667f*(txz[i3][i2+1][i1]-txz[i3][i2-2][i1])
			+0.0095703125000f*(txz[i3][i2+2][i1]-txz[i3][i2-3][i1])
			-0.0006975446429f*(txz[i3][i2+3][i1]-txz[i3][i2-4][i1]);
		diff3 = 1.1962890625000f*(tyz[i3][i2][i1]-tyz[i3-1][i2][i1])
			-0.0797526041667f*(tyz[i3+1][i2][i1]-tyz[i3-2][i2][i1])
			+0.0095703125000f*(tyz[i3+2][i2][i1]-tyz[i3-3][i2][i1])
			-0.0006975446429f*(tyz[i3+3][i2][i1]-tyz[i3-4][i2][i1]);
		diff4 = 1.1962890625000f*(txx[i3][i2][i1]-txx[i3][i2-1][i1])
			-0.0797526041667f*(txx[i3][i2+1][i1]-txx[i3][i2-2][i1])
			+0.0095703125000f*(txx[i3][i2+2][i1]-txx[i3][i2-3][i1])
			-0.0006975446429f*(txx[i3][i2+3][i1]-txx[i3][i2-4][i1]);
		diff5 = 1.1962890625000f*(txz[i3][i2][i1]-txz[i3][i2][i1-1])
			-0.0797526041667f*(txz[i3][i2][i1+1]-txz[i3][i2][i1-2])
			+0.0095703125000f*(txz[i3][i2][i1+2]-txz[i3][i2][i1-3])
			-0.0006975446429f*(txz[i3][i2][i1+3]-txz[i3][i2][i1-4]);
		diff6 = 1.1962890625000f*(txy[i3][i2][i1]-txy[i3-1][i2][i1])
			-0.0797526041667f*(txy[i3+1][i2][i1]-txy[i3-2][i2][i1])
			+0.0095703125000f*(txy[i3+2][i2][i1]-txy[i3-3][i2][i1])
			-0.0006975446429f*(txy[i3+3][i2][i1]-txy[i3-4][i2][i1]);
		diff7 = 1.1962890625000f*(txy[i3][i2][i1]-txy[i3][i2-1][i1])
			-0.0797526041667f*(txy[i3][i2+1][i1]-txy[i3][i2-2][i1])
			+0.0095703125000f*(txy[i3][i2+2][i1]-txy[i3][i2-3][i1])
			-0.0006975446429f*(txy[i3][i2+3][i1]-txy[i3][i2-4][i1]);
		diff8 = 1.1962890625000f*(tyz[i3][i2][i1]-tyz[i3][i2][i1-1])
			-0.0797526041667f*(tyz[i3][i2][i1+1]-tyz[i3][i2][i1-2])
			+0.0095703125000f*(tyz[i3][i2][i1+2]-tyz[i3][i2][i1-3])
			-0.0006975446429f*(tyz[i3][i2][i1+3]-tyz[i3][i2][i1-4]);
		diff9 = 1.1962890625000f*(tyy[i3][i2][i1]-tyy[i3-1][i2][i1])
			-0.0797526041667f*(tyy[i3+1][i2][i1]-tyy[i3-2][i2][i1])
			+0.0095703125000f*(tyy[i3+2][i2][i1]-tyy[i3-3][i2][i1])
			-0.0006975446429f*(tyy[i3+3][i2][i1]-tyy[i3-4][i2][i1]);
		uvz[i3][i2][i1]+=dt*rho[i3][i2][i1]*(_dx*diff2+_dz*diff1+_dy*diff3);
		uvx[i3][i2][i1]+=dt*rho[i3][i2][i1]*(_dx*diff4+_dz*diff5+_dy*diff6);
		uvy[i3][i2][i1]+=dt*rho[i3][i2][i1]*(_dx*diff7+_dz*diff8+_dy*diff9);
// 		sf_warning("dt=%g,uvz[i3][i2][i1]=%g,uvx[i3][i2][i1]=%g,uvy[i3][i2][i1]=%g",dt,uvz[i3][i2][i1],uvx[i3][i2][i1],uvy[i3][i2][i1]);
// 		if(diff1!=0 || diff2!=0 || diff3!=0)
// 		{
// // 		sf_warning("diff1=%g,diff2=%g,diff3=%g",diff1,diff2,diff3);
// // 		sf_warning("diff4=%g,diff5=%g,diff6=%g",diff4,diff5,diff6);
// // 		sf_warning("diff7=%g,diff8=%g,diff9=%g",diff7,diff8,diff9);
// 		sf_warning("i1=%d,i2=%d,i3=%d",i1,i2,i3);
// 		}
		
// 		if(uvz[i3][i2][i1]!=0) nnz++;
// 		if(uvx[i3][i2][i1]!=0) nnx++;
// 		if(uvy[i3][i2][i1]!=0) nny++;
// 		nsum++;
	}
// 	sf_warning("UV: nnz=%d,nnx=%d,nny=%d,nsum=%d",nnz,nnx,nny,nsum);
}


//#define C1 +0.800000   /* +4/5    */
//#define C2 -0.200000   /* -1/5    */
//#define C3 +0.038095   /* +4/105  */
//#define C4 -0.003571   /* -5/280  */

#define C1 +1.1962890625000f   /* +4/5    */
#define C2 -0.0797526041667f   /* -1/5    */
#define C3 +0.0095703125000f   /* +4/105  */
#define C4 -0.0006975446429f   /* -5/280  */

//#define Dx(a,ix,iz,s) (C4*(a[ix+4][iz] - a[ix-4][iz]) +		\
//		       C3*(a[ix+3][iz] - a[ix-3][iz]) +		\
//		       C2*(a[ix+2][iz] - a[ix-2][iz]) +		\
//		       C1*(a[ix+1][iz] - a[ix-1][iz])  )*s
//#define Dz(a,ix,iz,s) (C4*(a[ix][iz+4] - a[ix][iz-4]) +		\
//		       C3*(a[ix][iz+3] - a[ix][iz-3]) +		\
//		       C2*(a[ix][iz+2] - a[ix][iz-2]) +		\
//		       C1*(a[ix][iz+1] - a[ix][iz-1])  )*s

/*#define Dx(a,ix,iz,s) (C4*(a[ix+4][iz] - a[ix-3][iz]) +		\
		       C3*(a[ix+3][iz] - a[ix-2][iz]) +		\
		       C2*(a[ix+2][iz] - a[ix-1][iz]) +		\
		       C1*(a[ix+1][iz] - a[ix][iz])  )*s
#define Dz(a,ix,iz,s) (C4*(a[ix][iz+4] - a[ix][iz-3]) +		\
		       C3*(a[ix][iz+3] - a[ix][iz-2]) +		\
		       C2*(a[ix][iz+2] - a[ix][iz-1]) +		\
		       C1*(a[ix][iz+1] - a[ix][iz])  )*s  */

#define Dxp(a,ix,iz,s) (C4*(a[ix+3][iz] - a[ix-4][iz]) +		\
		       C3*(a[ix+2][iz] - a[ix-3][iz]) +		\
		       C2*(a[ix+1][iz] - a[ix-2][iz]) +		\
		       C1*(a[ix][iz] - a[ix-1][iz])  )*s
#define Dzp(a,ix,iz,s) (C4*(a[ix][iz+3] - a[ix][iz-4]) +		\
		       C3*(a[ix][iz+2] - a[ix][iz-3]) +		\
		       C2*(a[ix][iz+1] - a[ix][iz-2]) +		\
		       C1*(a[ix][iz] - a[ix][iz-1])  )*s


#define Dxs(a,ix,iz,s) (C4*(a[ix+4][iz] - a[ix-3][iz]) +		\
		       C3*(a[ix+3][iz] - a[ix-2][iz]) +		\
		       C2*(a[ix+2][iz] - a[ix-1][iz]) +		\
		       C1*(a[ix+1][iz] - a[ix][iz])  )*s
#define Dzs(a,ix,iz,s) (C4*(a[ix][iz+4] - a[ix][iz-3]) +		\
		       C3*(a[ix][iz+3] - a[ix][iz-2]) +		\
		       C2*(a[ix][iz+2] - a[ix][iz-1]) +		\
		       C1*(a[ix][iz+1] - a[ix][iz])  )*s
		       		       					       
void forward_txx_tyy_tzz_txz_txy_tyz(float ***uvx, float ***uvy, float ***uvz, float ***txx, float ***tyy, float ***tzz, float ***txz, float ***txy, float ***tyz, float ***vp, float ***vs)
/*< forward step: update txx, tzz, txz >*/
{
	int i1, i2, i3;
	float diff1, diff12, diff13, diff2, diff21, diff23, diff3, diff31,diff32;

// #ifdef _OPENMP
// #pragma omp parallel for default(none) 	\
// 	private(i1,i2,diff1,diff2,diff3,diff4)		\
// 	shared(uvx,uvz,txx,tzz,txz,vp,vs,nxpad,nzpad,dt,_dx,_dy,_dz)
// #endif
int nnz=0, nnx=0, nny=0, nsum=0;
	for(i3=3; i3<nypad-4; i3++)
	for(i2=3; i2<nxpad-4; i2++)
	for(i1=3; i1<nzpad-4; i1++)
	{
			
		diff1 = 1.1962890625000f*(uvz[i3][i2][i1+1]-uvz[i3][i2][i1])
			-0.0797526041667f*(uvz[i3][i2][i1+2]-uvz[i3][i2][i1-1])
			+0.0095703125000f*(uvz[i3][i2][i1+3]-uvz[i3][i2][i1-2])
			-0.0006975446429f*(uvz[i3][i2][i1+4]-uvz[i3][i2][i1-3]);
		diff12 = 1.1962890625000f*(uvz[i3][i2+1][i1]-uvz[i3][i2][i1])
			-0.0797526041667f*(uvz[i3][i2+2][i1]-uvz[i3][i2-1][i1])
			+0.0095703125000f*(uvz[i3][i2+3][i1]-uvz[i3][i2-2][i1])
			-0.0006975446429f*(uvz[i3][i2+4][i1]-uvz[i3][i2-3][i1]);
		diff13 = 1.1962890625000f*(uvz[i3+1][i2][i1]-uvz[i3][i2][i1])
			-0.0797526041667f*(uvz[i3+2][i2][i1]-uvz[i3-1][i2][i1])
			+0.0095703125000f*(uvz[i3+3][i2][i1]-uvz[i3-2][i2][i1])
			-0.0006975446429f*(uvz[i3+4][i2][i1]-uvz[i3-3][i2][i1]);

		diff2 = 1.1962890625000f*(uvx[i3][i2+1][i1]-uvx[i3][i2][i1])
			-0.0797526041667f*(uvx[i3][i2+2][i1]-uvx[i3][i2-1][i1])
			+0.0095703125000f*(uvx[i3][i2+3][i1]-uvx[i3][i2-2][i1])
			-0.0006975446429f*(uvx[i3][i2+4][i1]-uvx[i3][i2-3][i1]);
		diff21 = 1.1962890625000f*(uvx[i3][i2][i1+1]-uvx[i3][i2][i1])
			-0.0797526041667f*(uvx[i3][i2][i1+2]-uvx[i3][i2][i1-1])
			+0.0095703125000f*(uvx[i3][i2][i1+3]-uvx[i3][i2][i1-2])
			-0.0006975446429f*(uvx[i3][i2][i1+4]-uvx[i3][i2][i1-3]);
		diff23 = 1.1962890625000f*(uvx[i3+1][i2][i1]-uvx[i3][i2][i1])
			-0.0797526041667f*(uvx[i3+2][i2][i1]-uvx[i3-1][i2][i1])
			+0.0095703125000f*(uvx[i3+3][i2][i1]-uvx[i3-2][i2][i1])
			-0.0006975446429f*(uvx[i3+4][i2][i1]-uvx[i3-3][i2][i1]);
						
		diff3 = 1.1962890625000f*(uvy[i3+1][i2][i1]-uvy[i3][i2][i1])
			-0.0797526041667f*(uvy[i3+2][i2][i1]-uvy[i3-1][i2][i1])
			+0.0095703125000f*(uvy[i3+3][i2][i1]-uvy[i3-2][i2][i1])
			-0.0006975446429f*(uvy[i3+4][i2][i1]-uvy[i3-3][i2][i1]);	
		diff31 = 1.1962890625000f*(uvy[i3][i2][i1+1]-uvy[i3][i2][i1])
			-0.0797526041667f*(uvy[i3][i2][i1+2]-uvy[i3][i2][i1-1])
			+0.0095703125000f*(uvy[i3][i2][i1+3]-uvy[i3][i2][i1-2])
			-0.0006975446429f*(uvy[i3][i2][i1+4]-uvy[i3][i2][i1-3]);	
		diff32 = 1.1962890625000f*(uvy[i3][i2+1][i1]-uvy[i3][i2][i1])
			-0.0797526041667f*(uvy[i3][i2+2][i1]-uvy[i3][i2-1][i1])
			+0.0095703125000f*(uvy[i3][i2+3][i1]-uvy[i3][i2-2][i1])
			-0.0006975446429f*(uvy[i3][i2+4][i1]-uvy[i3][i2-3][i1]);				
			
	txx[i3][i2][i1]+=dt*(vp[i3][i2][i1]*_dx*diff2+(vp[i3][i2][i1]-2*vs[i3][i2][i1])*(_dy*diff3+_dz*diff1));
	tyy[i3][i2][i1]+=dt*(vp[i3][i2][i1]*_dy*diff3+(vp[i3][i2][i1]-2*vs[i3][i2][i1])*(_dx*diff2+_dz*diff1));
	tzz[i3][i2][i1]+=dt*(vp[i3][i2][i1]*_dz*diff1+(vp[i3][i2][i1]-2*vs[i3][i2][i1])*(_dx*diff2+_dy*diff3));	
		txy[i3][i2][i1]+=dt*vs[i3][i2][i1]*(_dy*diff23+_dx*diff32);
		txz[i3][i2][i1]+=dt*vs[i3][i2][i1]*(_dz*diff21+_dx*diff12);
		tyz[i3][i2][i1]+=dt*vs[i3][i2][i1]*(_dz*diff31+_dy*diff13);

// 		if(vp[i3][i2][i1]!=0) nnx++;
// 		if(vs[i3][i2][i1]!=0) nny++;

// 		if(txx[i3][i2][i1]!=0) nnx++;
// 		if(tyy[i3][i2][i1]!=0) nny++;
// 		if(tzz[i3][i2][i1]!=0) nnz++;
// 		nsum++;
	}
}


void apply_sponge(float ***a /*3-D matrix*/, float *bndr) 
/*< boundary decay (simple ABC but stable and works)>*/
{
    int i;
    int iz, iy, ix;
	
    /* top */
#ifdef _OPENMP
#pragma omp parallel default(shared) private(iz,ix,iy,i)
{
#endif

#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nb; iz++) {  
        for (ix=0; ix < nxpad; ix++) {
        for (iy=0; iy < nypad; iy++) {
	  a[iy][ix][iz] *= bndr[iz];
        }
        }
    }
    
    /* bottom */
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nb; iz++) {  
        for (ix=0; ix < nxpad; ix++) {
        for (iy=0; iy < nypad; iy++) {
	  a[iy][ix][nzpad-1-iz] *= bndr[iz];
        }
    }
    }
      
    /* left x*/
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nzpad; iz++) {  
        for (ix=0; ix < nb; ix++) {
        for (iy=0; iy < nypad; iy++) { 
	  a[iy][ix][iz] *= bndr[ix];
        }
        }
    }
    
    /* right x*/
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nzpad; iz++) {  
        for (ix=0; ix < nb; ix++) {
        for (iy=0; iy < nypad; iy++) {     
          a[iy][nxpad-1-ix][iz] *= bndr[ix];
        }
        }
    }
        
    /* left y*/
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nzpad; iz++) {  
       for (ix=0; ix < nxpad; ix++) {
        for (iy=0; iy < nb; iy++) { 
	  a[iy][ix][iz] *= bndr[iy];
        }
        }
    }
        
    /* right y*/
#ifdef _OPENMP
#pragma omp for
#endif
    for (iz=0; iz < nzpad; iz++) {  
       for (ix=0; ix < nxpad; ix++) {
        for (iy=0; iy < nb; iy++) {    
          a[nypad-1-iy][ix][iz] *= bndr[iy];
        }
        }
    }
    
#ifdef _OPENMP
}
#endif
}


void apply_free(float ***tzz, float ***txz, float ***tyz)
/*< apply free-surface boundary condition >*/
{
	int ix,iz,iy;

// #ifdef _OPENMP
// #pragma omp parallel for	    \
//     private(ix,iz)		    \
//     shared(bndr,u)
// #endif
	for(iy=0; iy<nypad; iy++)
	{
	for(ix=0; ix<nxpad; ix++)
	{
		for(iz=0;iz<nb;iz++){	// top ABC			
			tzz[iy][ix][iz]=0;
			txz[iy][ix][iz]=0;
			tyz[iy][ix][iz]=0;
			}
		}
		}
}

int main(int argc, char* argv[])
{
	bool verb, ifwfd;
	int jsnap, ft, it, ib, sx, sy, sz, ix, iy, iz, n4;
	float a, ct, *wlt, *bndr;
	float ***vp0, ***vs0, ***rho0, ***vp, ***vs, ***rho, ***uvx, ***uvy, ***uvz, ***tmp, ***tmp1, ***qpx, ***qsx, ***qpz, ***qsz, ***txx, ***tyy, ***tzz, ***txz, ***txy, ***tyz, ***Rvx, ***Rvy, ***Rvz;
	float mt[9];

	sf_file Fvp, Fvs, Frho, Fwavx,Fwavz,Fwavy,Frvx,Frvy,Frvz;
    
    	sf_init(argc,argv);
#ifdef _OPENMP
    	omp_init();
#endif

	Fvp = sf_input("in");/* p-wave veloctiy */
	Fvs = sf_input("vs");/* s-wave veloctiy */
	Frho = sf_input("rho");/* density */

	Frvz = sf_output("out");/* z-component of wavefield */	
	Frvx = sf_output("rvx");/* x-component of wavefield */	
	Frvy = sf_output("rvy");/* y-component of wavefield */	
		
    	if(!sf_getbool("verb",&verb)) verb=false;    /* verbosity */
    	if(!sf_getbool("ifwfd",&ifwfd)) ifwfd=false;    /* if output wavefield */
    	if (!sf_histint(Fvp,"n1",&nz)) sf_error("No n1= in input");/* veloctiy model: nz */
    	if (!sf_histint(Fvp,"n2",&nx)) sf_error("No n2= in input");/* veloctiy model: nx */
        if (!sf_histint(Fvp,"n3",&ny)) sf_error("No n3= in input");/* veloctiy model: ny */
    	if (!sf_histfloat(Fvp,"d1",&dz)) sf_error("No d1= in input");/* veloctiy model: dz */
    	if (!sf_histfloat(Fvp,"d2",&dx)) sf_error("No d2= in input");/* veloctiy model: dx */
        if (!sf_histfloat(Fvp,"d3",&dy)) sf_error("No d3= in input");/* veloctiy model: dy */
    	if (!sf_getint("nb",&nb)) nb=30; /* thickness of ABC layer */
    	if (!sf_getint("nt",&nt)) sf_error("nt required");/* number of time steps */
		if (jsnap>nt) sf_error("make sure jsnap<=nt");
    	if (!sf_getfloat("dt",&dt)) sf_error("dt required");/* time sampling interval */
    	if (!sf_getfloat("fm",&fm)) fm=20.0; /*dominant freq of Ricker wavelet */
   	if (!sf_getint("ft",&ft)) ft=0; /* first recorded time */
    	if (!sf_getint("jsnap",&jsnap)) jsnap=1;	/* interval for snapshots  */
      if(!sf_getfloat("ct",&ct)) ct=0.01;/*for absorbing boundary*/

	if(ifwfd)
	{
	Fwavx=sf_output("wavx");
	Fwavy=sf_output("wavy");
	Fwavz=sf_output("wavz");
	sf_putint(Fwavx,"d4",dt*jsnap);
	sf_putint(Fwavy,"d4",dt*jsnap);
	sf_putint(Fwavz,"d4",dt*jsnap);	

	sf_putint(Fwavx,"o4",0);
	sf_putint(Fwavy,"o4",0);
	sf_putint(Fwavz,"o4",0);	
	n4=ceilf(nt/(jsnap*1.0));
	sf_putint(Fwavx,"n4",n4);
	sf_putint(Fwavy,"n4",n4);
	sf_putint(Fwavz,"n4",n4);
	}else
	{
	Fwavx=NULL;
	Fwavy=NULL;
	Fwavz=NULL;
	}
	
	sf_putint(Frvx,"n1",nt);
	sf_putint(Frvx,"n2",nx);
	sf_putint(Frvx,"n3",ny);
	sf_putfloat(Frvx,"d1",dt);
	sf_putfloat(Frvx,"d2",dx);
	sf_putfloat(Frvx,"d3",dy);
	sf_putfloat(Frvx,"o1",0);
	sf_putfloat(Frvx,"o2",0);
	sf_putfloat(Frvx,"o3",0);
	sf_putint(Frvy,"n1",nt);
	sf_putint(Frvy,"n2",nx);
	sf_putint(Frvy,"n3",ny);
	sf_putfloat(Frvy,"d1",dt);
	sf_putfloat(Frvy,"d2",dx);
	sf_putfloat(Frvy,"d3",dy);
	sf_putfloat(Frvy,"o1",0);
	sf_putfloat(Frvy,"o2",0);
	sf_putfloat(Frvy,"o3",0);
	sf_putint(Frvz,"n1",nt);
	sf_putint(Frvz,"n2",nx);
	sf_putint(Frvz,"n3",ny);
	sf_putfloat(Frvz,"d1",dt);
	sf_putfloat(Frvz,"d2",dx);
	sf_putfloat(Frvz,"d3",dy);
	sf_putfloat(Frvz,"o1",0);
	sf_putfloat(Frvz,"o2",0);
	sf_putfloat(Frvz,"o3",0);

	nzpad=nz+2*nb;
	nxpad=nx+2*nb;
	nypad=ny+2*nb;		
	_dz=1.0/dz;
	_dx=1.0/dx;
	_dy=1.0/dy;
	sx=nxpad/2;
	sy=nypad/2;
	sz=nzpad/2;
	
    if (!sf_getfloats("M",mt,9)){mt[0]=1;mt[4]=1;mt[8]=1;}
    sf_warning("M11=%g",mt[0]);//Mxx
    sf_warning("M21=%g",mt[1]);//Myx
    sf_warning("M31=%g",mt[2]);//Mzx
    sf_warning("M12=%g",mt[3]);//Mxy
    sf_warning("M22=%g",mt[4]);//Myy
    sf_warning("M32=%g",mt[5]);//Mzy
    sf_warning("M13=%g",mt[6]);//Mxz
    sf_warning("M23=%g",mt[7]);//Myz
    sf_warning("M33=%g",mt[8]);//Mzz

    float Mxx,Myx,Mzx,Mxy,Myy,Mzy,Mxz,Myz,Mzz;
    Mxx=mt[0];
    Myx=mt[1];
    Mzx=mt[2];
    Mxy=mt[3];
    Myy=mt[4];
    Mzy=mt[5];
    Mxz=mt[6];
    Myz=mt[7];
    Mzz=mt[8];

    sf_warning("nz=%d,nx=%d,ny=%d,nzpad=%d,nxpad=%d,nypad=%d",nz,nx,ny,nzpad,nxpad,nypad);
    
	if (!sf_getint("sx",&sx)) sx=nxpad/2; 
	if (!sf_getint("sy",&sy)) sy=nypad/2; 
	if (!sf_getint("sz",&sz)) sz=nzpad/2; 
	
	if (!sf_getbool("free",&iffree)) iffree=false;  
	/*if free surface*/

	/* allocate memory for variables */
	wlt=sf_floatalloc(nt);
	bndr=sf_floatalloc(nb);
	vp0=sf_floatalloc3(nz,nx,ny); 
	vs0=sf_floatalloc3(nz,nx,ny); 
	rho0=sf_floatalloc3(nz,nx,ny);
	vp=sf_floatalloc3(nzpad, nxpad, nypad);
	vs=sf_floatalloc3(nzpad, nxpad, nypad);
	rho=sf_floatalloc3(nzpad, nxpad, nypad);
	uvx=sf_floatalloc3(nzpad, nxpad, nypad);
	uvy=sf_floatalloc3(nzpad, nxpad, nypad);
	uvz=sf_floatalloc3(nzpad, nxpad, nypad);
	qpx=sf_floatalloc3(nzpad, nxpad, nypad);
	qsx=sf_floatalloc3(nzpad, nxpad, nypad);
	qpz=sf_floatalloc3(nzpad, nxpad, nypad);
	qsz=sf_floatalloc3(nzpad, nxpad, nypad);
	tmp=sf_floatalloc3(nzpad, nxpad, nypad);
	Rvx=sf_floatalloc3(nt,nx,ny);
	Rvz=sf_floatalloc3(nt,nx,ny);
	Rvy=sf_floatalloc3(nt,nx,ny);
	tmp1=sf_floatalloc3(nzpad, nxpad, nypad);
	txx=sf_floatalloc3(nzpad, nxpad, nypad);
	tyy=sf_floatalloc3(nzpad, nxpad, nypad);
	tzz=sf_floatalloc3(nzpad, nxpad, nypad);
	txz=sf_floatalloc3(nzpad, nxpad, nypad);
	txy=sf_floatalloc3(nzpad, nxpad, nypad);
	tyz=sf_floatalloc3(nzpad, nxpad, nypad);
	/* initialization */
	for(it=0;it<nt;it++)
	{
		a=SF_PI*fm*(it*dt-1.0/fm);a*=a;
		wlt[it]=(1.0-2.0*a)*expf(-a);
	}
	for(ib=0;ib<nb;ib++)
	{
		a=ct*(nb-ib-1);
		bndr[ib]=expf(-a*a);
	}
	
	sf_floatread(vp0[0][0], nz*nx*ny, Fvp);
	sf_floatread(vs0[0][0], nz*nx*ny, Fvs);
	sf_floatread(rho0[0][0], nz*nx*ny, Frho);

	for(iy=0; iy<ny; iy++)
	for(ix=0; ix<nx; ix++)
	for(iz=0; iz<nz; iz++)
	{
		vp0[iy][ix][iz]=rho0[iy][ix][iz]*vp0[iy][ix][iz]*vp0[iy][ix][iz];
		vs0[iy][ix][iz]=rho0[iy][ix][iz]*vs0[iy][ix][iz]*vs0[iy][ix][iz];
		rho0[iy][ix][iz]=1.0/rho0[iy][ix][iz];
	}
	
	float h3;
	h3=dt/(dx*dy*dz);
	
	expand3d(vp, vp0);
	expand3d(vs, vs0);
	expand3d(rho, rho0);
	memset(uvx[0][0],0,nzpad*nxpad*nypad*sizeof(float));
	memset(uvy[0][0],0,nzpad*nxpad*nypad*sizeof(float));
	memset(uvz[0][0],0,nzpad*nxpad*nypad*sizeof(float));
	memset(txx[0][0],0,nzpad*nxpad*nypad*sizeof(float));
	memset(tyy[0][0],0,nzpad*nxpad*nypad*sizeof(float));
	memset(tzz[0][0],0,nzpad*nxpad*nypad*sizeof(float));
	memset(txz[0][0],0,nzpad*nxpad*nypad*sizeof(float));
	memset(txy[0][0],0,nzpad*nxpad*nypad*sizeof(float));
	memset(tyz[0][0],0,nzpad*nxpad*nypad*sizeof(float));
	memset(qpx[0][0],0,nzpad*nxpad*nypad*sizeof(float));
	memset(qsx[0][0],0,nzpad*nxpad*nypad*sizeof(float));
	memset(qpz[0][0],0,nzpad*nxpad*nypad*sizeof(float));
	memset(qsz[0][0],0,nzpad*nxpad*nypad*sizeof(float));
    memset(tmp[0][0],0,nzpad*nxpad*nypad*sizeof(float));
    memset(tmp1[0][0],0,nzpad*nxpad*nypad*sizeof(float));  
    
	for(it=0; it<nt; it++)
	{
// 		sf_warning("it=%d",it);
// 		txx[sy+nb][sx+nb][sz+nb]+=wlt[it];
// 		tzz[sy+nb][sx+nb][sz+nb]+=wlt[it];
// 		tyy[sy+nb][sx+nb][sz+nb]+=wlt[it];

// 		uvx[sy+nb][sx+nb][sz+nb]+=wlt[it];
// 		uvy[sy+nb][sx+nb][sz+nb]+=wlt[it];
// 		uvz[sy+nb][sx+nb][sz+nb]+=wlt[it];
		
		/*add moment tensor source on particle velocity*/
		/*x force*/
/*		uvx[sy+nb  ][sx+nb+1][sz+nb]+=Mxx*h3/dx*wlt[it];
		uvx[sy+nb  ][sx+nb  ][sz+nb]-=Mxx*h3/dx*wlt[it];
		uvx[sy+nb+1][sx+nb  ][sz+nb] =Mxy*h3/dy/4*wlt[it];
		uvx[sy+nb+1][sx+nb+1][sz+nb] =Mxy*h3/dy/4*wlt[it];
		uvx[sy+nb-1][sx+nb  ][sz+nb]-=Mxy*h3/dy/4*wlt[it];
		uvx[sy+nb-1][sx+nb+1][sz+nb]-=Mxy*h3/dy/4*wlt[it];		
		uvx[sy+nb  ][sx+nb  ][sz+nb+1] =Mxz*h3/dz/4*wlt[it];
		uvx[sy+nb  ][sx+nb+1][sz+nb+1] =Mxz*h3/dz/4*wlt[it];
		uvx[sy+nb  ][sx+nb  ][sz+nb-1]-=Mxz*h3/dz/4*wlt[it];
		uvx[sy+nb  ][sx+nb+1][sz+nb-1]-=Mxz*h3/dz/4*wlt[it];*/
		
		
		/*y force*/
/*		uvy[sy+nb+1][sx+nb  ][sz+nb  ]+=Myy*h3/dy*wlt[it];
		uvy[sy+nb  ][sx+nb  ][sz+nb  ]-=Myy*h3/dy*wlt[it];
		uvy[sy+nb  ][sx+nb+1][sz+nb  ] =Myx*h3/dx/4*wlt[it];
		uvy[sy+nb+1][sx+nb+1][sz+nb  ] =Myx*h3/dx/4*wlt[it];
		uvy[sy+nb  ][sx+nb-1][sz+nb  ]-=Myx*h3/dx/4*wlt[it];
		uvy[sy+nb+1][sx+nb-1][sz+nb  ]-=Myx*h3/dx/4*wlt[it];		
		uvy[sy+nb  ][sx+nb  ][sz+nb+1] =Myz*h3/dz/4*wlt[it];
		uvy[sy+nb+1][sx+nb  ][sz+nb+1] =Myz*h3/dz/4*wlt[it];
		uvy[sy+nb  ][sx+nb  ][sz+nb-1]-=Myz*h3/dz/4*wlt[it];
		uvy[sy+nb+1][sx+nb  ][sz+nb-1]-=Myz*h3/dz/4*wlt[it];*/
		
		
		/*z force*/
/*		uvz[sy+nb  ][sx+nb  ][sz+nb+1]+=Mzz*h3/dz*wlt[it];
		uvz[sy+nb  ][sx+nb  ][sz+nb  ]-=Mzz*h3/dz*wlt[it];
		uvz[sy+nb  ][sx+nb+1][sz+nb  ] =Mzx*h3/dx/4*wlt[it];
		uvz[sy+nb  ][sx+nb+1][sz+nb+1] =Mzx*h3/dx/4*wlt[it];
		uvz[sy+nb  ][sx+nb-1][sz+nb  ]-=Mzx*h3/dx/4*wlt[it];
		uvz[sy+nb  ][sx+nb-1][sz+nb+1]-=Mzx*h3/dx/4*wlt[it];			
		uvz[sy+nb+1][sx+nb  ][sz+nb  ] =Mzy*h3/dy/4*wlt[it];
		uvz[sy+nb+1][sx+nb  ][sz+nb+1] =Mzy*h3/dy/4*wlt[it];
		uvz[sy+nb-1][sx+nb  ][sz+nb  ]-=Mzy*h3/dy/4*wlt[it];
		uvz[sy+nb-1][sx+nb  ][sz+nb+1]-=Mzy*h3/dy/4*wlt[it];*/	
						
		/*add moment tensor source on stress*/
		txx[sy+nb  ][sx+nb  ][sz+nb] +=h3*Mxx*wlt[it];
		tyy[sy+nb  ][sx+nb  ][sz+nb] +=h3*Myy*wlt[it];
		tzz[sy+nb  ][sx+nb  ][sz+nb] +=h3*Mzz*wlt[it];
		
		txy[sy+nb  ][sx+nb  ][sz+nb] +=0.25*h3*Mxy*wlt[it];
		txy[sy+nb  ][sx+nb-1][sz+nb] +=0.25*h3*Mxy*wlt[it]; 
		txy[sy+nb-1][sx+nb  ][sz+nb] +=0.25*h3*Mxy*wlt[it];
		txy[sy+nb-1][sx+nb-1][sz+nb] +=0.25*h3*Mxy*wlt[it];
		
		tyz[sy+nb  ][sx+nb  ][sz+nb ] +=0.25*h3*Myz*wlt[it];
		tyz[sy+nb  ][sx+nb  ][sz+nb-1]+=0.25*h3*Myz*wlt[it]; 
		tyz[sy+nb-1][sx+nb  ][sz+nb ] +=0.25*h3*Myz*wlt[it];
		tyz[sy+nb-1][sx+nb  ][sz+nb-1]+=0.25*h3*Myz*wlt[it];
		
		txz[sy+nb  ][sx+nb  ][sz+nb] +=0.25*h3*Mxz*wlt[it];
		txz[sy+nb  ][sx+nb-1][sz+nb] +=0.25*h3*Mxz*wlt[it]; 
		txz[sy+nb  ][sx+nb  ][sz+nb-1] +=0.25*h3*Mxz*wlt[it];
		txz[sy+nb  ][sx+nb-1][sz+nb-1] +=0.25*h3*Mxz*wlt[it];
		
		forward_uvy_uvx_uvz(uvx, uvy, uvz, txx, tyy, tzz, txz, txy, tyz, rho);
		forward_txx_tyy_tzz_txz_txy_tyz(uvx, uvy, uvz, txx, tyy, tzz, txz, txy, tyz, vp, vs);

		apply_sponge(uvz, bndr);
		apply_sponge(uvx, bndr);
		apply_sponge(uvy, bndr);
		apply_sponge(txx, bndr);
		apply_sponge(tyy, bndr);
		apply_sponge(tzz, bndr);
		apply_sponge(txz, bndr);
		apply_sponge(txy, bndr);
		apply_sponge(tyz, bndr);
		
		if(iffree) /*be cautious, problematic*/
		{
		apply_free(tzz,txz,tyz);
		}
		
		if (it%jsnap==0 && ifwfd==1)
		{
			if (verb) sf_warning("writing wavefield at it=%d/%d",it,nt);
			n4++;
			window3d(vp0, uvx);
			sf_floatwrite(vp0[0][0], nz*nx*ny, Fwavx);
			window3d(vs0, uvz);
			sf_floatwrite(vs0[0][0], nz*nx*ny, Fwavz);
			window3d(vs0, uvy);
			sf_floatwrite(vs0[0][0], nz*nx*ny, Fwavy);	
		}
		
		for(iy=0;iy<ny;iy++)
		for(ix=0;ix<nx;ix++)
		{	
			Rvy[iy][ix][it]=uvy[iy+nb][ix+nb][nb+1];/*record at all spatial locations*/
			Rvx[iy][ix][it]=uvx[iy+nb][ix+nb][nb+1];/*record at all spatial locations*/
			Rvz[iy][ix][it]=uvz[iy+nb][ix+nb][nb+1];/*record at the surface*/
		}	
		if (verb) sf_warning("%d of %d;", it, nt);
	}
	
	sf_floatwrite(Rvx[0][0], nt*nx*ny, Frvx);
	sf_floatwrite(Rvz[0][0], nt*nx*ny, Frvz);
	sf_floatwrite(Rvy[0][0], nt*nx*ny, Frvy);
	
	free(wlt);
	free(bndr);
	free(**vp0);free(*vp0); free(vp0);
	free(**vs0); free(*vs0); free(vs0);
	free(**rho0);free(*rho0); free(rho0);
	free(**vp);free(*vp); free(vp);
	free(**vs);free(*vs); free(vs);
	free(**rho);free(*rho); free(rho);
	free(**uvx);free(*uvx); free(uvx);
	free(**uvz); free(*uvz); free(uvz);
	free(**txx);free(*txx); free(txx);
	free(**tzz);free(*tzz); free(tzz);
	free(**txz);free(*txz); free(txz);
	free(**qpx);free(*qpx); free(qpx);
	free(**qsx);free(*qsx); free(qsx);
	free(**qpz);free(*qpz); free(qpz);
	free(**qsz);free(*qsz); free(qsz);		
    exit(0);
}

