/* Anisotropic kernel routines. */

/*
  Copyright (C) 2006 The Board of Trustees of Stanford University
  
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

#include "anikernel.h"
#include "ani.h"
#include "arrayindex.h"

#ifndef _anikernel_h

struct shot_ker_ani_type
{
    float phi,f;
    float o_ep,d_ep,o_dl,d_dl,o_w_vp,d_w_vp;
    int m,n,n_ep,n_dl,n_w_vp;
    sf_complex *contablepp,*contablepm,*conapp,*conapm;
};
struct shot_ker_ani_implicit_type
{
    float o_ep,d_ep,o_dl,d_dl;
    int n_ep,n_dl;
    float ***fdcoe;
    int nfdcoe;
};
/*^*/

#endif

struct shot_ker_ani_type shot_ker_ani_init(float *phi,int nx,float dx,float dy,float dz){
    struct shot_ker_ani_type ani_par;
    sf_file anitable;
    float f,o_ep,d_ep,o_dl,d_dl,o_w_vp,d_w_vp;
    int n_ep,n_dl,n_w_vp,nangle,iangle;
    int m,n;
    float pphi,dkx,oangle,dangle;
    int stat,blocksize;

    pphi=fabsf(*phi);
    anitable=sf_input("Extable");
  
    sf_histfloat(anitable,"o5",&oangle);
    sf_histfloat(anitable,"d5",&dangle);
    sf_histint(anitable,"n5",&nangle);

    iangle=(pphi/SF_PI*180.0-oangle)/dangle; 
    if (iangle >=nangle) iangle=nangle-1;
  
    if (*phi <0 ){
	if ((float)(iangle)*dangle+oangle < 15.0) iangle=0; 
    } else{
	if ((float)(iangle)*dangle+oangle < 15.0) iangle=0;
    }


/*
  if (*phi <0 ){
  if ((float)(iangle)*dangle+oangle < 15.0) iangle=0; //iangle=(20.0-oangle)/dangle;
  }
  else{
  if ((float)(iangle)*dangle+oangle < 26.0) iangle=0;
  }
*/ 


 
    //if ((float)(iangle)*dangle+oangle > 25.0) iangle=(50.0-oangle)/dangle; 
    pphi=(float)(iangle)*dangle+oangle; 
    pphi=pphi/180.0*SF_PI; 
    if ( *phi<0 ) *phi=-pphi; else *phi=pphi; 

    sf_histfloat(anitable,"o2",&ani_par.o_dl);
    sf_histfloat(anitable,"d2",&ani_par.d_dl);
    sf_histint(anitable,"n2",&ani_par.n_dl);

    sf_histfloat(anitable,"o3",&ani_par.o_ep);
    sf_histfloat(anitable,"d3",&ani_par.d_ep);
    sf_histint(anitable,"n3",&ani_par.n_ep);

    sf_histfloat(anitable,"o4",&ani_par.o_w_vp);
    sf_histfloat(anitable,"d4",&ani_par.d_w_vp);
    sf_histint(anitable,"n4",&ani_par.n_w_vp);

    sf_histint(anitable,"n1",&ani_par.n); ani_par.n /= 2;

    ani_par.m=100; ani_par.f=1.0-10.0*10.0/(27.0*27.0);  
    ani_par.phi=pphi; 
    //printf("1111111\n");
    //printf("size of contable: nn=%d\n",ani_par.n*ani_par.n_dl*ani_par.n_ep*ani_par.n_w_vp);
    ani_par.contablepp=sf_complexalloc(ani_par.n*ani_par.n_dl*ani_par.n_ep*ani_par.n_w_vp);
    ani_par.contablepm=sf_complexalloc(ani_par.n*ani_par.n_dl*ani_par.n_ep*ani_par.n_w_vp);
    ani_par.conapp=sf_complexalloc(ani_par.n);  ani_par.conapm=sf_complexalloc(ani_par.n);
    //printf("22222222222\n");
    //printf("aniparnnnnnn=%d\n",ani_par.n);
    blocksize=(ani_par.n*2)*ani_par.n_dl*ani_par.n_ep*ani_par.n_w_vp;

    //printf("iangle=%d\n",iangle);
    sf_seek(anitable,iangle*blocksize*sizeof(sf_complex),SEEK_SET);
    sf_complexread(ani_par.contablepp,blocksize/2,anitable);
    sf_complexread(ani_par.contablepm,blocksize/2,anitable);
    return ani_par;
}


struct shot_ker_ani_type shot_ker_ani_init_impulse(float phi,int nx,float dx,float dy,float dz){
    struct shot_ker_ani_type ani_par;
    float f;
    float o_ep,d_ep,o_dl,d_dl,o_w_vp,d_w_vp;
    int n_ep,n_dl,n_w_vp;
    int m,n;
    float pphi,dkx;

    printf("begin here\n");
    pphi=fabsf(phi);
    f=1.0-10.0*10.0/(27.0*27.0);
    f=0.6913580;
    o_ep=0.03; d_ep=0.01; n_ep=40; o_dl=0.0; d_dl=0.01; n_dl=5; o_w_vp=0.003; d_w_vp=0.001; n_w_vp=200;
    n=5;  m=nx/2; m=150; 
    dkx=2.0*SF_PI/( (float)(nx)*dx );
    dkx=2.0*SF_PI/( (float)(2*m)*dx);

    //o_w_vp=0.041870;//o_w_vp=0.154369; 
    //n_w_vp=2; 
    o_ep=0.4; n_ep=2; o_dl=0.2; n_dl=2;

    //o_ep=0.05; d_ep=0.01; n_ep=30;
    //o_dl=0.0; d_dl=0.01; n_dl=2;
    o_w_vp=0.001; d_w_vp=0.001; n_w_vp=300;

//********************************
    o_w_vp=0.00710831; n_w_vp=224;
    o_ep=0.0-0.01; n_ep=18; d_ep=0.01;
    o_dl=0.0-0.01; n_dl=10;  d_dl=0.01;
//********************************

    ani_par.o_ep=o_ep; ani_par.d_ep=d_ep; ani_par.n_ep=n_ep;
    ani_par.o_dl=o_dl; ani_par.d_dl=d_dl; ani_par.n_dl=n_dl;
    ani_par.o_w_vp=o_w_vp; ani_par.d_w_vp=d_w_vp; ani_par.n_w_vp=n_w_vp;
  
    ani_par.n=n; ani_par.m=m; ani_par.f=f; ani_par.phi=pphi;
    ani_par.contablepp=sf_complexalloc(n*n_dl*n_ep*n_w_vp);
    ani_par.contablepm=sf_complexalloc(n*n_dl*n_ep*n_w_vp);
    ani_par.conapp=sf_complexalloc(n);  ani_par.conapm=sf_complexalloc(n);
    explicit_table(ani_par.contablepp,ani_par.contablepm,ani_par.phi,ani_par.f,dx,dkx,dz,m,n,o_ep,d_ep,n_ep,o_dl,d_dl,n_dl,o_w_vp,d_w_vp,n_w_vp);
    return ani_par;
}


void shot_ker_ani_release(struct shot_ker_ani_type *ani_par){
    free(ani_par->conapp); free(ani_par->conapm); free(ani_par->contablepp); free(ani_par->contablepm);
}



void inadd(sf_complex *app,sf_complex *p,float scale,int n){
    int i;
    for(i=0;i<n;i++)
	app[i]+=scale*p[i];
}

void  get_convlv_coe(float w_vp,float ep,float dl,struct shot_ker_ani_type *ani_par){
    int i_w_vp,i_ep,i_dl;
    float w_w_vp, w_ep,w_dl;
    int n,n_ep,n_dl,n_w_vp,i;
    float w00,w01,w10,w11,w000,w001,w010,w011,w100,w101,w110,w111;
    sf_complex *conapp,*conapm,*contablepp,*contablepm;
    sf_complex *p000,*p001,*p010,*p011,*p100,*p101,*p110,*p111;
    conapp=ani_par->conapp; conapm=ani_par->conapm; contablepp=ani_par->contablepp;  contablepm=ani_par->contablepm;
    n=ani_par->n; n_ep=ani_par->n_ep; n_dl=ani_par->n_dl; n_w_vp=ani_par->n_w_vp;
    i_w_vp=(int)((w_vp-ani_par->o_w_vp)/ani_par->d_w_vp);  w_w_vp=(w_vp-ani_par->o_w_vp)/ani_par->d_w_vp-(float)(i_w_vp);
    i_ep=(int)((ep-ani_par->o_ep)/ani_par->d_ep);          w_ep=(ep-ani_par->o_ep)/ani_par->d_ep - (float)(i_ep);
    i_dl=(int)((dl-ani_par->o_dl)/ani_par->d_dl);          w_dl=(dl-ani_par->o_dl)/ani_par->d_dl-(float)(i_dl);
//printf("i_w_vp=%d,i_ep=%d,i_dl=%d,w_w_vp=%f,w_ep=%f,w_dl=%f\n",i_w_vp,i_ep,i_dl,w_w_vp,w_ep,w_dl);
    if (i_w_vp==ani_par->n_w_vp-1){ i_w_vp=i_w_vp-1; w_w_vp=1.0;}
    if (i_ep==ani_par->n_ep-1){i_ep=i_ep-1;   w_ep=1.0;}
    if (i_dl==ani_par->n_dl-1){  i_dl=i_dl-1;  w_dl=1.0; }
    if (i_w_vp <0 || i_w_vp >ani_par->n_w_vp-1 ) {printf("bad thing\n"); printf("%f,%d\n",w_vp,i_w_vp);}
    if (i_ep <0 || i_ep >ani_par->n_ep-1 ) {printf("bad thing ep\n"); printf("%f,%d\n",ep,i_ep);}
    if (i_dl <0 || i_dl >ani_par->n_dl-1 ) printf("bad thing dl");
    w00=(1.0-w_dl)*(1.0-w_ep); w000=w00*(1.0-w_w_vp); w001=w00*w_w_vp;
    w01=(1.0-w_dl)*w_ep      ; w010=w01*(1.0-w_w_vp); w011=w01*w_w_vp;
    w10=w_dl*(1.0-w_ep)      ; w100=w10*(1.0-w_w_vp); w101=w10*w_w_vp;
    w11=w_dl*w_ep            ; w110=w11*(1.0-w_w_vp); w111=w11*w_w_vp;
//printf("%f,%f,%f,%f,%f,%f,%f,%f\n",w000,w001,w010,w011,w100,w101,w110,w111);
//printf("nnnnnn=%d\n",n);
    p000=i_dl*n+i_ep*n_dl*n+i_w_vp*n_ep*n_dl*n+contablepp; p100=p000+n; p010=p000+n_dl*n; p110=p010+n; 
    p001=p000+n_ep*n_dl*n; p101=p001+n; p011=p001+n_dl*n; p111=p011+n;
    p111=contablepp+(i_dl+1)*n+(i_ep+1)*n_dl*n+(i_w_vp+1)*n_ep*n_dl*n;
    vector_value_c(conapp,sf_cmplx(0.0,0.0),ani_par->n);
    inadd(conapp,p000,w000,n); inadd(conapp,p100,w100,n); inadd(conapp,p010,w010,n); inadd(conapp,p110,w110,n); 
    inadd(conapp,p001,w001,n); inadd(conapp,p101,w101,n); inadd(conapp,p011,w011,n); inadd(conapp,p111,w111,n);
//for(i=0;i<n;i++)
//   printf("here0 conapp[%d]=(%f,%f)",i,real(conapp[i]),__imag__ conapp[i]); 
//printf("\n");
 

    vector_value_c(conapm,sf_cmplx(0.0,0.0),ani_par->n);
    p000=i_dl*n+i_ep*n_dl*n+i_w_vp*n_ep*n_dl*n+contablepm; p100=p000+n; p010=p000+n_dl*n; p110=p010+n; 
    p001=p000+n_ep*n_dl*n; p101=p001+n; p011=p001+n_dl*n; p111=p011+n;
    inadd(conapm,p000,w000,n); inadd(conapm,p100,w100,n); inadd(conapm,p010,w010,n); inadd(conapm,p110,w110,n); 
    inadd(conapm,p001,w001,n); inadd(conapm,p101,w101,n); inadd(conapm,p011,w011,n); inadd(conapm,p111,w111,n);
/*
  conapp=contablepp(:,i_dl,i_ep,i_w_vp)*w000+contablepp(:,i_dl+1,i_ep,i_w_vp)*w100+contablepp(:,i_dl,i_ep+1,i_w_vp)*w010+contablepp(:,i_dl+1,i_ep+1,i_w_vp)*w110+contablepp(:,i_dl,i_ep,i_w_vp+1)*w001+contablepp(:,i_dl+1,i_ep,i_w_vp+1)*w101+contablepp(:,i_dl,i_ep+1,i_w_vp+1)*w011+contablepp(:,i_dl+1,i_ep+1,i_w_vp)*w111;
  conapm=contablepm(:,i_dl,i_ep,i_w_vp)*w000+contablepm(:,i_dl+1,i_ep,i_w_vp)*w100+contablepm(:,i_dl,i_ep+1,i_w_vp)*w010+contablepm(:,i_dl+1,i_ep+1,i_w_vp)*w110+contablepm(:,i_dl,i_ep,i_w_vp+1)*w001+contablepm(:,i_dl+1,i_ep,i_w_vp+1)*w101+contablepm(:,i_dl,i_ep+1,i_w_vp+1)*w011+contablepm(:,i_dl+1,i_ep+1,i_w_vp)*w111;
*/
}




void shot_ani_convlv3dtest(sf_complex *wld_z,sf_complex *conapp,sf_complex *conapm,int tnx,int tny,struct shot_ker_par_type ker_par, int signn)
{
    int nx,ny,tn;
    sf_complex  *tmpwld_z; //(:,:)
    sf_complex *conap,*conam;
    int  iy,ix,ic,i,icx,icy;
    tn=tnx*tny;
    nx=ker_par.nx; ny=ker_par.ny;
    tmpwld_z=ker_par.tmp_wld;   //allocate(tmpwld_z(nx,ny)) ker_par.tmp_wld is allocated in ker_init in kernal.c
    conap=sf_complexalloc(tn); conam=sf_complexalloc(tn);
    vector_cp_c(tmpwld_z,wld_z,nx*ny); //tmpwld_z=wld_z;
    vector_cp_c(conap,conapp,tn); vector_cp_c(conam,conapm,tn);
    vector_value_c(wld_z,sf_cmplx(0.0,0.0),nx*ny);
    if (signn==+1) {
	for(ic=0;ic<tn;ic++){
	    conap[ic]=conjg(conapp[ic]); conam[ic]=conjg(conapm[ic]);
	} //for(ic)
    }
    else{
	for(ic=0;ic<tn;ic++){
	    conap[ic]=(conapp[ic]); conam[ic]=(conapm[ic]);
	} //for(ic)
    }
    for(icy=1;icy<tny;icy++){
	for(icx=0; icx<tnx; icx++){
	    ic=icy*tnx+icx;
	    conap[ic]=0.5*conap[ic];
	    conam[ic]=0.5*conam[ic];
	}
    }

//  for (ic=0;ic<tn;ic++){
//    printf("ic=%d,(%f,%f)",ic,__real__ conap[ic],__imag__ conap[ic]);
//  }
    printf("begin convlv\n");

    for( iy=0;iy<ny;iy++){
	for(ix=0;ix<nx;ix++){
	    wld_z[i2(iy,ix,nx)]=conap[0]*tmpwld_z[i2(iy,ix,nx)]; 
                            
	    for(icy=1,icx=0; icy<tny;icy++){ 
		ic=icy*tnx+icx;    
		if (iy+icy <ny)    
		    wld_z[i2(iy,ix,nx)]+=conap[ic]*(tmpwld_z[i2(iy+icy,ix,nx)]);
		if (iy-icy >=0)
		    wld_z[i2(iy,ix,nx)]+=conap[ic]*(tmpwld_z[i2(iy-icy,ix,nx)]);
	    }

	    for(icy=0,icx=1;icx<tnx; icx++){
		ic=icy*tnx+icx;
		if (ix+icx <nx) wld_z[i2(iy,ix,nx)]+=conap[ic]*tmpwld_z[i2(iy,ix+icx,nx)];
		if (ix-ic >= 0 )wld_z[i2(iy,ix,nx)]+=conam[ic]*tmpwld_z[i2(iy,ix-icx,nx)];
	    }
	    for(icy=1;icy<tny; icy++){
		if (iy+icy <ny && iy-icy >=0 ){
		    for(icx=1;icx<tnx; icx++){
			ic=icy*tnx+icx;
			if (ix+icx <nx) wld_z[i2(iy,ix,nx)]+=conap[ic]*(tmpwld_z[i2(iy+icy,ix+icx,nx)]+tmpwld_z[i2(iy-icy,ix+icx,nx)]);
			if (ix-icx >=0) wld_z[i2(iy,ix,nx)]+=conam[ic]*(tmpwld_z[i2(iy+icy,ix-icx,nx)]+tmpwld_z[i2(iy-icy,ix-icx,nx)]);
		    }
		}
		if (iy+icy<ny && iy-icy <0){
		    for(icx=1;icx<tnx; icx++){
			ic=icy*tnx+icx;
			if (ix+icx <nx) wld_z[i2(iy,ix,nx)]+=conap[ic]*(tmpwld_z[i2(iy+icy,ix+icx,nx)]);
			if (ix-icx >=0) wld_z[i2(iy,ix,nx)]+=conam[ic]*(tmpwld_z[i2(iy+icy,ix-icx,nx)]);
		    }
		}
		if (iy+icy >=ny && iy-icy >=0){
		    for(icx=1;icx<tnx; icx++){
			ic=icy*tnx+icx;
			if (ix+icx <nx) wld_z[i2(iy,ix,nx)]+=conap[ic]*(tmpwld_z[i2(iy-icy,ix+icx,nx)]);
			if (ix-icx >=0) wld_z[i2(iy,ix,nx)]+=conam[ic]*(tmpwld_z[i2(iy-icy,ix-icx,nx)]);
		    }
		}

	    } 
     
	}// for(ix) 
    }// for(iy)
    free(conap); free(conam);
  
}



void shot_phase_shift_ani(sf_complex *wld_z,float w_v0_z,float ep,float dl,float phi,float f,struct shot_ker_par_type ker_par,int signn)
{// wld_z(ker_par.ny,ker_par.nx)
    int nx,ny,weitmx1,weitmx2,weitmy1,weitmy2;
    float fkx,fky,dkx,dky,dz;
  
    int ix,iy,ipi,isin,ikz2;
    float kx,ky,kz,kz2,phsft,pio2,wpi,tphsft,wsin,sin_phsft,cos_phsft;
    sf_complex kz1;
    sf_complex *kza;
    pio2=SF_PI*2.0;
    nx=ker_par.nx; ny=ker_par.ny; 
    fkx=ker_par.fkx; fky=ker_par.fky; 
    dkx=ker_par.dkx; dky=ker_par.dky;
    dz=ker_par.dz;
  
    kza=sf_complexalloc(nx*ny);
    for(iy=0;iy<ny; iy++){
	ky=fky+iy*dky;
	for(ix=0;ix<nx;ix++){
	    kx=fkx+ix*dkx;
	    kza[i2(iy,ix,nx)]=kzani3d(phi,ep,dl,w_v0_z,f,kx,ky);

	    if ( fabsf(aimag( kza[i2(iy,ix,nx)]  ))<0.00000000000001){
		if ( fabsf(aimag( kza[i2(iy,ix,nx)]  )) <0.00000000001*fabsf(real( kza[i2(iy,ix,nx)]  )) )
		    kza[i2(iy,ix,nx)]=sf_cmplx(fabsf(real(kza[i2(iy,ix,nx)])),0.0);    
		else
		    kza[i2(iy,ix,nx)]=sf_cmplx(0.0,-fabsf(aimag(kza[i2(iy,ix,nx)])));
	    }
	    else
		kza[i2(iy,ix,nx)]=sf_cmplx(0.0,-fabsf(aimag(kza[i2(iy,ix,nx)])));
	}
    }  

    for (iy=0;iy<ny;iy++){
	for(ix=0;ix<nx;ix++){
	    if (fabsf(__imag__ kza[i2(iy,ix,nx)])<0.000001){
		if ( fabsf(__imag__ kza[i2(iy,ix+1,nx)])>0.000001 && fabsf(__imag__ kza[i2(iy,ix-1,nx)])>0.000001   )
		    kza[i2(iy,ix,nx)]=sf_cmplx( 0.0,1.0);
		if ( fabsf(__imag__ kza[i2(iy+1,ix,nx)])>0.000001 && fabsf(__imag__ kza[i2(iy-1,ix,nx)])>0.000001   )
		    kza[i2(iy,ix,nx)]=sf_cmplx( 0.0,1.0);

	    }
	}
    }

    iy=ny/2;
    for(ix=nx/2;ix<nx;ix++)
	if (fabsf(__imag__ kza[i2(iy,ix,nx)])>0.001 && fabsf(__real__ kza[i2(iy,ix,nx)])<0.00000001) break;
    weitmx2=ix;
    for(ix=0;ix<nx/2;ix++)
	if (fabsf(__imag__ kza[i2(iy,ix,nx)])>0.001 && fabsf(__real__ kza[i2(iy,ix,nx)])<0.00000001) break;
    weitmx1=ix; 
    ix=nx/2;
    for(iy=ny/2;iy<ny;iy++)
	if (fabsf(__imag__ kza[i2(iy,ix,nx)])>0.001 && fabsf(__real__ kza[i2(iy,ix,nx)])<0.00000001) break;
    weitmy2=iy;
    for(iy=0;iy<ny/2;iy++)
	if (fabsf(__imag__ kza[i2(iy,ix,nx)])>0.001 && fabsf(__real__ kza[i2(iy,ix,nx)])<0.00000001) break;
    weitmy1=iy;
    for (iy=0;iy<ny;iy++){
	for(ix=0;ix<nx;ix++){
	    if (iy<weitmy1 || iy>weitmy2)
		kza[i2(iy,ix,nx)]=sf_cmplx( 0.0,1.0);
	    if (ix<weitmx1 || ix>weitmx2)
		kza[i2(iy,ix,nx)]=sf_cmplx( 0.0,1.0);
	}
    }



    for(iy=0;iy<ny;iy++){ //do iy=1,ny
	if (ker_par.is2d) 
	    ky=0.0;
	else
	    ky=fky+iy*dky;
	for(ix=0;ix<nx;ix++){  //do ix=1,nx
	    kx=fkx+ix*dkx;
	    //kz1=kzani3d(phi,ep,dl,w_v0_z,f,kx,ky);
	    kz1=kza[i2(iy,ix,nx)];
	    if ( fabsf(aimag(kz1)) <0.0000001*fabsf(real(kz1)) ){
		phsft=-signn*dz*fabsf(real(kz1));
		wld_z[i2(iy,ix,nx)]=wld_z[i2(iy,ix,nx)]*sf_cmplx(cos(phsft),sin(phsft));
	    }
	    else{
		wld_z[i2(iy,ix,nx)]=sf_cmplx(0.0,0.0);
		//phsft=dz*fabsf(aimag(kz1));  
		//wld_z[i2(iy,ix,nx)]=wld_z[i2(iy,ix,nx)]*exp(-phsft);
	    }
      
	} 
    }


    free(kza);
}


void shot_phase_shift_vti(sf_complex *wld_z,float w_v0_z,float ep,float dl,struct shot_ker_par_type ker_par,int signn)
{// wld_z(ker_par.ny,ker_par.nx)
    int nx,ny,weitmx1,weitmx2,weitmy1,weitmy2;
    float fkx,fky,dkx,dky,dz;
  
    int ix,iy,ipi,isin,ikz2;
    float kx,ky,kz,kz2,phsft,pio2,wpi,tphsft,wsin,sin_phsft,cos_phsft;
    sf_complex kz1;
    sf_complex *kza;
    pio2=SF_PI*2.0;
    nx=ker_par.nx; ny=ker_par.ny; 
    fkx=ker_par.fkx; fky=ker_par.fky; 
    dkx=ker_par.dkx; dky=ker_par.dky;
    dz=ker_par.dz;
    iy=0;
    for (ix=0;ix<nx;ix++){
	kx=fkx+ix*dkx;
	kx=kx/w_v0_z;
	kz2=(1-kx*kx*(1.0+2.0*ep))/(1.0-kx*kx*2.0*(ep-dl));
	if (kz2>=0){
	    phsft=-signn*dz*w_v0_z*(sqrt(kz2));
	    wld_z[i2(iy,ix,nx)]=wld_z[i2(iy,ix,nx)]*sf_cmplx(cos(phsft),sin(phsft));
	}
	else{
	    if (kz2<-1) kz2=-1;
	    phsft=dz*w_v0_z*(sqrt(-kz2));
	    //wld_z[i2(iy,ix,nx)]=wld_z[i2(iy,ix,nx)]*exp(-phsft);
	    wld_z[i2(iy,ix,nx)]=0.0;
	}
    }
  
}





