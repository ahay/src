/* Kernel routines. */

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

#include "kernel.h"
#include "arrayindex.h"
#include "fft.h"
#include "mymath.h"
#include "anikernel.h"
#include "ani.h"

#ifndef _kernel_h

struct shot_ker_par_type
{
    int op,nx,ny,nref;
    float dx,dy,dz,fkx,dkx,fky,dky;
    float ax,ay,bx,by,axnoa,aynoa,bxnob,bynob,tx,ty,fda,fdb;
    int is2d;
    sf_complex *ux,*vx,*uy,*vy,*ra,*rb,*rc;
    float *cax,*cbx,*cay,*cby;
    float *axvector,*bxvector,*avector,*bvector;
    sf_complex *tmp_fft, *add_wld,*tmp_wld;
    float *v0_z,*v_zw,*v_zw0;
    float *sin,*cos,*sqrt;
    float dsin,dqrt;
    int Nsin;
    float rotate;
};
/*^*/

#endif

struct shot_ker_par_type shot_ker_p_init(int op,int nx,int ny,float dx,float dy,float dz,
					 float a,float b,float trick,int is2d,float rotate)
/*< initialize >*/
{
    struct shot_ker_par_type ker_par;


    int i;
    ker_par.nref=3;
    ker_par.op=op;
    ker_par.fda=a; ker_par.fdb=b;
    ker_par.nx=nx; ker_par.ny=ny;
    ker_par.dx=dx; ker_par.dy=dy; ker_par.dz=dz;
    ker_par.fkx=-SF_PI/ker_par.dx; ker_par.fky=-SF_PI/ker_par.dy;  //change tilt 
    ker_par.dkx=2.0*SF_PI/(nx*ker_par.dx); ker_par.dky=2.0*SF_PI/(ny*ker_par.dy);  //change in tilt
    ker_par.ax=a*dz/(2.0*ker_par.dx*ker_par.dx); ker_par.ay=a*dz/(2.0*ker_par.dy*ker_par.dy);
    ker_par.bx=b/(ker_par.dx*ker_par.dx); ker_par.by=b/(ker_par.dy*ker_par.dy);
    ker_par.axnoa=dz/(2.0*ker_par.dx*ker_par.dx); ker_par.aynoa=dz/(2.0*ker_par.dy*ker_par.dy);
    ker_par.bxnob=1.0/(ker_par.dx*ker_par.dx); ker_par.bynob=1.0/(ker_par.dy*ker_par.dy);
    ker_par.tx=trick;     ker_par.ty=trick;
    ker_par.is2d=is2d;
    ker_par.axvector=sf_floatalloc(nx); ker_par.bxvector=sf_floatalloc(nx);
    ker_par.avector=sf_floatalloc(nx); ker_par.bvector=sf_floatalloc(nx);
    ker_par.ux=sf_complexalloc(nx);  ker_par.vx=sf_complexalloc(nx); 
    ker_par.cax=sf_floatalloc(nx); ker_par.cbx=sf_floatalloc(nx);
    ker_par.uy=sf_complexalloc(ny);  ker_par.vy=sf_complexalloc(ny); 
    ker_par.cay=sf_floatalloc(ny); ker_par.cby=sf_floatalloc(ny);
    ker_par.ra=sf_complexalloc(nx);  ker_par.rb=sf_complexalloc(nx); ker_par.rc=sf_complexalloc(nx);
    ker_par.tmp_fft=sf_complexalloc(2*nx); ker_par.add_wld=sf_complexalloc(nx*ny); ker_par.tmp_wld=sf_complexalloc(nx*ny);
    //printf ("nx in initial kernel:nx=%d\n",nx);
    ker_par.v_zw=sf_floatalloc(nx*ny); ker_par.v_zw0=sf_floatalloc(nx*ny);
    ker_par.v0_z=sf_floatalloc(ker_par.nref);
    ker_par.rotate=rotate;

    ker_par.Nsin=2000;
    ker_par.sin=sf_floatalloc(ker_par.Nsin); 
    ker_par.cos=sf_floatalloc(ker_par.Nsin); 

    ker_par.sqrt=sf_floatalloc(ker_par.Nsin);
    ker_par.dsin=SF_PI*2.0/(float)(ker_par.Nsin);
    ker_par.dqrt=1.0/(float)(ker_par.Nsin);

    for (i=0;i<ker_par.Nsin;i++){
	ker_par.sin[i]=sin(ker_par.dsin*i);
	ker_par.cos[i]=cos(ker_par.dsin*i);
	ker_par.sqrt[i]=sqrt(ker_par.dqrt*i);
    }
  
    return ker_par;
}



void shot_ker_release(struct shot_ker_par_type *ker_par)
/*< release >*/
{
    free(ker_par->ux); 
    free(ker_par->vx); free(ker_par->cax); free(ker_par->cbx);
    free(ker_par->uy); free(ker_par->vy); free(ker_par->cay); free(ker_par->cby);
    free(ker_par->ra); free(ker_par->rb); free(ker_par->rc); 
    free(ker_par->tmp_fft);free(ker_par->add_wld); free(ker_par->tmp_wld);
    free(ker_par->v_zw); free(ker_par->v_zw0); free(ker_par->v0_z);
    free(ker_par->sin); 
    free(ker_par->cos); 
    free(ker_par->sqrt);
    free(ker_par->axvector); free(ker_par->bxvector); free(ker_par->avector); free(ker_par->bvector);
}

void tilt_nx_dkx_change_iz(struct shot_ker_par_type *ker_par, int nx_iz, float dkx_iz)
{
    ker_par->dkx=dkx_iz; ker_par->nx=nx_iz;
}

void shot_ker_depth_onestep(sf_complex *wld_z,float *v_z,float w,int signn, struct shot_ker_par_type ker_par)
{
    int op,nref,nx,ny,is2d;
    float wvx,dref;
    sf_complex  *tmp_fft,*add_wld,*tmp_wld;
    int iref,ix,iy,ivx;
    float w_v0_z,maxvz,minvz;
    float *v0_z,*v_zw,*v_zw0;
    nref=ker_par.nref;  nx=ker_par.nx;  ny=ker_par.ny;  is2d=ker_par.is2d;  op=ker_par.op; //! op:1:ssf 2:fd 3:ffd
    v0_z=ker_par.v0_z; v_zw=ker_par.v_zw; v_zw0=ker_par.v_zw0;
//old tmp_fft=ker_par.tmp_fft; 
    add_wld=ker_par.add_wld; tmp_wld=ker_par.tmp_wld;
    tmp_fft=sf_complexalloc(nx*2);
//old !v0_z(1)=minval(v_z);//!v0_z(nref)=maxval(v_z);//!dref=(v0_z(nref)-v0_z(1))/(nref-1)



    maxvz=maxval(v_z,nx*ny); minvz=minval(v_z,nx*ny);  dref=(maxvz-minvz)/(nref+1);
    if ( fabs(maxvz-minvz)<50.0 ) { ///?????????????wait
	nref=1;  v0_z[0]=0.5*(maxvz+minvz);
    }
    else{
	v0_z[0]=minvz;  v0_z[nref-1]=maxvz;  dref=(v0_z[nref-1]-v0_z[0])/(float)(nref-1); v0_z[1]=(minvz+maxvz)*0.5;  
	//old v0_z(2)=(minval(v_z)+maxval(v_z))*0.5; //!  do iref=1,nref //!    v0_z(iref)=iref*dref //!  end do
    }


//old!v0_z(2)=0.5*v0_z(1)+0.5*v0_z(3);//!v0_z(1)=0.5*v0_z(1)+0.5*v0_z(2);//!v0_z(3)=0.5*v0_z(3)+0.5*v0_z(2)
//old  v0_z[0]=0.5*(v0_z[0]+v0_z[1]); v0_z[2]=0.5*(v0_z[1]+v0_z[2]); dref=v0_z[1]-v0_z[0];
//old nref=1;  v0_z[0]=minvz;
//old printf("nx=%d,ny=%d\n",nx,ny); 



    for (iy=0;iy<ny;iy++)
	for(ix=0;ix<nx;ix++)
	    v_zw[i2(iy,ix,nx)]=v_z[i2(iy,ix,nx)]/w;

    if (op == 2){
   
	w_v0_z=0.0;
	shot_ker_ssf(wld_z,w_v0_z,w,v_z,ker_par,signn);
	shot_ker_fd(wld_z,v_zw,v_zw,ker_par,signn);
      

   
	if (is2d){
	    vector_value_c(tmp_fft,sf_cmplx(0.0,0.0),nx);
	    rowc(wld_z,tmp_fft,0,nx,1);  //tmp_fft=wld_z(1,:);
	    kissfft(tmp_fft,nx,+1);
	    rowc(wld_z,tmp_fft,0,nx,0);  //wld_z(1,:)=tmp_fft(:);
	}
	w_v0_z=w/(0.75*minvz+0.25*maxvz);
	li_filter_new(wld_z,w_v0_z,ker_par,signn); 
	if (is2d){
	    rowc(wld_z,tmp_fft,0,nx,1); //tmp_fft=tmp_wld(:,1);
	    kissfft(tmp_fft,nx,-1);
	    rowc(wld_z,tmp_fft,0,nx,0);  //tmp_wld(:,1)=tmp_fft(:);
	}
   

    }
    else{
	if (is2d){
	    //printf("nx in fft nx=%d\n",nx);
	    vector_value_c(tmp_fft,sf_cmplx(0.0,0.0),nx);
	    //printf("in before iso max0=%f\n",maxvalc(wld_z,nx));
	    rowc(wld_z,tmp_fft,0,nx,1);  //tmp_fft=wld_z(1,:);
	    kissfft(tmp_fft,nx,+1);
	    rowc(wld_z,tmp_fft,0,nx,0);  //wld_z(1,:)=tmp_fft(:);
	    //printf("in iso max0=%f\n",maxvalc(wld_z,nx));
	}
	else{
	    //printf ("bbbbb iam 3d");
	    for(iy=0;iy<ny;iy++){
		rowc(wld_z,tmp_fft,iy,nx,1);  //tmp_fft=wld_z(1,:);
		kissfft(tmp_fft,nx,+1);
		rowc(wld_z,tmp_fft,iy,nx,0);  //wld_z(1,:)=tmp_fft(:);
	    }
	    if (ny!=1){
		for(ix=0;ix<nx;ix++){
		    colc(wld_z,tmp_fft,ix,nx,ny,1);
		    kissfft(tmp_fft,ny,+1);
		    colc(wld_z,tmp_fft,ix,nx,ny,0);   // wait tmp_fftx tmp_ffty
		}
	    }
  
	}
  
  
	for( iy=0;iy<ny;iy++)
	    for(ix=0;ix<nx;ix++)
		add_wld[i2(iy,ix,nx)]=sf_cmplx(0.0,0.0);   //add_wld=0.0;

	for( iref=0;iref<nref;iref++){
	    matrix_equ_c(tmp_wld,wld_z,ny,nx);//tmp_wld=wld_z;
	    w_v0_z=w/v0_z[iref];
	    shot_phase_shift(tmp_wld,w_v0_z,ker_par,signn);
	    //shot_phase_shift_vti(tmp_wld,w_v0_z,0.4,0.2,ker_par,signn);
    

	    if (is2d){
		rowc(tmp_wld,tmp_fft,0,nx,1); //tmp_fft=tmp_wld(:,1);
		kissfft(tmp_fft,nx,-1);
		rowc(tmp_wld,tmp_fft,0,nx,0);  //tmp_wld(:,1)=tmp_fft(:);
	    }
	    else{
		if (ny!=1){
		    for(ix=0;ix<nx;ix++){
			colc(tmp_wld,tmp_fft,ix,nx,ny,1);
			kissfft(tmp_fft,ny,-1);
			colc(tmp_wld,tmp_fft,ix,nx,ny,0);   // wait tmp_fftx tmp_ffty
		    }
		}
		for(iy=0;iy<ny;iy++){
		    rowc(tmp_wld,tmp_fft,iy,nx,1);  //tmp_fft=wld_z(1,:);
		    kissfft(tmp_fft,nx,-1);
		    rowc(tmp_wld,tmp_fft,iy,nx,0);  //wld_z(1,:)=tmp_fft(:);
		}

	    }
	    shot_ker_ssf(tmp_wld,w_v0_z,w,v_z,ker_par,signn);
    
	    if (op==3){ 
		for( iy=0;iy<ny;iy++)
		    for(ix=0;ix<nx;ix++)
			v_zw0[i2(iy,ix,nx)]=v_zw[i2(iy,ix,nx)]-v0_z[iref]/w;
		shot_ker_fd(tmp_wld,v_zw0,v_zw,ker_par,signn);
	    }
    
	    if ( nref==1) vector_cp_c(add_wld,tmp_wld,ny*nx); //add_wld=tmp_wld;
	    else{
      
		for(iy=0;iy<ny;iy++){      
		    for(ix=0;ix<nx;ix++){  
			ivx=(int)((v_z[i2(iy,ix,nx)]-v0_z[0])/dref); //wait ???
			wvx=(v_z[i2(iy,ix,nx)]-(v0_z[0]+(float)(ivx)*dref))/dref;   //wait ???
			if (v_z[i2(iy,ix,nx)] >= v0_z[nref-1]){
			    if (iref==nref-1) add_wld[i2(iy,ix,nx)]+=tmp_wld[i2(iy,ix,nx)];
			}
			else{
			    if (v_z[i2(iy,ix,nx)] <= v0_z[0]){ 
				if (iref==0)  add_wld[i2(iy,ix,nx)]+=tmp_wld[i2(iy,ix,nx)];
			    }
			    else{
				if (ivx==iref) add_wld[i2(iy,ix,nx)] += tmp_wld[i2(iy,ix,nx)]*(1-wvx);
				if (ivx+1==iref) add_wld[i2(iy,ix,nx)] += tmp_wld[i2(iy,ix,nx)]*wvx;
			    }
			}
		    }
		}
      
	    }
    
	}
	vector_cp_c(wld_z,add_wld,ny*nx);  //wld_z=add_wld;
    }
    free(tmp_fft);
}



void li_filter(sf_complex *wld_z,float w_v0_z,struct shot_ker_par_type ker_par,int signn)
{
    int nx;
    float fkx,dkx,dz;
    int ix,iy;
    float kx,kz2,phsft,kx2;
    float a,b;
    nx=ker_par.nx; 
    fkx=ker_par.fkx;
    dkx=ker_par.dkx; 
    dz=ker_par.dz; a=ker_par.fda; b=ker_par.fdb;
    for(iy=0;iy<1;iy++){
	for(ix=0;ix<nx;ix++){  
	    kx=fkx+ix*dkx;
	    kx2=(kx*kx)/(w_v0_z*w_v0_z);
	    if (kx2 <=1) {
		kz2=1.0-kx2;
		phsft=-signn*dz*w_v0_z*(sqrt(kz2)-(1-(a*kx2)/(1-b*kx2) ));
		wld_z[i2(iy,ix,nx)]=wld_z[i2(iy,ix,nx)]*sf_cmplx(cos(phsft),sin(phsft));
	    }
	    else{
		if (kx2>2) kx2=2.0;
		kz2=fabs(kx2- 1);
		phsft=dz*w_v0_z*(sqrt(kz2));
		wld_z[i2(iy,ix,nx)]=0.0;
	    }
	}
    }
}


void li_filter_new(sf_complex *wld_z,float w_v0_z,struct shot_ker_par_type ker_par,int signn)
/*< LI filter >*/
{
    int nx;
    float fkx,dkx,dz;
    int ix,iy;
    float kx,kx2;
    float a,b;
    float deltx2,dx,dx2,aterm,bterm,tr;
    sf_complex leftterm,rightterm;

    nx=ker_par.nx; 
    fkx=ker_par.fkx; 
    dkx=ker_par.dkx; 
    dz=ker_par.dz; a=ker_par.fda; b=ker_par.fdb;
    dx=ker_par.dx; dx2=dx*dx; tr=ker_par.tx;
    for(iy=0;iy<1;iy++){
	for(ix=0;ix<nx;ix++){  
	    kx=fkx+ix*dkx;
	    kx2=(kx*kx)/(w_v0_z*w_v0_z);
	    if (kx2 <=1) {
		deltx2=2.0*(cos(kx*dx)-1.0)/dx2;
		bterm=tr*dx2+b/(w_v0_z*w_v0_z);
		aterm=-a*dz/(2.0*w_v0_z);
		leftterm=1.0+sf_cmplx(bterm,-aterm)*deltx2;
		rightterm=1.0+sf_cmplx(bterm,aterm)*deltx2; 
		wld_z[i2(iy,ix,nx)]=wld_z[i2(iy,ix,nx)]*(rightterm/leftterm);

	    } else{
		if (kx2>2) kx2=2.0;
		//wld_z[i2(iy,ix,nx)]=wld_z[i2(iy,ix,nx)]*exp(-phsft);
		wld_z[i2(iy,ix,nx)]=0.0;  
	    }
	}
    }
}

void shot_phase_shift(sf_complex *wld_z,float w_v0_z,struct shot_ker_par_type ker_par,int signn)
/*< phase shift >*/
{
    int nx,ny;
    float fkx,fky,dkx,dky,dz;
    int ix,iy;
    float kx,ky,kz2,phsft;
    nx=ker_par.nx; ny=ker_par.ny; 
    fkx=ker_par.fkx; fky=ker_par.fky; 
    dkx=ker_par.dkx; dky=ker_par.dky;
    dz=ker_par.dz;

    for(iy=0;iy<ny;iy++){ 
	if (ker_par.is2d) 
	    ky=0.0;
	else
	    ky=fky+iy*dky;
	if (ny==1) ky=0.0;
 
	for(ix=0;ix<nx;ix++){  
	    kx=fkx+ix*dkx;
	    kz2=(kx*kx+ky*ky)/(w_v0_z*w_v0_z);
	    if (kz2 <=1.0) {
		kz2=1.0-kz2;
		phsft=-signn*dz*w_v0_z*(sqrt(kz2));
		wld_z[i2(iy,ix,nx)]=wld_z[i2(iy,ix,nx)]*sf_cmplx(cos(phsft),sin(phsft));
	    }
	    else{
		if (kz2>2) kz2=2.0;
		kz2=kz2-1;
		phsft=dz*w_v0_z*(sqrt(kz2));
		wld_z[i2(iy,ix,nx)]=wld_z[i2(iy,ix,nx)]*exp(-phsft);
	    }     
	} 
    }
}


void shot_ker_ssf(sf_complex *wld_z,float w_v0_z,float w,float *v_z, struct shot_ker_par_type ker_par,int signn)
/*< SSF >*/
{   
    int nx,ny,ix,iy;
    float phsft,dz;

    nx=ker_par.nx; ny=ker_par.ny; dz=ker_par.dz;
  
    for(iy=0;iy<ny;iy++){
	for(ix=0;ix<nx;ix++){
	    phsft=-signn*dz*(w/v_z[i2(iy,ix,nx)]-w_v0_z);
	    wld_z[i2(iy,ix,nx)]=wld_z[i2(iy,ix,nx)]*sf_cmplx(cos(phsft),sin(phsft));
	}
    }
}

void shot_ker_fd(sf_complex *wld_z,float *ca,float *cb,struct shot_ker_par_type ker_par,int signn)
/*< FD >*/
{

    int nx,ny,iy;
    float tx,ax,bx;
    sf_complex *ux,*vx;
    float *cax,*cbx;

    ax=-signn*ker_par.ax; 
    bx=ker_par.bx; 
    nx=ker_par.nx; ny=ker_par.ny;
    tx=ker_par.tx;

    ux=ker_par.ux; vx=ker_par.vx; cax=ker_par.cax; cbx=ker_par.cbx;

    for(iy=0;iy<ny;iy++){
	rowc(wld_z,ux,iy,nx,1);  //ux=wld_z(iy,:)
	rowf(ca,cax,iy,nx,1);   //  cax=ca(iy,:)
	rowf(cb,cbx,iy,nx,1);   //  cbx=cb(iy,:)
	fd(ux,vx,tx,ax,bx,cax,cbx,nx,ker_par); //,rax,rbx,rcx);
	rowc(wld_z,vx,iy,nx,0);//  wld_z(iy,:,:)=vx
    }
}


void fd(sf_complex *u,sf_complex *v,float tr,float a,float b,float *ca,float *cb,int nx,
	struct shot_ker_par_type ker_par)
/*< finite differences >*/
{ 
    sf_complex *ra,*rb,*rc;
    sf_complex  url,urr,ror,rol,rorr,rorl;
    int  ix;

    ra=ker_par.ra; rb=ker_par.rb; rc=ker_par.rc;
    for(ix=0;ix<nx;ix++){
	if (ix==0)    url=sf_cmplx(0.0,0.0);
	else          url=u[ix-1];
	if (ix==nx-1) urr=sf_cmplx(0.0,0.0);
	else          urr=u[ix+1];
	ror=sf_cmplx(tr+b*cb[ix]*cb[ix],a*ca[ix]);
	if (ix==nx-1) rorr=ror;
	else rorr=sf_cmplx(tr+b*cb[ix]*cb[ix+1],a*ca[ix]);
	if (ix==0) rorl=ror;
	else rorl=sf_cmplx(tr+b*cb[ix]*cb[ix-1],a*ca[ix]);
	v[ix]=u[ix]-2.0*ror*u[ix]+rorr*urr+rorl*url;
    

	rol=conjf(ror);
	ra[ix]=1.0-2.0*rol;
	if (ix <nx-1) rb[ix] =sf_cmplx(tr+b*cb[ix]*cb[ix+1],-a*ca[ix]);
	if (ix >0) rc[ix-1]=sf_cmplx(tr+b*cb[ix]*cb[ix-1],-a*ca[ix]);
    }
    thryj(ra,rb,rc,v,nx,1);
}




void shot_ker_depth_onestep_ani_phaseshift(sf_complex *wld_z,float *v_z,float w,float ep,float dl,float phi, float f,int signn, struct shot_ker_par_type ker_par){
    int op,nref,nx,ny,is2d;

    float wvx,dref;
    sf_complex  *tmp_fft,*add_wld,*tmp_wld;
    int iref,ix,iy,ivx;
    float w_v0_z,maxvz,minvz;
    float *v0_z,*v_zw,*v_zw0;
    nref=ker_par.nref;  nx=ker_par.nx;  ny=ker_par.ny;  is2d=ker_par.is2d;  op=ker_par.op; //! op:1:ssf 2:fd 3:ffd
    v0_z=ker_par.v0_z; v_zw=ker_par.v_zw; v_zw0=ker_par.v_zw0;
    tmp_fft=ker_par.tmp_fft; add_wld=ker_par.add_wld; tmp_wld=ker_par.tmp_wld;

    maxvz=maxval(v_z,nx*ny); minvz=minval(v_z,nx*ny);  dref=(maxvz-minvz)/(nref+1);

    if ( fabs(maxvz-minvz)<50.0 ) { ///?????????????wait
	nref=1;  v0_z[0]=0.5*(maxvz+minvz);
    }
    else{
	v0_z[0]=minvz;  v0_z[nref-1]=maxvz;  dref=(v0_z[nref-1]-v0_z[0])/(float)(nref-1); v0_z[1]=(minvz+maxvz)*0.5;  
    }

    for (iy=0;iy<ny;iy++)
	for(ix=0;ix<nx;ix++)
	    v_zw[i2(iy,ix,nx)]=v_z[i2(iy,ix,nx)]/w;
    if (op == 2){
	w_v0_z=0.0;
	shot_ker_ssf(wld_z,w_v0_z,w,v_z,ker_par,signn);
	shot_ker_fd(wld_z,v_zw,v_zw,ker_par,signn);
    }
    else{
	if (is2d){
	    rowc(wld_z,tmp_fft,0,nx,1);  //tmp_fft=wld_z(1,:);
	    kissfft(tmp_fft,nx,+1);
	    rowc(wld_z,tmp_fft,0,nx,0);  //wld_z(1,:)=tmp_fft(:);
	}
	else{
	    for(iy=0;iy<ny;iy++){
		rowc(wld_z,tmp_fft,iy,nx,1);  //tmp_fft=wld_z(1,:);
		kissfft(tmp_fft,nx,+1);
		rowc(wld_z,tmp_fft,iy,nx,0);  //wld_z(1,:)=tmp_fft(:);
	    }
	    for(ix=0;ix<nx;ix++){
		colc(wld_z,tmp_fft,ix,nx,ny,1);
		kissfft(tmp_fft,ny,+1);
		colc(wld_z,tmp_fft,ix,nx,ny,0);   // wait tmp_fftx tmp_ffty
	    }
  
	}

	for( iy=0;iy<ny;iy++)
	    for(ix=0;ix<nx;ix++)
		add_wld[i2(iy,ix,nx)]=sf_cmplx(0.0,0.0);   //add_wld=0.0;

	for( iref=0;iref<nref;iref++){
	    matrix_equ_c(tmp_wld,wld_z,ny,nx);//tmp_wld=wld_z;
	    w_v0_z=w/v0_z[iref];
	    shot_phase_shift_ani(tmp_wld,w_v0_z,ep,dl,phi,f,ker_par,signn);
	    if (is2d){
		rowc(tmp_wld,tmp_fft,0,nx,1); //tmp_fft=tmp_wld(:,1);
		kissfft(tmp_fft,nx,-1);
		rowc(tmp_wld,tmp_fft,0,nx,0);  //tmp_wld(:,1)=tmp_fft(:);
	    }
	    else{
		for(ix=0;ix<nx;ix++){
		    colc(tmp_wld,tmp_fft,ix,nx,ny,1);
		    kissfft(tmp_fft,ny,-1);
		    colc(tmp_wld,tmp_fft,ix,nx,ny,0);   // wait tmp_fftx tmp_ffty
		}
		for(iy=0;iy<ny;iy++){
		    rowc(tmp_wld,tmp_fft,iy,nx,1);  //tmp_fft=wld_z(1,:);
		    kissfft(tmp_fft,nx,-1);
		    rowc(tmp_wld,tmp_fft,iy,nx,0);  //wld_z(1,:)=tmp_fft(:);
		}

	    }
	    //shot_ker_ssf(tmp_wld,w_v0_z,w,v_z,ker_par,signn);
	    if (op==3){ 
		for( iy=0;iy<ny;iy++)
		    for(ix=0;ix<nx;ix++)
			v_zw0[i2(iy,ix,nx)]=v_zw[i2(iy,ix,nx)]-v0_z[iref]/w;
		shot_ker_fd(tmp_wld,v_zw0,v_zw,ker_par,signn);
	    }
	    if ( nref==1) vector_cp_c(add_wld,tmp_wld,ny*nx); //add_wld=tmp_wld;
	    else{
		for(iy=0;iy<ny;iy++){      
		    for(ix=0;ix<nx;ix++){  
			ivx=(int)((v_z[i2(iy,ix,nx)]-v0_z[0])/dref); //wait ???
			wvx=(v_z[i2(iy,ix,nx)]-(v0_z[0]+(float)(ivx)*dref))/dref;   //wait ???
			if (v_z[i2(iy,ix,nx)] >= v0_z[nref-1]){
			    if (iref==nref-1) add_wld[i2(iy,ix,nx)]+=tmp_wld[i2(iy,ix,nx)];
			}
			else{
			    if (v_z[i2(iy,ix,nx)] <= v0_z[0]){ 
				if (iref==0)  add_wld[i2(iy,ix,nx)]+=tmp_wld[i2(iy,ix,nx)];
			    }
			    else{
				if (ivx==iref) add_wld[i2(iy,ix,nx)] += tmp_wld[i2(iy,ix,nx)]*(1-wvx);
				if (ivx+1==iref) add_wld[i2(iy,ix,nx)] += tmp_wld[i2(iy,ix,nx)]*wvx;
			    }
			}
		    }
		}
	    }
    
	}

	vector_cp_c(wld_z,add_wld,ny*nx);  //wld_z=add_wld;
    }

}




void shot_ker_depth_onestep_vti_phaseshift(sf_complex *wld_z,float w_v0_z,float ep,float dl,int signn, struct shot_ker_par_type ker_par){
    int nx;
    sf_complex  *tmp_fft;
    nx=ker_par.nx;  //! op:1:ssf 2:fd 3:ffd
    tmp_fft=ker_par.tmp_fft; 

    rowc(wld_z,tmp_fft,0,nx,1);  //tmp_fft=wld_z(1,:);
    kissfft(tmp_fft,nx,+1);
    rowc(wld_z,tmp_fft,0,nx,0);  //wld_z(1,:)=tmp_fft(:);

    shot_phase_shift(wld_z,w_v0_z,ker_par,signn);

    rowc(wld_z,tmp_fft,0,nx,1); //tmp_fft=tmp_wld(:,1);
    kissfft(tmp_fft,nx,-1);
    rowc(wld_z,tmp_fft,0,nx,0);  //tmp_wld(:,1)=tmp_fft(:);

}

void shot_ker_depth_onestep_impulse(sf_complex *wld_z,float *v_z,float w,int signn, float ep,float dl,
				    struct shot_ker_par_type ker_par)
/*< onestep impulse >*/
{
    int op,nref,nx,ny,is2d;
    sf_complex  *tmp_fft;
    int ix,iy;
    float w_v0_z,maxvz,minvz;
    float *v0_z,*v_zw;
    nref=ker_par.nref;  nx=ker_par.nx;  ny=ker_par.ny;  is2d=ker_par.is2d;  op=ker_par.op; //! op:1:ssf 2:fd 3:ffd
    v0_z=ker_par.v0_z; v_zw=ker_par.v_zw;
    tmp_fft=sf_complexalloc(nx*2);
    maxvz=maxval(v_z,nx*ny); minvz=minval(v_z,nx*ny); 
    if ( fabs(maxvz-minvz)<50.0 ) { ///?????????????wait
	nref=1;  v0_z[0]=0.5*(maxvz+minvz);
    }
    else{
	v0_z[0]=minvz;  v0_z[nref-1]=maxvz;  v0_z[1]=(minvz+maxvz)*0.5;  
    }

    for (iy=0;iy<ny;iy++)
	for(ix=0;ix<nx;ix++)
	    v_zw[i2(iy,ix,nx)]=v_z[i2(iy,ix,nx)]/w;

    if (op == 2){
	w_v0_z=0.0;
	shot_ker_ssf(wld_z,w_v0_z,w,v_z,ker_par,signn);
	shot_ker_fd(wld_z,v_zw,v_zw,ker_par,signn);
    }

    else{
	if (is2d){
	    vector_value_c(tmp_fft,sf_cmplx(0.0,0.0),nx);
	    rowc(wld_z,tmp_fft,0,nx,1);  //tmp_fft=wld_z(1,:);
	    kissfft(tmp_fft,nx,+1);
	    rowc(wld_z,tmp_fft,0,nx,0);  //wld_z(1,:)=tmp_fft(:);
	}
	w_v0_z=w/v0_z[0];
	shot_phase_shift_vti(wld_z,w_v0_z,ep,dl,ker_par,signn);
	if (is2d){
	    rowc(wld_z,tmp_fft,0,nx,1); //tmp_fft=tmp_wld(:,1);
	    kissfft(tmp_fft,nx,-1);
	    rowc(wld_z,tmp_fft,0,nx,0);  //tmp_wld(:,1)=tmp_fft(:);
	}
    }
    free(tmp_fft);
}




void shot_ker_depth_onestep_fd2(sf_complex *wld_z,float *v_z,float w,int signn, 
				float a1,float b1,float a2,float b2,struct shot_ker_par_type ker_par)
/*< onestep fd2 >*/
{
    int op,nref,nx,ny,is2d;

    sf_complex  *tmp_fft;
    int ix,iy;
    float w_v0_z,maxvz,minvz;
    float *v0_z,*v_zw;
    nref=ker_par.nref;  nx=ker_par.nx;  ny=ker_par.ny;  is2d=ker_par.is2d;  op=ker_par.op; //! op:1:ssf 2:fd 3:ffd
    v0_z=ker_par.v0_z; v_zw=ker_par.v_zw;
    tmp_fft=sf_complexalloc(nx*2);

    maxvz=maxval(v_z,nx*ny); minvz=minval(v_z,nx*ny);  
    if ( fabs(maxvz-minvz)<50.0 ) { ///?????????????wait
	nref=1;  v0_z[0]=0.5*(maxvz+minvz);
    }
    else{
	v0_z[0]=minvz;  v0_z[nref-1]=maxvz;  v0_z[1]=(minvz+maxvz)*0.5;  
    }
    for (iy=0;iy<ny;iy++)
	for(ix=0;ix<nx;ix++)
	    v_zw[i2(iy,ix,nx)]=v_z[i2(iy,ix,nx)]/w;

    if (op == 2){
	if (is2d){
	    vector_value_c(tmp_fft,sf_cmplx(0.0,0.0),nx);
	    rowc(wld_z,tmp_fft,0,nx,1);  //tmp_fft=wld_z(1,:);
	    kissfft(tmp_fft,nx,+1);
	    rowc(wld_z,tmp_fft,0,nx,0);  //wld_z(1,:)=tmp_fft(:);
	}
	w_v0_z=w/(0.75*minvz+0.25*maxvz);
	li_filter_fd2(wld_z,w_v0_z,a1,b1,a2,b2,ker_par,signn); 
	if (is2d){
	    rowc(wld_z,tmp_fft,0,nx,1); //tmp_fft=tmp_wld(:,1);
	    kissfft(tmp_fft,nx,-1);
	    rowc(wld_z,tmp_fft,0,nx,0);  //tmp_wld(:,1)=tmp_fft(:);
	}
	w_v0_z=0.0;
	shot_ker_ssf(wld_z,w_v0_z,w,v_z,ker_par,signn);
   
	ker_par.ax=a1*ker_par.dz/(2.0*ker_par.dx*ker_par.dx); ker_par.ay=a1*ker_par.dz/(2.0*ker_par.dy*ker_par.dy);
	ker_par.bx=b1/(ker_par.dx*ker_par.dx);        ker_par.by=b1/(ker_par.dy*ker_par.dy);
	shot_ker_fd(wld_z,v_zw,v_zw,ker_par,signn);
	ker_par.ax=a2*ker_par.dz/(2.0*ker_par.dx*ker_par.dx); ker_par.ay=a2*ker_par.dz/(2.0*ker_par.dy*ker_par.dy);
	ker_par.bx=b2/(ker_par.dx*ker_par.dx);        ker_par.by=b2/(ker_par.dy*ker_par.dy);
	shot_ker_fd(wld_z,v_zw,v_zw,ker_par,signn);
   
   
    }
    free(tmp_fft);
}


void li_filter_fd2(sf_complex *wld_z,float w_v0_z,float a1,float b1,float a2,float b2,
		   struct shot_ker_par_type ker_par,int signn)
/*< Li filter >*/
{
    int nx;
    float fkx,dkx,dz;
    int ix,iy;
    float kx,kz2,phsft,kx2;
    nx=ker_par.nx; 
    fkx=ker_par.fkx; 
    dkx=ker_par.dkx; 
    dz=ker_par.dz;
    for(iy=0;iy<1;iy++){
	for(ix=0;ix<nx;ix++){  //do ix=1,nx
	    kx=fkx+ix*dkx;
	    kx2=(kx*kx)/(w_v0_z*w_v0_z);
	    if (kx2 <=1) {
		if (kx2> sinf(SF_PI*80/180)*sinf(SF_PI*80/180)){
		    kz2=1.0-kx2;
		    phsft=-signn*dz*w_v0_z*(sqrtf(kz2)-  (1-(a1*kx2)/(1-b1*kx2)-(a2*kx2)/(1-b2*kx2) )   );
		    wld_z[i2(iy,ix,nx)]=wld_z[i2(iy,ix,nx)]*sf_cmplx(cosf(phsft),sinf(phsft));
		}
	    }
	    else{
		if (kx2>2) kx2=2.0;
		kz2=fabsf(kx2- 1);
		phsft=dz*w_v0_z*(sqrtf(kz2));
		wld_z[i2(iy,ix,nx)]=0.0;
	    }


	}
    }

}

void shot_ani_convlv(sf_complex *wld_z,float *v_z,float *ep_z,float *dl_z,float w,struct shot_ker_par_type ker_par,struct shot_ker_ani_type ani_par,int signn)
{
    float w_vp,ep,dl;
    int nx,ny;
    sf_complex  *tmpwld_z; //(:,:)
    sf_complex *conapp,*conapm;
    int  iy,ix,ic;
    nx=ker_par.nx; ny=ker_par.ny;
    tmpwld_z=ker_par.tmp_wld;   //allocate(tmpwld_z(nx,ny)) ker_par.tmp_wld is allocated in ker_init in kernal.c
    conapp=ani_par.conapp; conapm=ani_par.conapm; 
    // both  ani_par.conapp ani_par.conapm allocated in ker_ani_init() released in ker_ani_release

    vector_cp_c(tmpwld_z,wld_z,nx*ny); //tmpwld_z=wld_z;
    vector_value_c(wld_z,sf_cmplx(0.0,0.0),nx*ny);
    for( iy=0;iy<ny;iy++){
	for(ix=0;ix<nx;ix++){
	    w_vp=w/v_z[i2(iy,ix,nx)]; ep=ep_z[i2(iy,ix,nx)]; dl=dl_z[i2(iy,ix,nx)];
	    if (ep == 0.0 && dl==0.0) wld_z[i2(iy,ix,nx)]=tmpwld_z[i2(iy,ix,nx)];
	    else{
		//printf("%d,%f,%f,%f\n",ix,w_vp,ep,dl);
		get_convlv_coe(w_vp,ep,dl,&ani_par);
         
		if (signn==+1) {
		    for(ic=0;ic<ani_par.n;ic++){
			conapp[ic]=conjf(conapp[ic]); conapm[ic]=conjf(conapm[ic]);
		    } //for(ic)
		}
		//for (ic=0;ic<ani_par.n;ic++){
		//  printf("conapp[%d]=(%f,%f),conapm[%d]=(%f,%f)\n",ic,__real__ conapp[ic],__imag__ conapp[ic],ic,__real__ conapm[ic],__imag__ conapm[ic]);
		//}       
 
		wld_z[i2(iy,ix,nx)]=conapp[0]*tmpwld_z[i2(iy,ix,nx)];
		for(ic=1;ic<=ani_par.n-1;ic++){
		    if (ix+ic <nx ) wld_z[i2(iy,ix,nx)]+=conapp[ic]*tmpwld_z[i2(iy,ix+ic,nx)];
		    if (ix-ic >= 0 )wld_z[i2(iy,ix,nx)]+=conapm[ic]*tmpwld_z[i2(iy,ix-ic,nx)]; 
		}// for(ic)
        

	    }// if
	}// for(ix) 
    }// for(iy)

}

void shot_ani_convlv3dtest(sf_complex *wld_z,sf_complex *conapp,sf_complex *conapm,int tnx,int tny,
			   struct shot_ker_par_type ker_par, int signn)
{
    int nx,ny,tn;
    sf_complex  *tmpwld_z; //(:,:)
    sf_complex *conap,*conam;
    int  iy,ix,ic,icx,icy;
    tn=tnx*tny;
    nx=ker_par.nx; ny=ker_par.ny;
    tmpwld_z=ker_par.tmp_wld;   //allocate(tmpwld_z(nx,ny)) ker_par.tmp_wld is allocated in ker_init in kernal.c
    conap=sf_complexalloc(tn); conam=sf_complexalloc(tn);
    vector_cp_c(tmpwld_z,wld_z,nx*ny); //tmpwld_z=wld_z;
    vector_cp_c(conap,conapp,tn); vector_cp_c(conam,conapm,tn);
    vector_value_c(wld_z,sf_cmplx(0.0,0.0),nx*ny);
    if (signn==+1) {
	for(ic=0;ic<tn;ic++){
	    conap[ic]=conjf(conapp[ic]); conam[ic]=conjf(conapm[ic]);
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

void shot_phase_shift_ani(sf_complex *wld_z,float w_v0_z,float ep,float dl,float phi,float f,
			  struct shot_ker_par_type ker_par,int signn)
/*< phase shift anisotropic >*/
{
    int nx,ny,weitmx1,weitmx2,weitmy1,weitmy2;
    float fkx,fky,dkx,dky,dz;
  
    int ix,iy;
    float kx,ky,phsft;
    sf_complex kz1;
    sf_complex *kza;
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

	    if ( fabsf(cimagf( kza[i2(iy,ix,nx)]  ))<0.00000000000001){
		if ( fabsf(cimagf( kza[i2(iy,ix,nx)]  )) <0.00000000001*fabsf(crealf( kza[i2(iy,ix,nx)]  )) )
		    kza[i2(iy,ix,nx)]=sf_cmplx(fabsf(crealf(kza[i2(iy,ix,nx)])),0.0);    
		else
		    kza[i2(iy,ix,nx)]=sf_cmplx(0.0,-fabsf(cimagf(kza[i2(iy,ix,nx)])));
	    }
	    else
		kza[i2(iy,ix,nx)]=sf_cmplx(0.0,-fabsf(cimagf(kza[i2(iy,ix,nx)])));
	}
    }  

    for (iy=0;iy<ny;iy++){
	for(ix=0;ix<nx;ix++){
	    if (fabsf(cimagf(kza[i2(iy,ix,nx)]))<0.000001){
		if ( fabsf(cimagf(kza[i2(iy,ix+1,nx)]))>0.000001 && fabsf(cimagf(kza[i2(iy,ix-1,nx)]))>0.000001   )
		    kza[i2(iy,ix,nx)]=sf_cmplx( 0.0,1.0);
		if ( fabsf(cimagf(kza[i2(iy+1,ix,nx)]))>0.000001 && fabsf(cimagf(kza[i2(iy-1,ix,nx)]))>0.000001   )
		    kza[i2(iy,ix,nx)]=sf_cmplx( 0.0,1.0);

	    }
	}
    }

    iy=ny/2;
    for(ix=nx/2;ix<nx;ix++)
	if (fabsf(cimagf(kza[i2(iy,ix,nx)]))>0.001 && fabsf(crealf(kza[i2(iy,ix,nx)]))<0.00000001) break;
    weitmx2=ix;
    for(ix=0;ix<nx/2;ix++)
	if (fabsf(cimagf(kza[i2(iy,ix,nx)]))>0.001 && fabsf(crealf(kza[i2(iy,ix,nx)]))<0.00000001) break;
    weitmx1=ix; 
    ix=nx/2;
    for(iy=ny/2;iy<ny;iy++)
	if (fabsf(cimagf(kza[i2(iy,ix,nx)]))>0.001 && fabsf(crealf(kza[i2(iy,ix,nx)]))<0.00000001) break;
    weitmy2=iy;
    for(iy=0;iy<ny/2;iy++)
	if (fabsf(cimagf(kza[i2(iy,ix,nx)]))>0.001 && fabsf(crealf(kza[i2(iy,ix,nx)]))<0.00000001) break;
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
	    if ( fabsf(cimagf(kz1)) <0.0000001*fabsf(crealf(kz1)) ){
		phsft=-signn*dz*fabsf(crealf(kz1));
		wld_z[i2(iy,ix,nx)]=wld_z[i2(iy,ix,nx)]*sf_cmplx(cos(phsft),sin(phsft));
	    }
	    else{
		wld_z[i2(iy,ix,nx)]=sf_cmplx(0.0,0.0);
		//phsft=dz*fabsf(cimagf(kz1));  
		//wld_z[i2(iy,ix,nx)]=wld_z[i2(iy,ix,nx)]*exp(-phsft);
	    }
      
	} 
    }


    free(kza);
}

void shot_phase_shift_vti(sf_complex *wld_z,float w_v0_z,float ep,float dl,
			  struct shot_ker_par_type ker_par,int signn)
/*< phase shift VTI >*/
{
    int nx;
    float fkx,dkx,dz; 
    int ix,iy;
    float kx,kz2,phsft;

    nx=ker_par.nx; 
    fkx=ker_par.fkx; 
    dkx=ker_par.dkx; 
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





