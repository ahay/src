/* Impulse response for plane-wave migration in tilted coordinates */
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

/*
  Guojian Shan
  email: shan@sep.stanford.edu
  This program runs plane migration in VTI media in tilted coordinates.
  Input File: Data=data.H Velocity=velocity.H Wavelet=wavelet.H Eps=ep.H(\epsilon) Dlt=dl.H(\delta)
  Input parameters: nhx=  mig_nz=  mig_dz= image_min_x=  image_max_x=
  Output File: Image_hx=image_h.H (horizontal ODCIG) Image_hz=image_v.H (vertical ODCIG)
  ---------------------------------------------------------------------------------------------
  run: anitiltplanemig.x Data=celldata.H Velocity=cellvelocity.H Wavelet=wavelet_exxon1.H 
  Eps=cellep.H Dlt=celldl.H Extable=table.H
  Image_hx=imagehx.H Image_hz=hz.H 
  nhx=1 mig_nz=1000 mig_dz=4 image_min_x=-1000 image_max_x=9000.0 
  ----------------------------------------------------------------------------------------------
*/
#include <rsf.h>

#include "data.h"
#include "image.h"
#include "kernel.h"
#include "rotate.h"
#include "velocity.h"
#include "arrayindex.h"

static void mig_grab_par(int *op,int *mig_nz,float *mig_dz,float *fd_a,float *fd_b,float *trick);
static void veltranspose(float *vel,int ny,int nx,int nz);

int main(int argc, char **argv)
{
    struct shot_data_par_type   data_par;
    struct shot_plane_par_type   plane_par;
    struct shot_image_par_type  image_par;
    struct shot_ker_par_type    ker_par;
    struct shot_rotate_par_type rotate_par;
    sf_file velocity,wave;
    float mig_dx,mig_dy,mig_dz,mig_min_x,mig_min_y,fd_a,fd_b,trick,w;
    int op,mig_nz,mig_nx,mig_ny,tmp_mig_nx;
    int mig_nx_rotate,mig_nz_rotate,mig_ny_rotate;
    float rotate_angle;
    float epvalue,dlvalue;

    int i_sy,i_sx,iw,iz;
    float *vel,*ep,*dl,*image_z,*image,*image_rotate,*v_z;
    sf_complex *wld_s,*wld_r,*wld_s_z0,*wld_r_z0,*wld_s_tilt,*wld_r_tilt;
    float *image_1hx,**image_rotate_allhx,**image_rotate_allhz;
    int ixz,ihx;
    int data_nsy, ix;
    int mig_size,image_size;
    float dkx; int bx,ex,nx_iz;
    float ttmprotate=1.0;
    float t;
    bool vti;
    float a1,b1,a2,b2;

    sf_init(argc,argv);
    if (!sf_getfloat("t",&t))    t=0.8;
    if (!sf_getfloat("ep",&epvalue)) epvalue=0.4;
    if (!sf_getfloat("dl",&dlvalue)) dlvalue=0.2;
    if (!sf_getbool("vti",&vti)) vti=true;
    a1=0.040315157; b1=0.873981642;
    a2=0.457289566; b2=0.222691983;

    a1=0.659558; b1=0.790347;
    a2=0.030433; b2=1.692681;

    velocity_init(&velocity);
    //velocity_ani_init(&eps, &dlt);
    data_par=shot_data_init();
    sf_warning("nsy=%d  nsx=%d",data_par.nsy,data_par.nsx); 
    shot_data_wavelet(data_par);
    image_par=shot_image_grab_par();

    mig_dx=data_par.dhx; mig_dy=data_par.dhy;
    mig_grab_par(&op,&mig_nz,&mig_dz,&fd_a,&fd_b,&trick);
    shot_image_init(mig_dx,mig_dy,mig_dz,mig_nz,&image_par); 
  
    wave = sf_output("wave");
    sf_putfloat(wave,"d1",mig_dx);
    sf_putfloat(wave,"d2",mig_dy);
    sf_putfloat(wave,"d3",mig_dz);
  
    data_nsy=data_par.nsy; 
    /* data_nsx=data_par.nsx; */

    for(i_sy=0;i_sy< data_nsy;i_sy++){
	for(i_sx=0;i_sx< 1;i_sx++){
	    shot_data_plane_grab_mig_par(i_sx,i_sy,&mig_min_x,&mig_min_y,&tmp_mig_nx,&mig_ny,data_par,&plane_par );
	    mig_nx=tmp_mig_nx;
	    sf_warning("nx---------------%d",mig_nx); 
	    rotate_angle=asin( (i_sx*data_par.dsx+data_par.osx)*1500.0);
      
	    rotate_angle=-50.0/180.0*SF_PI;
	    rotate_angle=0.0/180.0*SF_PI;

	    rotate_par=shot_rotate_init(rotate_angle,tmp_mig_nx,mig_nz,mig_dx,mig_dz);
	    sf_warning("nx=%d,nz=%d,ny=%d",rotate_par.nx_rotate,rotate_par.nz_rotate,rotate_par.ny_rotate);
	    mig_nx_rotate=rotate_par.nx_rotate; mig_ny_rotate=rotate_par.ny_rotate;
	    mig_nz_rotate=rotate_par.nz_rotate;
      
	    mig_size=mig_nx_rotate*mig_ny_rotate; image_size=mig_size*image_par.nhy*image_par.nhx*mig_nz_rotate;

	    ker_par=shot_ker_p_init(op,mig_nx_rotate,mig_ny_rotate,mig_dx,mig_dy,mig_dz,fd_a,fd_b,trick,1,ttmprotate);
	    vel=sf_floatalloc(mig_size*mig_nz_rotate); 
	    wld_s=sf_complexalloc(SF_MAX(mig_size,mig_nx*mig_ny_rotate)); 
	    wld_r=sf_complexalloc(SF_MAX(mig_size,mig_nx*mig_ny_rotate));
	    ep=sf_floatalloc(mig_size*mig_nz_rotate); dl=sf_floatalloc(mig_size*mig_nz_rotate);
	    sf_warning("nz_s_rotate=%d",rotate_par.nz_s_rotate);
	    //if (0.0!=rotate_angle){
	    wld_s_z0=sf_complexalloc(mig_size*rotate_par.nz_s_rotate); wld_r_z0=sf_complexalloc(mig_size*rotate_par.nz_s_rotate);
	    //}
	    image=sf_floatalloc(image_size); vector_value_f(image,0.0,image_size);
      
	    sf_putint(wave,"n2",mig_nx_rotate);
	    sf_putfloat(wave,"o2",0.0);
	    sf_putint(wave,"n1",mig_nz_rotate);

	    if (rotate_par.rotate_angle>=0.0) 
		sf_putfloat(wave,"o2",mig_min_x); 
	    else 
		sf_putfloat(wave,"o2",mig_min_x+(mig_nx-1)*mig_dx);

	    velocity_read(velocity,vel,wave,rotate_par.rotate_angle,rotate_par.ix_rotate_center); 
	    veltranspose(vel,mig_ny_rotate,mig_nx_rotate,mig_nz_rotate);
      
	    for(iw=0;iw<data_par.nw; iw++){ 
		w=iw*data_par.dw+data_par.ow;
		sf_warning("iw=%d  w=%f;",iw,w);
           
		shot_plane_source(mig_min_x,mig_min_y,mig_nx,mig_ny,iw,i_sx,i_sy,wld_s,data_par,plane_par);
		shot_plane_receiver(mig_min_x,mig_min_y,mig_nx,mig_ny,iw,i_sx,i_sy,wld_r,data_par,plane_par);
		shot_rotate_data_rotate(wld_r,wld_r_z0,&rotate_par);
		shot_rotate_data_rotate(wld_s,wld_s_z0,&rotate_par);
		vector_value_c(wld_s,sf_cmplx(0.0,0.0),mig_size); 
		vector_value_c(wld_r,sf_cmplx(0.0,0.0),mig_size); //wld_r=0.0  wld_s=0.0
        
		for(iz=0;iz<mig_nz_rotate;iz++){
		    image_z=image+mig_size*image_par.nhy*image_par.nhx*iz;
		    if (iz <rotate_par.nz_s_rotate){
			for(ix=0;ix<mig_nx_rotate;ix++) wld_s[ix]+=wld_s_z0[i2(iz,ix,mig_nx_rotate)];  //wld_s+=wld_s_z0(iz,:,:);
			for(ix=0;ix<mig_nx_rotate;ix++) wld_r[ix]+=wld_r_z0[i2(iz,ix,mig_nx_rotate)];  //wld_r+=wld_r_z0(iz,:,:);
		    }
		    if ( 0.0!=rotate_angle){
			shot_rotate_get_bx_ex_iz(&dkx,&bx,&ex,iz,20,rotate_par);
			nx_iz=ex-bx+1;  ker_par.dkx=dkx; ker_par.nx=nx_iz;
			v_z=vel+iz*mig_size+bx; wld_s_tilt=wld_s+bx; wld_r_tilt=wld_r+bx; image_z+=image_par.nhy*image_par.nhx*bx;
			/* ep_z=ep+iz*mig_size+bx; dl_z=dl+iz*mig_size+bx; */
		    }
		    else{
			v_z=vel+iz*mig_size; wld_s_tilt=wld_s; wld_r_tilt=wld_r; /* ep_z=ep+iz*mig_size; dl_z=dl+iz*mig_size; */
			nx_iz=mig_nx_rotate;
		    }
		    if (vti)
			shot_ker_depth_onestep_impulse(wld_s_tilt,v_z,w,-1,epvalue,dlvalue,ker_par);
		    else
			shot_ker_depth_onestep_fd2(wld_s_tilt,v_z,w,-1,a1,b1,a2,b2,ker_par); 
		    //shot_ker_depth_onestep(wld_s_tilt,v_z,w,-1,ker_par);
		    //shot_ker_depth_onestep(wld_r_tilt,v_z,w,+1,ker_par); 
		    vector_value_c(wld_r,sf_cmplx(cosf(w*t),sinf(w*t)),mig_nx_rotate);
		    //printf("max=%f\n",maxvalc(wld_r,mig_nx_rotate));
		    shot_image_cross_tilt(wld_s_tilt,wld_r_tilt,image_z,mig_nx_rotate,mig_ny_rotate,image_par,nx_iz);
          
		}
	    }
	    sf_warning(".");
	    shot_ker_release(&ker_par);
	    free(wld_s); free(wld_r);free(vel); free(ep); free(dl); 
	    if (0.0!=rotate_angle){
		free(wld_s_z0); free(wld_r_z0); 
	    }

	    image_rotate=sf_floatalloc(mig_nx*mig_nz);  
	    image_1hx=sf_floatalloc(mig_size*mig_nz_rotate);
	    image_rotate_allhx=sf_floatalloc2(image_par.nhx,mig_nx*mig_nz);
     
	    for (ihx=0;ihx<image_par.nhx;ihx++){
		//printf("ihx=%d\n",ihx);
		for(ixz=0;ixz<mig_size*mig_nz_rotate;ixz++) image_1hx[ixz]=image[i2(ixz,ihx,image_par.nhx)];
		shot_rotate_image_rotate(image_1hx,image_rotate,rotate_par); 
		for(ixz=0;ixz<mig_nx*mig_nz;ixz++) image_rotate_allhx[ixz][ihx]=image_rotate[ixz]*1.0;
	    }
      
       
 
	    free(image); //rotate image 
	    //stat=auxputch ("rotate_angle","f",&rotate_angle,"Image_hx");

	    if (rotate_angle !=0.0){
		image_rotate_allhz=sf_floatalloc2(mig_nx*mig_nz, image_par.nhx);
		shot_image_decompose(*image_rotate_allhx, *image_rotate_allhz,mig_nx,mig_nz,rotate_angle,image_par);
		sf_warning("ishot=%d,maxhx=%f,maxhz=%f",
			   i_sx,
			   maxval(*image_rotate_allhx,mig_nx*mig_nz*image_par.nhx),
			   maxval(*image_rotate_allhz,mig_nx*mig_nz*image_par.nhx));
        
		shot_image_stack(*image_rotate_allhz,mig_min_x,mig_min_y,mig_nx,mig_ny,image_par.imagez,image_par);
		free(*image_rotate_allhz);
		free(image_rotate_allhz);
	    }


	    shot_image_stack(*image_rotate_allhx,mig_min_x,mig_min_y,mig_nx,mig_ny,image_par.imagex,image_par); 
	    free(image_rotate); // rite image
      
	    free(image_1hx); 
	    free(*image_rotate_allhx);
	    free(image_rotate_allhx);
	    shot_rotate_release(&rotate_par);
	}
    }
    shot_data_release(&data_par);
    shot_image_release(&image_par);
     
    exit(0);
}

static void mig_grab_par(int *op,int *mig_nz,float *mig_dz,float *fd_a,float *fd_b,float *trick)
{
    if (!sf_getint  ("operator",op))    *op=1;
    if (!sf_getfloat("mig_dz",mig_dz))  *mig_dz=4.0;
    if (!sf_getint  ("mig_nz",mig_nz))  *mig_nz=1000;
    if (!sf_getfloat("fd_a",fd_a))      *fd_a=0.5;
    if (!sf_getfloat("fd_b",fd_b))      *fd_b=0.0;
    if (!sf_getfloat("fd_trick",trick)) *trick=1.0/6.0;
}

static void veltranspose(float *vel,int ny,int nx,int nz)
{
    float *tmpvel=sf_floatalloc(ny*nx*nz);
    int nvel[2],ntmpvel[2];
    int i,iz,iy,ix;
    d3(ny,nx,nvel); d3(nx,nz,ntmpvel);
    for(i=0;i<ny*nx*nz;i++)
	tmpvel[i]=vel[i];
    for(iz=0;iz<nz;iz++)
	for(iy=0;iy<ny;iy++)
	    for(ix=0;ix<nx;ix++)
		vel[i3(iz,iy,ix,nvel)]=tmpvel[i3(iy,ix,iz,ntmpvel)];
    free(tmpvel); 
}
