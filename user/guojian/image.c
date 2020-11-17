/* Image type routines */

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
#include <rsf_su.h>

#include "image.h"
#include "arrayindex.h"

#ifndef _image_h

struct shot_image_par_type
{
    float min_x,max_x,min_y,max_y;
    float dx,dy,dz;
    float ohx,ohy;
    int nx,ny,nz,nhx,nhy;
    float *sincx_table; 
    sf_file imagex,imagez;
};
/*^*/

#endif

struct shot_image_par_type shot_image_grab_par(void)
/*< initialize >*/
{
    struct shot_image_par_type image_par;
    float min_x,max_x,min_y,max_y;
    int nhx,nhy,ihx;
    if (!sf_getint("nhx",&nhx) ) nhx=1;
    if (!sf_getint("nhy",&nhy) ) nhy=1;
    if (!sf_getfloat("image_min_x",&min_x)) min_x=1000.0;
    if (!sf_getfloat("image_max_x",&max_x)) max_x=10000.0;
    if (!sf_getfloat("image_min_y",&min_y)) min_y=0.0;
    if (!sf_getfloat("image_max_y",&max_y)) max_y=0.0;
    image_par.nhx=nhx; image_par.nhy=nhy;
    image_par.min_x=min_x; image_par.max_x=max_x;
    image_par.min_y=min_y; image_par.max_y=max_y;
    image_par.sincx_table=sf_floatalloc(8*2001);
    for(ihx=0;ihx<2001;ihx++) mksinc((float)(ihx)*0.0005,8,image_par.sincx_table+ihx*8); 
    sf_warning("Image: minx=%f  maxx=%f\n",image_par.min_x,image_par.max_x);
    return image_par; 
}

void shot_image_init(float dx,float dy,float dz,int nz,struct shot_image_par_type *image_par)
/*< initialize image >*/
{
    sf_file imagex, imagez;
    float *tmpimg;
    float ohx,ohy;
    int i; 

    tmpimg=sf_floatalloc(nz);
 
    vector_value_f(tmpimg,0.0,nz);
    image_par->dx=dx; image_par->dy=dy; image_par->dz=dz;
    image_par->nx=(image_par->max_x-image_par->min_x)/(image_par->dx)+1;
    image_par->ny=(image_par->max_y-image_par->min_y)/(image_par->dy)+1;
    image_par->nz=nz;
    if (!sf_getfloat("imageohx",&ohx)) ohx=-image_par->dx*image_par->nhx/2;
    if (!sf_getfloat("imageohy",&ohy)) ohy=-image_par->dy*image_par->nhy/2;
    image_par->ohx=ohx;   image_par->ohy=ohy;
    
    imagex = sf_output("Image_hx");
    imagez = sf_output("Image_hz");

    sf_putfloat(imagex,"o1",0.0);              sf_putfloat(imagex,"d1",image_par->dz); sf_putint(imagex,"n1",image_par->nz);
    sf_putfloat(imagex,"o2",image_par->ohx);   sf_putfloat(imagex,"d2",image_par->dx); sf_putint(imagex,"n2",image_par->nhx);
    sf_putfloat(imagex,"o3",image_par->ohy);   sf_putfloat(imagex,"d3",image_par->dy); sf_putint(imagex,"n3",image_par->nhy);
    sf_putfloat(imagex,"o4",image_par->min_x); sf_putfloat(imagex,"d4",image_par->dx); sf_putint(imagex,"n4",image_par->nx);
    sf_putfloat(imagex,"o5",image_par->min_y); sf_putfloat(imagex,"d5",image_par->dy); sf_putint(imagex,"n5",image_par->ny);

    sf_putfloat(imagez,"o1",0.0);              sf_putfloat(imagez,"d1",image_par->dz); sf_putint(imagez,"n1",image_par->nz);
    sf_putfloat(imagez,"o2",image_par->ohx);   sf_putfloat(imagez,"d2",image_par->dx); sf_putint(imagez,"n2",image_par->nhx);
    sf_putfloat(imagez,"o3",image_par->ohy);   sf_putfloat(imagez,"d3",image_par->dy); sf_putint(imagez,"n3",image_par->nhy);
    sf_putfloat(imagez,"o4",image_par->min_x); sf_putfloat(imagez,"d4",image_par->dx); sf_putint(imagez,"n4",image_par->nx);
    sf_putfloat(imagez,"o5",image_par->min_y); sf_putfloat(imagez,"d5",image_par->dy); sf_putint(imagez,"n5",image_par->ny);

    for(i=0;i<image_par->ny*image_par->nx*image_par->nhy*image_par->nhx;i++) {
	sf_floatwrite(tmpimg,image_par->nz,imagex);
	sf_floatwrite(tmpimg,image_par->nz,imagez);
    }
    free(tmpimg);

    image_par->imagex = imagex;
    image_par->imagez = imagez;
}

void shot_image_sinc_table(struct shot_image_par_type *image_par){
    int ihx;
    float *sincx_table=image_par->sincx_table;
    image_par->sincx_table=sf_floatalloc(8*2001);
    for(ihx=0;ihx<2001;ihx++) mksinc((float)(ihx)*0.0005,8,sincx_table+ihx*8);
}

void shot_image_cross(sf_complex *wld_s,sf_complex *wld_r,float *image_z,int nx,int ny,struct shot_image_par_type image_par){
    sf_complex ws,wr;
    int ihx,ihy,ix,iy,iihx,iihy,iohx,iohy;
    int shift_ix_s,shift_iy_s,shift_ix_r,shift_iy_r;
    int nimg[3];
    d4(nx,image_par.nhy,image_par.nhx,nimg);
    iohy=(int)(image_par.ohy/image_par.dy); iohx=(int)(image_par.ohx/image_par.dx);
    for(iy=0;iy<ny;iy++){
	for(ix=0;ix<nx;ix++){
	    for(ihy=0; ihy<image_par.nhy;ihy++){
		//hy=ihy*image_par.dy+image_par.ohy;   iihy=(int)(hy/image_par.dy);
		iihy=ihy+iohy;
		for(ihx=0;ihx<image_par.nhx;ihx++){    
		    //hx=ihx*image_par.dx+image_par.ohx;  iihx=(int)(hx/image_par.dx);
		    iihx=ihx+iohx;
		    shift_ix_s=ix-iihx;  shift_iy_s=iy-iihy; shift_ix_r=ix+iihx; shift_iy_r=iy+iihy;
		    if (shift_ix_s <0) shift_ix_s=0; 
            if (shift_ix_s >=nx) shift_ix_s=nx-1;
		    if (shift_ix_r <0) shift_ix_r=0; 
            if (shift_ix_r >=nx) shift_ix_r=nx-1;
		    if (shift_iy_s <0) shift_iy_s=0; 
            if (shift_iy_s >=ny) shift_iy_s=ny-1;
		    if (shift_iy_r <0) shift_iy_r=0; 
            if (shift_iy_r >=ny) shift_iy_r=ny-1;
		    ws=wld_s[i2(shift_iy_s,shift_ix_s,nx)]; wr=wld_r[i2(shift_iy_r,shift_ix_r,nx)];  
		    image_z[i4(iy,ix,ihy,ihx,nimg)]+=crealf( conjf(ws)*wr);
            
		}
	    }
	}
    }

}


void shot_image_cross_tilt(sf_complex *wld_s,sf_complex *wld_r,float *image_z,int nx,int ny,
			   struct shot_image_par_type image_par,int nx_iz)
/*< cross correlation >*/
{
    sf_complex ws,wr;
    int ihx,ihy,ix,iy,iihx,iihy,iohx,iohy;
    int shift_ix_s,shift_iy_s,shift_ix_r,shift_iy_r;
    int nimg[3];
    d4(nx,image_par.nhy,image_par.nhx,nimg);
    iohy=(int)(image_par.ohy/image_par.dy); iohx=(int)(image_par.ohx/image_par.dx);
    for(iy=0;iy<ny;iy++){
	for(ix=0;ix<nx_iz;ix++){
	    for(ihy=0; ihy<image_par.nhy;ihy++){
		//hy=ihy*image_par.dy+image_par.ohy;   iihy=(int)(hy/image_par.dy);
		iihy=ihy+iohy;
		for(ihx=0;ihx<image_par.nhx;ihx++){
		    //hx=ihx*image_par.dx+image_par.ohx;  iihx=(int)(hx/image_par.dx);
		    iihx=ihx+iohx;
		    shift_ix_s=ix-iihx;  shift_iy_s=iy-iihy; shift_ix_r=ix+iihx; shift_iy_r=iy+iihy;
		    if (shift_ix_s <0) shift_ix_s=0; 
            if (shift_ix_s >=nx_iz) shift_ix_s=nx_iz-1;
		    if (shift_ix_r <0) shift_ix_r=0; 
            if (shift_ix_r >=nx_iz) shift_ix_r=nx_iz-1;
		    if (shift_iy_s <0) shift_iy_s=0; 
            if (shift_iy_s >=ny) shift_iy_s=ny-1;
		    if (shift_iy_r <0) shift_iy_r=0; 
            if (shift_iy_r >=ny) shift_iy_r=ny-1;
		    ws=wld_s[i2(shift_iy_s,shift_ix_s,nx)]; wr=wld_r[i2(shift_iy_r,shift_ix_r,nx)]; // wait for nx
/*
  de=0.0;
          
  for (ismooth=1;ismooth<=nsmooth;ismooth++){
  ixsmooth=shift_ix_s-nsmooth/2+ismooth;
  if (ixsmooth<0) ixsmooth=0;
  if (ixsmooth>=nx_iz) ixsmooth=nx_iz-1;
  de=de+wld_s[i2(shift_iy_s,ixsmooth,nx)]*conjg(wld_s[i2(shift_iy_s,ixsmooth,nx)]);
  }
  de=de/nsmooth;
*/ 
         
		    //image_z[i4(iy,ix,ihy,ihx,nimg)]+=real( conjg(ws)*wr/(conjg(ws)*ws+eps)   );
		    //image_z[i4(iy,ix,ihy,ihx,nimg)]+=real( conjg(ws)*wr/(de)   );
		    image_z[i4(iy,ix,ihy,ihx,nimg)]+=crealf( conjf(ws)*wr);
		    //image_z[i4(iy,ix,ihy,ihx,nimg)]+=real( conjg(ws)*ws);
       
		}
	    }
	}
    }

}

void shot_image_stack(float *img,float mig_min_x,float mig_min_y,int mig_nx,int mig_ny,sf_file tag,
		      struct shot_image_par_type image_par)
/*< stack image >*/
{
    float *tmp_img;
    float img_min_x,img_max_x,img_min_y,img_dx,img_dy;
    float mig_max_x,rite_min_x,rite_max_x;
    int img_ny,img_nx,rite_nx;
    int n_hy_hx,img_nhx,img_nhy,nz;
    int ntmp_img[2],nimg[3];
    int rite_ix_b,array_ix_b;
    int block_size,block_i;
    int iy,rite_iy,ix,i_hy_hx,iz;

    img_min_x=image_par.min_x; img_max_x=image_par.max_x; img_min_y=image_par.min_y;
    img_dx=image_par.dx; img_dy=image_par.dy;             img_ny=image_par.ny; img_nx=image_par.nx;
    img_nhx=image_par.nhx; img_nhy=image_par.nhy;  nz=image_par.nz;

    mig_max_x=((float)(mig_nx-1))*img_dx+mig_min_x;

    rite_min_x=SF_MAX(img_min_x,mig_min_x); 
    rite_max_x=SF_MIN(img_max_x,mig_max_x);
    rite_nx=(rite_max_x-rite_min_x)/img_dx+1;

    n_hy_hx= img_nhy*img_nhx;  d3(n_hy_hx,nz,ntmp_img); d4(mig_ny,mig_nx,n_hy_hx,nimg);

    if (rite_nx >=1){
	rite_ix_b= (rite_min_x-img_min_x)/img_dx;
	array_ix_b=(rite_min_x-mig_min_x)/img_dx;
	tmp_img=sf_floatalloc(rite_nx*n_hy_hx*nz);
	block_size=4*n_hy_hx*nz;
	for(iy=0;iy<mig_ny;iy++){
	    vector_value_f(tmp_img,0.0,rite_nx*n_hy_hx*nz);
	    rite_iy= iy+(mig_min_y-img_min_y)/img_dy;
	    if (rite_iy < img_ny && rite_iy>= 0){
		block_i=rite_iy*img_nx+rite_ix_b;
		sf_seek(tag,block_i*block_size,SEEK_SET);
		sf_floatread(tmp_img,rite_nx*n_hy_hx*nz,tag);
		for(ix=0;ix<rite_nx;ix++){
		    for(i_hy_hx=0;i_hy_hx<n_hy_hx;i_hy_hx++){
			for(iz=0;iz<nz;iz++){
			    tmp_img[i3(ix,i_hy_hx,iz,ntmp_img)]+=img[i4(iz,iy,array_ix_b+ix,i_hy_hx,nimg)]; 
			}
		    }
		} 
		sf_seek(tag,block_i*block_size,SEEK_SET);
		sf_floatwrite(tmp_img,rite_nx*n_hy_hx*nz,tag);
	    }
	}
	free(tmp_img);
    }
}




void shot_image_decompose(float *image_hx,float *image_hz,int nx,int nz,float rotate_angle,struct shot_image_par_type image_par)
/*< decompose >*/
{

    int ix,iz;
    int iblock,pointmove;
    float *tmp_image_hx,*tmp_image_hz;

    float tmpmax;

    for(iz=54;iz<55;iz++){
	for(ix=635;ix<636;ix++){
//for(iz=0;iz<nz;iz++){
//  for(ix=0;ix<nx;ix++){
	    iblock=iz*nx+ix;
	    pointmove=image_par.nhx*image_par.nhy*iblock;
	    tmp_image_hx=image_hx+pointmove;  tmp_image_hz=image_hz+pointmove;
	    shot_image_decompose_1imagepoint(tmp_image_hx,tmp_image_hz,nx,nz,rotate_angle,image_par); 
	    tmpmax=maxval(tmp_image_hz,image_par.nhx);
	    if (tmpmax>10000000.0) sf_warning("ix=%d,iz=%d,max=%f",ix,iz,tmpmax);
	}
    }
}



void shot_image_decompose_1imagepoint(float *image_hx,float *image_hz,int nx,int nz,float rotate_angle,
				      struct shot_image_par_type image_par)
/*< one image point >*/
{
    float *image_dip,*sincx_table;
    int i_ihz,i_ihx,ihdip,ihx,ihz,ii_hz,ii_hx,iwx;
    float hx,hz,whx,whz,ohx,hdip,ddx;
    float sincz[8],sincx[8];
    float *i0sinc;
    float dwx;
    int nit=8;

    dwx=0.0005;

    sincx_table=image_par.sincx_table;
    image_dip=sf_floatalloc(image_par.nhx*nit);
    //wait ohx
    ohx=-image_par.dx*image_par.nhx/2;
    ddx=image_par.dx/(float)(nit);
    vector_value_f(image_dip,0.0,image_par.nhx*nit);
    for(ihdip=0;ihdip<image_par.nhx*nit;ihdip++){
	hdip=(int)(ihdip)*ddx+ohx;
	ihx=(int)((hdip-ohx)/image_par.dx);   whx=(hdip-(ohx+(float)(ihx)*image_par.dx))/image_par.dx;
	iwx=(int)(whx/dwx);                 //wwx=(whx- ((float)(iwx)*dwx))/dwx;
	i0sinc=sincx_table+8*iwx;
	for(i_ihx=0;i_ihx<8;i_ihx++)
	    sincx[i_ihx]=i0sinc[i_ihx];
	//sincx[i_ihx-1]=(1-wwx)*sincx_table[i2(iwx,i_ihx-1,8)]+wwx*sincx_table[i2(iwx+1,i_ihx-1,8)];
	for(i_ihx=1;i_ihx<=8;i_ihx++){
	    ii_hx=i_ihx+ihx-4;
	    if (ii_hx>=0 && ii_hx < image_par.nhx )  image_dip[ihdip]+=sincx[i_ihx-1]*image_hx[ii_hx];
	} //for(i_ihx)
    }//for(ihdip)
    printf("maxhdip=%f\n",maxval(image_dip,image_par.nhx*nit));
    vector_value_f(image_hx,0.0,image_par.nhx*image_par.nhy);
    vector_value_f(image_hz,0.0,image_par.nhx*image_par.nhy);
    for(ihdip=0;ihdip<image_par.nhx*nit;ihdip++){
	hdip=ohx+(float)(ihdip)*ddx;
	if (hdip==0){
	    hx=0.0;   hz=0.0;
	}
	else{
	    hx=hdip/cos(rotate_angle);  hz=hdip/sin(rotate_angle);
	}
	ihx=(int)((hx-ohx)/image_par.dx);  ihz=(int)((hz-ohx)/image_par.dx);
    
	if (ihx >= 0 && ihx<image_par.nhx ){ 
	    whx=( hx-(ohx+(float)(ihx)*image_par.dx) )/image_par.dx;  
	    iwx=(int)(whx/dwx);   //wwx=(whx- ((float)(iwx)*dwx))/dwx;
	    i0sinc=sincx_table+8*iwx;
	    for(i_ihx=0;i_ihx<8;i_ihx++) 
		sincx[i_ihx]=i0sinc[i_ihx];
	    //sincx[i_ihx-1]=(1-wwx)*sincx_table[i2(iwx,i_ihx-1,8)]+wwx*sincx_table[i2(iwx+1,i_ihx-1,8)];
	    //sincx(:)=(1-wwx)*sincx_table(:,iwx)+wwx*sincx_table(:,iwx+1);
	    for(i_ihx=1;i_ihx<=8;i_ihx++){
		ii_hx=i_ihx+ihx-4;
		if (ii_hx >=0 && ii_hx < image_par.nhx )  image_hx[ii_hx]+= sincx[i_ihx-1]*image_dip[ihdip];
	    }//for(i_ihx)
	}
	if (ihz >=0 && ihz <image_par.nhx ){
	    whz=( hz-(ohx+(float)(ihz)*image_par.dx) )/image_par.dx;
	    iwx=(int)(whz/dwx);   //wwx=(whz-((float)(iwx)*dwx))/dwx;
	    i0sinc=sincx_table+8*iwx;
	    for(i_ihx=0;i_ihx<8;i_ihx++) 
		sincz[i_ihx]=i0sinc[i_ihx];
	    //sincz[i_ihx-1]=(1-wwx)*sincx_table[i2(iwx,i_ihx-1,8)]+wwx*sincx_table[i2(iwx+1,i_ihx-1,8)];
	    for(i_ihz=1;i_ihz<=8;i_ihz++){
		ii_hz=i_ihz+ihz-4;
		if (ii_hz >0 && ii_hz <= image_par.nhx ) image_hz[ii_hz]+= sincz[i_ihz-1]*image_dip[ihdip];
		if (ihdip==125) printf("i_ihz=%d,imagehz=%f,sinc=%f,imagedip=%f\n",i_ihz,image_hz[ii_hz],sincz[i_ihz-1],image_dip[ihdip]);
	    }//for(i_ihz)
	}
	printf("iz=%d,maxhz=%f\n",ihdip,maxval(image_hz,image_par.nhx));
    }// for(ihdip)

    free(image_dip);
}


void shot_image_release(struct shot_image_par_type *image_par)
/*< release >*/
{
    free(image_par->sincx_table);
}
