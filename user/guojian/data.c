/* Shot data routines. */
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

#include "data.h"
#include "arrayindex.h"

#ifndef _data_h

/* begin data */
struct shot_data_par_type
{
    float dhx,dhy,dw,dsx,dsy;
    float ohx,ohy,ow,osx,osy;
    int nhx,nhy,nw,nsx,nsy;
    sf_complex *wavelet;
    sf_complex *tmp_wld;
    sf_file data;
};
/* begin plane */
struct shot_plane_par_type
{
    float ow,dw,ox,oy,dx,dy,opx,dpx,opy,dpy;
    int   npx,npy,nx,ny,nw;
    float osx,osy,dsx,dsy;
    int   nsx,nsy;
    float ohx,ohy,dhx,dhy;
    int   nhx,nhy;
};
/* end plane */
/*^*/

#endif

struct shot_data_par_type shot_data_init(void)
/*< initialize >*/
{
    struct shot_data_par_type data_par;
    sf_file data;
    int ndims, n[SF_MAX_DIM];
   
    data=sf_input("Data");
    ndims = sf_filedims(data,n);

    printf("nndims=%d\n",ndims);
    if (ndims<5){
	data_par.nsy=1; data_par.dsy=10.0; data_par.osy=0.0;
    }
    if (ndims==5)
    {
	data_par.nsy=n[4]; sf_histfloat(data,"d5",&data_par.dsy); sf_histfloat(data,"o5",&data_par.osy);
    }
    if (ndims<4){
	data_par.nsx=1; data_par.dsx=10.0; data_par.osx=0.0;
    }
    else{
	data_par.nsx=n[3]; sf_histfloat(data,"d4",&data_par.dsx); sf_histfloat(data,"o4",&data_par.osx);
    }
    data_par.nw= n[2]; sf_histfloat(data,"d3",&data_par.dw);  sf_histfloat(data,"o3",&data_par.ow);
    data_par.nhy=n[1]; sf_histfloat(data,"d2",&data_par.dhy); sf_histfloat(data,"o2",&data_par.ohy);
    data_par.nhx=n[0]; sf_histfloat(data,"d1",&data_par.dhx); sf_histfloat(data,"o1",&data_par.ohx);
  
    if (data_par.nsx==0) data_par.nsx=1;
    data_par.wavelet=sf_complexalloc(data_par.nw);
    data_par.tmp_wld=sf_complexalloc(data_par.nhx*data_par.nhy);

    data_par.data = data;
    return data_par;

}

void shot_data_release(struct shot_data_par_type *data_par)
/*< release >*/
{
    free(data_par->wavelet); free(data_par->tmp_wld);
}


void shot_data_wavelet(struct shot_data_par_type data_par)
/*< data wavelet should run after shot_data_init >*/
{
    sf_complex *wavelet=data_par.wavelet;
    sf_file wavlt;
    int n1;

    wavlt=sf_input("Wavelet");
    if (!sf_histint(wavlt,"n1",&n1) || n1 != data_par.nw ) 
	sf_error("wavelet should have the same frequency  as data");
    sf_complexread(wavelet,n1,wavlt);

    sf_fileclose(wavlt);
}

// data_1shot_init in f90
void shot_data_grab_mig_par(int i_sx,int i_sy,float *mig_x_min,float *mig_y_min, int *nx,int *ny,struct shot_data_par_type data_par )
{
    float sx,sy;
    sx=((float)i_sx)*data_par.dsx+data_par.osx; sy=((float) i_sy)*data_par.dsy+data_par.osy; 
    *mig_x_min=sx+data_par.ohx;                  *mig_y_min=sy+data_par.ohy;
    *nx=data_par.nhx;                           *ny=data_par.nhy;
}

void data_source(float min_x,float min_y,int nx,int ny,int iw,int i_sx,int i_sy,sf_complex *wld,struct shot_data_par_type data_par)
{
    float sx,sy;
    int ix,iy;
    float data_dsx, data_osx, data_dsy, data_osy, data_dhx, data_dhy;
    sf_complex * wavelet=data_par.wavelet;
    data_dsx=data_par.dsx; data_osx=data_par.osx; data_dsy=data_par.dsy; data_osy=data_par.osy;
    data_dhx=data_par.dhx; data_dhy=data_par.dhy;

    sx=i_sx*data_dsx+data_osx;               sy=i_sy*data_dsy+data_osy;
    ix=(sx-min_x)/data_dhx; iy=(sy-min_y)/data_dhy;
    if (ix <0 || iy <0 || ix>=nx || iy>=ny) sf_error("source point is outside of wavefield area");
    vector_value_c(wld,sf_cmplx(0.0,0.0),nx*ny);  //wld=0.0
    wld[i2(iy,ix,nx)]=wavelet[iw]; if (ix+1 <nx)   wld[i2(iy,ix+1,nx)]=0.5*wavelet[iw]; if (ix-1 >=0  ) wld[i2(iy,ix-1,nx)]=0.5*wavelet[iw];
}



void data_receiver(float min_x,float min_y,int nx,int ny,int iw,int i_sx,int i_sy,sf_complex *wld,struct shot_data_par_type data_par)
{
    sf_complex *tmp_wld=data_par.tmp_wld;
    float sx,sy,max_x,max_y,min_x_d,min_y_d,max_x_d,max_y_d;
    float max_x_read,max_y_read,min_x_read,min_y_read;
    int  nx_read,ny_read,i_ox_read,i_oy_read,i_ox_wld,i_oy_wld;
    int n_block,blocksize;
    int i_s,ix,iy,ix_wld,iy_wld,ix_read,iy_read;
    float data_dhx,data_dhy,data_ohx,data_ohy,data_osx,data_osy,data_dsx,data_dsy;
    int data_nhx,data_nhy,data_nsx,data_nw;
    data_dhx=data_par.dhx; data_ohx=data_par.ohx; data_nhx=data_par.nhx; 
    data_dhy=data_par.dhy; data_ohy=data_par.ohy; data_nhy=data_par.nhy; 
    data_osx=data_par.osx; data_dsx=data_par.dsx; data_nsx=data_par.nsx; 
    data_osy=data_par.osy; data_dsy=data_par.dsy;
    data_nw=data_par.nw;


    vector_value_c(wld,sf_cmplx(0.0,0.0),nx*ny);  //wld=0.0
    sx=i_sx*data_dsx+data_osx;                    sy=i_sy*data_dsy+data_osy; 
    max_x=min_x+(nx-1)*data_dhx;                  max_y=min_y+(ny-1)*data_dhy;
    min_x_d=sx+data_ohx;                          min_y_d=sy+data_ohy;         
    max_x_d=min_x_d+data_dhx*(data_nhx-1);        max_y_d=min_y_d+data_dhy*(data_nhy-1);
    max_x_read=SF_MIN(max_x_d,max_x);                max_y_read=SF_MIN(max_y_d,max_y);
    min_x_read=SF_MAX(min_x_d,min_x);                min_y_read=SF_MAX(min_y_d,min_y);
    nx_read=(max_x_read-min_x_read)/data_dhx+1;   ny_read=(max_y_read-min_y_read)/data_dhy+1;
    i_ox_read=(min_x_read-min_x_d)/data_dhx;      i_oy_read=(min_y_read-min_y_d)/data_dhy;
    i_ox_wld=(min_x_read-min_x)/data_dhx;         i_oy_wld=(min_y_read-min_y)/data_dhy;
    i_s=i_sy*data_nsx+i_sx;                       n_block=i_s*data_nw+iw;                 
    blocksize=data_nhy*data_nhx*8;
    //printf("blocksize=%d  n_block=%d\n",blocksize,n_block);
    sf_seek(data_par.data,n_block*blocksize,SEEK_SET);
    sf_complexread(tmp_wld,blocksize,data_par.data);
    for(iy=0;iy<ny_read;iy++){
	iy_wld=iy+i_oy_wld;iy_read=iy+i_oy_read;
	for(ix=0;ix<nx_read;ix++){
	    ix_wld=ix+i_ox_wld;ix_read=ix+i_ox_read; 
	    wld[i2(iy_wld,ix_wld,nx)]=tmp_wld[i2(iy_read,ix_read,data_nhx)];
	}
    }
    //wld(i_ox_wld:i_ox_wld+nx_read-1,i_oy_wld:i_oy_wld+ny_read-1)=tmp_wld(i_ox_read:i_ox_read+nx_read-1,i_oy_read:i_oy_read+ny_read-1);
}

void shot_data_plane_grab_mig_par(int i_sx,int i_sy,float *mig_x_min,float *mig_y_min,int *nx,int *ny,
				  struct shot_data_par_type data_par,struct shot_plane_par_type *plane_par )
/*< get mig par >*/
{
    float plane_osx,plane_dsx,plane_osy,plane_dsy;
    int   plane_nsx,plane_nsy;

    *mig_x_min=data_par.ohx;  *mig_y_min=data_par.ohy;  *nx=data_par.nhx;   *ny=data_par.nhy;

    if (!sf_histfloat(data_par.data,"plane_os_x", &plane_osx)) plane_osx=0.0;
    if (!sf_histfloat(data_par.data,"plane_ds_x", &plane_dsx)) plane_dsx=1.0;
    if (!sf_histint(data_par.data,  "plane_ns_x", &plane_nsx)) plane_nsx=1;
    if (!sf_histfloat(data_par.data,"plane_os_y", &plane_osy)) plane_osy=0.0;
    if (!sf_histfloat(data_par.data,"plane_ds_y", &plane_dsy)) plane_dsy=1.0;
    if (!sf_histint(data_par.data,  "plane_ns_y", &plane_nsy)) plane_nsy=1;
  
    plane_par->osx=plane_osx; plane_par->dsx=plane_dsx; plane_par->nsx=plane_nsx;
    plane_par->osy=plane_osy; plane_par->dsy=plane_dsy; plane_par->nsy=plane_nsy;
    plane_par->ohx=data_par.ohx; plane_par->dhx=data_par.dhx; plane_par->nhx=data_par.nhx;
    plane_par->ohy=data_par.ohy; plane_par->dhy=data_par.dhy; plane_par->nhy=data_par.nhy;
}

void shot_plane_source(float min_x,float min_y,int nx,int ny,int iw,int i_px,int i_py,
		       sf_complex *wld,struct shot_data_par_type data_par,struct shot_plane_par_type plane_par)
/*< plane-wave source >*/
{
    float x_dt,t_x;
    int iy,ix,isx,isy;
    float plane_s_x,px,w;
    float plane_dx,plane_dw,plane_ow,plane_opx,plane_dpx,plane_osx,plane_dsx;
    int plane_nsx,plane_nsy;
    sf_complex * wavelet=data_par.wavelet;
    plane_dw=data_par.dw; plane_ow=data_par.ow; 
    plane_dpx=data_par.dsx; plane_opx=data_par.osx; 
    plane_nsx=plane_par.nsx; plane_dsx=plane_par.dsx; plane_osx=plane_par.osx;
    plane_nsy=plane_par.nsy; 
    plane_dx= data_par.dhx;
 
    w=(float)(iw)*plane_dw+plane_ow;
    px=(float)(i_px)*plane_dpx+plane_opx;
    vector_value_c(wld,sf_cmplx(0.0,0.0),nx*ny);  //wld=0.0
    x_dt=plane_dsx*px;   //plane_par.dsx is dsx of point source of original dataset
    for(isy=0;isy<plane_nsy;isy++){  // plane_nsy  nsy of point source of original dataset
	iy=isy;
	for(isx=0;isx<plane_nsx;isx++){  // plane_nsx nsx of point source
	    plane_s_x=(float)(isx)*plane_dsx+plane_osx;
	    ix=(plane_s_x-min_x)/fabs(plane_dx);      t_x=(float)(isx)*x_dt; //t_x=(float)(isx-isx0)*x_dt;
	    if (ix>=0 && ix <nx)
		wld[i2(iy,ix,nx)]+= wavelet[iw]*sf_cmplx(cosf(w*t_x),sinf(w*t_x));
	    if (ix-1>=0 && ix-1<nx)
		wld[i2(iy,ix-1,nx)]+= 0.5*wavelet[iw]*sf_cmplx(cosf(w*t_x),sinf(w*t_x));
	    if (ix+1>=0 && ix+1<nx)
		wld[i2(iy,ix+1,nx)]+= 0.5*wavelet[iw]*sf_cmplx(cosf(w*t_x),sinf(w*t_x));
	}
    }



}


void shot_plane_source_3d(float min_x,float min_y,int nx,int ny,int iw,int i_px,int i_py,sf_complex *wld,struct shot_data_par_type data_par,struct shot_plane_par_type plane_par)
{
    float x_dt,y_dt,t_x,t_y;
    int iy,ix,isx,isy;
    float plane_s_x,plane_s_y,px,py,w;
    float plane_dx,plane_dy,plane_dw,plane_ow,plane_opx,plane_dpx,plane_opy,plane_dpy,plane_osx,plane_dsx,plane_osy,plane_dsy;
    int plane_nsx,plane_nsy;
    sf_complex * wavelet=data_par.wavelet;
    sf_complex phaseshift;
    plane_dw=data_par.dw; plane_ow=data_par.ow; 
    plane_dx= data_par.dhx;  plane_dy=data_par.dhy;
    plane_dpx=data_par.dsx; plane_opx=data_par.osx; 
    plane_dpy=data_par.dsy; plane_opy=data_par.osy; 
    plane_nsx=plane_par.nsx; plane_dsx=plane_par.dsx; plane_osx=plane_par.osx;
    plane_nsy=plane_par.nsy; plane_dsy=plane_par.dsy; plane_osy=plane_par.osy;

    w=(float)(iw)*plane_dw+plane_ow;
    px=(float)(i_px)*plane_dpx+plane_opx;     py=(float)(i_py)*plane_dpy+plane_opy;
    //printf("px=%f,py=%f,iw=%d\n",px,py,iw);
    vector_value_c(wld,sf_cmplx(0.0,0.0),nx*ny);  
    //plane_par.dsx is dsx of point source of original dataset
    x_dt=plane_dsx*px;     y_dt=plane_dsy*py;
    //printf("plane_nsy=%d,plane_nsx=%d\n",plane_nsy,plane_nsx);
  
    for(isy=0;isy<plane_nsy;isy++){  // plane_nsy  nsy of point source of original dataset
	plane_s_y=(float)(isy)*plane_dsy+plane_osy;
	iy=(plane_s_y-min_y)/fabs(plane_dy); t_y=(float)(isy)*y_dt;
	//printf("-----------------------iysource=%d,ny=%d\n",iy,ny);
	for(isx=0;isx<plane_nsx;isx++){  // plane_nsx nsx of point source
	    plane_s_x=(float)(isx)*plane_dsx+plane_osx;
	    ix=(plane_s_x-min_x)/fabs(plane_dx);       t_x=(float)(isx)*x_dt;
	    phaseshift=sf_cmplx(cosf(w*t_x),sinf(w*t_x))*sf_cmplx(cosf(w*t_y),sinf(w*t_y));
	    wld[i2(iy,ix,nx)]+=wavelet[iw]*phaseshift;
	    if (ix+1 < nx){
		wld[i2(iy,ix+1,nx)]+=0.5*wavelet[iw]*phaseshift;
		if (iy+1 <ny) wld[i2(iy+1,ix+1,nx)]+=0.35*wavelet[iw]*phaseshift;
		if (iy-1 >=0) wld[i2(iy-1,ix+1,nx)]+=0.35*wavelet[iw]*phaseshift;
	    }
	    if (ix-1 >=0){
		wld[i2(iy,ix-1,nx)]+=0.5*wavelet[iw]*phaseshift;
		if (iy+1 <ny) wld[i2(iy+1,ix-1,nx)]+=0.35*wavelet[iw]*phaseshift; 
		if (iy-1 >=0) wld[i2(iy-1,ix-1,nx)]+=0.35*wavelet[iw]*phaseshift; 
	    }
	    if (iy+1 <ny) wld[i2(iy+1,ix,nx)]+=0.5*wavelet[iw]*phaseshift; 
	    if (iy-1 >=0) wld[i2(iy-1,ix,nx)]+=0.5*wavelet[iw]*phaseshift;
	}
    }
}





void shot_plane_receiver(float min_x,float min_y,int nx,int ny,int iw,int i_sx,int i_sy,
			 sf_complex *wld,struct shot_data_par_type data_par,struct shot_plane_par_type plane_par)
/*< plane-wave receiver >*/
{  

    float  min_x_d,min_y_d,max_x_d,max_y_d,max_x,max_y;
    float max_x_read,max_y_read,min_x_read,min_y_read;
    int nx_read,ny_read,i_ox_read,i_oy_read,i_ox_wld,i_oy_wld,n_block,blocksize,i_s;
    float plane_dx,plane_dy,plane_ox,plane_oy;
    int plane_nx,plane_ny;
    int ix,iy;
    sf_complex *tmp_wld=data_par.tmp_wld;

    plane_dx=data_par.dhx;  plane_dy=data_par.dhy; plane_ox=data_par.ohx; plane_oy=data_par.ohy;
    plane_nx=data_par.nhx;  plane_ny=data_par.nhy;


    i_s=i_sy*data_par.nsx+i_sx;
    max_x=min_x+(float)(nx-1)*fabs(plane_dx);              max_y=min_y+(float)(ny-1)*fabs(plane_dy);
    min_x_d=plane_ox;                                      min_y_d=plane_oy;
    max_x_d=min_x_d+fabs(plane_dx)*(plane_nx-1);           max_y_d=min_y_d+fabs(plane_dy)*(plane_ny-1);
    max_x_read=SF_MIN(max_x_d,max_x);                             max_y_read=SF_MIN(max_y_d,max_y);
    min_x_read=SF_MAX(min_x_d,min_x);                             min_y_read=SF_MAX(min_y_d,min_y);
    nx_read=(max_x_read-min_x_read)/fabs(plane_dx)+1;      ny_read=(max_y_read-min_y_read)/fabs(plane_dy)+1;
    i_ox_read=(min_x_read-min_x_d)/fabs(plane_dx);         i_oy_read=(min_y_read-min_y_d)/fabs(plane_dy);
    i_ox_wld=(min_x_read-min_x)/fabs(plane_dx);            i_oy_wld=(min_y_read-min_y)/fabs(plane_dy);
    n_block=i_s*data_par.nw+iw;        blocksize=plane_ny*plane_nx;
//printf("plane_ny=%d, plane_nx=%d,nblock=%d,i_s=%d\n",plane_ny,plane_nx,n_block,i_s);
    sf_seek(data_par.data,n_block*blocksize*sizeof(sf_complex),SEEK_SET);
    sf_complexread(tmp_wld,blocksize,data_par.data);
    for(iy=0;iy<ny_read;iy++)
	for(ix=0;ix<nx_read;ix++)
	    wld[i2(i_oy_wld+iy,i_ox_wld+ix,nx)]=tmp_wld[i2(i_oy_read+iy,i_ox_read+ix,plane_nx)];

} // shot_plane_receiver




/*
  sf_complex * shot_data_wavelet_read(filehead *wavlt)
  {
  int stat;
  sf_complex *wavelet;
  stat=file_open_input("Wavelet",wavlt);
  printf("n is equal to %d",wavlt->n[0]);
  wavelet=allocatec(wavlt->n[0]);
  if ( sreed("Wavelet",wavelet,wavlt->n[0]*8) != wavlt->n[0]*8 ) sf_error("err wavelet");
  return wavelet;
  }

  void shot_data_clear(sf_complex *wavelet){
  free(wavelet);
  }

  void shot_data_source_modeling(float sx,float min_x,int nx,sf_complex *wld,sf_complex *wavelet, int iw, float dhx)
  {

  int ix,iy;
  int ny;
  int iix,iiy;
  ny=1;iy=0;  ix=(sx-min_x)/dhx;
  if (ix <0 || ix>=nx ) sf_error("err create data_source");

  for ( iiy=0;iiy < ny; iiy++){
  for (iix=0;iix<nx; iix++){
  wld[i2(iiy,iix,nx)]=sf_cmplx(0.0,0.0);
  }
  }
  wld[i2(iy,ix,nx)]=wavelet[iw];
  if (ix+1 <=nx) wld[i2(iy,ix+1,nx)]=0.5*wavelet[iw];
  if (ix-1 >0  ) wld[i2(iy,ix-1,nx)]=0.5*wavelet[iw];
  }
*/
