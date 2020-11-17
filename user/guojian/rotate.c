/* Rotate routines. */
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

#include "rotate.h"
#include "fft.h"
#include "arrayindex.h"

#ifndef _rotate_h

struct shot_rotate_par_type
{
    float rotate_angle; 
    int ix_rotate_center;
    int nx,nz,nx_rotate,nz_rotate,ny_rotate,nz_left,nz_right;
    int nz_s_rotate;
    float dx,dz;
    sf_complex *wld_before_rotate, *wld_after_rotate;
    float *sinc;
    int Nsinc;
    /* add in 3d */
    float dy,oyout_azimuth,oxout_tilt,azimuth;
    int ny,nyout_azimuth,nxout_azimuth,nxout_tilt,nzout_tilt;
    float sinazimuth,cosazimuth,sintilt,costilt;
    float ox,oy;
    int nx_s_rotate;
    int *bx_s_rotate;  
    /* end add in 3d */
};
/*^*/

#endif

struct shot_rotate_par_type shot_rotate_init(float rotate_angle,int nx,int nz,float dx, float dz)
/*< initialize >*/
{
    struct shot_rotate_par_type rotate_par;
    float x,z,wx;
    int i;
    rotate_par.nx=nx; rotate_par.dx=dx; rotate_par.dz=dz; rotate_par.nz=nz;
    rotate_par.rotate_angle=rotate_angle; 
    if( 0.0==rotate_angle){
	rotate_par.nx_rotate=kissfftn(nx); rotate_par.nz_rotate=nz; rotate_par.ny_rotate=1;
	rotate_par.nz_s_rotate=1;  rotate_par.ix_rotate_center=0;
    }
    else{
	rotate_par.nz_rotate=((nz-1)*dz*cos(fabs(rotate_angle))+(nx-1)*dx*sin(fabs(rotate_angle)))/dz+1;
	rotate_par.ix_rotate_center=((nz-1)*dz*sin(fabs(rotate_angle))-0.0)/dx;
	//printf("rororotate=%d, nx=%d,center=%d\n",tmpn,nx,rotate_par.ix_rotate_center);
	rotate_par.nx_rotate=kissfftn( ((nx-1)*dx*cosf(fabsf(rotate_angle)))/dx+rotate_par.ix_rotate_center+1);
	rotate_par.nz_s_rotate=((nx-1)*dx*sinf(fabsf(rotate_angle)))/dz+1;
	rotate_par.ny_rotate=1;
    
    
	x=(float)(rotate_par.ix_rotate_center)*dx; z=x*fabs(cos(rotate_par.rotate_angle)/sin(rotate_par.rotate_angle));
	rotate_par.nz_left=(int)(z/dz)+1; rotate_par.nz_right=rotate_par.nz_rotate-rotate_par.nz_left;
    }
    rotate_par.Nsinc=2000;
    rotate_par.sinc=sf_floatalloc(8*(rotate_par.Nsinc+1));
  
    for (i=0;i<=rotate_par.Nsinc;i++){
	wx=(float)(i)*1.0/(float)rotate_par.Nsinc;
	mksinc(wx,8,rotate_par.sinc+8*i); 
    }/*
       for (i=0;i<rotate_par.Nsinc;i++){
       printf("%d\n",i);
       for(j=0;j<8;j++) printf(" %f ",*(rotate_par.sinc+8*i+j));
       printf("\n");
       }*/
    return rotate_par;
}

void shot_rotate_release(struct shot_rotate_par_type *rotate_par)
/*< release >*/
{
    free(rotate_par->sinc); 
}

void shot_rotate_get_bx_ex_iz(float *dkx,int *bx,int *ex, int iz,int boundary_nx, struct shot_rotate_par_type rotate_par)
/*< rotate >*/
{

    float dx,dz,rotate_angle;
    int nx_rotate,ix_rotate_center;
    int nz_left,nz_right;
    int tmp_nx,nx_iz;
    float x,z,tan_rotate_angle;
    dx=rotate_par.dx; dz=rotate_par.dz; rotate_angle=rotate_par.rotate_angle;
    nx_rotate=rotate_par.nx_rotate;     ix_rotate_center=rotate_par.ix_rotate_center;
    nz_left=rotate_par.nz_left; nz_right=rotate_par.nz_right;
    tan_rotate_angle=sin(fabs(rotate_angle))/cos(fabs(rotate_angle));
    if (iz < nz_left) {
	z=(float)(iz)*dz; 
	x=z*tan_rotate_angle;
	tmp_nx=(int)(x/dx)+boundary_nx;
	*bx=SF_MAX(0,ix_rotate_center-tmp_nx);
    }
    else{
	z=(float)(iz-nz_left)*dz; x=z/tan_rotate_angle;
	tmp_nx=(int)(x/dx)-boundary_nx;
	*bx=SF_MAX(0,tmp_nx);
    }
    if (iz <nz_right){
	z=(float)(iz)*dz; x=z/tan_rotate_angle;
	tmp_nx=(int)(x/dx)+boundary_nx;
	*ex=SF_MIN(nx_rotate-1,ix_rotate_center+tmp_nx); 
    }
    else{
	z=(float)(iz-nz_right)*dz; x=z*tan_rotate_angle;
	tmp_nx=(int)(x/dx)-boundary_nx;
	*ex=SF_MIN(nx_rotate-1,nx_rotate-tmp_nx);
    }
    nx_iz=kissfftn(*ex-*bx+1); *ex=*bx+nx_iz-1;
    if ( *ex > nx_rotate-1 ){
	//*bx=0; *ex=nx_rotate-1; nx_iz=nx_rotate;
	*ex=nx_rotate-1; *bx=*ex+1-nx_iz; 
	if (*bx<0){*ex=nx_rotate-1; *bx=0; nx_iz=nx_rotate;}
    }
    *dkx=2.0*SF_PI/((float)(nx_iz)*dx);
}


void shot_rotate_data_rotate(sf_complex * z0_s,sf_complex * z0_s_rotate,struct shot_rotate_par_type *rotate_par)
/*< rotation >*/
{
    float rotate_angle,cos_rotate_angle,sin_rotate_angle,dx,dz;
    int nx,nx_rotate,nz_s_rotate,ix_rotate_center;

    int  ix,ix_rotate,iz_rotate,i_ix_rotate,i_iz_rotate;
    float x,x_rotate,z_rotate,wx_rotate,wz_rotate;
    float  sincx[16],sincz[16];
    sf_complex *z0_s_smooth;
    int nd; float ddx;
    int ix_fin_rotate,iz_fin_rotate;

    rotate_angle=rotate_par->rotate_angle;
    dx=rotate_par->dx; dz=rotate_par->dz; nx=rotate_par->nx; 
    nx_rotate=rotate_par->nx_rotate;
    nz_s_rotate=rotate_par->nz_s_rotate;
    ix_rotate_center=rotate_par->ix_rotate_center;
    cos_rotate_angle=cos(fabs(rotate_angle)); sin_rotate_angle=sin(fabs(rotate_angle));
  
    nd=16; 
    ddx=dx/(float)(nd);
    z0_s_smooth=sf_complexalloc(nd*nx);

    vector_value_c(z0_s_smooth,sf_cmplx(0.0,0.0),nd*nx);
    for(ix=0;ix<nd*nx;ix++){
	wx_rotate=(float)(ix%nd)/(float)(nd);
	ix_rotate=ix/nd;
	//printf("ix=%d,ixrotate=%d\n",ix,ix_rotate);
    
	//if (wx_rotate < 0) wx_rotate=1+wx_rotate;
	//iwx=(int)(wx_rotate/dwx); wwx=wx_rotate/dwx-(float)iwx; i0sinc=rotate_par->sinc+8*iwx; i1sinc=i0sinc+8; 
	//for(i_ix_rotate=0;i_ix_rotate<8;i_ix_rotate++) 
	//  sincx[i_ix_rotate]=(1-wwx)*i0sinc[i_ix_rotate]+wwx*i1sinc[i_ix_rotate];
	wx_rotate=0.0; 
	mksinc(wx_rotate,16,sincx);
	//printf("ix=%d\n",ix);
	//for(i_ix_rotate=0;i_ix_rotate<16;i_ix_rotate++)
	//  printf("%f  ",sincx[i_ix_rotate]);
	//printf("\n");
	//printf("ix=%d\n",ix);
	for(i_ix_rotate=1;i_ix_rotate<=16;i_ix_rotate++){
	    if (ix_rotate-8+i_ix_rotate >=0 && ix_rotate-8+i_ix_rotate <nx ){ 
		z0_s_smooth[ix]+=z0_s[ix_rotate-8+i_ix_rotate]*sincx[i_ix_rotate-1];
		//printf("ix=%d,r=%f  ",ix_rotate-8+i_ix_rotate, sincx[i_ix_rotate-1]);
	    }
	}  //for(i_ix_rotate=0;i_ix_rotate<8;i_ix_rotate++)
	//printf("\n");
    }// for(ix=0;ix<nd*nx;ix++)
    //printf("max1=%f\n",maxvalc(z0_s_smooth,nd*nx));
    vector_value_c(z0_s_rotate,sf_cmplx(0.0,0.0),nx_rotate*nz_s_rotate);
    for(ix=0;ix<nx*nd;ix++){
	if (rotate_angle >=0 )  x=(float)(ix)*ddx;
	else x=(float)(nx*nd-1-ix)*ddx;
	x_rotate=x*cos_rotate_angle; z_rotate=x*sin_rotate_angle;
	ix_rotate=(int)((x_rotate-0.0)/dx);iz_rotate=(int)((z_rotate-0.0)/dz);
	wx_rotate=(x_rotate-(float)(ix_rotate)*dx)/dx;
	wz_rotate=(z_rotate-(float)(iz_rotate)*dz)/dz;
	//printf("ix=%d,xr=%f,zr=%f,ixr=%d,izr=%d\n",ix,wx_rotate,wz_rotate,ix_rotate,iz_rotate);
    
	//if (wx_rotate < 0) wx_rotate=1+wx_rotate; if (wz_rotate < 0) wz_rotate=1+wz_rotate;
	//iwx=(int)(wx_rotate/dwx); wwx=wx_rotate/dwx-(float)iwx; i0sinc=rotate_par->sinc+8*iwx; i1sinc=i0sinc+8; 
	//for(i_ix_rotate=0;i_ix_rotate<8;i_ix_rotate++) 
	//  sincx[i_ix_rotate]=(1-wwx)*i0sinc[i_ix_rotate]+wwx*i1sinc[i_ix_rotate];
	//iwx=(int)(wz_rotate/dwx); wwx=wz_rotate/dwx-(float)iwx; i0sinc=rotate_par->sinc+8*iwx; i1sinc=i0sinc+8;
	//for(i_ix_rotate=0;i_ix_rotate<8;i_ix_rotate++) 
	//  sincz[i_ix_rotate]=(1-wwx)*i0sinc[i_ix_rotate]+wwx*i1sinc[i_ix_rotate];   
 
	mksinc(wx_rotate,16,sincx); 
	mksinc(wz_rotate,16,sincz);
	for( i_ix_rotate=1; i_ix_rotate<=16;i_ix_rotate++){
	    ix_fin_rotate=ix_rotate-8+i_ix_rotate+ix_rotate_center; 
	    if (ix_fin_rotate >=0 && ix_fin_rotate <nx_rotate){
		for(i_iz_rotate=1;i_iz_rotate<=16;i_iz_rotate++){
		    iz_fin_rotate=iz_rotate-8+i_iz_rotate;
		    if (iz_fin_rotate >=0 && iz_fin_rotate <nz_s_rotate) 
			z0_s_rotate[i2(iz_fin_rotate,ix_fin_rotate,nx_rotate)]+=
			    z0_s_smooth[ix]*sincx[i_ix_rotate-1]*sincz[i_iz_rotate-1];
		}//for(i_iz_rotate=1;i_iz_rotate<=8;i_iz_rotate++)
	    }
	}//for( i_ix_rotate=1; i_ix_rotate<=8;i_ix_rotate++)

    }//for(ix=0;ix<nx*nd;ix++)
    //printf("max1=%f\n",maxvalc(z0_s_rotate,rotate_par->nz_s_rotate*rotate_par->nx_rotate));

  
    free(z0_s_smooth);
} // shot_rotate_data_rotate

void shot_rotate_image_rotate(float *image,float *image_rotate,struct shot_rotate_par_type rotate_par)
/*< rotate >*/
{
    float rotate_angle,cos_rotate_angle,sin_rotate_angle,dx,dz;
    int ix_rotate_center, nx, nz, nx_rotate, nz_rotate;

    int ix,iz,ix_rotate,iz_rotate,i_ix_rotate,i_iz_rotate,ii_ix_rotate,ii_iz_rotate;
    float x,z,x_rotate,z_rotate,wx_rotate,wz_rotate;
    float sincx[16],sincz[16];
    float *i0sinc; 
    int iwx; 
    float dwx;

    rotate_angle=rotate_par.rotate_angle; dx=rotate_par.dx; dz=rotate_par.dz;
    ix_rotate_center=rotate_par.ix_rotate_center; nx=rotate_par.nx; nz=rotate_par.nz;
    nx_rotate=rotate_par.nx_rotate; nz_rotate=rotate_par.nz_rotate;
    cos_rotate_angle=cos(fabs(rotate_angle)); sin_rotate_angle=sin(fabs(rotate_angle));
    vector_value_f(image_rotate,0.0,nx*nz);

    dwx=1.0/(float)(rotate_par.Nsinc);
    for(iz=0;iz<nz;iz++){
	z=(float)(iz)*dz;
	for(ix=0;ix<nx;ix++){
	    if (rotate_angle >=0 ) x=ix*dx; 
	    else x=(nx-1-ix)*dx;
	    x_rotate=x*cos_rotate_angle-z*sin_rotate_angle;   z_rotate=x*sin_rotate_angle+z*cos_rotate_angle;
	    ix_rotate=(int)((x_rotate-0.0)/dx);               iz_rotate=(int)((z_rotate-0.0)/dz);
	    wx_rotate=(x_rotate-(float)(ix_rotate)*dx)/dx;    wz_rotate=(z_rotate-(float)(iz_rotate)*dz)/dz; 
	    /* 
	       if (wx_rotate < 0) wx_rotate=1+wx_rotate; if (wz_rotate < 0) wz_rotate=1+wz_rotate;
	       iwx=(int)(wx_rotate/dwx); wwx=wx_rotate/dwx-(float)iwx; i0sinc=rotate_par.sinc+8*iwx; i1sinc=i0sinc+8; 
	       for(ii_ix_rotate=0;ii_ix_rotate<8;ii_ix_rotate++) 
	       sincx[ii_ix_rotate]=(1-wwx)*i0sinc[ii_ix_rotate]+wwx*i1sinc[ii_ix_rotate];
	       iwx=(int)(wz_rotate/dwx); wwx=wz_rotate/dwx-(float)iwx; i0sinc=rotate_par.sinc+8*iwx; i1sinc=i0sinc+8;
	       for(ii_ix_rotate=0;ii_ix_rotate<8;ii_ix_rotate++) 
	       sincz[ii_ix_rotate]=(1-wwx)*i0sinc[ii_ix_rotate]+wwx*i1sinc[ii_ix_rotate];
	    */
	    if (wx_rotate < 0) wx_rotate=1+wx_rotate; 
        if (wz_rotate < 0) wz_rotate=1+wz_rotate;
	    iwx=(int)(wx_rotate/dwx); i0sinc=rotate_par.sinc+8*iwx;
	    for(ii_ix_rotate=0;ii_ix_rotate<8;ii_ix_rotate++) sincx[ii_ix_rotate]=i0sinc[ii_ix_rotate];
	    iwx=(int)(wz_rotate/dwx); i0sinc=rotate_par.sinc+8*iwx;
	    for(ii_ix_rotate=0;ii_ix_rotate<8;ii_ix_rotate++) sincz[ii_ix_rotate]=i0sinc[ii_ix_rotate];

	    for(ii_ix_rotate=1;ii_ix_rotate<=8; ii_ix_rotate++){
		i_ix_rotate=ix_rotate-4+ii_ix_rotate+ix_rotate_center;
		if (i_ix_rotate >=0 && i_ix_rotate <nx_rotate){ 
		    for(ii_iz_rotate=1;ii_iz_rotate<=8;ii_iz_rotate++){
			i_iz_rotate=iz_rotate-4+ii_iz_rotate;
			if (i_iz_rotate >=0 && i_iz_rotate <nz_rotate) 
			    image_rotate[i2(iz,ix,nx)]+=
				image[i2(i_iz_rotate,i_ix_rotate,nx_rotate)]*sincx[ii_ix_rotate-1]*sincz[ii_iz_rotate-1];
		    }// for(ii_iz_rotate)
		}
	    } // for(ii_ix_rotate)


	    //sinc_mksinc(wx_rotate,16,sincx);             sinc_mksinc(wz_rotate,16,sincz);
	    //for(ii_ix_rotate=1;ii_ix_rotate<=16; ii_ix_rotate++){
	    //   i_ix_rotate=ix_rotate-8+ii_ix_rotate+ix_rotate_center;
	    //  if (i_ix_rotate >=0 && i_ix_rotate <nx_rotate){ 
	    //    for(ii_iz_rotate=1;ii_iz_rotate<=16;ii_iz_rotate++){
	    //      i_iz_rotate=iz_rotate-8+ii_iz_rotate;
	    //      if (i_iz_rotate >=0 && i_iz_rotate <nz_rotate) 
	    //        image_rotate[i2(iz,ix,nx)]+=
	    //         image[i2(i_iz_rotate,i_ix_rotate,nx_rotate)]*sincx[ii_ix_rotate-1]*sincz[ii_iz_rotate-1];
	    //    }// for(ii_iz_rotate)
	    //  }
	    //} // for(ii_ix_rotate)  


	}//for(ix)
    } //for(iz)

}




void shot_rotate_data_tilt_2d(sf_complex *wldin, sf_complex **wldout,struct shot_rotate_par_type rotate_par){
    int nxin,nxout,nzout,nd;
    float *sincx,*sincz;
    sf_complex *wldsmooth;
    int *bx_s_rotate;
    float oxout,costilt,sintilt,tilt,dx,dz,dwx,ddx;
    int ix,ix1,iz1,iwx,iwz,jz,jx,ixout,izout,bxout,exout,nx;
    float x,x1,z1,wx,wz;
    nxin=rotate_par.nxout_azimuth;
    nxout=rotate_par.nx_s_rotate;   nzout=rotate_par.nz_s_rotate; oxout=rotate_par.oxout_tilt;
    costilt=rotate_par.costilt; sintilt=rotate_par.sintilt; tilt=rotate_par.rotate_angle;
    dx=rotate_par.dx; dz=rotate_par.dz;
    bx_s_rotate=rotate_par.bx_s_rotate;
    dwx=1.0/(float)(rotate_par.Nsinc);
    nx=nxin;

    nd=16;
    ddx=dx/(float)(nd);
    wldsmooth=sf_complexalloc(nd*nx);
    vector_value_c(wldsmooth,sf_cmplx(0.0,0.0),nd*nx);
    vector_value_c(*wldout,sf_cmplx(0.0,0.0),nxout*nzout);
    for(ix=0;ix<=nd*(nx-1);ix++){
	wx=(float)(ix%nd)/(float)(nd);
	ix1=ix/nd;
	wx=0.0;
	iwx=(int)(wx/dwx); 
	sincx=rotate_par.sinc+8*iwx;
	for(jx=1;jx<=8;jx++){
	    if(ix1-4+jx>=0 && ix1-4+jx <nx){
		wldsmooth[ix]+=wldin[ix1-4+jx]*sincx[jx-1];
	    }
	}
    }

    for(ix=0;ix<=(nx-1)*nd;ix++){
	if (tilt >=0) x=(float)(ix)*ddx;
	else  x=(float)((nx-1)*nd-ix)*ddx;
	x1=x*costilt; z1=x*sintilt;
	ix1=(int)((x1-oxout)/dx); iz1=(int)((z1-0.0)/dz);
	wx=(x1-((float)(ix1)*dx+oxout))/dx;
	wz=(z1-(float)(iz1)*dz)/dz;
	iwx=(int)(wx/dwx); iwz=(int)(wz/dwx);
	sincx=rotate_par.sinc+8*iwx; sincz=rotate_par.sinc+8*iwz;
 
	for( jz=1; jz<=8;jz++){
	    izout=iz1-4+jz; 
	    if (izout >=0 && izout <nzout){
		bxout=bx_s_rotate[izout]; exout=bxout+rotate_par.nx_s_rotate-1;
		for(jx=1;jx<=8;jx++){
		    ixout=ix1-4+jx;
		    if (ixout >=bxout && ixout <=exout) {
			wldout[izout][ixout-bxout]+=wldsmooth[ix]*sincz[jz-1]*sincx[jx-1];
		    }
		    else{
			if (0.0!=rotate_par.rotate_angle){
			    printf("report bug in x z rotation\n");
			    printf("ixout=%d,bxout=%d,exout=%d\n",ixout,bxout,exout);
			}
		    }
		}
	    }
	}
    }

    //printf("outsub=%f\n",maxvalc(*wldout,rotate_par.nz_s_rotate*rotate_par.nx_s_rotate )); 

/*
  for(ix=0;ix<nx*nd;ix++){
  if (tilt >=0) x=(float)(ix)*ddx;
  else  x=(float)(nx*nd-1-ix)*ddx;
  x1=x*costilt-dz*sintilt; z1=x*sintilt+dz*costilt;
  ix1=(int)((x1-oxout)/dx); iz1=(int)((z1-0.0)/dz);
  wx=(x1-((float)(ix1)*dx+oxout))/dx;
  wz=(z1-(float)(iz1)*dz)/dz;
  iwx=(int)(wx/dwx); iwz=(int)(wz/dwx);
  sincx=rotate_par.sinc+8*iwx; sincz=rotate_par.sinc+8*iwz; 
  for( jz=1; jz<=8;jz++){
  izout=iz1-4+jz; 
  if (izout >=0 && izout <nzout){
  bxout=bx_s_rotate[izout]; exout=bxout+rotate_par.nx_s_rotate-1;
  for(jx=1;jx<=8;jx++){
  ixout=ix1-4+jx;
  if (ixout >=bxout && ixout <=exout) {
  wldout[izout][ixout-bxout]+=0.5*wldsmooth[ix]*sincz[jz-1]*sincx[jx-1];
  }
  else{
  printf("report bug in x z rotation\n");
  }
  }
  }
  }
  }

  for(ix=16;ix<nx*nd;ix++){
  if (tilt >=0) x=(float)(ix)*ddx;
  else  x=(float)(nx*nd-1-ix)*ddx;
  x1=x*costilt+dz*sintilt; z1=x*sintilt-dz*costilt;
  ix1=(int)((x1-oxout)/dx); iz1=(int)((z1-0.0)/dz);
  wx=(x1-((float)(ix1)*dx+oxout))/dx;
  wz=(z1-(float)(iz1)*dz)/dz;
  iwx=(int)(wx/dwx); iwz=(int)(wz/dwx);
  sincx=rotate_par.sinc+8*iwx; sincz=rotate_par.sinc+8*iwz; 
  for( jz=1; jz<=8;jz++){
  izout=iz1-4+jz; 
  if (izout >=0 && izout <nzout){
  bxout=bx_s_rotate[izout]; exout=bxout+rotate_par.nx_s_rotate-1;
  for(jx=1;jx<=8;jx++){
  ixout=ix1-4+jx;
  if (ixout >=bxout && ixout <=exout) {
  wldout[izout][ixout-bxout]+=0.5*wldsmooth[ix]*sincz[jz-1]*sincx[jx-1];
  }
  else{
  printf("report bug in x z rotation\n");
  }
  }
  }
  }
  }
*/

    free(wldsmooth); 

}

void shot_rotate_velocity_tilt_2d(float **velin2d,float **velout2d,sf_file velocityin,
				  struct shot_rotate_par_type rotate_par)
{
    int nxin,nzin,nxout,nzout;
    float dzout,dzin,dxout,dxin,oxout,oxin,ozin;
    int ixin,izin,ixout,izout;
    float zout,xout,xin,zin,wx,wz,w00,w01,w10,w11;
    float costilt,sintilt;

    costilt=rotate_par.costilt; sintilt=rotate_par.sintilt;
    sf_histint(velocityin,"n2",&nxin); sf_histfloat(velocityin,"d2",&dxin); sf_histfloat(velocityin,"o2",&oxin);
    sf_histint(velocityin,"n1",&nzin); sf_histfloat(velocityin,"d1",&dzin); sf_histfloat(velocityin,"o1",&ozin);
  
    nxout=rotate_par.nx_rotate; dxout=rotate_par.dx; oxout=rotate_par.oxout_tilt;
    nzout=rotate_par.nz_rotate; dzout=rotate_par.dz; 
    for(izout=0;izout<nzout;izout++){
	zout=(float)(izout)*dzout;
	for(ixout=0;ixout<nxout;ixout++){
	    xout=(float)(ixout)*dxout+oxout;
	    if (rotate_par.rotate_angle >0 ){
		xin= costilt*xout+sintilt*zout+rotate_par.ox;
	    }
	    else{
		xin=(rotate_par.ox+(rotate_par.nx-1)*rotate_par.dx)-(costilt*xout+sintilt*zout);  
	    }
	    zin=-sintilt*xout+costilt*zout;
	    ixin=(int)((xin-oxin)/dxin); izin=(int)((zin-ozin)/dzin);
	    wx=(xin-(oxin+(float)(ixin)*dxin))/dxin;  wz=(zin-(ozin+(float)(izin)*dzin))/dzin;
	    if (ixin <0){ ixin=0; wx=0.5;}  if (ixin>=nxin-1){ixin=nxin-2; wx=0.5;}
	    if (izin <0){ izin=0; wz=0.5;}  if (izin>=nzin-1){izin=nzin-2; wz=0.5;}
	    w00=(1.0-wz)*(1.0-wx);  w01=(1.0-wz)*wx; w10=wz*(1.0-wx); w11=wz*wx;
	    velout2d[izout][ixout]=velin2d[izin][ixin]*w00+velin2d[izin][ixin+1]*w01+velin2d[izin+1][ixin]*w10+velin2d[izin+1][ixin+1]*w11;
	}
    }  

}

void shot_rotate_image_tilt_2d(float **image,float **image_rotate,struct shot_rotate_par_type rotate_par)
{
    float rotate_angle,cos_rotate_angle,sin_rotate_angle,dx,dz;
    int nx, nz, nx_rotate, nz_rotate;
    int ix,iz,ix_rotate,iz_rotate,ix_rotate_shift,iz_rotate_shift,ix_shift,iz_shift;
    float x,z,x_rotate,z_rotate,wx_rotate,wz_rotate;
    float *sincx,*sincz;
    int iwx,iwz; 
    float dwx;

    rotate_angle=rotate_par.rotate_angle; dx=rotate_par.dx; dz=rotate_par.dz;
    nx=rotate_par.nx; nz=rotate_par.nz;
    nx_rotate=rotate_par.nx_rotate; nz_rotate=rotate_par.nz_rotate;
    cos_rotate_angle=rotate_par.costilt; sin_rotate_angle=rotate_par.sintilt;
    vector_value_f(*image,0.0,nx*nz);

    dwx=1.0/(float)(rotate_par.Nsinc);
    for(iz=0;iz<nz;iz++){
	z=(float)(iz)*dz;
	for(ix=0;ix<nx;ix++){
	    if (rotate_angle >=0 ) x=(float)(ix)*dx; else x=(float)(nx-1-ix)*dx;
	    x_rotate=x*cos_rotate_angle-z*sin_rotate_angle;   z_rotate=x*sin_rotate_angle+z*cos_rotate_angle;
	    ix_rotate=(int)((x_rotate-rotate_par.oxout_tilt)/dx);               iz_rotate=(int)((z_rotate-0.0)/dz);
	    wx_rotate=(x_rotate-((float)(ix_rotate)*dx+rotate_par.oxout_tilt))/dx;    wz_rotate=(z_rotate-(float)(iz_rotate)*dz)/dz;
	    if (ix_rotate==0 && wx_rotate <0) wx_rotate=1+wx_rotate;
	    iwx=(int)(wx_rotate/dwx); iwz=(int)(wz_rotate/dwx);
	    //printf("iwx=%d,iwz=%d,wx_rotate=%f,wz_rotate=%f,dwx=%f,ix_rotate=%d\n",iwx,iwz,wx_rotate,wz_rotate,dwx,ix_rotate); 
	    sincx=rotate_par.sinc+8*iwx; sincz=rotate_par.sinc+8*iwz;
      
	    for(ix_shift=1;ix_shift<=8; ix_shift++){
		ix_rotate_shift=ix_rotate-4+ix_shift;
		if (ix_rotate_shift >=0 && ix_rotate_shift <nx_rotate){ 
		    for(iz_shift=1;iz_shift<=8;iz_shift++){
			iz_rotate_shift=iz_rotate-4+iz_shift;
			if (iz_rotate_shift >=0 && iz_rotate_shift <nz_rotate){ 
			    image[iz][ix]+=image_rotate[iz_rotate_shift][ix_rotate_shift]*sincx[ix_shift-1]*sincz[iz_shift-1];
			}
		    }// for(ii_iz_rotate)
		}
	    } // for(ii_ix_rotate)  
	}//for(ix)
    } //for(iz)
}




