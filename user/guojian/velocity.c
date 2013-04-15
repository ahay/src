/* Reading velocity information */

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

void velocity_init(sf_file *velocity)
/*< open velocity file >*/
{
   *velocity=sf_input("Velocity");
}

void velocity_ani_init(sf_file *eps,sf_file *dlt)
/*< open anisotropy parameter file ep and dlt if neccesaray >*/
{
   *eps=sf_input("Eps");
   *dlt=sf_input("Dlt");
}

void velocity_read(sf_file velocity,float *vel,sf_file wave,float theta,int ix_center,char *tag)
/*< read velocity or anisotropy parameter in tilted coordinates >*/
{
  float *tmpvel;
  int nvel[2],ntmpvel[2];
  float mig_dx,mig_dz,mig_nx,mig_nz,mig_ox,mig_oz;
  float vel_dx,vel_dz,vel_nx,vel_nz,vel_ox,vel_oz;
  int iy,ix,iz,n_block,blocksize; 
  float tmp_x,tmp_z,tmp_wx,tmp_wz,dtmp_x,dtmp_z;
  int tmp_ix,tmp_iz;
 
  mig_dx=wave.d[1]; mig_dz=wave.d[0]; mig_nx=wave.n[1]; mig_nz=wave.n[0]; mig_ox=wave.o[1]; mig_oz=wave.o[0];
  vel_dx=velocity.d[1]; vel_dz=velocity.d[0]; vel_nx=velocity.n[1]; vel_nz=velocity.n[0];
  vel_ox=velocity.o[1]; vel_oz=velocity.o[0];

  d3(mig_nx,mig_nz,nvel);
  d3(vel_nx,vel_nz,ntmpvel);

  tmpvel=allocatef(vel_nz*vel_nx);
  if (theta >=0)
  {
    for(iy=0;iy<1;iy++){    //do iy=1,ny
      n_block=iy*vel_nx;
      if (sseek_block(tag,n_block,vel_nz*4,0)!= n_block) seperr("velocity seek wrong");
      blocksize=vel_nz*vel_nx*4;
      if(sreed(tag,tmpvel,blocksize)!=blocksize) seperr("velocity read wrong");

      for(ix=0;ix<mig_nx;ix++){
        dtmp_x=ix*mig_dx-(ix_center-1.0)*mig_dx;
        for(iz=0;iz<mig_nz;iz++){
          dtmp_z=iz*mig_dz;

          tmp_x=mig_ox+dtmp_x*cos(theta)+dtmp_z*sin(theta);
          tmp_x=(tmp_x-vel_ox)/vel_dx;
          tmp_ix=(int)(tmp_x);  tmp_wx=tmp_x-tmp_ix;
          if (tmp_x <0){
            tmp_ix=0;  tmp_wx=0.5;
          }
          if (tmp_ix >= vel_nx-1 ){
            tmp_ix=vel_nx-2;  tmp_wx=0.5;
          }

          tmp_z=mig_oz-dtmp_x*sin(theta)+dtmp_z*cos(theta);
          tmp_z=(tmp_z- vel_oz)/vel_dz;
          tmp_iz=(int)(tmp_z); tmp_wz=tmp_z-tmp_iz;
          if (tmp_z < 0){
            tmp_iz=0;   tmp_wz=0.5;
          }
          if (tmp_iz >= vel_nz-1 ){ 
            tmp_iz=vel_nz-2; tmp_wz=0.5;
          }
         
          vel[i3(iy,ix,iz,nvel)]=tmpvel[i3(0,tmp_ix,tmp_iz,ntmpvel)]*(1-tmp_wz)*(1-tmp_wx)+ 
                                 tmpvel[i3(0,tmp_ix,tmp_iz+1,ntmpvel)]*tmp_wz*(1-tmp_wx)  + 
                                 tmpvel[i3(0,tmp_ix+1,tmp_iz,ntmpvel)]*(1-tmp_wz)*tmp_wx  + 
                                 tmpvel[i3(0,tmp_ix+1,tmp_iz+1,ntmpvel)]*tmp_wz*tmp_wx;    
        }
      }
    }
  }
  else
  
  {
    theta=-theta;
    for(iy=0;iy<1;iy++){
      n_block=iy*vel_nx;
      if (sseek_block(tag,n_block,vel_nz*4,0)!= n_block) seperr("velocity seek wrong");
      blocksize=vel_nz*vel_nx*4;
      if(sreed(tag,tmpvel,blocksize)!=blocksize) seperr("velocity read wrong");
        
      for(ix=mig_nx-1;ix>=0;ix--){
        dtmp_x=(mig_nx-1-ix)*mig_dx-(ix_center-1)*mig_dx;
        for(iz=0;iz<mig_nz;iz++){
          dtmp_z=iz*mig_dz;

          tmp_x=mig_ox-dtmp_x*cos(theta)-dtmp_z*sin(theta);
          tmp_x=(tmp_x-vel_ox)/vel_dx;
          tmp_ix=(int)(tmp_x);  tmp_wx=tmp_x-tmp_ix;
          if (tmp_x <0){ 
            tmp_ix=0;  tmp_wx=0.5;
          }
          if (tmp_ix >= vel_nx-1 ){
            tmp_ix=vel_nx-2;  tmp_wx=0.5;
          }
          tmp_z=mig_oz-dtmp_x*sin(theta)+dtmp_z*cos(theta);
          tmp_z=(tmp_z-vel_oz)/vel_dz;
          tmp_iz=(int)(tmp_z); tmp_wz=tmp_z-tmp_iz;
          if (tmp_z < 0){
            tmp_iz=0;   tmp_wz=0.5;
          }
          if (tmp_iz >= vel_nz-1 ){
            tmp_iz=vel_nz-2; tmp_wz=0.5;
          }
          vel[i3(iy,(mig_nx-1)-ix,iz,nvel)]=tmpvel[i3(0,tmp_ix,tmp_iz,ntmpvel)]*(1-tmp_wz)*(1-tmp_wx)+ 
                                      tmpvel[i3(0,tmp_ix,tmp_iz+1,ntmpvel)]*tmp_wz*(1-tmp_wx)  + 
                                      tmpvel[i3(0,tmp_ix+1,tmp_iz,ntmpvel)]*(1-tmp_wz)*tmp_wx  + 
                                      tmpvel[i3(0,tmp_ix+1,tmp_iz+1,ntmpvel)]*tmp_wz*tmp_wx;
        }
      }        
    }  
 
  }
 
  free(tmpvel);
}

/* find max and min value of velocity or anisotropy parameters */
void velocity_minmax(char *tag,sf_file velocity,float *minvel,float *maxvel){
  int n,idim;
  float *vel;
  int ndim=velocity.ndims;
  for(n=1,idim=0;idim<ndim; idim++){
    n*=velocity.n[idim];
  }
  vel=allocatef(n);
  if(sreed(tag,vel,n*4)!=n*4) seperr("velocity read wrong");
  *maxvel=maxval(vel,n);
  *minvel=minval(vel,n);
  free(vel);
}

