/* calculate gradient by adjoint state method*/
/*
  Copyright (C) 2019 
  coded by carrot @ China University of Petroleum 
  (Email:1017014258@qq.com)
  
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
#include <time.h>
#include "acoustic1D.h"

float* adjoint_gradient(float *vel_ini,float *src,float dx,int nx,float dt, int nt, int sx, int rx,float *rec,int flag)
/*<calculate gradient by adjoint state method>*/
{


  float *rec_ini=sf_floatalloc(nt);float **rec_all_ini=sf_floatalloc2(nt,nx);
     // sf_warning("nx=%d, dx=%g, nt=%d, dt=%g, sx=%d, rx=%d",nx,dx,nt,dt,sx,rx);

   // sf_warning("vel_ini[0]=%f",vel_ini[0]);
  //  sf_warning("src[1000]=%f",src[1000]);
     //   sf_warning("rec[1499]=%f",rec[1499]);


  acoustic1D(vel_ini,src,dx,nx,dt,nt,sx,rx,rec_ini, rec_all_ini);
  float *rec_res=sf_floatalloc(nt);for (int it = 0; it < nt; ++it) rec_res[nt-1-it]=rec[it]-rec_ini[it];//reverse source


  float *rec_adjoint=sf_floatalloc(nt);float **rec_all_adjoint=sf_floatalloc2(nt,nx);
  acoustic1D(vel_ini,rec_res,dx,nx,dt,nt,sx,rx,rec_adjoint, rec_all_adjoint);

  float *gradient=sf_floatalloc(nx);
  for (int ix = 0; ix < nx; ++ix)
  {
    gradient[ix]=0.0;
    for (int it = 1; it < nt-1; ++it) gradient[ix]=gradient[ix]+(rec_all_ini[ix][it+1]+rec_all_ini[ix][it-1]-2*rec_all_ini[ix][it])/dt/dt*rec_all_adjoint[ix][nt-1-it]; 
    gradient[ix]=2*gradient[ix]/(vel_ini[ix]*vel_ini[ix]*vel_ini[ix]);
  }


  float error_sum=0.0;for (int it = 0; it < nt; ++it)  error_sum=error_sum+(rec_ini[it]-rec[it])*(rec_ini[it]-rec[it]); 
  if (flag==1) sf_warning("error_sum=%f\n",error_sum);



  free(rec_ini);
  free(rec_res);
  free(rec_adjoint);
  return gradient;





}