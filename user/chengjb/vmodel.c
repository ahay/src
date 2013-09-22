/*************************************************************************
 * sub-programs to prepare velocity models
 *
 * *************************************************************************/
/*
   Copyright (C) 2012 Tongji University, Shanghai, China 
   Authors: Jiubing Cheng, Wei Kang and Tengfei Wang
     
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
#include "_cjb.h"

void vmodelboundary3d(float ***v, int nx, int ny, int nz, int nxpad, int nypad, int nzpad, int bd)
/*< filling padded boundaries of 3D model: nz, nx, ny are fastest, mid, and slowest dimension >*/
{
     int i, j, k;

     float vmin=9999999999.0, vmax=-9999999999.0;

     for(k=bd;k<nzpad-bd;k++){
       for(i=bd;i<nxpad-bd;i++){
          for(j=0;j<bd;j++)           v[j][i][k] = v[bd][i][k];
          for(j=nypad-bd;j<nypad;j++) v[j][i][k] = v[nypad-bd-1][i][k];
       }
       for(j=0;j<nypad;j++){
          for(i=0;i<bd;i++)           v[j][i][k] = v[j][bd][k];
          for(i=nxpad-bd;i<nxpad;i++) v[j][i][k] = v[j][nxpad-bd-1][k];
       }
     }
     for(i=0;i<nxpad;i++)
     for(j=0;j<nypad;j++){
       for(k=0;k<bd;k++)              v[j][i][k] = v[j][i][bd];
       for(k=nzpad-bd;k<nzpad;k++)    v[j][i][k] = v[j][i][nzpad-bd-1];
       for(k=0;k<nzpad;k++){
          if(v[j][i][k]>vmax) vmax=v[j][i][k];
          if(v[j][i][k]<vmin) vmin=v[j][i][k];
       }
     }
     sf_warning("model parameter minimum and maximum value: %f %f\n",vmin,vmax);
}

void vmodelboundary2d(float **v, int nx, int nz, int nxpad, int nzpad, int bd)
/*< filling padded boundaries of 2D model: nz is fastest dimension >*/
{
     int i, k;

     float vmin=9999999999.0, vmax=-9999999999.0;

     for(k=bd;k<nzpad-bd;k++){
       for(i=bd;i<nxpad-bd;i++)    v[i][k] = v[bd][k];
       for(i=nxpad-bd;i<nxpad;i++) v[i][k] = v[nxpad-bd-1][k];
     }
     for(i=0;i<nxpad;i++){
       for(k=0;k<bd;k++)              v[i][k] = v[i][bd];
       for(k=nzpad-bd;k<nzpad;k++)    v[i][k] = v[i][nzpad-bd-1];
       for(k=0;k<nzpad;k++){
          if(v[i][k]>vmax) vmax=v[i][k];
          if(v[i][k]<vmin) vmin=v[i][k];
       }
     }
     sf_warning("model parameter minimum and maximum value: %f %f\n",vmin,vmax);
}
