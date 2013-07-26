/*  Source load */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

typedef struct {
    int srange;
    float alpha;
    int decay;
    float trunc;
}spara; /*source parameter*/
/*^*/

void explsourcet(float **vtxx/*@out@*/,
		 float *vwavlet, 
		 int vit, float vdt,
                 int vsx, int vsdepth, 
		 int vnx, int vnz,
		 spara *vps/*decay parameters*/)
/*<explosive source for stress txx>*/ 
{
    float phi = 0.0;
    int cent = vps->srange/2;
    int ix, iz;

    if (vps->decay ==1){
	for (ix=0; ix<2*cent; ix++)
	    for (iz=0; iz<2*cent; iz++) {
		phi = exp( -1*vps->alpha*vps->alpha*((ix-cent)*(ix-cent)+(iz-cent)*(iz-cent)) );
		vtxx[vsx-cent+ix][vsdepth-cent+iz] += vwavlet[vit]*phi*vdt;
	    }
    } else {
	vtxx[vsx][vsdepth] += vwavlet[vit]*vdt;
    } 
}

void explsourcet1(float *vtxx/*@out@*/,
		  float *vwavlet, float dt,
		 int it, int vsdepth, int vnt, 
		 spara *vps/*decay parameters*/)
/*<1D explosive source for stress txx>*/ 
{
    float phi = 0.0;
    int cent = vps->srange/2;
    int ix;

    if (vps->decay ==1 && it <vnt-1){
	for (ix=0; ix<2*cent; ix++){
	    phi = exp(-1*vps->alpha*vps->alpha*(ix-cent)*(ix-cent));
	    vtxx[vsdepth-cent+ix] += 0.5*(vwavlet[it]+vwavlet[it+1])*phi*dt;
	}
    } else {
	vtxx[vsdepth] += 0.5*(vwavlet[it]+vwavlet[it+1])*dt;
    } 
}

void distrsrc(float **vtxx, float **dissrc, int nx, int nz, int shiftx, int shiftz)
/*<2D distribute source>*/
{
    int ix, iz;
    for (ix=0; ix<nx; ix++) {
	for (iz=0; iz<nz; iz++) {
	    vtxx[ix+shiftx][iz+shiftz] += dissrc[ix][iz];
	}
    }
}



