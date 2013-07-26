/*  source for sglfd RTM */
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
#include "srcm.h"

#ifndef _srcm_h

#define SRCRANGE 10
#define SRCALPHA 0.5
#define SRCTRUNC 100
#define SRCDECAY false
/*^*/

typedef struct Srcpar{
    float *wavelet;
    int nt;
    float dt;
    int range;
    float alpha;
    bool decay;
    float trunc;
} *srcpar; /*source parameter*/
/*^*/

#endif

srcpar createsrc(void)
/*<Create source struct>*/
{
    srcpar srcp;

    srcp = (srcpar)sf_alloc(1, sizeof(*srcp));
    srcp->wavelet = NULL;
    srcp->nt = 0;
    srcp->dt = 0.0;
    srcp->range = SRCRANGE;
    srcp->alpha = SRCALPHA;
    srcp->decay = SRCDECAY;
    srcp->trunc = SRCTRUNC;
    return srcp;
}

void loadsrc(srcpar srcp, sf_file srcfile)
/*<allocate source wavelet>*/
{
    if (srcp->nt == 0) sf_error("Need nt in srcpar!");
    if (srcp->wavelet != NULL) sf_error("source has been loaded!");
    srcp->wavelet = sf_floatalloc(srcp->nt);
    sf_floatread(srcp->wavelet, srcp->nt, srcfile);
}

void freesrc(srcpar srcp)
/*<Free allocated storage>*/
{
    free(srcp->wavelet);
    free(srcp);
}


void explsourcet(float **vtxx/*@out@*/,
		 float *vwavlet, 
		 int vit, float vdt,
                 int vsx, int vsdepth, 
		 int vnx, int vnz,
		 srcpar vps/*decay parameters*/)
/*<explosive source for stress txx>*/ 
{
    float phi = 0.0;
    int cent = vps->range/2;
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
		 srcpar vps/*decay parameters*/)
/*<1D explosive source for stress txx>*/ 
{
    float phi = 0.0;
    int cent = vps->range/2;
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



