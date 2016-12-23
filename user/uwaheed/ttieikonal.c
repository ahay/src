/* Fast sweeping TTI eikonal solver (2D) */
/*
  Copyright (C) 2016 Princeton University
  
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


#include <math.h>
#include <rsf.h>
/*^*/

#include "ttieikonal.h"


bool sf_init_fast_sweep (float *t,
                         int n2, int n1,
                         float o2, float o1,
                         float d2, float d1,
                         int shoty, int shotz)
/*< initialize >*/
{
    int i, n12;

    if (NULL == t)
        return false;

    if (shoty < 0 || shoty >= n2 ||
        shotz < 0 || shotz >= n1)
        return false;

    n12 = n1 * n2;

    for (i = 0; i < n12; i++)
        t[i] = SF_HUGE;

    t[shoty * n1 + shotz] = 0.0;

    return true;
}

bool isCausalRoot(float root, float uz, float uy, 
	              int sz, int sy,
		          float a, float b, float c,
                  float dz, float dy){

	float dtdz, dtdy, vgz, vgy;
  
	dtdz = (root-uz)/(sz*dz);
	dtdy = (root-uy)/(sy*dy);
	vgz = b*dtdz + c*dtdy;
	vgy = a*dtdy + c*dtdz;


	if(vgz*fabs(dtdz)*sz>0.0 && vgy*fabs(dtdy)*sy>0.0) return true;
	else return false;
}


static void sf_fast_sweep_2d_stencil (float *t, float *vz, float *vx,
                                      float *theta, float *rhs,
                                      int i, int j, int n1, int n2,
                                      float d1, float d2,
                                      int shotz, int shoty) {


    float vzc = vz[j * n1 + i], vxc = vx[j * n1 + i], thetac = theta[j*n1+i], rhsc = rhs[j*n1+i];
    float tij = t[j * n1 + i];
    int sy, sz;
    float dz = d1, dz2 = dz*dz, dx = d2, dx2 = dx*dx;
    float vzc2 = vzc*vzc, vxc2 = vxc*vxc;
    float ct = cos(thetac), st = sin(thetac);
    float ct2 = ct*ct, st2 = st*st;    
    float a = vxc2*ct2 + vzc2*st2, b = vxc2*st2 + vzc2*ct2, c = (vxc2-vzc2)*st*ct; 
    float tymin, tzmin, ty, tz;
    float ttemp = tij, num, den, root;



    if(i==shotz && j==shoty) return;

    if(i==0 || i==n1-1){
        if(i==0){
            tzmin = t[j * n1 + 1]; sz = -1; 
        }
        else{
            tzmin = t[j * n1 + n1-2]; sz = 1; 
        }
    }
    else{
        if(t[j*n1+i-1] < t[j*n1+i+1]){
            tzmin = t[j * n1 + i - 1]; sz = 1;
        }
		else{
            tzmin = t[j * n1 + i + 1]; sz = -1; 
        }

    }

    if(j==0 || j==n2-1){
        if(j==0){
            tymin = t[n1 + i]; sy = -1; 
        }
        else{
            tymin = t[(n2-2)*n1 + i]; sy = 1; 
        }
    }
    else{
        if( t[(j-1)*n1+i] < t[(j+1)*n1+i] ){
            tymin = t[(j-1)*n1+i]; sy = 1;
        }
		else{
            tymin = t[(j+1) * n1 + i ]; sy = -1; 
        }

    }

    if (tymin == SF_HUGE && tzmin == SF_HUGE)
        return;

    ty = tymin + dx*sqrt(b*rhsc/(a*b-c*c)); 
    tz = tzmin + dz*sqrt(a*rhsc/(a*b-c*c));  




    if (tymin == SF_HUGE) {
        ttemp = tz;
    } else if (tzmin == SF_HUGE) {
        ttemp = ty;
    } else {


        num = a*dz2*tymin + dx*sz*(b*dx*sz*tzmin + c*dz*sy*(tymin+tzmin)) + dx*dz*sqrt(a*(dz2*rhsc 
              - b*(tymin-tzmin)*(tymin-tzmin)) + sz*(2*c*dx*dz*rhsc*sy + b*dx2*rhsc*sz + c*c*sz*(tymin-tzmin)
              *(tymin-tzmin)));

        den = a*dz2 + dx*sz*(2*c*dz*sy + b*dx*sz);


        root = num/den;


	    if(isCausalRoot(root,tzmin,tymin,sz,sy,a,b,c,dz,dx)) ttemp = root;
        else ttemp = ty<tz?ty:tz;

    }

    if (ttemp <= tij)
        t[j * n1 + i] = ttemp;

}



void sf_run_fast_sweep (float *t, float *vz, float *vx,
                        float *theta, float *rhs,
                        int niter,
                        int n2, int n1,
                        float o2, float o1,
                        float d2, float d1,
                        int shoty, int shotz)
/*< run sweeps >*/
{
    int i, j ,k = 0, l = 0;

    for (l = 0; l < niter; l++) {

        for (j = n2 - 1; j >= 0; j--) {
            for (i = 0; i < n1; i++) {
                sf_fast_sweep_2d_stencil (t, vz, vx, theta, rhs,
                                          i, j, n1, n2, d1, d2, shotz, shoty);
            }
        }

        for (j = n2 - 1; j >= 0; j--) {
            for (i = n1 - 1; i >= 0; i--) {
                sf_fast_sweep_2d_stencil (t, vz, vx, theta, rhs,
                                          i, j, n1, n2, d1, d2, shotz, shoty);
            }
        }

        for (j = 0; j < n2; j++) {
            for (i = 0; i < n1; i++) {
                sf_fast_sweep_2d_stencil (t, vz, vx, theta, rhs,
                                          i, j, n1, n2, d1, d2, shotz, shoty);
            }
        }

        for (j = 0; j < n2; j++) {
            for (i = n1 - 1; i >= 0; i--) {
                sf_fast_sweep_2d_stencil (t, vz, vx, theta, rhs,
                                          i, j, n1, n2, d1, d2, shotz, shoty);
            }
        }

    }


}


