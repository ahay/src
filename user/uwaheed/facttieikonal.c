/* Fast sweeping factored TTI eikonal solver (2D) */
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

bool sf_init_fast_sweep (float *tau,
                         int n2, int n1,
                         float o2, float o1,
                         float d2, float d1,
                         int shoty, int shotz, int fac)
/*< initialize >*/
{
    int i, n12;

    if (NULL == tau)
        return false;

    if (shoty < 0 || shoty >= n2 ||
        shotz < 0 || shotz >= n1)
        return false;

    n12 = n1 * n2;

    for (i = 0; i < n12; i++)
        tau[i] = SF_HUGE;

    if(fac==0) tau[shoty * n1 + shotz] = 0.0;
    else if(fac==1) tau[shoty * n1 + shotz] = 1.0;
    else sf_error("Choose fac=0(Additive) or fac=1(Multiplicative) factorization");

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


static void sf_fast_sweep_2d_stencil_add (float *tau, float *T0,
                                      float *py0, float *pz0,
                                      float *vz, float *vx,
                                      float *theta, float *rhs,
                                      int i, int j, int n1, int n2,
                                      float d1, float d2,
                                      int shotz, int shoty) {

    /* A variable with c at the end means its value of current grid point.
       For example: vzc is vz at (i,j) or the current grid point */

    /* Sign variable is chosen assuming that the information is spread along 
       the positive y and z dimensions */

    float vzc = vz[j * n1 + i], vxc = vx[j * n1 + i], thetac = theta[j*n1+i], rhsc = rhs[j*n1+i];
    float T0c = T0[j * n1 + i];
    float tauij = tau[j * n1 + i], root;
    int sy=0, sz=0;
    float dz = d1, dz2 = dz*dz, dy = d2, dy2 = dy*dy;
    float vzc2 = vzc*vzc, vxc2 = vxc*vxc;
    float ct = cos(thetac), st = sin(thetac);
    float ct2 = ct*ct, st2 = st*st;    
    float a = vxc2*ct2 + vzc2*st2, b = vxc2*st2 + vzc2*ct2, c = (vxc2-vzc2)*st*ct; 
    float py0c = py0[j * n1 + i], pz0c=pz0[j * n1 + i];
    float tauy, tauz, T0y, T0z, tauC, tauCy, tauCz;
    float num, den, t1, t2;

    if(i==shotz && j==shoty) return;

    if(i==0 || i==n1-1){
        if(i==0){
            tauz = tau[j * n1 + 1]; T0z = T0[j * n1 + 1]; sz = -1; 
        }
        else{
            tauz = tau[j * n1 + n1-2]; T0z = T0[j * n1 + n1-2]; sz = 1; 
        }
    }
    else{
        if((tau[j*n1+i-1]+T0[j*n1+i-1]) < (tau[j*n1+i+1]+T0[j*n1+i+1])){
            tauz = tau[j * n1 + i - 1]; T0z = T0[j * n1 + i - 1]; sz = 1;
        }
		else{
            tauz = tau[j * n1 + i + 1]; T0z = T0[j * n1 + i + 1]; sz = -1; 
        }

    }

    if(j==0 || j==n2-1){
        if(j==0){
            tauy = tau[n1 + i]; T0y = T0[n1 + i]; sy = -1; 
        }
        else{
            tauy = tau[(n2-2)*n1 + i]; T0y = T0[(n2-2)*n1 + i]; sy = 1; 
        }
    }
    else{
        if((tau[(j-1)*n1+i]+T0[(j-1)*n1+i]) < (tau[(j+1)*n1+i]+T0[(j+1)*n1+i])){
            tauy = tau[(j-1)*n1+i]; T0y = T0[(j-1)*n1+i]; sy = 1;
        }
		else{
            tauy = tau[(j+1) * n1 + i ]; T0y = T0[(j+1) * n1 + i ]; sy = -1; 
        }

    }
    

    if (SF_HUGE == tauy && SF_HUGE == tauz)
        return;


    /* Took absolute of py0c and pz0c because point C is in the upwind direction */

	tauCy = tauy - dy*fabs(py0c) + dy*sqrtf((b*rhsc)/(a*b-c*c));
	tauCz = tauz - dz*fabs(pz0c) + dz*sqrtf((a*rhsc)/(a*b-c*c));



    /* The condition used below is different than the one in the isotropic eikonal solver 
       and it holds real importance in getting the correct solution. */

    if (tauy == SF_HUGE) {
        tauC = tauCz;
    } else if (tauz == SF_HUGE) {
        tauC = tauCy;
    } else {

        t1 = a*dz2*tauy + dy2*sz*(-dz*(c*py0c + b*pz0c) + b*sz*tauz) 
              + dy*dz*sy*(-a*dz*py0c + c*(-dz*pz0c + sz*(tauy+tauz)));

        t2 = 2*c*dy*dz*rhsc*sy*sz + b*dy2*rhsc + c*c*(dz*pz0c*sy - dy*py0c*sz 
             + sy*sz*tauy -sy*sz*tauz)*(dz*pz0c*sy - dy*py0c*sz + sy*sz*tauy -sy*sz*tauz) 
             + a*(dz2*rhsc - b*(dz*pz0c*sy - dy*py0c*sz + sy*sz*tauy - sy*sz*tauz)
             *(dz*pz0c*sy - dy*py0c*sz + sy*sz*tauy - sy*sz*tauz));
              

        num = t1 + dy2*dz2*sqrt(t2/(dy2*dz2));

        den = a*dz2 + dy*sz*(2*c*dz*sy + b*dy*sz);


        root = num/den;


	    if(isCausalRoot(root+T0c,tauz+T0z,tauy+T0y,sz,sy,a,b,c,dz,dy)) tauC=root;
        else {
            if(tauCy+T0y < tauCz+T0z) tauC=tauCy;
            else tauC = tauCz;
        }

    }

    if (tauC+T0c < tauij+T0c)
        tau[j * n1 + i] = tauC;
}

static void sf_fast_sweep_2d_stencil_mul (float *tau, float *T0,
                                      float *py0, float *pz0,
                                      float *vz, float *vx,
                                      float *theta, float *rhs,
                                      int i, int j, int n1, int n2,
                                      float d1, float d2,
                                      int shotz, int shoty) {

    /* A variable with c at the end means its value of current grid point.
       For example: vzc is vz at (i,j) or the current grid point */

    /* Sign variable is chosen assuming that the information is spread along 
       the positive y and z dimensions */

    float vzc = vz[j * n1 + i], vxc = vx[j * n1 + i], thetac = theta[j*n1+i], rhsc = rhs[j*n1+i];
    float T0c = T0[j * n1 + i], T0c2 = T0c*T0c;
    float tauij = tau[j * n1 + i], root;
    int sy=0, sz=0;
    float dz = d1, dz2 = dz*dz, dy = d2, dy2 = dy*dy;
    float vzc2 = vzc*vzc, vxc2 = vxc*vxc;
    float ct = cos(thetac), st = sin(thetac);
    float ct2 = ct*ct, st2 = st*st;    
    float a = vxc2*ct2 + vzc2*st2, b = vxc2*st2 + vzc2*ct2, c = (vxc2-vzc2)*st*ct; 
    float py0c = py0[j * n1 + i], pz0c=pz0[j * n1 + i];
    float py0c2 = py0c*py0c, pz0c2 = pz0c*pz0c;
    float tauy, tauz, T0y, T0z, tauC, tauCy, tauCz;
    float num, den, t1, t2;

    if(i==shotz && j==shoty) return;

    if(i==0 || i==n1-1){
        if(i==0){
            tauz = tau[j * n1 + 1]; T0z = T0[j * n1 + 1]; sz = -1; 
        }
        else{
            tauz = tau[j * n1 + n1-2]; T0z = T0[j * n1 + n1-2]; sz = 1; 
        }
    }
    else{
        if((tau[j*n1+i-1]*T0[j*n1+i-1]) < (tau[j*n1+i+1]*T0[j*n1+i+1])){
            tauz = tau[j * n1 + i - 1]; T0z = T0[j * n1 + i - 1]; sz = 1;
        }
		else{
            tauz = tau[j * n1 + i + 1]; T0z = T0[j * n1 + i + 1]; sz = -1; 
        }

    }

    if(j==0 || j==n2-1){
        if(j==0){
            tauy = tau[n1 + i]; T0y = T0[n1 + i]; sy = -1; 
        }
        else{
            tauy = tau[(n2-2)*n1 + i]; T0y = T0[(n2-2)*n1 + i]; sy = 1; 
        }
    }
    else{
        if((tau[(j-1)*n1+i]*T0[(j-1)*n1+i]) < (tau[(j+1)*n1+i]*T0[(j+1)*n1+i])){
            tauy = tau[(j-1)*n1+i]; T0y = T0[(j-1)*n1+i]; sy = 1;
        }
		else{
            tauy = tau[(j+1) * n1 + i ]; T0y = T0[(j+1) * n1 + i ]; sy = -1; 
        }

    }
    

    if (SF_HUGE == tauy && SF_HUGE == tauz)
        return;


    /* Took absolute of py0c and pz0c because point C is in the upwind direction */

    tauCy = (T0c*tauy + dy*sqrtf(b*rhsc/(a*b-c*c)))/(T0c + fabs(py0c)*dy);
    tauCz = (T0c*tauz + dz*sqrtf(a*rhsc/(a*b-c*c)))/(T0c + fabs(pz0c)*dz);


    // The condition used below is different than the one in the isotropic eikonal solver 
    // and it holds real importance in getting the correct solution.

    if (tauy == SF_HUGE) {
        tauC = tauCz;
    } else if (tauz == SF_HUGE) {
        tauC = tauCy;
    } else {

        t1 = a*dz2*T0c2*tauy + dy2*sz*T0c*(c*dz*py0c + b*dz*pz0c + b*sz*T0c)*tauz 
             + dy*dz*sy*T0c*(a*dz*py0c*tauy + c*(dz*pz0c*tauy + sz*T0c*(tauy+tauz)));

        t2 = 2*c*dy*dz*rhsc*(dy*py0c+sy*T0c)*(dz*pz0c+sz*T0c) 
             + b*dy2*rhsc*(dz*pz0c + sz*T0c)*(dz*pz0c + sz*T0c) + c*c*T0c*T0c*(dz*pz0c*sy*tauy 
             + sy*sz*T0c*(tauy-tauz) - dy*py0c*sz*tauz)*(dz*pz0c*sy*tauy + sy*sz*T0c*(tauy-tauz) 
             - dy*py0c*sz*tauz) + a*(-T0c2*(dz2*(-rhsc+b*pz0c2*tauy*tauy) 
             + 2*b*dz*pz0c*sz*T0c*tauy*(tauy-tauz) + b*T0c2*(tauy-tauz)*(tauy-tauz)) 
             + 2*dy*py0c*sy*T0c*(dz2*rhsc + b*dz*pz0c*sz*T0c*tauy*tauz + b*T0c2*(tauy-tauz)*tauz) 
             + dy2*py0c2*(dz2*rhsc - b*T0c2*tauz*tauz));
           

        num = t1 + dy2*dz2*sqrtf(t2/(dy2*dz2));

        den = a*dz2*(dy*py0c + sy*T0c)*(dy*py0c + sy*T0c) + dy*(dz*pz0c + sz*T0c)*(2*c*dz*(dy*py0c 
              + sy*T0c) + b*dy*(dz*pz0c + sz*T0c));


        root = num/den;



	    if(isCausalRoot(root*T0c,tauz*T0z,tauy*T0y,sz,sy,a,b,c,dz,dy)) tauC=root;
        else {
            if(tauCy*T0y < tauCz*T0z) tauC=tauCy;
            else tauC = tauCz;
        }

    }

    if (tauC*T0c < tauij*T0c)
        tau[j * n1 + i] = tauC;
}


void sf_run_fast_sweep (float *tau, float *T0, 
                        float *py0, float *pz0,
                        float *vz, float *vx,
                        float *theta, float *rhs,
                        int niter,
                        int n2, int n1,
                        float o2, float o1,
                        float d2, float d1,
                        int shoty, int shotz, int fac)
/*< run sweeps >*/
{
    int i, j , l = 0;

    if(fac==0){

        for (l = 0; l < niter; l++) {

            for (j = n2 - 1; j >= 0; j--) {
                for (i = 0; i < n1; i++) {
                    sf_fast_sweep_2d_stencil_add (tau, T0, py0, pz0, vz, vx,
                                              theta, rhs, i, j, n1, n2, d1, d2, shotz, shoty);
                }
            }
            for (j = n2 - 1; j >= 0; j--) {
                for (i = n1 - 1; i >= 0; i--) {
                    sf_fast_sweep_2d_stencil_add (tau, T0, py0, pz0, vz, vx,
                                              theta, rhs, i, j, n1, n2, d1, d2, shotz, shoty);
                }
            }
            for (j = 0; j < n2; j++) {
                for (i = 0; i < n1; i++) {
                    sf_fast_sweep_2d_stencil_add (tau, T0, py0, pz0, vz, vx,
                                          theta, rhs, i, j, n1, n2, d1, d2, shotz, shoty);
                }
            }


            for (j = 0; j < n2; j++) {
                for (i = n1 - 1; i >= 0; i--) {
                    sf_fast_sweep_2d_stencil_add (tau, T0, py0, pz0, vz, vx,
                                              theta, rhs, i, j, n1, n2, d1, d2, shotz, shoty);
                }
            }

        }

    }
    else if(fac==1){

        for (l = 0; l < niter; l++) {

            for (j = n2 - 1; j >= 0; j--) {
                for (i = 0; i < n1; i++) {
                    sf_fast_sweep_2d_stencil_mul (tau, T0, py0, pz0, vz, vx,
                                              theta, rhs, i, j, n1, n2, d1, d2, shotz, shoty);
                }
            }
            for (j = n2 - 1; j >= 0; j--) {
                for (i = n1 - 1; i >= 0; i--) {
                    sf_fast_sweep_2d_stencil_mul (tau, T0, py0, pz0, vz, vx,
                                              theta, rhs, i, j, n1, n2, d1, d2, shotz, shoty);
                }
            }
            for (j = 0; j < n2; j++) {
                for (i = 0; i < n1; i++) {
                    sf_fast_sweep_2d_stencil_mul (tau, T0, py0, pz0, vz, vx,
                                          theta, rhs, i, j, n1, n2, d1, d2, shotz, shoty);
                }
            }


            for (j = 0; j < n2; j++) {
                for (i = n1 - 1; i >= 0; i--) {
                    sf_fast_sweep_2d_stencil_mul (tau, T0, py0, pz0, vz, vx,
                                              theta, rhs, i, j, n1, n2, d1, d2, shotz, shoty);
                }
            }

        }

    }
    else sf_error("Choose fac=0(Additive) or fac=1(Multiplicative) factorization");
}


