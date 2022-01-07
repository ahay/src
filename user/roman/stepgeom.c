/*************************************************************************

Copyright Rice University, 2009.
All rights reserved.

Permission is hereby granted, free of charge, to any person obtaining a
copy of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, provided that the above copyright notice(s) and this
permission notice appear in all copies of the Software and that both the
above copyright notice(s) and this permission notice appear in supporting
documentation.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT OF THIRD PARTY
RIGHTS. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR HOLDERS INCLUDED IN THIS
NOTICE BE LIABLE FOR ANY CLAIM, OR ANY SPECIAL INDIRECT OR CONSEQUENTIAL
DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR
PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE OF
THIS SOFTWARE.

Except as contained in this notice, the name of a copyright holder shall
not be used in advertising or otherwise to promote the sale, use or other
dealings in this Software without prior written authorization of the
copyright holder.

**************************************************************************/

/* Modified for distribution with Madagascar */
#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

#ifdef IWAVE_USE_MPI
#include <mpi.h>
#endif

/* target_p [nz * nx], nx-fast
   source_p
   v = velocity model
*/
/*
const int NBxx = 100;
*/

#define SIGN(X) (((X) > 0) ? (1) : ((X) == 0 ? 0 : (-1)))

#define NB 100
/*^*/
#define SQR(X) ((X)*(X))

typedef struct AgeomB {
/*^*/
    int len;
/*^*/
    int ix1[NB];
/*^*/
    int ix2[NB];
/*^*/
    int iz1[NB];
/*^*/
    int iz2[NB];
/*^*/
    float e[NB];
/*^*/
    float r1[NB];
/*^*/
} geomB;
/*^*/

float norm1(float *wvf, int nz, int nx)
{
/*< norma-L1 >*/
    int ix, iz, ioff;
    float sum = 0.f;

    for (ix=1;ix<nx-1;ix++) {
	for (iz=1;iz<nz-1;iz++) {
	    
	    ioff=iz+ix*nz;
	    sum += fabs(wvf[ioff]);
	}
    }
    return sum / (nx*nz);
}

/* arrB[nx*nz] - array of B-s [nx*nz]
   v[nx*nz] - velo squared(!)
 */

void construct_Boundaries(geomB     *arrB, 
			  float *v,
			  float vmax,
			  float dt, float dz, float dx,
                          int nz, int nx)
{
/*< geoemtrical step forward >*/
    int max_lenB = 4 * (floor(fmaxf(vmax*dt/dx, vmax*dt/dz)) + 1);
    int max_lenX = floor(vmax*dt/dx) + 1;
    int max_lenZ = floor(vmax*dt/dz) + 1;
    int max_lenX2 = max_lenX + max_lenX + 1;
    int max_lenZ2 = max_lenZ + max_lenZ + 1;

    int iz;
    int ix;
    int ioff; 

    int lenB;
    geomB *pB;
    float  vel, vel1, dl;
    int ix1, iz1, ix2, iz2, ioff1, ioff2;
    float TT[NB*NB], Len[NB*NB];

    assert(NB > max_lenB);

    for (ix=1;ix<nx-1;ix++) {
	for (iz=1;iz<nz-1;iz++) {
	    
	    ioff=iz+ix*nz;
	    
	    /* Dejkstra construct TT-traveltime-matrix */
	    vel = sqrtf(v[ioff]);
	    vel1 = 1.f / vel;
	    /* local coords centered at (ix, iz) */
	    for (ix1 = 0; ix1 < max_lenX2; ix1++) {
		for (iz1 = 0; iz1 < max_lenZ2; iz1++) {

		    ioff1 = iz1 + ix1 * NB;

		    /* local domain is of size max_lexX2 x max_lenZ2 */
		    /* local center is at max_lenX,max_lenZ */
		    dl = sqrtf((SQR( (ix1 - max_lenX) * dx ) + SQR( (iz1 - max_lenZ)*dz)));

		    *(TT  + ioff1) = vel1 * dl;

		    *(Len + ioff1) = dl;
		}
	    }
	    /* construct boundary B -list of pairs with e-factor */
	    pB = arrB + ioff;

	    lenB = 0;

	    for (ix1 = 0; ix1 < max_lenX2; ix1++) {
		for (iz1 = 0; iz1 < max_lenZ2; iz1++) {

		    ioff1 = iz1 + ix1 * NB;

		    /* (ix2, iz2) is a prev point on the ray */
		    ix2 = ix1 + SIGN(max_lenX - ix1); /* (int)(ix1 - 0.5f * dx * (float)(ix1 - max_lenX)); */
		    iz2 = iz1 + SIGN(max_lenZ - iz1); /* (int)(iz1 - 0.5f * dz * (float)(iz1 - max_lenZ)); */
		    assert(abs(ix2 - ix1) + abs(iz2 - iz1) > 0 || (ix2 == max_lenX && iz2==max_lenZ));

		    ioff2 = iz2 + ix2 * NB;

		    if (TT[ioff1] > 1e-6 + dt && TT[ioff2]+ 1e-6f < dt) {
			/* add to boundary */
			pB->ix1[lenB] = ix1 - max_lenX;
			pB->iz1[lenB] = iz1 - max_lenZ;

			pB->ix2[lenB] = ix2 - max_lenX;
			pB->iz2[lenB] = iz2 - max_lenZ;
			
			assert (fabs(TT[ioff1] - TT[ioff2]) > 1e-4);

			pB->e [lenB] =  (dt - TT[ioff2])/(TT[ioff1] - TT[ioff2]);
			assert (pB->e[lenB] >= 0.f && pB->e[lenB] <= 1.f);

			pB->r1 [lenB] = 1.f / (pB->e[lenB] * Len[ioff1] + (1.f - pB->e[lenB]) * Len[ioff2]);
		
			lenB++;

			assert(lenB < NB);
		    }		    
		}
	    }

	    pB->len = lenB;

	    /* return; const vel => same B */

	}  /* end iz */
    }  /* end ix */
}


void step_forward_g(geomB * arrB,
		    float * tgt_p, 
		    float * src_p, 
		    float * v,
		    int nz, int nx,
		    float dz, float dx, float dt,
                    int izsrc, int ixsrc)
{
/*< geometrical step forward >*/
    int iz;
    int ix;
    int ioff; 
    int lenB;
    
    geomB *pB;
    float e, sum, r1, sum_r1, prod; /* l2, l02 */ 
    int ix1, iz1, ix2, iz2, k, ioff1, ioff2;
    int nxnz = nx * nz;
    const float dx2 = dx*dx, dz2 = dz*dz;

    for (ix=1;ix<nx-1;ix++) {
	for (iz=1;iz<nz-1;iz++) {

	    ioff=iz+ix*nz;

	    /* const velocity */
	    /* if (!(fabs(vel - v[ioff]) < 1e-6f))
	       sf_warning("v0=%g v %d = %g",vel, ioff, v[ioff]); */

	    pB = arrB + ioff; /* 1 + nz;  + ioff; */

	    // l2 = SQR((ix - ixsrc)*dx) + SQR((iz - izsrc)*dz); 
		   
	    lenB = 0;
	    sum = 0.f;
	    sum_r1 = 0.f;
	    for (k = 0; k < pB->len; k++) {

		e    = pB->e[k];
		assert (e >= 0.f && e <= 1.f);

		r1    = pB->r1[k];

		ix1 = ix + pB->ix1[k];
		iz1 = iz + pB->iz1[k];
		ioff1=iz1+ix1*nz;

		ix2 = ix + pB->ix2[k];
		iz2 = iz + pB->iz2[k];
		ioff2=iz2+ix2*nz;
		
		if (ioff1 >= 0 && ioff2 >= 0 && ioff1 < nxnz && ioff2 < nxnz) {

		    // l02 = SQR((ix1 - ixsrc)*dx) + SQR((iz1 - izsrc)*dz);
		    // if (l2 > l02) {
		    
		    prod = (ix - ix1)*(ix1 - ixsrc)*dx2 + (iz - iz1)*(iz1 - izsrc)*dz2;
		    if (prod > 0.f) {
			sum += r1 * (e * src_p[ioff1] + (1-e) * src_p[ioff2]);
			sum_r1 += r1;
			lenB++;
		    }
		}
	    } /* ik */

	    if (lenB > 0) {
		tgt_p[ioff] = sum / sum_r1; /* / (float)lenB; */
	    }

	    } /* iz */
	} /* ix */
    }

/* here dz *= k and dx *= k, so that Vmax dt < CFL max(dx,dz)
   target_p [nz * nx], nx-fast
   source_p
   v = velocity model
 */

int step_forward_k(float * tgt_p, 
		   float * src_p, 
		   float * v,
		   int nz, int nx,
		   float rz, 
		   float rx, 
		   float s, 
		   int k,
		   float dt, float CFL, float dx, float dz)
/*< step forward >*/
{

    int iz;
    int ix;
    int ioff; 
    float two=2.0;

    // rz=wi.dt*wi.dt/(wi.dz*wi.dz);
    // s =2.0*(rz+rx);

    float 
	k2 = k * k,
	k_rz = rz / k2, 
	k_rx = rx / k2, 
	k_s = 2.f * (k_rz + k_rx);

    int 
	k_nz = k * nz, 
	kp1 = k + 1;

    float
	src_iz_plus_1 ,
	src_iz_minus_1 ,
	src_ix_iz,
	src_ix_plus_1  ,
	src_ix_minus_1,
	tgt_ix_iz,
	vel;
    int
	ioffkzp,
	ioffkzm,
	ioffkxp,
	ioffkxm,
	num_cfl = 0;

    float *tgt = (float*)malloc(nx*nz*sizeof(float));
    memcpy(tgt, tgt_p,nx*nz*sizeof(float));

    for (ix=kp1;ix<nx-kp1;ix++) {
	for (iz=kp1;iz<nz-kp1;iz++) {

	    ioff=iz+ix*nz;

	    vel = sqrtf(v[ioff]);
	    /*
	      if (sqrtf(v[ioff])*dt<CFL*fmaxf(dx,dz)) {

		tgt_p[ioff]=two*src_p[ioff]-tgt_p[ioff]
		    +v[ioff]*(rz*(src_p[ioff+1]+src_p[ioff-1]) +
			      rx*(src_p[ioff+nz]+src_p[ioff-nz]) - 
			      s*src_p[ioff]);    
		num_cfl++;
		continue;
		} */
	    k = floor(vel * dt / (CFL*fmaxf(dx,dz))) + 1;
	    k2 = k * k;
	    k_rz = rz / k2; 
	    k_rx = rx / k2;
	    k_s = 2.f * (k_rz + k_rx);

	    k_nz = k * nz;
	    kp1 = k + 1;

	    /* src_ix_iz = src_p[ioff];
	    tgt_ix_iz = tgt_p[ioff];

	    src_iz_plus_1  = src_p[ioff + k];
	    src_iz_minus_1 = src_p[ioff - k];

	    src_ix_plus_1  = src_p[ioff + k_nz];
	    src_ix_minus_1 = src_p[ioff - k_nz];
	    */
	    /*
	    src_ix_iz = 0.25f *(2.f*src_p[ioff] + src_p[ioff+1] + src_p[ioff+nz]);

	    src_iz_plus_1  = 0.5 * (src_p[ioff + k1] + src_p[ioff + k]);
	    src_iz_minus_1 = 0.5 * (src_p[ioff - k1] + src_p[ioff - k]);

	    src_ix_plus_1  = 0.5f * (src_p[ioff + k_nz] + src_p[ioff + k1_nz]);
	    src_ix_minus_1 = 0.5f * (src_p[ioff - k_nz] + src_p[ioff - k1_nz]);
	    */
	    ioffkzp = ioff + k;
	    ioffkzm = ioff - k;
	    ioffkxp = ioff + k_nz;
	    ioffkxm = ioff - k_nz;

	    tgt_ix_iz = 0.125f *(4.f*tgt[ioff]+tgt[ioff+1]+tgt[ioff+nz]+tgt[ioff-1]+tgt[ioff-nz]);

	    src_ix_iz = 0.125f *(4.f*src_p[ioff]+src_p[ioff+1]+src_p[ioff+nz]+src_p[ioff-1] + src_p[ioff-nz]);

	    src_ix_plus_1 = 0.125f *(4.f*src_p[ioffkxp]+src_p[ioffkxp+1]+src_p[ioffkxp+nz]+src_p[ioffkxp-1] + src_p[ioffkxp-nz]);
	    src_ix_minus_1 = 0.125f *(4.f*src_p[ioffkxm]+src_p[ioffkxm+1]+src_p[ioffkxm+nz]+src_p[ioffkxm-1] + src_p[ioffkxm-nz]);

	    src_iz_plus_1 = 0.125f *(4.f*src_p[ioffkzp]+src_p[ioffkzp+1]+src_p[ioffkzp+nz]+src_p[ioffkzp-1] + src_p[ioffkzp-nz]);
	    src_iz_minus_1 = 0.125f *(4.f*src_p[ioffkzm]+src_p[ioffkzm+1]+src_p[ioffkzm+nz]+src_p[ioffkzm-1] + src_p[ioffkzm-nz]);

	    vel = 0.125f *(4.f*v[ioff]+v[ioff+1]+v[ioff+nz]+v[ioff-1]+v[ioff-nz]);


	    tgt_p[ioff]=two*src_ix_iz - tgt_ix_iz
/*
		+v[ioff]*(rz*(src_p[ioff+1]+src_p[ioff-1]) +
			  rx*(src_p[ioff+nz]+src_p[ioff-nz]) - 
			  s*src_p[ioff]);    */
		+vel*(k_rz*(src_iz_plus_1 + src_iz_minus_1) +
			  k_rx*(src_ix_plus_1 + src_ix_minus_1) -
			  k_s * src_ix_iz);
	}
    }

    free(tgt);

    return num_cfl;
}


/* target_p [nz * nx], nx-fast
   source_p
   v = velocity model
 */

void step_forward(float * tgt_p, 
		  float * src_p, 
		  float * v,
		  int nz, int nx,
		  float rz, float rx, float s) 
/*< step forward in k-scale >*/
{

    int iz;
    int ix;
    int ioff; 
    float two=2.0;
    
    for (ix=1;ix<nx-1;ix++) {
	for (iz=1;iz<nz-1;iz++) {
	    ioff=iz+ix*nz;
	    tgt_p[ioff]=two*src_p[ioff]-tgt_p[ioff]
		+v[ioff]*(rz*(src_p[ioff+1]+src_p[ioff-1]) +
			  rx*(src_p[ioff+nz]+src_p[ioff-nz]) - 
			  s*src_p[ioff]);    
	}
    }
}

