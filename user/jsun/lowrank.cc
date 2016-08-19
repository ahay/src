// Lowrank decomposition for 3-D isotropic wave propagation. 
//   Copyright (C) 2016 University of Texas at Austin
//  
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License, or
//   (at your option) any later version.
//  
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//  
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#include <rsf.hh>
#include "vecmatop.hh"
#include "serialize.hh"
/* head files aumatically produced from C programs */
extern "C"{
#include "rtmutil.h"
}

using namespace std;

static std::valarray<float> vs;
static std::valarray<float> ks;
static float dt,eps;
static int m, n, npk, seed;
static vector<int> lidx, ridx;
static CpxNumMat mid;
static vector<int> ms, ns, js;

int sample(vector<int>& rs, vector<int>& cs, CpxNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,cpx(0.0f,0.0f));
    for(int a=0; a<nr; a++) {
	for(int b=0; b<nc; b++) {
	    float phase = vs[rs[a]]*ks[cs[b]]*dt;
	    res(a,b) = cpx(cos(phase),sin(phase));
	}
    }
    return 0;
}

void lowrank_init(int jump,
                  int seed_,
                  int npk_,
                  float eps_,
                  float dt_,
                  float *vel,
                  fdm3d fdm,
                  dft3d dft)
/*< initialize lowrank >*/
{   
    seed = seed_;
    npk  = npk_;
    eps  = eps_;
    dt   = dt_;
    m = fdm->nzpad*fdm->nxpad*fdm->nypad;
    n = dft->nkz*dft->nkx*dft->nky;
    ms.resize(3); ms[0] = fdm->nzpad; ms[1] = fdm->nxpad; ms[2] = fdm->nypad;
    ns.resize(3); ns[0] = dft->nkz;   ns[1] = dft->nkx;   ns[2] = dft->nky;
    js.resize(3); js[0] = jump;       js[1] = jump;       js[2] = jump;

    vs.resize(m);
    for (int i=0; i<m; i++) vs[i] = vel[i];
    
    ks.resize(n);
    for (int iy=0; iy < dft->nky; iy++) {
	float ky = dft->oky+iy*dft->dky;
	for (int ix=0; ix < dft->nkx; ix++) {
	    float kx = dft->okx+ix*dft->dkx;
	    for (int iz=0; iz < dft->nkz; iz++) {
		float kz = dft->okz+iz*dft->dkz;
		ks[iz+dft->nkz*(ix+dft->nkx*iy)] = 2*SF_PI*sqrt(kx*kx+kz*kz+ky*ky);
	    }
	}
    }
}

int lowrank_rank()
/*< perform lowrank decomposition >*/
{
    int nrank;

    srand48(seed);

    //iC( lowrank(m,n,sample,eps,npk,lidx,ridx,mid) );
    iC( lowrank(ms,ns,js,sample,eps,npk,lidx,ridx,mid) );

    nrank = mid.n();
    return nrank;
}

void lowrank_mat(sf_complex **lt, sf_complex **rt)
/*< output lowrank matrices >*/
{
    int m2=mid.m();
    int n2=mid.n();

    vector<int> midx(m), nidx(n);
    for (int k=0; k < m; k++) 
	midx[k] = k;
    for (int k=0; k < n; k++) 
	nidx[k] = k;    

    CpxNumMat lmat(m,m2);
    iC ( sample(midx,lidx,lmat) );

    CpxNumMat lmat2(m,n2);
    iC( zgemm(1.0, lmat, mid, 0.0, lmat2) );

    cpx *ldat = lmat2.data();

    for (int i=0; i < n2; i++)
        for (int k=0; k < m; k++)
            lt[i][k] = sf_cmplx(real(ldat[i*m+k]),imag(ldat[i*m+k]));

    CpxNumMat rmat(n2,n);
    iC ( sample(ridx,nidx,rmat) );
    cpx *rdat = rmat.data();

    for (int k=0; k < n; k++)
        for (int i=0; i < n2; i++)
            rt[k][i] = sf_cmplx(real(rdat[k*n2+i]),imag(rdat[k*n2+i]));

}
