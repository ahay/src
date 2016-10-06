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

static std::valarray<float> vs, vx, vz, q, t;
static std::valarray<float> ks, kzs, kxs, kys;
static float dt,eps;
static int m, n, npk, seed;
static vector<int> lidx, ridx;
static ZpxNumMat mid;
static vector<int> ms, ns, js;

static int (*sample)(vector<int>& rs, vector<int>& cs, ZpxNumMat& res);

static int sample_iso(vector<int>& rs, vector<int>& cs, ZpxNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,zpx(0.0,0.0));
    for(int a=0; a<nr; a++) {
	for(int b=0; b<nc; b++) {
	    double phase = vs[rs[a]]*ks[cs[b]]*dt;
	    res(a,b) = zpx(cos(phase),sin(phase));
	}
    }
    return 0;
}

static int sample_tti(vector<int>& rs, vector<int>& cs, ZpxNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,zpx(0.0,0.0));
    for(int a=0; a<nr; a++) {
	int i=rs[a];
	double wx = vx[i]*vx[i];
	double wz = vz[i]*vz[i];
	double qq = q[i];
	double tt = t[i];
	double c = cos(tt);
	double s = sin(tt);
	
	for(int b=0; b<nc; b++) {
	    int j = cs[b];
	    double r;
	    double x0 = kxs[j];
	    double z0 = kzs[j];
	    // rotation of coordinates
	    double x = x0*c+z0*s;
	    double z = z0*c-x0*s;

            z = wz*z*z;
            x = wx*x*x;
            r = x+z;
            r = r+sqrt(r*r-qq*x*z);
            r = sqrt(0.5*r);

            res(a,b) = zpx(cos(r*dt),sin(r*dt));
	}
    }
    return 0;
}


void lowrank_init(int jump,
                  int seed_,
                  int npk_,
                  float eps_,
                  float dt_,
                  int media,
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

    // media
    switch(media) {
        case 1:
            sample = &sample_tti;
            kzs.resize(n); kxs.resize(n); kys.resize(n);
            for (int iy=0; iy < dft->nky; iy++) {
                float ky = dft->oky+iy*dft->dky;
                for (int ix=0; ix < dft->nkx; ix++) {
                    float kx = dft->okx+ix*dft->dkx;
                    for (int iz=0; iz < dft->nkz; iz++) {
                        float kz = dft->okz+iz*dft->dkz;
                        kzs[iz+dft->nkz*(ix+dft->nkx*iy)] = 2*SF_PI*kz;
                        kxs[iz+dft->nkz*(ix+dft->nkx*iy)] = 2*SF_PI*kx;
                        kys[iz+dft->nkz*(ix+dft->nkx*iy)] = 2*SF_PI*ky;
                    }
                }
            }
            break;
        default:
            sample = &sample_iso;
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
            break;
    }
}

void lowrank_iso(float *vel)
/*< iso model setup >*/
{
    vs.resize(m);
    for (int i=0; i<m; i++) vs[i] = vel[i];
}

void lowrank_tti(float *velx, float *velz, float *eta, float *theta)
/*< iso model setup >*/
{
    vx.resize(m);
    vz.resize(m);
    q.resize(m);
    t.resize(m);
    for (int i=0; i<m; i++) {
        vx[i] = velx[i];
        vz[i] = velz[i];
        q[i]  = eta[i];
        t[i]  = theta[i]*SF_PI/180.; /* from degrees to radians */
    }
}

int lowrank_rank()
/*< perform lowrank decomposition >*/
{
    int nrank;

    srand48(seed);

    if(ms[2]>1) iC( ddlowrank(ms,ns,js,sample,(double)eps,npk,lidx,ridx,mid) );
    else iC( ddlowrank(m,n,sample,(double)eps,npk,lidx,ridx,mid) );

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

    ZpxNumMat lmat(m,m2);
    iC ( sample(midx,lidx,lmat) );

    ZpxNumMat lmat2(m,n2);
    iC( zzgemm(1.0, lmat, mid, 0.0, lmat2) );

    zpx *ldat = lmat2.data();

    for (int i=0; i < n2; i++)
        for (int k=0; k < m; k++)
            lt[i][k] = sf_cmplx(real(ldat[i*m+k]),imag(ldat[i*m+k]));

    ZpxNumMat rmat(n2,n);
    iC ( sample(ridx,nidx,rmat) );
    zpx *rdat = rmat.data();

    for (int k=0; k < n; k++)
        for (int i=0; i < n2; i++)
            rt[k][i] = sf_cmplx(real(rdat[k*n2+i]),imag(rdat[k*n2+i]));

}
