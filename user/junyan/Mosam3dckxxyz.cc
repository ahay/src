// Lowrank decomposition for 3-D transversely isotropic wave propagation. 
// Copyright (C)2014 Institute of Geology and Geophysics, Chinese Academy of Sciences (Jun Yan) 
//					 2009 University of Texas at Austin
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
#include <time.h>

#include <rsf.hh>

#include "vecmatop.hh"
#include "serialize.hh"

using namespace std;

static std::valarray<float>  vs;
static std::valarray<float> ks;
static std::valarray<double> kx, ky, kz;
static float dt;

int sample(vector<int>& rs, vector<int>& cs, FltNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,0.0f);
    for(int a=0; a<nr; a++) {
	for(int b=0; b<nc; b++) {
       res(a,b) = -1.0*kx[cs[b]]*kx[cs[b]]*ky[cs[b]]*kz[cs[b]]/(ks[cs[b]]*ks[cs[b]])*2.0*(cos(vs[rs[a]]*ks[cs[b]]*dt)-1.0)
					/(vs[rs[a]]*vs[rs[a]]*dt*dt*ks[cs[b]]*ks[cs[b]]);
	}
    }
    return 0;
}

int main(int argc, char** argv)
{   
    sf_init(argc,argv); // Initialize RSF

    iRSF par(0);
    int seed;

    par.get("seed",seed,time(NULL)); // seed for random number generator
    srand48(seed);

    float eps;
    par.get("eps",eps,1.e-4); // tolerance

    int npk;
    par.get("npk",npk,20); // maximum rank

    par.get("dt",dt); // time step

    iRSF vel;

    int nz,nx,ny;
    vel.get("n1",nz);
    vel.get("n2",nx);
    vel.get("n3",ny);
    int m = nx*ny*nz;
    std::valarray<float> vels(m);
    vel >> vels;
    vs.resize(m);
    vs = vels;
    
    iRSF fft("fft");

    int nkz,nkx,nky;
    fft.get("n1",nkz);
    fft.get("n2",nkx);
    fft.get("n3",nky);

    float dkz,dkx,dky;
    fft.get("d1",dkz);
    fft.get("d2",dkx);
    fft.get("d3",dky);
    
    float kz0,kx0,ky0;
    fft.get("o1",kz0);
    fft.get("o2",kx0);
    fft.get("o3",ky0);

    float kx1,ky1,kz1;

    int n = nkx*nky*nkz;
    std::valarray<float> k(n);
    kx.resize(n);
    ky.resize(n);
    kz.resize(n);
    for (int iy=0; iy < nky; iy++) {
	ky1 = ky0+iy*dky;
	for (int ix=0; ix < nkx; ix++) {
	    kx1 = kx0+ix*dkx;
	    for (int iz=0; iz < nkz; iz++) {
		kz1 = kz0+iz*dkz;
		k[iz+nkz*(ix+nkx*iy)] = 2*SF_PI*sqrt(kx1*kx1+kz1*kz1+ky1*ky1);
		kx[iz+nkz*(ix+nkx*iy)] = 2*SF_PI*kx1;
		ky[iz+nkz*(ix+nkx*iy)] = 2*SF_PI*ky1;
		kz[iz+nkz*(ix+nkx*iy)] = 2*SF_PI*kz1;
	    }
	}
    }
    ks.resize(n);
    ks = k;

    vector<int> lidx, ridx;
    FltNumMat mid;

    iC( lowrank(m,n,sample,eps,npk,lidx,ridx,mid) );

    int m2=mid.m();
    int n2=mid.n();

    vector<int> midx(m), nidx(n);
    for (int k=0; k < m; k++) 
	midx[k] = k;
    for (int k=0; k < n; k++) 
	nidx[k] = k;    

    FltNumMat lmat(m,m2);
    iC ( sample(midx,lidx,lmat) );

    FltNumMat lmat2(m,n2);
    iC( dgemm(1.0, lmat, mid, 0.0, lmat2) );

    float *ldat = lmat2.data();

    std::valarray<float> ldata(m*n2);
    for (int k=0; k < m*n2; k++) 
	ldata[k] = ldat[k];
    oRSF left("left");
    left.put("n1",m);
    left.put("n2",n2);
    left.put("n3",1);
    left << ldata;

    FltNumMat rmat(n2,n);
    iC ( sample(ridx,nidx,rmat) );
    float *rdat = rmat.data();

    std::valarray<float> rdata(n2*n);    
    for (int k=0; k < n2*n; k++) 
	rdata[k] = rdat[k];
    oRSF right;
    right.put("n1",n2);
    right.put("n2",n);
    right.put("n3",1);
    right << rdata;

    return 0;
}
