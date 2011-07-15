// Lowrank decomposition for 2-D isotropic wave propagation (pseudo-depth equation). 

//   Copyright (C) 2011 KAUST
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
static std::valarray<double> ks;
static float dt;

static std::valarray<float> vxxs,vxzs,vzzs;
static std::valarray<double> kxxs,kxzs,kzzs;

int sample(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,0.0);
    double w;
    for(int a=0; a<nr; a++) {
	for(int b=0; b<nc; b++) {
	    //res(a,b) = 2*(cos(vs[rs[a]]*ks[cs[b]]*dt)-1); 
	    w=2*SF_PI*sqrt(vxxs[rs[a]]*kxxs[cs[b]]
		          +vxzs[rs[a]]*kxzs[cs[b]]
		          +vzzs[rs[a]]*kzzs[cs[b]]);
	    // res(a,b) = 2*(cos(w*dt)-1); 
	    res(a,b) = 2*cos(w*dt); 
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

    float v0;
    par.get("v0",v0,1.0); // scaling velocity

    iRSF vel;

    int nz,nx;
    vel.get("n1",nz);
    vel.get("n2",nx);
    int m = nx*nz;
    std::valarray<float> vels(m);
    vel >> vels;
    vs.resize(m);
    vs = vels;
    
    iRSF sig("sig"); // sigma
    std::valarray<float> ss(m);
    sig >> ss;

    iRSF fft("fft");

    int nkz,nkx;
    fft.get("n1",nkz);
    fft.get("n2",nkx);

    float dkz,dkx;
    fft.get("d1",dkz);
    fft.get("d2",dkx);
    
    float kz0,kx0;
    fft.get("o1",kz0);
    fft.get("o2",kx0);

    float kx, kz;

    int n = nkx*nkz;  
    std::valarray<double> k(n);
    std::valarray<double> kxx(n);
    std::valarray<double> kxz(n);
    std::valarray<double> kzz(n);
    int i;
    for (int ix=0; ix < nkx; ix++) {
	kx = kx0 + ix*dkx;
	for (int iz=0; iz < nkz; iz++) {
	    kz = kz0 + iz*dkz;
	    i = iz+ix*nkz;
	    k[i] = 2*SF_PI*hypot(kx,kz);
	    kxx[i]=kx*kx;
	    kxz[i]=kx*kz;
	    kzz[i]=kz*kz;
	}
    }
    ks.resize(n);
    ks = k;
    kxxs.resize(n); kxxs=kxx;
    kxzs.resize(n); kxzs=kxz;
    kzzs.resize(n); kzzs=kzz;

    std::valarray<float> vxx(m);
    std::valarray<float> vxz(m);
    std::valarray<float> vzz(m);
    float vsq;
    for (int ix=0; ix<nx; ix++) {
	for (int iz=0; iz<nz; iz++) {
	    i = iz+ix*nz;
	    vsq = vs[i]*vs[i];
	    vxx[i] = vsq;
	    vxz[i] = 2*v0*ss[i]*vsq;
	    vzz[i] = v0*v0*(vsq*ss[i]*ss[i] + 1);
	}
    }
    vxxs.resize(m); vxxs=vxx;
    vxzs.resize(m); vxzs=vxz;
    vzzs.resize(m); vzzs=vzz;

    vector<int> lidx, ridx;
    DblNumMat mid;

    iC( lowrank(m,n,sample,(double) eps,npk,lidx,ridx,mid) );

    int m2=mid.m();
    int n2=mid.n();
    double *dmid = mid.data();

    std::valarray<float> fmid(m2*n2);
    for (int k=0; k < m2*n2; k++) {
	fmid[k] = dmid[k];
    }

    oRSF middle;
    middle.put("n1",m2);
    middle.put("n2",n2);
    middle << fmid;

    vector<int> midx(m), nidx(n);
    for (int k=0; k < m; k++) 
	midx[k] = k;
    for (int k=0; k < n; k++) 
	nidx[k] = k;    

    DblNumMat lmat(m,m2);
    iC ( sample(midx,lidx,lmat) );
    double *ldat = lmat.data();

    std::valarray<float> ldata(m*m2);
    for (int k=0; k < m*m2; k++) 
	ldata[k] = ldat[k];
    oRSF left("left");
    left.put("n1",m);
    left.put("n2",m2);
    left << ldata;

    DblNumMat rmat(n2,n);
    iC ( sample(ridx,nidx,rmat) );
    double *rdat = rmat.data();

    std::valarray<float> rdata(n2*n);    
    for (int k=0; k < n2*n; k++) 
	rdata[k] = rdat[k];
    oRSF right("right");
    right.put("n1",n2);
    right.put("n2",n);
    right << rdata;

    return 0;
}
