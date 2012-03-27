// Lowrank decomposition for 2-D prestack exploding reflector in V(z)

//   Copyright (C) 2010 University of Texas at Austin
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
static float dt;
static int nh;

int sample(vector<int>& rs, vector<int>& cs, FltNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,0.0f);
    for(int a=0; a<nr; a++) {
	for(int b=0; b<nc; b++) {
	    res(a,b) = 2*(cos(vs[rs[a]/nh]*ks[cs[b]]*dt)-1); 
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

    float tol;
    par.get("tol",tol,1.e-4); // tolerance

    int npk;
    par.get("npk",npk,20); // maximum rank

    par.get("dt",dt); // time step

    iRSF vel;

    int nz,nx; // velocity comes transposed
    vel.get("n1",nx);
    vel.get("n2",nz);
    
    int m = nx*nz;
    std::valarray<float> vels(m);
    vel >> vels;    
    vs.resize(m);
    vs = vels;

    par.get("nh",nh); /* number of offsets */
    m *= nh;
    
    iRSF fft("fft");

    int nkh,nkx,nkz;
    fft.get("n1",nkh);
    fft.get("n2",nkx);
    fft.get("n3",nkz);

    int n = nkh*nkx*nkz;

    float dkh,dkx,dkz;
    fft.get("d1",dkh);
    fft.get("d2",dkx);
    fft.get("d3",dkz);
    
    float kh0,kx0,kz0;
    fft.get("o1",kh0); if (nkh==1) kh0=0.0;
    fft.get("o2",kx0);
    fft.get("o3",kz0);

    float kh, kx, kz, kzmin;

    std::valarray<float> k(n);
    for (int ix=0; ix < nkx; ix++) {
	kx = kx0+ix*dkx;
	kx *= kx;
	for (int ih=0; ih < nkh; ih++) {
	    kh = kh0+ih*dkh;
	    kh *= kh;

	    kzmin = sqrt(kh*kx);

	    for (int iz=0; iz < nkz; iz++) {
		kz = kz0+iz*dkz;
		kz *= kz;

		kz = SF_MAX(kz,kzmin);
		k[ih+nkh*(ix+iz*nkx)] = SF_PI*sqrt(kx*kh/kz+kz+kx+kh);
	    }
	}
    }
    ks.resize(n);
    ks = k;

    vector<int> lidx, ridx;
    FltNumMat mid;

    iC( lowrank(m,n,sample,tol,npk,lidx,ridx,mid) );

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
    right << rdata;

    return 0;
}
