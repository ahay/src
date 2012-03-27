// Lowrank decomposition for prestack exploding reflector in v(z). 
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
static float dt, kx, kh, dkz;

int sample(vector<int>& rs, vector<int>& cs, FltNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,0.0f);
    for(int a=0; a<nr; a++) {
	for(int b=0; b<nc; b++) {
	    float kz = cs[b]? cs[b]*dkz: dkz;
	    float x = (kz*kz+kx*kx)/kz;
	    float h = (kz*kz+kh*kh)/kz;

	    res(a,b) = 2*(cos(SF_PI*vs[rs[a]]*sqrt(x*h)*dt)-1); 
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
    par.get("npk",npk,5); // maximum rank

    par.get("dt",dt); // time step

    iRSF vel;

    // read velocity
    int nz;
    vel.get("n1",nz);
    std::valarray<float> vels(nz);
    vel >> vels;
    vs.resize(nz);
    vs = vels;
    
    iRSF fft("fft");

    int nkz,nkh,nkx;
    fft.get("n1",nkz);
    fft.get("n2",nkh);
    fft.get("n3",nkx);

    float dkh,dkx;
    fft.get("d1",dkz);
    fft.get("d2",dkh);
    fft.get("d3",dkx);

    vector<int> lidx, ridx;
    FltNumMat mid;

    int nk = nkx*nkh;

    oRSF left("left");
    left.put("n1",nz);
    left.put("n2",npk);
    left.put("n3",nk);

    oRSF right("right");
    right.put("n1",npk);
    right.put("n2",nkz);
    right.put("n3",nk);

    oRSF size("size");
    size.put("n1",2);
    size.put("n2",nk);
    size.type(SF_INT);   

    vector<int> midx(nz), nidx(nkz);
    for (int k=0; k < nz; k++) 
	midx[k] = k;
    for (int k=0; k < nkz; k++) 
	nidx[k] = k;    

    std::valarray<float> ldata(nz*npk);
    std::valarray<float> rdata(npk*nkz);  
    std::valarray<int>   nsize(2); 

    for (int ix=0; ix < nkx; ix++) {
	kx = ix*dkx;
	for (int ih=0; ih < nkh; ih++) {
	    kh = ih*dkh;

	    iC( lowrank(nz,nkz,sample,eps,npk,lidx,ridx,mid) );

	    int m2=mid.m();
	    int n2=mid.n();

	    FltNumMat lmat(nz,m2);   
	    iC ( sample(midx,lidx,lmat) );

	    FltNumMat lmat2(nz,n2);
	    iC( dgemm(1.0, lmat, mid, 0.0, lmat2) );

	    float *ldat = lmat2.data();

	    for (int k=0; k < nz*n2; k++) 
		ldata[k] = ldat[k];
	    left << ldata;

	    FltNumMat rmat(n2,nkz);
	    iC ( sample(ridx,nidx,rmat) );
	    float *rdat = rmat.data();

	    for (int iz=0; iz < nkz; iz++) {
		for (int j=0; j < n2; j++) {
		    rdata[j+iz*npk] = rdat[j+iz*n2];
		}
	    }
	    right << rdata;

	    nsize[0] = m2;
	    nsize[1] = n2;
	    size << nsize;
	}
    }

    return 0;
}
