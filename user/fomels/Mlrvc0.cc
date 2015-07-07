// Lowrank decomposition for zero-offset time migration

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

static std::valarray<float>  v;
static std::valarray<sf_complex>  d;
static float w, k, v0;
static int nt,nx;
static int nw,nk;
static float w0,k0,wmin;
static float dw,dk;

static int sample(vector<int>& rs, vector<int>& cs, CpxNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,cpx(0.0f,0.0f));
    for(int a=0; a<nr; a++) { /* space index */
	int i = rs[a];
	float v2 = v[i]; 
	v2 = v2*v2-v0*v0;

	for(int b=0; b<nc; b++) { /* wavenumber index */
	    int j = cs[b];

	    int iw = j%nw; j /= nw; w = w0+iw*dw; 	    
	    int ik = j%nk;          k = k0+ik*dk; 
	    
	    float phi = k * k * v2 * SF_PI * 0.125 / SF_MAX(w,wmin);

	    res(a,b) = cpx(cos(phi),-sin(phi)) * cpx(crealf(d[j]),cimagf(d[j]));
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

    iRSF vel("vel");
    vel.get("n1",nt);
    vel.get("n2",nx);

    int m = nt*nx;

    std::valarray<float> vels(m);
    vel >> vels;    
    v.resize(m);
    v = vels;

    par.get("v0",v0,0.0); // minimum velocity
    
    iRSF fft;

    fft.get("n1",nw);
    fft.get("n2",nk);

    int n = nw*nk;

    std::valarray<sf_complex> dats(n);
    fft >> dats;    
    d.resize(n);
    d = dats;

    fft.get("d1",dw);
    fft.get("d2",dk);
    
    fft.get("o1",w0); 
    fft.get("o2",k0);

    par.get("fmin",wmin,dw); // minimum frequency

    vector<int> lidx, ridx;
    CpxNumMat mid;

    iC( lowrank(m,n,sample,tol,npk,lidx,ridx,mid) );

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

    std::valarray<sf_complex> ldata(m*n2);
    for (int k=0; k < m*n2; k++) 
	ldata[k] = sf_cmplx(real(ldat[k]),imag(ldat[k]));

    oRSF left("left");
    left.put("n1",m);
    left.put("n2",n2);
    left << ldata;

    CpxNumMat rmat(n2,n);
    iC ( sample(ridx,nidx,rmat) );
    cpx *rdat = rmat.data();
    
    std::valarray<sf_complex> rdata(n2*n);    
    for (int k=0; k < n2*n; k++) 
	rdata[k] = sf_cmplx(real(rdat[k]),imag(rdat[k]));

    oRSF right;
    right.put("n1",n2);
    right.put("n2",n);
    right << rdata;
    
    return 0;
}
