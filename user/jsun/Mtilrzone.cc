// Lowrank decomposition for 2-D anisotropic wave propagation using exact phase velocity (2 step time marching). 
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
#include <assert.h>

#include <rsf.hh>

#include "vecmatop.hh"
#include "serialize.hh"

using namespace std;

static std::valarray<float> C11,C13,C33,C55,Theta;
static std::valarray<float> kx, kz;
static float dt;
static int mode;
static bool approx;

static int sample(vector<int>& rs, vector<int>& cs, FltNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,0.0f);
    for(int a=0; a<nr; a++) {
	int i=rs[a];
	float c11 = C11[i];
	float c13 = C13[i];
	float c33 = C33[i];
	float c55 = C55[i];
	float c = cosf(Theta[i]);
	float s = sinf(Theta[i]);
	for(int b=0; b<nc; b++) {
	    int j=cs[b];
	    float x0 = kx[j];
	    float z0 = kz[j];
	    // rotation of coordinates
	    float x = x0*c+z0*s;
	    float z = z0*c-x0*s;
	    float r;
	    if (approx) {
		float n11= x*x; // n1=x
		float n33= z*z; // n3=z
		float e  = c11*n11 + c33*n33;
		float qv = (powf((c13+c55),2) + c55*(c33-c55))/(c11*(c33-c55));
		float qh = (powf((c13+c55),2) + c55*(c11-c55))/(c33*(c11-c55));
//		float qm = qh*n11 + qv*n33;
		float qm;
		if(x==0 && z==0) qm=qv;
		else qm = qh*(x/hypotf(x,z))*(x/hypotf(x,z)) + qv*(z/hypotf(x,z))*(z/hypotf(x,z));
		float sv = (c11*(c33-c11)*(qh-1)*(qv-1)*(qv-1))/(2*(c33*(1-qh)+c11*(qv-1))*(c33*(qv-qh)+c11*(qh+qv-qv*qv-1)));
		float sh = (c33*(c33-c11)*(qh-1)*(qh-1)*(qv-1))/(2*(c33*(1-qh)+c11*(qv-1))*(c11*(qh-qv)+c33*(qh+qv-qh*qh-1)));
//		float sm = sh*n11 + sv*n33;
		float sm;
		if(x==0 && z==0) sm=0.5;
		else sm=sh*(x/hypotf(x,z))*(x/hypotf(x,z)) + sv*(z/hypotf(x,z))*(z/hypotf(x,z));
//		if (sm==0) sf_warning("sm-0! sh=%g, sv=%g, n11=%g, n33=%g, theta=%g, x0=%g, z0=%g, c=%g, s=%g",sh,sv,n11,n33,Theta[i],x0,z0,c,s);
		r  = sqrt(e*(1-sm) + sm*sqrt(e*e + 2*(qm-1)*c11*c33*n11*n33/sm));
		if (mode==1) sf_error("No approximation for Sv wave!");
	    } else {
		float p1 = (c11+c55)*x*x + (c33+c55)*z*z;
		float p2 = (c11-c55)*x*x - (c33-c55)*z*z;
		float p3 = 4*(c13+c55)*(c13+c55)*x*x*z*z;
		r  = sqrt(0.5*p1 + 0.5*sqrt(p2*p2+p3));
		if (mode==1) r = sqrt(0.5*p1 - 0.5*sqrt(p2*p2+p3));
	    }
	    res(a,b) = 2*(cos(r*dt)-1); 
	}
    }
    return 0;
}

int main(int argc, char** argv)
{   
    sf_init(argc,argv); // Initialize RSF

    iRSF par(0);

    par.get("mode",mode,0); // wave mode (0=p wave, 1=Sv wave)

    int seed;
    par.get("seed",seed,time(NULL)); // seed for random number generator
    srand48(seed);

    float eps;
    par.get("eps",eps,1.e-4); // tolerance

    int npk;
    par.get("npk",npk,20); // maximum rank

    par.get("dt",dt); // time step

    par.get("approx",approx,false); // whether to use zone's approximation
    if (approx) sf_warning("Using zone's approximation!");

    iRSF c11, c13("c13"), c33("c33"), c55("c55"), theta("theta");

    int nz,nx;
    c11.get("n1",nz);
    c11.get("n2",nx);
    int m = nx*nz;

    C11.resize(m);
    C13.resize(m);
    C33.resize(m);
    C55.resize(m);
    Theta.resize(m);

    c11 >> C11;
    c13 >> C13;
    c33 >> C33;
    c55 >> C55;
    theta >> Theta;

    /* from degrees to radians */
    for (int im=0; im < m; im++) {
	Theta[im] *= SF_PI/180.;
    }
    
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

    int n = nkx*nkz;
    kx.resize(n);
    kz.resize(n);
    int i = 0;
    for (int ix=0; ix < nkx; ix++) {
	for (int iz=0; iz < nkz; iz++) {
	    kx[i] = 2*SF_PI*(kx0+ix*dkx);
	    kz[i] = 2*SF_PI*(kz0+iz*dkz);
	    i++;
	}
    }

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

    exit(0);
}
