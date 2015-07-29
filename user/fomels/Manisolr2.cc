// Lowrank decomposition for 2-D anisotropic wave propagation. 
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

static std::valarray<float>  vx, vz, q, t, xtap, ktap;
static std::valarray<double> kx, kz;
static double dt;
bool ktaper, xtaper;

static int sample(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,0.0);
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
	    double x0 = kx[j];
	    double z0 = kz[j];
	    // rotation of coordinates
	    double x = x0*c+z0*s;
	    double z = z0*c-x0*s;

	    z = wz*z*z;
	    x = wx*x*x;
	    double r = x+z;
	    r = r+sqrt(r*r-qq*x*z);
	    r = sqrt(0.5*r);
	    r = 2*(cos(r*dt)-1);

	    if (xtaper) r *= xtap[i];
	    if (ktaper) r *= ktap[j];

	    res(a,b) = r;
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

    iRSF velz, velx("velx"), eta("eta"), theta("theta");

    int nz,nx;
    velz.get("n1",nz);
    velz.get("n2",nx);
    int m = nx*nz;

    vx.resize(m);
    vz.resize(m);
    q.resize(m);
    t.resize(m);

    velx >> vx;
    velz >> vz;
    eta >> q;
    theta >> t;

    /* from eta to q */
    for (int im=0; im < m; im++) {
	q[im] = 8*q[im]/(1.0+2*q[im]);
    }

    /* fram degrees to radians */
    for (int im=0; im < m; im++) {
	t[im] *= SF_PI/180.;
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

    par.get("xtaper",xtaper,false); // if taper in x
    par.get("ktaper",ktaper,false); // if taper in k

    if (xtaper) {
	iRSF fxtap("xtap");
	xtap.resize(m);
	fxtap >> xtap;
    }

    if (ktaper) {
	iRSF fktap("ktap");
	ktap.resize(n);
	fktap >> ktap;
    }

    vector<int> lidx, ridx;
    DblNumMat mid;

    iC( ddlowrank(m,n,sample,eps,npk,lidx,ridx,mid) );

    int m2=mid.m();
    int n2=mid.n();

    vector<int> midx(m), nidx(n);
    for (int k=0; k < m; k++) 
	midx[k] = k;
    for (int k=0; k < n; k++) 
	nidx[k] = k;    

    DblNumMat lmat(m,m2);
    iC ( sample(midx,lidx,lmat) );

    DblNumMat lmat2(m,n2);
    iC( ddgemm(1.0, lmat, mid, 0.0, lmat2) );

    double *ldat = lmat2.data();
    std::valarray<float> ldata(m*n2);
    for (int k=0; k < m*n2; k++) 
	ldata[k] = ldat[k];
    oRSF left("left");
    left.put("n1",m);
    left.put("n2",n2);
    left << ldata;

    DblNumMat rmat(n2,n);
    iC ( sample(ridx,nidx,rmat) );

    double *rdat = rmat.data();
    std::valarray<float> rdata(n2*n);    
    for (int k=0; k < n2*n; k++) 
	rdata[k] = rdat[k];
    oRSF right;
    right.put("n1",n2);
    right.put("n2",n);
    right << rdata;

    exit(0);
}
