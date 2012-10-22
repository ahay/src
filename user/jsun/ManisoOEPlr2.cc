// Lowrank decomposition for 2-D anisotropic wave propagation using exact P phase velocity in homogeneous orthorhombic media. 
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

// static std::valarray<float>  vx, vz, vs, q, t;
static std::valarray<float> C11,C12,C13,C22,C23,C33,C44,C55,C66,eps2,del2,f2;
static std::valarray<double> kx, kz;
static double dt;

static int sample(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,0.0);
//    double itta,f,p1,p2,p3,r;
    for(int a=0; a<nr; a++) {
	int i=rs[a];
		double ep = eps2[i];
		double de = del2[i];
		double f   = f2[i];
		double vp0 = sqrt(C33[i])*1000; // convert from km/s to m/s
	//	double tt = t[i];
	//	double c = cos(tt);
	//	double s = sin(tt);

	for(int b=0; b<nc; b++) {
	    double x0 = kx[cs[b]];
	    double z0 = kz[cs[b]];
	    // rotation of coordinates
	    //	    double x = x0*c+z0*s;
	    //	    double z = z0*c-x0*s;
	    //	    z = wz*z*z;
	    //	    x = wx*x*x;

	    double p1 = 1+ep*x0*x0-f/2;
	    double p2 = 1+2*ep*x0*x0/f;
	    double p3 = 2*(ep-de)*4*x0*x0*z0*z0/f;
	    double r  = vp0*sqrt(p1+sqrt(p2*p2-p3)*f/2);
	    // double r  = vp0*sqrt(p1-sqrt(p2*p2-p3)*f/2); //SV wave exact phase velocity

	    res(a,b) = 2*(cos(r*dt)-1); 
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

    iRSF c11, c12("c12"), c13("c13"), c22("c22"), c23("c23"), c33("c33"), c44("c44"), c55("c55"), c66("c66") ;

    int nz,nx;
    c11.get("n1",nz);
    c11.get("n2",nx);
    int m = nx*nz;

    C11.resize(m);
    C12.resize(m);
    C13.resize(m);
    C22.resize(m);
    C23.resize(m);
    C33.resize(m);
    C44.resize(m);
    C55.resize(m);
    C66.resize(m);
    eps2.resize(m);
    del2.resize(m);
    f2.resize(m);

    c11 >> C11;
    c12 >> C12;
    c13 >> C13;
    c22 >> C22;
    c23 >> C23;
    c33 >> C33;
    c44 >> C44;
    c55 >> C55;
    c66 >> C66;

    /* Calculate epsilon2 */
    for (int im=0; im < m; im++) {
      eps2[im] = (C11[im]-C33[im])/(2*C33[im]);
    }
    
    /* Calculate delta2 */
    for (int im=0; im < m; im++) {
      del2[im] = ((C13[im]+C55[im])*(C13[im]+C55[im])-(C33[im]-C55[im])*(C33[im]-C55[im]))/(2*C33[im]*(C33[im]-C55[im]));
    }
    
    /* Calculate f2 */
    for (int im=0; im < m; im++) {
      f2[im] = 1-C55[im]/C33[im];
    }
    
    /* ******************************************************************************************************* */

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
