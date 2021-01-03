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

static std::valarray<float>  vx, vz, q, t;
static std::valarray<float> kx, kz;
static float dt,ct,cb,cl,cr;
static int nz,nx,nbt,nbb,nbl,nbr;

static int sample(vector<int>& rs, vector<int>& cs, CpxNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,cpx(0.0f,0.0f));
    for(int a=0; a<nr; a++) {
	int i=rs[a];
	float wx = vx[i]*vx[i];
	float wz = vz[i]*vz[i];
	float qq = q[i];
	float tt = t[i];
	float c = cos(tt);
	float s = sin(tt);
	
	for(int b=0; b<nc; b++) {
	    float x0 = kx[cs[b]];
	    float z0 = kz[cs[b]];
	    // rotation of coordinates
	    float x = x0*c+z0*s;
	    float z = z0*c-x0*s;
	    int iz  = rs[a] % nz;
	    int ix  = (int) rs[a]/nz;
	    z = wz*z*z;
	    x = wx*x*x;
	    float r = x+z;
	    r = r+sqrt(r*r-qq*x*z);
	    r = sqrt(0.5*r)*dt;
	    float hypk = hypot(z0,x0);
	    float phf=1.;
	    if (iz < nbt)
		phf *= exp(-powf(ct*(nbt-iz)*(z0/hypk),2));
	    else if (iz > nz-1-nbb)
		phf *= exp(-powf(cb*(iz-nz+1+nbb)*(z0/hypk),2));
	    if (ix < nbl)
		phf *= exp(-powf(cl*(nbl-ix)*(x0/hypk),2));
	    else if (ix > nx-1-nbr)
		phf *= exp(-powf(cr*(ix-nx+1+nbr)*(x0/hypk),2));
	    res(a,b) = cpx(cos(r),sin(r))*phf;
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

    par.get("nbt",nbt,0);
    par.get("nbb",nbb,0);
    par.get("nbl",nbl,0);
    par.get("nbr",nbr,0);

    par.get("ct",ct,0.0);
    par.get("cb",cb,0.0);
    par.get("cl",cl,0.0);
    par.get("cr",cr,0.0);


    iRSF velz, velx("velx"), eta("eta"), theta("theta");

//    int nz,nx;
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

    vector<int> lidx, ridx;
    CpxNumMat mid;

    iC( lowrank(m,n,sample,eps,npk,lidx,ridx,mid) );

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
    left.type(SF_COMPLEX);
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
    right.type(SF_COMPLEX);
    right.put("n1",n2);
    right.put("n2",n);
    right << rdata;

    exit(0);
}
