// Lowrank decomposition for 2-D SV-wave propagation. 
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

static std::valarray<float>  vp, vs, e, d, t;
static std::valarray<float>  vx, q, c11, c13, c33, c55;
static std::valarray<double> kx, kz;
static double dt;
bool exact, half, pwave;

static int sample(vector<int>& rs, vector<int>& cs, CpxNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,cpx(0.0f,0.0f));
    for(int a=0; a<nr; a++) {
	int i=rs[a];
	double wx = vx[i]*vx[i];
	double wz = vp[i]*vp[i];
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
	    r = r-sqrt(r*r-qq*x*z);
	    r = 0.5*sqrt(0.5*r)*dt;
		res(a,b)=cpx(cos(r), sin(r));
	}
    }
    return 0;
}

static int sample1(vector<int>& rs, vector<int>& cs, CpxNumMat& res)
{
	int nr = rs.size();
	int nc = cs.size();
	res.resize(nr, nc);
	setvalue(res, cpx(0.0f,0.0f));
	for(int a=0; a<nr; a++){
		int i=rs[a];
		double c15 = c11[i]+c55[i];
		double c35 = c33[i]+c55[i];
		double ic15 = c11[i]-c55[i];
		double ic35 = c33[i]-c55[i];
		double wc13 = c13[i];
		double tt = t[i];
		double c = cos(tt);
		double s = sin(tt);

		for(int b=0; b<nc; b++){
			int j = cs[b];
			double x0 = kx[j];
			double z0 = kz[j];
			//rotation of coordinates
			double x=x0*c+z0*s;
			double z=z0*c-x0*s;

			x=x*x;
			z=z*z;
			double r=ic15*x-ic35*z;
			if(pwave){
				r=c15*x+c35*z + sqrt(r*r +4*wc13*x*z);
			}else{
				r=c15*x+c35*z - sqrt(r*r + 4.0*wc13*x*z);
			}
			if(half){
				r = 0.5*sqrt(0.5*r)*dt;
			}else{
				r = sqrt(0.5*r)*dt;
			}
			res(a,b) = cpx(cos(r), sin(r));
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

	par.get("exact",exact,false); // if y, use exact SV-wave velocity; if n, use approximate velocity

	if(exact){
		par.get("half",half,false); // if y, do decomposition for zero-offset migration.
		par.get("pwave",pwave,true); // if y, p-wave phase velocity.
	}

    iRSF vp0, vs0("vs0"), epsilon("epsilon"), delta("delta"), theta("theta");

    int nz,nx;
    vp0.get("n1",nz);
    vp0.get("n2",nx);
    int m = nx*nz;

    vp.resize(m);
    vs.resize(m);
	e.resize(m);
	d.resize(m);
	t.resize(m);

	vp0 >> vp;
	vs0 >> vs;
	epsilon >> e;
	delta >> d;
	theta >> t;

	if(exact){
		c11.resize(m);
		c33.resize(m);
		c55.resize(m);
		c13.resize(m);

		for(int im=0; im<m; im++){
			c33[im] = vp[im]*vp[im];
			c55[im] = vs[im]*vs[im];
			c11[im] = c33[im]*(2.0*e[im]+1);
			float tmp = c33[im]-c55[im];
			c13[im] = tmp*tmp+2.0*d[im]*tmp;
		}
	}else{
		vx.resize(m);
		q.resize(m);

		for(int im=0; im<m; im++){
			vx[im] = vp[im]*sqrt(2.0*e[im]+1);
			q[im] = 8.0*(e[im]-d[im])/(1.0+2.0*e[im]);
		}
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

	if(exact){
		iC( lowrank(m,n,sample1,eps,npk,lidx,ridx,mid) );
	}else{
		iC( lowrank(m,n,sample,eps,npk,lidx,ridx,mid) );
	}

    int m2=mid.m();
    int n2=mid.n();

    vector<int> midx(m), nidx(n);
    for (int k=0; k < m; k++) 
	midx[k] = k;
    for (int k=0; k < n; k++) 
	nidx[k] = k;    

    CpxNumMat lmat(m,m2);
	if(exact){
		iC ( sample1(midx,lidx,lmat) );
	}else{
		iC ( sample(midx,lidx,lmat) );
	}

    CpxNumMat lmat2(m,n2);
    iC( zgemm(1.0, lmat, mid, 0.0, lmat2) );

    cpx *ldat = lmat2.data();
    std::valarray<sf_complex> ldata(m*n2);
    for (int k=0; k < m*n2; k++) 
	ldata[k] = sf_cmplx(real(ldat[k]), imag(ldat[k]));
    oRSF left("left");
	left.type(SF_COMPLEX);
    left.put("n1",m);
    left.put("n2",n2);
    left << ldata;

    CpxNumMat rmat(n2,n);
	if(exact){
		iC ( sample1(ridx,nidx,rmat) );
	}else{
		iC ( sample(ridx,nidx,rmat) );
	}

    cpx *rdat = rmat.data();
    std::valarray<sf_complex> rdata(n2*n);    
    for (int k=0; k < n2*n; k++) 
	rdata[k] = sf_cmplx(real(rdat[k]), imag(rdat[k]));
    oRSF right;
	right.type(SF_COMPLEX);
    right.put("n1",n2);
    right.put("n2",n);
    right << rdata;

    exit(0);
}
