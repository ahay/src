// Complex lowrank decomposition for 2-D isotropic wave propagation. 
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

static std::valarray<float> vs,qs;
static std::valarray<float> ksz,ksx;
static float ct,cb,cl,cr;
static int nkzs,nz,nx,nbt,nbb,nbl,nbr;
static float dt,c0,w0;
static bool rev;
static int mode,sign;

int sample(vector<int>& rs, vector<int>& cs, CpxNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,cpx(0.0f,0.0f));
    for(int a=0; a<nr; a++) {
	for(int b=0; b<nc; b++) {
	    int ikz = cs[b] % nkzs;
	    int ikx = (int) cs[b]/nkzs;
	    int iz  = rs[a] % nz;
	    int ix  = (int) rs[a]/nz;
	    float hypk = hypot(ksz[ikz],ksx[ikx]);
	    cpx phf;
	    if (mode == 0) { /*viscoacoustic*/
		float gamma = atanf(1./qs[rs[a]])/SF_PI;
		float eta = -powf(c0,2.*gamma)*powf(w0,-2.*gamma)*cosf(SF_PI*gamma);
		float tao = -powf(c0,2.*gamma-1.)*powf(w0,-2.*gamma)*sinf(SF_PI*gamma);
		float p1  = tao*powf(vs[rs[a]],2)*powf(hypk,2.*gamma+1.);
		float p2  = -powf(p1,2) - 4*eta*powf(vs[rs[a]],2)*powf(hypk,2.*gamma+2.);
		if (p2 < 0) sf_warning("square root is imaginary!!!");
		sf_complex phr = sf_cmplx(p1*dt/2.,0);
		sf_complex phi = sf_cmplx(0,0);
		sf_complex phase = sf_cmplx(0,0);
#ifdef SF_HAS_COMPLEX_H
		phi = (sign==0)? sf_cmplx(0,1)*csqrtf(sf_cmplx(p2,0))*dt/2. : sf_complex(0,1)*(-1*csqrtf(sf_complex(p2,0))*dt/2.);
		phase = phr + phi;
		if (rev) phase*=-1;
#else
		if (sign==0)
		    phi = sf_cmul(sf_cmplx(0,1),sf_crmul(csqrtf(sf_cmplx(p2,0)),dt/2.));
		else
		    phi = sf_cmul(sf_cmplx(0,1),sf_crmul(csqrtf(sf_cmplx(p2,0)),-1*dt/2.));
		phase = sf_cadd(phr,phi);
		if (rev) phase=sf_crmul(phase,-1);
#endif
		phf = cpx(crealf(cexpf(phase)),cimagf(cexpf(phase)));
	    }
	    else if (mode == 1) { /*loss dominated*/
		float gamma = atanf(1./qs[rs[a]])/SF_PI;
		float tao = -powf(c0,2.*gamma-1.)*powf(w0,-2.*gamma)*sinf(SF_PI*gamma);
		float p1  = tao*powf(vs[rs[a]],2)*powf(hypk,2.*gamma+1.);
		float p2  = -powf(p1,2) + 4*powf(vs[rs[a]],2)*powf(hypk,2);
		if (p2 < 0) sf_warning("square root is imaginary!!!");
		sf_complex phr = sf_cmplx(p1*dt/2.,0);
		sf_complex phi = sf_cmplx(0,0);
		sf_complex phase = sf_cmplx(0,0);
#ifdef SF_HAS_COMPLEX_H
		phi = (sign==0)? sf_cmplx(0,1)*csqrtf(sf_cmplx(p2,0))*dt/2. : sf_complex(0,1)*(-1*csqrtf(sf_complex(p2,0))*dt/2.);
		phase = phr + phi;
		if (rev) phase*=-1;
#else
		if (sign==0)
		    phi = sf_cmul(sf_cmplx(0,1),sf_crmul(csqrtf(sf_cmplx(p2,0)),dt/2.));
		else
		    phi = sf_cmul(sf_cmplx(0,1),sf_crmul(csqrtf(sf_cmplx(p2,0)),-1*dt/2.));
		phase = sf_cadd(phr,phi);
		if (rev) phase=sf_crmul(phase,-1);
#endif
		phf = cpx(crealf(cexpf(phase)),cimagf(cexpf(phase)));
	    }
	    else if (mode == 2) { /*dispersion-dominated*/
		float gamma = atanf(1./qs[rs[a]])/SF_PI;
		float eta = -powf(c0,2.*gamma)*powf(w0,-2.*gamma)*cosf(SF_PI*gamma);
		float phase = sqrtf(-eta*powf(vs[rs[a]],2)*powf(hypk,2.*gamma+2.))*dt;
		if (rev) phase*=-1;
		phf = cpx(cosf(phase),sinf(phase));
	    }
	    else { /*acoustic*/
		float phase = vs[rs[a]]*hypk*dt; 
		if (rev) phase*=-1;
		phf = cpx(cosf(phase),sinf(phase)); 
	    }
	    /* absorbing boundary */
	    if (iz < nbt)
		phf *= exp(-powf(ct*(nbt-iz)*abs(ksz[ikz]/hypk),2));
	    else if (iz > nz-1-nbb)
		phf *= exp(-powf(cb*(iz-nz+1+nbb)*abs(ksz[ikz]/hypk),2));
	    if (ix < nbl)
		phf *= exp(-powf(cl*(nbl-ix)*abs(ksx[ikx]/hypk),2));
	    else if (ix > nx-1-nbr)
		phf *= exp(-powf(cr*(ix-nx+1+nbr)*abs(ksx[ikx]/hypk),2));

	    res(a,b) = phf;
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
    par.get("c0",c0); // reference velocity
    par.get("w0",w0); // reference frequency
    
    par.get("rev",rev,false); // reverse propagation
    par.get("mode",mode,0); // mode of propagation: 0 is viscoacoustic (default); 1 is loss-dominated; 2 is dispersion dominated; 3 is acoustic
    par.get("sign",sign,0); // sign of solution: 0 is positive, 1 is negative
    
    par.get("nbt",nbt,0);
    par.get("nbb",nbb,0);
    par.get("nbl",nbl,0);
    par.get("nbr",nbr,0);

    par.get("ct",ct,0.0);
    par.get("cb",cb,0.0);
    par.get("cl",cl,0.0);
    par.get("cr",cr,0.0);

    iRSF vel, q("q");

    vel.get("n1",nz);
    vel.get("n2",nx);
    int m = nx*nz;

    vs.resize(m);
    qs.resize(m);

    vel >> vs;
    q >> qs;

    iRSF fft("fft");

    int nkz,nkx;
    fft.get("n1",nkz);
    fft.get("n2",nkx);
    nkzs = nkz;

    float dkz,dkx;
    fft.get("d1",dkz);
    fft.get("d2",dkx);
    
    float kz0,kx0;
    fft.get("o1",kz0);
    fft.get("o2",kx0);

    int n = nkx*nkz;

    std::valarray<float> kz(nkz),kx(nkx);
    for (int iz=0; iz < nkz; iz++)
	kz[iz] = 2*SF_PI*(kz0+iz*dkz);
    for (int ix=0; ix < nkx; ix++)
	kx[ix] = 2*SF_PI*(kx0+ix*dkx);

    ksz.resize(nkz);
    ksz = kz;
    ksx.resize(nkx);
    ksx = kx;

    vector<int> lidx, ridx;
    CpxNumMat mid;

    iC( lowrank(m,n,sample,eps,npk,lidx,ridx,mid) );

    int n2=mid.n();
    int m2=mid.m();

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
    for (int k=0; k < m*n2; k++) {
	ldata[k] = sf_cmplx(real(ldat[k]),imag(ldat[k]));
    }

    oRSF left("left");
    left.type(SF_COMPLEX);
    left.put("n1",m);
    left.put("n2",n2);
    left << ldata;

    CpxNumMat rmat(n2,n);
    iC ( sample(ridx,nidx,rmat) );

    cpx *rdat = rmat.data();
    std::valarray<sf_complex> rdata(n2*n);    
    for (int k=0; k < n2*n; k++) {
	rdata[k] = sf_cmplx(real(rdat[k]),imag(rdat[k]));
    }
    oRSF right;
    right.type(SF_COMPLEX);
    right.put("n1",n2);
    right.put("n2",n);
    right << rdata;

    exit(0);
}
