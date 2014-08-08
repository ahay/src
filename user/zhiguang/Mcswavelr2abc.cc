// Lowrank decomposition for 2-D P-SV wave propagation. 
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

static std::valarray<float>  vp, vs, e, d, t, p;
static std::valarray<float>  f, s;
static std::valarray<float> kx, kz;
static float dt, ct, cb, cl, cr;
static int nz, nx, nbt, nbb, nbl, nbr;
bool exact, half, pwave;

static int samplep(vector<int>& rs, vector<int>& cs, CpxNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,cpx(0.0f,0.0f));
    for(int a=0; a<nr; a++) {
	int i=rs[a];
	float pp = vp[i];
	float ee = e[i];
	float dd = d[i];
	float tt = t[i];
	float c = cos(tt);
	float s = sin(tt);
	
	for(int b=0; b<nc; b++) {
	    int j = cs[b];
	    float x0 = kx[j];
	    float z0 = kz[j];
	    // rotation of coordinates
	    float x = x0*c+z0*s;
	    float z = z0*c-x0*s;
		float k = sqrt(x*x+z*z);
			int iz=i%nz;
			int ix=i/nz;
			float hypk=hypot(z0, x0);
			float phf=1.0;
			if(iz < nbt)
				phf *= exp(-powf(ct*(nbt-iz)*(z0/hypk),2));
			else if(iz > nz-1-nbb)
				phf *= exp(-powf(cb*(iz-nz+1+nbb)*(z0/hypk),2));
			if(ix < nbl)
				phf *= exp(-powf(cl*(nbl-ix)*(x0/hypk),2));
			else if(ix>nx-1-nbr)
				phf *= exp(-powf(cr*(ix-nx+1+nbr)*(x0/hypk),2));

		x=x/k;
		z=z/k;
		x=x*x;
		z=z*z;

		float r = pp*(1+dd*x*z+ee*x*x)*k*dt;
		if(!half)
			res(a,b)=cpx(cos(r), sin(r))*phf;
		if(half)
			res(a,b)=cpx(cos(0.5*r), sin(0.5*r))*phf;
	}
    }
    return 0;
}

static int samplesv(vector<int>& rs, vector<int>& cs, CpxNumMat& res)
{
	int nr = rs.size();
	int nc = cs.size();
	res.resize(nr, nc);
	setvalue(res, cpx(0.0f,0.0f));

	for(int a=0; a<nr; a++){
		int i=rs[a];
		float vss=vs[i];
		float ss=s[i];
		float tt=t[i];
		float c=cos(tt);
		float s=sin(tt);

		for(int b=0; b<nc; b++){
			int j=cs[b];
			float x0=kx[j];
			float z0=kz[j];
			// rotation of coordinates
			float x=x0*c+z0*s;
			float z=z0*c-x0*s;
			float k=sqrt(x*x+z*z);
			int iz=i%nz;
			int ix=i/nz;
			float hypk=hypot(z0, x0);
			float phf=1.0;
			if(iz < nbt)
				phf *= exp(-powf(ct*(nbt-iz)*(z0/hypk),2));
			else if(iz > nz-1-nbb)
				phf *= exp(-powf(cb*(iz-nz+1+nbb)*(z0/hypk),2));
			if(ix < nbl)
				phf *= exp(-powf(cl*(nbl-ix)*(x0/hypk),2));
			else if(ix>nx-1-nbr)
				phf *= exp(-powf(cr*(ix-nx+1+nbr)*(x0/hypk),2));

			x=x/k;
			z=z/k;
			x=x*x;
			z=z*z;

			float r=vss*(1+ss*x*z)*k*dt;
			if(!half)
				res(a,b) = cpx(cos(r), sin(r))*phf;
			if(half)
				res(a,b) = cpx(cos(0.5*r), sin(0.5*r))*phf;
		}
	}
	return 0;
}

static int sample(vector<int>& rs, vector<int>& cs, CpxNumMat& res)
{
	int nr = rs.size();
	int nc = cs.size();
	res.resize(nr, nc);
	setvalue(res, cpx(0.0f,0.0f));
	for(int a=0; a<nr; a++){
		int i=rs[a];
		float ff = f[i];
		float pp = vp[i];
		float ee = e[i];
		float dd = d[i];
		float tt = t[i];
		float c = cos(tt);
		float s = sin(tt);

		for(int b=0; b<nc; b++){
			int j = cs[b];
			float x0 = kx[j];
			float z0 = kz[j];
			//rotation of coordinates
			float x=x0*c+z0*s;
			float z=z0*c-x0*s;
			float k=sqrt(x*x+z*z);
			int iz=i%nz;
			int ix=i/nz;
			float hypk=hypot(z0, x0);
			float phf=1.0;
			if(iz < nbt)
				phf *= exp(-powf(ct*(nbt-iz)*(z0/hypk),2));
			else if(iz > nz-1-nbb)
				phf *= exp(-powf(cb*(iz-nz+1+nbb)*(z0/hypk),2));
			if(ix < nbl)
				phf *= exp(-powf(cl*(nbl-ix)*(x0/hypk),2));
			else if(ix>nx-1-nbr)
				phf *= exp(-powf(cr*(ix-nx+1+nbr)*(x0/hypk),2));

			x=x/k;
			z=z/k;
			x=x*x;
			z=z*z;

			float r = 2.0*x/ff;
			r = max(0., 1.0+2*r*(2.0*dd*z-(z-x)*ee)+r*r*ee*ee);
			if(pwave){
				if(!half){
					r=max(0., 1.0+ee*x-(1.0-sqrt(r))*ff/2.0);
					r=pp*sqrt(r)*k*dt;
					res(a,b)=cpx(cos(r), sin(r))*phf;
				}else{
					r=1.0+ee*x-(1.0-sqrt(r))*ff/2.0;
					r=0.5*pp*sqrt(r)*k*dt;
					res(a,b)=cpx(cos(r), sin(r))*phf;
				}
			}else{
				if(!half){
					r=max(0., 1.0+ee*x-(1.0+sqrt(r))*ff/2.0);
					r=pp*sqrt(r)*k*dt;
					res(a,b)=cpx(cos(r), sin(r))*phf;
				}else{
					r=1.0+ee*x-(1.0+sqrt(r))*ff/2.0;
					r=0.5*pp*sqrt(r)*k*dt;
					res(a,b)=cpx(cos(r), sin(r))*phf;
				}
			}
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
    par.get("eps",eps,1.e-6); // tolerance

    int npk;
    par.get("npk",npk,30); // maximum rank

    par.get("dt",dt); // time step

	par.get("nbt",nbt,0);
	par.get("nbb",nbb,0);
	par.get("nbl",nbl,0);
	par.get("nbr",nbr,0);

	par.get("ct",ct,0.0);
	par.get("cb",cb,0.0);
	par.get("cl",cl,0.0);
	par.get("cr",cr,0.0);

	par.get("exact",exact,true); 
	// if y, use exact P-SV phase velocities; if n, use Thomsen's weak-anisotropy approximations.
	par.get("half", half, false);
	// if y, do half velocity approximation for zero-offset migration.
	par.get("pwave",pwave, true);
	// if y, yield left and right matrices for P-wave; else for SV-wave.

    iRSF vp0, vs0("vs0"), epsilon("epsilon"), delta("delta"), theta("theta"), phi("phi");

    vp0.get("n1",nz);
    vp0.get("n2",nx);
    int m = nx*nz;

    vp.resize(m);
    vs.resize(m);
	e.resize(m);
	d.resize(m);
	t.resize(m);
	p.resize(m);

	vp0 >> vp;
	vs0 >> vs;
	epsilon >> e;
	delta >> d;
	theta >> t;
	phi >> p;

	if(exact){
		f.resize(m);
		for(int im=0; im<m; im++)
			f[im]=1.0-vs[im]*vs[im]/vp[im]/vp[im];
	}

	if(!exact && !pwave){
		s.resize(m);
		for(int im=0; im<m; im++)
			s[im]=powf(vp[im]/vs[im], 2)*(e[im]-d[im]);
	}

    /* fram degrees to radians */
    for (int im=0; im < m; im++) {
	t[im] *= cos(p[im]);
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
		iC( lowrank(m,n,sample,eps,npk,lidx,ridx,mid) );
	}else{
		if(pwave){ 
			iC( lowrank(m,n,samplep,eps,npk,lidx,ridx,mid) );
		}else{
			iC( lowrank(m,n,samplesv,eps,npk,lidx,ridx,mid) );
		}
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
		iC ( sample(midx,lidx,lmat) );
	}else{
		if(pwave){
			iC ( samplep(midx,lidx,lmat) );
		}else{
			iC ( samplesv(midx,lidx,lmat) );
		}
	}

    CpxNumMat lmat2(m,n2);
    iC( zgemm(1.0, lmat, mid, 0.0, lmat2) );

    cpx *ldat = lmat2.data();
    std::valarray<sf_complex> ldata(m*n2);
    for (int k=0; k < m*n2; k++) 
	ldata[k] = sf_cmplx(real(ldat[k]), imag(ldat[k]));
    oRSF left;
	left.type(SF_COMPLEX);
    left.put("n1",m);
    left.put("n2",n2);
    left << ldata;

    CpxNumMat rmat(n2,n);
	if(exact){
		iC ( sample(ridx,nidx,rmat) );
	}else{
		if(pwave){
			iC ( samplep(ridx,nidx,rmat) );
		}else{
			iC ( samplesv(ridx,nidx,rmat) );
		}
	}

    cpx *rdat = rmat.data();
    std::valarray<sf_complex> rdata(n2*n);    
    for (int k=0; k < n2*n; k++) 
	rdata[k] = sf_cmplx(real(rdat[k]), imag(rdat[k]));
    oRSF right("right");
	right.type(SF_COMPLEX);
    right.put("n1",n2);
    right.put("n2",n);
    right << rdata;

    exit(0);
}
