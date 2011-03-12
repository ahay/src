// Lowrank decomposition for 2-D prestack exploding reflector in V(x,z)
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
static float dt, kh, kx, kz, kzmin;
static int nh,nx,nz;
static int nkh,nkx,nkz;
static float kh0,kx0,kz0;
static float dkh,dkx,dkz;

int sample(vector<int>& rs, vector<int>& cs, DblNumMat& res)
// assuming dx=dh
{
    int ix, ih, iz;

    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,0.0);
    for(int a=0; a<nr; a++) { /* space index */
	int i = rs[a];

	ih = i%nh; i /= nh; 
	ix = i%nx; i /= nx; 
	iz = i%nz; 

	int is = SF_MIN(SF_MAX(0,ix-ih),nx-1);
	int ir = SF_MIN(SF_MAX(0,ix+ih),nx-1);
 
	float vs = v[is+iz*nx]; 
	float vr = v[ir+iz*nx]; 

	float vp2 = 0.5*(vs+vr);
	vp2 *= vp2;

	float vsr = vs*vr;
	vs *= vs;
	vr *= vr;

	float vp = vs+vr;
	float vm = vs-vr;
    
	for(int b=0; b<nc; b++) { /* wavenumber index */
	    int i = cs[b];

	    ih = i%nkh; i /= nkh; kh = kh0+ih*dkh; 	    
	    ix = i%nkx; i /= nkx; kx = kx0+ix*dkx; 
	    iz = i%nkz; i /= nkz; kz = kz0+iz*dkz; 
	    
	    float kp = kx+kh; 
	    float km = kx-kh;
	    float kxh = kx*kh;
	    
	    kp *= kp;
	    km *= km;
	    kz *= kz;
	    kh *= kh;
	    kx *= kx;
	    
	    kzmin = sqrt(kh*kx);
	    /* kz = SF_MAX(kz,kzmin); */
	    
	    float phi = (kh + kz)*(kx + kz)*vr*vs/
		(-kxh*vm + kz*vp + 
		 sqrt(fabs(kz*(2*(kh + kx + 2*kz)*vr*vs - km*vs*vs - kp*vr*vr)))); // exact
	    
	    phi = (kh + kz)*(kx + kz)*vs*vr/(4*kz*vp2); // v(z)-like

	    phi = (vsr*((kh + kx + 2*kxh)*vr + (kh + kx - 2*kxh)*vs + 2*(kh + kx + 2*kz)*vsr + 
			vsr*sqrt(((kh + kx - 2*kxh)/vr - (kh + kx + 2*kxh)/vs)*((kh + kx - 2*kxh)/vr - (kh + kx + 2*kxh)/vs)*vsr*
				 (vr + vs + 2*vsr) + ((kh + kx + 2*kxh)*vr + (kh + kx - 2*kxh)*vs + 
						      2*(kh + kx + 2*kz)*vsr)*((kh + kx + 2*kxh)*vr + (kh + kx - 2*kxh)*vs + 
									       2*(kh + kx + 2*kz)*vsr)/(vr*vs))))/(8.*(vr + vs + 2*vsr)); // tariq's
	    
	    phi = ((kh + kx + kz+sqrt(4*kx*kh+(kh+kx+kz)*(kh+kx+kz)))*vs*vr)/(8.*vp2); // tariq's v(z)-like

	    res(a,b) = 2*(cos(2*SF_PI*sqrt(phi)*dt)-1); 
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
    vel.get("n1",nx);
    vel.get("n2",nz);

    int nxz = nx*nz;

    std::valarray<float> vels(nxz);
    vel >> vels;    
    v.resize(nxz);
    v = vels;

    par.get("nh",nh); /* number of offsets */
    int m = nxz*nh;
    
    iRSF fft("fft");

    fft.get("n1",nkh);
    fft.get("n2",nkx);
    fft.get("n3",nkz);

    int n = nkh*nkx*nkz;

    fft.get("d1",dkh);
    fft.get("d2",dkx);
    fft.get("d3",dkz);
    
    fft.get("o1",kh0); if (nkh==1) kh0=0.0;
    fft.get("o2",kx0);
    fft.get("o3",kz0);

    vector<int> lidx, ridx;
    DblNumMat mid;

    iC( lowrank(m,n,sample,(double) tol,npk,lidx,ridx,mid) );

    int m2=mid.m();
    int n2=mid.n();
    double *dmid = mid.data();

    std::valarray<float> fmid(m2*n2);
    for (int k=0; k < m2*n2; k++) {
	fmid[k] = dmid[k];
    }

    oRSF middle;
    middle.put("n1",m2);
    middle.put("n2",n2);
    middle << fmid;

    vector<int> midx(m), nidx(n);
    for (int k=0; k < m; k++) 
	midx[k] = k;
    for (int k=0; k < n; k++) 
	nidx[k] = k;    

    DblNumMat lmat(m,m2);
    iC ( sample(midx,lidx,lmat) );
    double *ldat = lmat.data();

    std::valarray<float> ldata(m*m2);
    for (int k=0; k < m*m2; k++) 
	ldata[k] = ldat[k];
    oRSF left("left");
    left.put("n1",m);
    left.put("n2",m2);
    left << ldata;

    DblNumMat rmat(n2,n);
    iC ( sample(ridx,nidx,rmat) );
    double *rdat = rmat.data();

    std::valarray<float> rdata(n2*n);    
    for (int k=0; k < n2*n; k++) 
	rdata[k] = rdat[k];
    oRSF right("right");
    right.put("n1",n2);
    right.put("n2",n);
    right << rdata;

    return 0;
}
