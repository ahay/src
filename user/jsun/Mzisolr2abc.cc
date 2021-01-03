// Complex lowrank decomposition for 2-D isotropic wave propagation with absorbing boundaries. 
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

static std::valarray<float> vs;
static std::valarray<double> kx,kz;
static float dt,ct,cb,cl,cr;
static int nz,nx,nbt,nbb,nbl,nbr,abc;
static bool rev;

int sample(vector<int>& rs, vector<int>& cs, ZpxNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,zpx(0.0,0.0));
    for(int a=0; a<nr; a++) {
        int i = rs[a];
        int iz  = i%nz;
        int ix  = (int) i/nz;
        double c0 = vs[i];
	for(int b=0; b<nc; b++) {
	    int j = cs[b];
	    double k = hypot(kz[j],kx[j]);
            double phase = c0*k*dt; 
            zpx phf = zpx(cos(phase),sin(phase)); 
	    //if (rev) phase*=-1;
            //double phf = 1;
 
            /* absorbing boundary */
            if (abc==0) {
                if (iz < nbt)
                    phf *= exp(-pow(ct*(nbt-iz)*abs(kz[j]/k),2));
                else if (iz > nz-1-nbb)
                    phf *= exp(-pow(cb*(iz-nz+1+nbb)*abs(kz[j]/k),2));
                if (ix < nbl)
                    phf *= exp(-pow(cl*(nbl-ix)*abs(kx[j]/k),2));
                else if (ix > nx-1-nbr)
                    phf *= exp(-pow(cr*(ix-nx+1+nbr)*abs(kx[j]/k),2));
            } else {
                if (iz < nbt)
                    phf *= exp(-pow(ct*(nbt-iz),2));
                else if (iz > nz-1-nbb)
                    phf *= exp(-pow(cb*(iz-nz+1+nbb),2));
                if (ix < nbl)
                    phf *= exp(-pow(cl*(nbl-ix),2));
                else if (ix > nx-1-nbr)
                    phf *= exp(-pow(cr*(ix-nx+1+nbr),2));
            }
	    res(a,b) = (rev) ? conj(phf) : phf;
	    //res(a,b) = zpx(cos(phase),sin(phase))*phf; 
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

    par.get("abc",abc,0); /*absorbing mode: 0-> direction dependent; 1-> direction independent.*/

    par.get("rev",rev,false);

    iRSF vel;

    vel.get("n1",nz);
    vel.get("n2",nx);
    int m = nx*nz;
    std::valarray<float> vels(m);
    vel >> vels;
    vs.resize(m);
    vs = vels;

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
    for (int ikx=0; ikx < nkx; ikx++) {
	for (int ikz=0; ikz < nkz; ikz++) {
	    int ik = ikz+ikx*nkz;
	    kx[ik] = 2*SF_PI*(kx0+ikx*dkx);
	    kz[ik] = 2*SF_PI*(kz0+ikz*dkz);
	}
    }

    vector<int> lidx, ridx;
    ZpxNumMat mid;

    iC( ddlowrank(m,n,sample,(double)eps,npk,lidx,ridx,mid) );

    int n2=mid.n();
    int m2=mid.m();

    vector<int> midx(m), nidx(n);
    for (int k=0; k < m; k++) 
	midx[k] = k;
    for (int k=0; k < n; k++) 
	nidx[k] = k;    

    ZpxNumMat lmat(m,m2);
    iC ( sample(midx,lidx,lmat) );

    ZpxNumMat lmat2(m,n2);
    iC( zzgemm(1.0, lmat, mid, 0.0, lmat2) );

    zpx *ldat = lmat2.data();
    std::valarray<sf_complex> ldata(m*n2);
    for (int k=0; k < m*n2; k++) {
	ldata[k] = sf_cmplx(real(ldat[k]),imag(ldat[k]));
    }

    oRSF left("left");
    left.type(SF_COMPLEX);
    left.put("n1",m);
    left.put("n2",n2);
    left << ldata;

    ZpxNumMat rmat(n2,n);
    iC ( sample(ridx,nidx,rmat) );

    zpx *rdat = rmat.data();
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
