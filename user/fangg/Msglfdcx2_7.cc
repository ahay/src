// 2D Lowrank FD coefficient of d/dx on staggered grid (optimized)

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

//FltNumVec vs; //c
//FltNumVec ks; //k
static std::valarray<float> vs;
static std::valarray<double> ks;
static std::valarray<double> kx, kz;
static float pi=SF_PI;
static float dt;

static float sinc(float x)
{
    if (fabs(x)<=SF_EPS) return 1.0;
    return sinf(x)/(x+SF_EPS);
}

int samplex(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,0.0);
    for(int a=0; a<nr; a++) {
	for(int b=0; b<nc; b++) {
	   res(a,b) = 2.0*pi*kx[cs[b]]*sinc(pi*vs[rs[a]]*ks[cs[b]]*dt);
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
    par.get("npk",npk,50); // maximum rank

    par.get("dt",dt); // time step

    int size;
    par.get("size",size,6); // stencil length 
    iRSF velf;
    oRSF outm;  // FD coefficient of d/dx
    oRSF fsx("sx");
    oRSF fsz("sz"); 
    float dx, dz;
    
    float wavnumcut;
    par.get("wavnumcut",wavnumcut,1.0); // wavenumber cut percentile

    int nz, nx;
    velf.get("n1",nz);
    velf.get("d1",dz);
    velf.get("n2",nx);
    velf.get("d2",dx);
    float dkz, dkx;
    dkx = 1.0/(dx*nx);
    dkz = 1.0/(dz*nz);

    int nxz = nx*nz;
    vs.resize(nxz);
    ks.resize(nxz);
    velf >> vs;

    int m = nxz;
    int n = nxz;

    int COUNT= 0;
    kx.resize(nxz);
    kz.resize(nxz);
    float kx1, kz1;
    float kx0 = -dkx*nx/2.0; 
    float kz0 = -dkz*nz/2.0; 
    //float a = nx/4.0*dkx; 
    //float b = nz/4.0*dkz; 
    float a = nx/3.0*dkx*wavnumcut; 
    float b = nz/3.0*dkz*wavnumcut;
    
    int i=0;
    float dkxz=dkx+dkz;
    int SMK=0;
    for (int ix=0; ix < nx; ix++) {
        kx1 = kx0+ix*dkx; 
        for (int iz=0; iz < nz; iz++) {
            kz1 = kz0+iz*dkz; 
            ks[i] = sqrtf(kx1*kx1+kz1*kz1);
	    if (((kx1/a)*(kx1/a)+(kz1/b)*(kz1/b))<1.0) COUNT++;
            if (ks[i] < (dkxz+0.00001)) SMK++;
            kx[i] = kx1;
            kz[i] = kz1;
            i++;
        }
    }

    vector<int> ksc(COUNT), smallk(SMK);
    sf_warning("COUNT=%d",COUNT);
    sf_warning("SMK=%d",SMK);
    int nk=0, mk=0; 
    i=0;
    for (int ix=0; ix < nx; ix++) {
        kx1 = kx0+ix*dkx; 
        for (int iz=0; iz < nz; iz++) {
            kz1 = kz0+iz*dkz; 
            if (((kx1/a)*(kx1/a)+(kz1/b)*(kz1/b))<1.0) {
               ksc[nk] = i;
               nk++;
            }
            if (ks[i] < (dkxz+0.00001)){ 
               smallk[mk] = i;
               mk++;
            }
            i++;
        }
    }
    sf_warning("nk=%d",nk);
    
    /*Low rank decomposition*/   

    vector<int> cidx, ridx;
    DblNumMat mid;
    iC( ddlowrank(m,n,samplex,eps,npk,cidx,ridx,mid) );
    
    sf_warning("cidx.size=%d", cidx.size());
    sf_warning("ridx.size=%d", ridx.size());

    DblNumMat M1(m,cidx.size());
    vector<int> rs(m);
    for(int k=0; k<m; k++) rs[k]=k;
    iC( samplex(rs,cidx,M1) );
    DblNumMat M2(ridx.size(),n);
    vector<int> cs(n);
    for(int k=0; k<n; k++) cs[k]=k;
    iC( samplex(ridx,cs,M2) );


    /*FD coefficient*/

    /* d/dx */
    int len=0;
    len = (size+1)/2; //x: -len ,... -1, 1, ... len;
    std::valarray<float> xtmp(2*len);
    std::valarray<float> ztmp(len);
    
    for (int ix = 0; ix<len; ix++) {
	xtmp[ix] = -(len-ix)+0.5;  // - (2l-1)/2
	xtmp[ix+len] = ix+0.5;     //   (2l-1)/2
	ztmp[ix] = ix;
    }
    cerr<<"len  = " <<len <<endl;
    cerr<<"xtmp = "; for (int ix=0; ix<2*len; ix++) cerr<<xtmp[ix]<<", "; cerr<<endl;
    cerr<<"ztmp = "; for (int ix=0; ix<len; ix++) cerr<<ztmp[ix]<<", ";   cerr<<endl;

    int gdc = 0;
    for (int ix=0; ix<2*len; ix++) {
	for (int iz=0; iz<len; iz++) { 
	    if ((ztmp[iz]>0 && xtmp[ix]<1.5 && xtmp[ix]>-1.5 )|| (xtmp[ix]>0 && ztmp[iz] == 0)) {
		if ( (xtmp[ix]*dx*xtmp[ix]*dx+ztmp[iz]*dz*ztmp[iz]*dz) < dx*dx*(xtmp[0]*xtmp[0]+0.0001)) {gdc++;
		}
	    }
	    if ( (xtmp[ix]==1.5||xtmp[ix]==-1.5) && ztmp[iz]==1) {gdc++;}
	}
    }

    nk =0;
    DblNumMat sx(gdc,1), sz(gdc,1);
    for (int ix=0; ix<2*len; ix++) {
	for (int iz=0; iz<len; iz++) { 
	    if ((ztmp[iz]>0 && xtmp[ix]<1.5 && xtmp[ix]>-1.5)|| (xtmp[ix]>0 && ztmp[iz] == 0)) {
		if ( (xtmp[ix]*dx*xtmp[ix]*dx+ztmp[iz]*dz*ztmp[iz]*dz) < dx*dx*(xtmp[0]*xtmp[0]+0.0001)) {
		    sz._data[nk] = ztmp[iz];
		    sx._data[nk] = xtmp[ix];
		    nk++;
		} 
	    }
	    if ( (xtmp[ix]==1.5||xtmp[ix]==-1.5) && ztmp[iz]==1) {
		sz._data[nk] = ztmp[iz];
		sx._data[nk] = xtmp[ix];
		nk++;
	    }
	}
    }
      
    cerr<<"[x,z]="; for(int k=0; k<sx._m; k++) cerr<<"["<<sx._data[k]<<","<<sz._data[k]<<"] ";    
    cerr<<endl;
    cerr<<"gdc "<<gdc<<" "<<endl;
    cerr<<"nk "<<nk<<" "<<endl;

    DblNumMat kxtmp(1,nxz); for(int k=0; k<nxz; k++) kxtmp._data[k]=kx[k];
    DblNumMat kztmp(1,nxz); for(int k=0; k<nxz; k++) kztmp._data[k]=kz[k];
    DblNumMat kxtmpc(1,COUNT); for(int k=0; k<COUNT; k++) kxtmpc._data[k]=kx[ksc[k]];
    DblNumMat kztmpc(1,COUNT); for(int k=0; k<COUNT; k++) kztmpc._data[k]=kz[ksc[k]];
    int LEN = sx._m;
    DblNumMat Bx(LEN,nxz), B(LEN,nxz);
    DblNumMat Bxc(LEN,COUNT), Bc(LEN,COUNT);
    iC(ddgemm(dz,sz,kztmp,0.0,B));
    iC(ddgemm(dx,sx,kxtmp,0.0,Bx));
    iC(ddgemm(dz,sz,kztmpc,0.0,Bc));
    iC(ddgemm(dx,sx,kxtmpc,0.0,Bxc));
    for(int k=0; k<B._m*B._n; k++) B._data[k]=sin(2.0*pi*(B._data[k]+Bx._data[k]));
    for(int k=0; k<Bc._m*Bc._n; k++) Bc._data[k]=sin(2.0*pi*(Bc._data[k]+Bxc._data[k]));
    DblNumMat IBc(COUNT,LEN);    iC( ddpinv(Bc, 1e-16, IBc) );
      
    DblNumMat coef(ridx.size(),LEN);
    DblNumMat M2c;
    iC( samplex(ridx,ksc,M2c) );
    iC(ddgemm(1.0,M2c,IBc,0.0,coef));
 
    DblNumMat G(nxz,LEN), tmpG(mid._m,LEN);
    iC(ddgemm(1.0,mid,coef,0.0,tmpG));
    iC(ddgemm(1.0,M1,tmpG,0.0,G));

    Bc.resize(LEN,SMK);

    for(int k=0; k<LEN; k++) {
       for (int j=0; j<SMK; j++) {
           Bc(k,j) =B(k,smallk[j]);
       }
    }
    DblNumMat tmpB(nxz,SMK), maxB(nxz,1);
    iC(ddgemm(1.0,G,Bc,0.0,tmpB));
    float tmpmax;
    for (int k=0; k<nxz; k++) {
        tmpmax=-9999999.0;
        for (int j=0; j<SMK; j++) {
            if (fabs(tmpB(k,j)) > tmpmax) tmpmax=fabs(tmpB(k,j));
        }
        maxB._data[k]=tmpmax;
    }
    i=0;
    for(int k=0; k<nxz; k++) if (maxB._data[k]>2.0) i++;
    sf_warning("i=%d",i);
    if (i>0) {
	for(int k=0; k<nxz; k++) { maxB._data[k]=2.0/fabs(maxB._data[k]); }
	for (int x=0; x<nxz; x++){
	    for (int k=0; k<LEN; k++){
		G(x,k) = G(x,k)*maxB._data[x];
	    }
	} 
    }

    std::valarray<float> fMlr(nxz*LEN);
    double *ldat = G.data();
    for (int k=0; k < nxz*LEN; k++) {
        fMlr[k] = ldat[k];
    } 
    outm.put("n1",nxz);
    outm.put("n2",LEN);
    outm << fMlr;

    fsx.put("n1", LEN);
    fsx.put("n2", 1);
    fsz.put("n1", LEN);
    fsz.put("n2", 1);

    std::valarray<float> fs(LEN);
    ldat = sx.data();
    for (int k=0; k < LEN; k++) {
         fs[k] = (ldat[k]+0.5);
    } 
    fsx << fs;
    ldat = sz.data();
    for (int k=0; k < LEN; k++) {
        fs[k] = ldat[k];
    } 
    fsz << fs;
    
   return 0;
}


