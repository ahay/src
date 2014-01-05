//  1D SG Lowrank FD and 1D FD coefficient(4-th 8-th 16-th 32-th) 

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
static std::valarray<double> ks;
static float pi=SF_PI;
static float dt,dx;

static float sinc(float x)
{
    if (fabs(x)<=SF_EPS) return 1.0;
    return sinf(x)/(x+SF_EPS);
}

int sample(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,0.0);
    for(int a=0; a<nr; a++) {
	for(int b=0; b<nc; b++) {
	    res(a,b) = 2.0*pi*ks[cs[b]]*sinc(pi*vs[rs[a]]*fabs(ks[cs[b]])*dt);
	}
    }
    return 0;
}

int fdx16(vector<int>& cs, DblNumMat& res)
{
    int nc = cs.size();
    res.resize(nc,1);  
    setvalue(res,0.0);
    float c1 = 1.2597508/dx;
    float c2 = -0.1284399/dx;
    float c3 = 0.0387632/dx;
    float c4 = -0.014999/dx;
    float c5 = 0.00610876/dx;
    float c6 = -0.0023589/dx;
    float c7 = 0.0007711/dx;
    float c8 = -0.0001545/dx;
    float v, k;

    for(int b=0; b<nc; b++) {
	v = vs[0];
	k = ks[cs[b]]*2*pi; 
	res(b,0) = c1*sin(k*dx/2.0)+c2*sin(k*dx*3.0/2.0)+c3*sin(k*dx*5.0/2.0)+c4*sin(k*dx*7.0/2.0)+c5*sin(k*dx*9.0/2.0)+
	    c6*sin(k*dx*11.0/2.0)+c7*sin(k*dx*13.0/2.0)+c8*sin(k*dx*15.0/2.0);
	res(b,0) = asin(dt*v*res(b,0))/(0.5*k*dt*v);// exp(ikx)-exp(-ikx) = 2*sin(kx)
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

    iRSF velf;
    oRSF outm;
    oRSF Mfdfile("Mfd");
    oRSF Mlrfile("Mlr");
    //oRSF outm, Mexactfile("Mexact"), Mlrfile("Mlr");

    int N;
    velf.get("n1",N);
    float dk;
    velf.get("d1",dx);
    dk = 1.0/(dx*N);

    vs.resize(N);
    ks.resize(N);
    velf >> vs;
    
    int SIZE;
    par.get("size",SIZE,16); // stencil length 
    //outm.put("n2",N);
    int m = N;
    int n = N;

    int count = 0;
    float CUT = N/3.0*dk;
    for (int k=0; k < N; k++) {
	ks[k] = -dk*N/2.0+k*dk;
        if (fabs(ks[k]) < CUT) count++;
    }
    vector<int> ksc(count);
    int nk = 0;
    for (int k=0; k < N/2; k++) {
	ks[k] = k*dk;
	if (fabs(ks[k]) < CUT) {
           ksc[nk] = k;
           nk++;
        }
    }
    for (int k=N/2; k < N; k++) {
	ks[k] = (-N+k)*dk;
        if (fabs(ks[k]) < CUT) { 
           ksc[nk] = k;
           nk++;
        }
    }

    /*Lowrank FD coefficients*/
    vector<int> cidx, ridx;
    DblNumMat mid;
    iC( ddlowrank(m,n,sample,(double) eps,npk,cidx,ridx,mid) );
    
    DblNumMat M1(m,cidx.size());
    vector<int> rs(m);
    for(int k=0; k<m; k++) rs[k]=k;
    iC( sample(rs,cidx,M1) );
    DblNumMat M2(ridx.size(),n);
    vector<int> cs(n);
    for(int k=0; k<n; k++) cs[k]=k;
    iC( sample(ridx,cs,M2) );

    float dk2=dk*2*pi;
    
    /* d/dx */
    int len=0;
    len = (SIZE+1)/2; //x: -len ,... -1, 1, ... len;
    std::valarray<float> xtmp(len);
    for (int ix = 1; ix<=len; ix++) {
	xtmp[ix-1] = (2.0*ix-1.0)/2.0;     //   (2l-1)/2
    }
    
    cerr<<"len  = " <<len <<endl;
    cerr<<"xtmp = "; for (int ix=0; ix<len; ix++) cerr<<xtmp[ix]<<", "; cerr<<endl;
    cerr<<"dx = "<<dx<<endl;
    cerr<<"COUNT = "<<count<<endl;
    
    DblNumMat s(len,1); 
    for(int k=0; k<len; k++) 
	s._data[k]=xtmp[k];
    
    DblNumMat ktmp(1,N); for(int k=0; k<N; k++) ktmp._data[k]=ks[k];
    DblNumMat ktmpc(1,count); for(int k=0; k<count; k++) ktmpc._data[k]=ks[ksc[k]];
    DblNumMat B(len,N), Bc(len,count);
    iC(ddgemm(2*pi*dx,s,ktmp,0.0,B));
    iC(ddgemm(2*pi*dx,s,ktmpc,0.0,Bc));
    for(int k=0; k<B._m*B._n; k++) B._data[k]=sin(B._data[k]);
    for(int k=0; k<Bc._m*Bc._n; k++) Bc._data[k]=sin(Bc._data[k]);
    DblNumMat IB(N,len);    iC( ddpinv(B, 1e-16, IB) );
    DblNumMat IBc(count,len);    iC( ddpinv(Bc, 1e-16, IBc) );
    DblNumMat coef(ridx.size(),len);
    DblNumMat M2c;
    iC( sample(ridx,ksc,M2c) );
    
    iC(ddgemm(1.0,M2c,IBc,0.0,coef));
        
    DblNumMat G(N,len), tmpG(mid._m,len);
    iC(ddgemm(1.0,mid,coef,0.0,tmpG));
    iC(ddgemm(1.0,M1,tmpG,0.0,G));

    /*Bc.resize(len,1);
    for(int k=0; k<len; k++) Bc._data[k]=1.0;//B(k,0);
    DblNumMat tmpB(N,1);
    iC(ddgemm(1.0,G,Bc,0.0,tmpB));
    for(int k=0; k<N; k++) tmpB._data[k]=2.0/tmpB._data[k];
    for (int x=0; x<N; x++){
        for (int k=0; k<len; k++){
            G(x,k) = G(x,k)*tmpB._data[x];
        }
	}*/
   
    /*   Bc.resize(2*len,1);
    for(int k=0; k<2*len; k++) Bc._data[k]=1.0; //B(k,0); 
    DblNumMat tmpB(N,1);
    iC(ddgemm(1.0,G,Bc,0.0,tmpB));
    for(int k=0; k<N; k++) {cerr<<"tmpB["<<k<<"]="<<tmpB._data[k]<<endl; tmpB._data[k]=2.0/tmpB._data[k];}
    for (int x=0; x<N; x++){
        for (int k=0; k<2*len; k++){
            G(x,k) = G(x,k)*tmpB._data[x];
        }
	}*/
   
    outm.put("n1",N);
    outm.put("n2",len);
    
    std::valarray<float> fG(N*len);
    double *ldat = G.data();
    fG.resize(N*len);
    for (int k=0; k < len*N; k++) {
        fG[k] = ldat[k];
    } 
    outm << fG;

    /*dispersion velocity*/
    
    Mlrfile.put("n1",N);
    Mlrfile.put("n2",1);
    Mlrfile.put("d2",dk2);
    Mlrfile.put("o2",0);
    Mfdfile.put("n1",N);
    Mfdfile.put("n2",1);
    Mfdfile.put("d2",dk2);
    Mfdfile.put("o2",0);

    DblNumMat Mlr(N,N);
    iC(ddgemm(1.0,G,B,0.0,Mlr));
    
    std::valarray<float> fMlr(N);
    //ldat = Mlr.data();
    double tg; 
    fMlr[0] = 1.0; 
    for (int k=1; k < N; k++) {
        tg = 0.0;
        for (int l=1; l <= len; l++) {
            tg += G(0,l-1)*sin((2*l-1)*dx*(ks[k])*pi);
        }
        fMlr[k] = asin( dt*vs[0]*tg/2.0 )/((ks[k])*pi*dt*vs[0]); 
    }
    Mlrfile << fMlr;

    
    Mlr.resize(N,1);
    iC( fdx16(cs,Mlr) );
    ldat = Mlr.data();
    for (int k=0; k < N; k++) {
        fMlr[k] = ldat[k];
    } 
    Mfdfile << fMlr;
    return 0;

}
