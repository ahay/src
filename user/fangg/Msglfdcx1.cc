// 1D Lowrank FD coefficient of d/dx on staggered grid

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
	    res(a,b) = 2.0*SF_PI*ks[cs[b]]*sinc(SF_PI*vs[rs[a]]*fabs(ks[cs[b]])*dt);
	}
    }
    return 0;
}

int fdx16(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
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
    for(int a=0; a<nr; a++) {
	for(int b=0; b<nc; b++) {
           v = vs[rs[a]];
           k = ks[cs[b]]*2*SF_PI*v; 
	   res(a,b) = c1*sin(k*dx/2.0)+c2*sin(k*dx*3.0/2.0)+c3*sin(k*dx*5.0/2.0)+c4*sin(k*dx*7.0/2.0)+c5*sin(k*dx*9.0/2.0)+
	              c6*sin(k*dx*11.0/2.0)+c7*sin(k*dx*13.0/2.0)+c8*sin(k*dx*15.0/2.0);
	   res(a,b) = res(a,b)*2.0;// exp(ikx)-exp(-ikx) = 2*sin(kx)
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

    iRSF velf;
    oRSF outm, Mexactfile("Mexact"), Mlrfile("Mlr"), Mappfile("Mapp"), Mfdfile("Mfd");
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
    par.get("size",SIZE,6); // stencil length 
    outm.put("n1",N);
    outm.put("n2",SIZE);
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
    DblNumMat Mexact(N,N);
    iC( sample(rs,cs,Mexact) );
    //for(int k=0;k<m;k++) cerr<<"rs"<<k<<"="<<rs[k]<<"  "<<vs[rs[k]]<<"cs"<<k<<"="<<cs[k]<<"  "<<ks[cs[k]]<<endl;
    

    float dk2=dk/2.0;
    Mexactfile.put("n1",N);
    Mexactfile.put("n2",N);
    Mexactfile.put("d2",dk2);
    Mexactfile.put("o2",0);
    std::valarray<float> fMlr(N*N);
    double *ldat = Mexact.data();
    for (int k=0; k < N*N; k++) {
        fMlr[k] = ldat[k];
    } 
    Mexactfile << fMlr;
    
    Mlrfile.put("n1",N);
    Mlrfile.put("n2",N);
    Mlrfile.put("d2",dk2);
    Mlrfile.put("o2",0);
    Mappfile.put("n1",N);
    Mappfile.put("n2",N);
    Mappfile.put("d2",dk2);
    Mappfile.put("o2",0);
    Mfdfile.put("n1",N);
    Mfdfile.put("n2",N);
    Mfdfile.put("d2",dk2);
    Mfdfile.put("o2",0);
    
    DblNumMat Mlr(N,N);
    iC( fdx16(rs,cs,Mlr) );
    ldat = Mlr.data();
    for (int k=0; k < N*N; k++) {
        fMlr[k] = ldat[k];
    } 
    Mfdfile << fMlr;

    cerr<<mid._m<<" ";
    cerr<<endl;
    cerr<<mid._n<<" ";
    cerr<<endl;
    cerr<<M2._m<<" ";
    cerr<<endl;
    cerr<<M2._n<<" ";
    cerr<<endl;

    DblNumMat tmpM(mid._m,M2._n);
    iC(ddgemm(1.0,mid,M2,0.0,tmpM));
    iC(ddgemm(1.0,M1,tmpM,0.0,Mlr));
    ldat = Mlr.data();
    for (int k=0; k < N*N; k++) {
        fMlr[k] = ldat[k];
    } 
    Mlrfile << fMlr;

    /*FD coefficient*/

    /* d/dx */
    int len=0;
    len = (SIZE+1)/2; //x: -len ,... -1, 1, ... len;
    std::valarray<float> xtmp(2*len);
    for (int ix = 1; ix<=2*len; ix++) {
	xtmp[ix-1] = (2.0*ix-1.0)/2.0;     //   (2l-1)/2
    }
    cerr<<"len  = " <<len <<endl;
    cerr<<"xtmp = "; for (int ix=0; ix<2*len; ix++) cerr<<xtmp[ix]<<", "; cerr<<endl;
    cerr<<"dx = "<<dx<<endl;
    cerr<<"COUNT = "<<count<<endl;
    
    DblNumMat s(2*len,1); for(int k=0; k<2*len; k++) s._data[k]=xtmp[k];  
    DblNumMat ktmp(1,N); for(int k=0; k<N; k++) ktmp._data[k]=ks[k];
    DblNumMat ktmpc(1,count); for(int k=0; k<count; k++) ktmpc._data[k]=ks[ksc[k]];
    DblNumMat B(SIZE,N), Bc(SIZE,count);
    iC(ddgemm(2*SF_PI*dx,s,ktmp,0.0,B));
    iC(ddgemm(2*SF_PI*dx,s,ktmpc,0.0,Bc));
    for(int k=0; k<B._m*B._n; k++) B._data[k]=sin(B._data[k]);
    for(int k=0; k<Bc._m*Bc._n; k++) Bc._data[k]=sin(Bc._data[k]);
    DblNumMat IB(N,SIZE);    iC( ddpinv(B, 1e-16, IB) );
    DblNumMat IBc(count,SIZE);    iC( ddpinv(Bc, 1e-16, IBc) );
    DblNumMat coef(ridx.size(),SIZE);
    DblNumMat M2c;
    iC( sample(ridx,ksc,M2c) );
    
    iC(ddgemm(1.0,M2c,IBc,0.0,coef));
    
    DblNumMat G(N,2*len), tmpG(mid._m,2*len);
    iC(ddgemm(1.0,mid,coef,0.0,tmpG));
    iC(ddgemm(1.0,M1,tmpG,0.0,G));
    
    /*
    Bc.resize(2*len,1);
    for(int k=0; k<2*len; k++) Bc._data[k]=1.0; //B(k,0); 
    DblNumMat tmpB(N,1);
    iC(ddgemm(1.0,G,Bc,0.0,tmpB));
    for(int k=0; k<N; k++) {cerr<<"tmpB["<<k<<"]="<<tmpB._data[k]<<endl; tmpB._data[k]=2.0/tmpB._data[k];}
    for (int x=0; x<N; x++){
        for (int k=0; k<2*len; k++){
            G(x,k) = G(x,k)*tmpB._data[x];
        }
    }
    */

    iC(ddgemm(1.0,G,B,0.0,Mlr));
    ldat = Mlr.data();
    for (int k=0; k < N*N; k++) {
        fMlr[k] = ldat[k];
    } 
    Mappfile << fMlr;

    ldat = G.data();
    fMlr.resize(N*SIZE);
    for (int k=0; k < SIZE*N; k++) {
        fMlr[k] = ldat[k];
    } 
    outm << fMlr;
  //  */
    return 0;
}
