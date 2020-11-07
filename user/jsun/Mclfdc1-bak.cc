// 1D 10th-order Lowrank Onestep FD coefficient

//   Copyright (C) 2015 University of Texas at Austin
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
static float dt,dx;

int sample(vector<int>& rs, vector<int>& cs, ZpxNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,zpx(0.,0.));
    for(int a=0; a<nr; a++) {
	for(int b=0; b<nc; b++) {
        double phase = 2*SF_PI*vs[rs[a]]*fabs(ks[cs[b]])*dt;
        res(a,b) = zpx(cos(phase),sin(phase));
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
    oRSF outm, Mexactfile("Mexact"), Mlrfile("Mlr"), Mappfile("Mapp");

    int N;
    velf.get("n1",N);
    float dk;
    velf.get("d1",dx);
    dk = 1.0/(dx*N);

    vs.resize(N);
    ks.resize(N);
    velf >> vs;
    
    int SIZE=11;
    outm.put("n1",N);
    outm.put("n2",SIZE);
    outm.type(SF_COMPLEX);
    int m = N;
    int n = N;

    int count = 0;
    float CUT = N/3.0*dk;
    float k0 = -dk*N/2.0;
    vector<int> ksc(N);
    for (int k=0; k < N; k++) {
	ks[k] = k0+k*dk;
        if (abs(ks[k]) < CUT) {
            ksc[count] = k;
            count++;
        }
    }
    ksc.resize(count);

    /*
    int nk=0;
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
    */
    
    vector<int> cidx, ridx;
    ZpxNumMat mid;
    iC( ddlowrank(m,n,sample,(double) eps,npk,cidx,ridx,mid) );

    ZpxNumMat M1(m,cidx.size());
    vector<int> rs(m);
    for(int k=0; k<m; k++) rs[k]=k;
    iC( sample(rs,cidx,M1) );
    ZpxNumMat M2(ridx.size(),n);
    vector<int> cs(n);
    for(int k=0; k<n; k++) cs[k]=k;
    iC( sample(ridx,cs,M2) );
    ZpxNumMat Mexact(N,N);
    iC( sample(rs,cs,Mexact) );

    Mexactfile.put("n1",N);
    Mexactfile.put("n2",N);
    Mexactfile.put("d2",dk);
    Mexactfile.put("o2",k0);
    Mexactfile.type(SF_COMPLEX);
    std::valarray<sf_complex> fMlr(N*N);
    zpx *ldat = Mexact.data();
    for (int k=0; k < N*N; k++) {
        fMlr[k] = sf_cmplx(real(ldat[k]),imag(ldat[k]));
    } 
    Mexactfile << fMlr;

    Mlrfile.put("n1",N);
    Mlrfile.put("n2",N);
    Mlrfile.put("d2",dk);
    Mlrfile.put("o2",k0);
    Mlrfile.type(SF_COMPLEX);
    Mappfile.put("n1",N);
    Mappfile.put("n2",N);
    Mappfile.put("d2",dk);
    Mappfile.put("o2",k0);
    Mappfile.type(SF_COMPLEX);

    ZpxNumMat Mlr(N,N);

    cerr<<mid._m<<" ";
    cerr<<endl;
    cerr<<mid._n<<" ";
    cerr<<endl;
    cerr<<M2._m<<" ";
    cerr<<endl;
    cerr<<M2._n<<" ";
    cerr<<endl;
    ZpxNumMat tmpM(mid._m,M2._n);
    iC(zzgemm(1.0,mid,M2,0.0,tmpM));
    iC(zzgemm(1.0,M1,tmpM,0.0,Mlr));
    ldat = Mlr.data();
    for (int k=0; k < N*N; k++) {
        fMlr[k] = sf_cmplx(real(ldat[k]),imag(ldat[k]));
    } 
    Mlrfile << fMlr;

    float stmp[] = {-5,-4,-3,-2,-1,0,1,2,3,4,5};
    DblNumMat s(SIZE,1); for(int k=0; k<SIZE; k++) s._data[k]=stmp[k];    
    DblNumMat ktmp(1,N); for(int k=0; k<N; k++) ktmp._data[k]=ks[k];
    DblNumMat ktmpc(1,count); for(int k=0; k<count; k++) ktmpc._data[k]=ks[ksc[k]];
    DblNumMat rB(SIZE,N), rBc(SIZE,count);
    iC(ddgemm(2*SF_PI*dx,s,ktmp,0.0,rB));
    iC(ddgemm(2*SF_PI*dx,s,ktmpc,0.0,rBc));
    ZpxNumMat B(SIZE,N), Bc(SIZE,count);
    for(int k=0; k<B._m*B._n; k++) B._data[k]=zpx(cos(rB._data[k]),sin(rB._data[k]));
    for(int k=0; k<Bc._m*Bc._n; k++) Bc._data[k]=zpx(cos(rBc._data[k]),sin(rBc._data[k]));
    ZpxNumMat IB(N,SIZE);    iC( ddpinv(B, 1e-16, IB) );
    ZpxNumMat IBc(count,SIZE);    iC( ddpinv(Bc, 1e-16, IBc) );
    ZpxNumMat coef(ridx.size(),SIZE);
    ZpxNumMat M2c;
    iC( sample(ridx,ksc,M2c) );
    
    iC(zzgemm(1.0,M2c,IBc,0.0,coef));

    ZpxNumMat G(N,SIZE), tmpG(mid._m,SIZE);
    iC(zzgemm(1.0,mid,coef,0.0,tmpG));
    iC(zzgemm(1.0,M1,tmpG,0.0,G));

    ///*
    Bc.resize(SIZE,1);
    //for(int k=0; k<SIZE; k++) Bc._data[k]=B(k,N/2);
    for(int k=0; k<SIZE; k++) Bc._data[k]=zpx(1,0);
    ZpxNumMat tmpB(N,1);
    iC(zzgemm(1.0,G,Bc,0.0,tmpB));
    for(int x=0; x<N; x++) { 
        cerr<<abs(tmpB._data[x])<<" ";
        tmpB._data[x]=1.0/abs(tmpB._data[x]);
        cerr<<tmpB._data[x]<<" ";
        cerr<<endl;
        }
    for (int x=0; x<N; x++){
        for (int k=0; k<SIZE; k++){
            G(x,k) = G(x,k)*tmpB._data[x];
        }
    }
    //*/
      
    iC(zzgemm(1.0,G,B,0.0,Mlr));
    ldat = Mlr.data();
    for (int k=0; k < N*N; k++) {
        fMlr[k] = sf_cmplx(real(ldat[k]),imag(ldat[k]));
    } 
    Mappfile << fMlr;

    ldat = G.data();
    fMlr.resize(N*SIZE);
    for (int k=0; k < SIZE*N; k++) {
        fMlr[k] = sf_cmplx(real(ldat[k]),imag(ldat[k]));
    } 
    outm << fMlr;

    exit(0);
}

