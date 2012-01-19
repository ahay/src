// 1D 10th-order Lowrank FD coefficient

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

static DblNumMat matrix;
DblNumVec vs; //c
DblNumVec ks; //k


int sample(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,0.0);
    for(int a=0; a<nr; a++) {
	for(int b=0; b<nc; b++) {
        res(a,b) = cos(2*pi*vs(rs[a])*ks(cs[b])*dt);
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

    iRSF in;
    oRSF out;

    int N;
    in.get("n1",N);
    float dx, dk;
    in.get("d1",dx,1.0);
    dk = 1.0/(dx*N);

    vs.resize(N);
    ks.resize(N);
    in >> vs;
    
    for (int k=0; k < N; k++) {
	ks[k] = -dk*N/2.0+(k-1)*dk;
    }
    
    int m = vs.m();
    int n = ks.m();
    
    vector<int> cidx, ridx;
    DblNumMat mid;

    iC( lowrank(m,n,&sample,(double) eps,npk,cidx,ridx,mid) );
    for(int k=0; k<cidx.size(); k++)
       cerr<<cidx[k]<<" ";
    cerr<<endl;
    for(int k=0; k<ridx.size(); k++)
       cerr<<ridx[k]<<" ";
    cerr<<endl;

    DblNumMat M1;
    vector<int> rs(m);
    for(int k=0; k<m; k++) rs[k]=k;
    iC( sample(rs,cidx,M1) );
    DblNumMat M2;
    vector<int> cs(n);
    for(int k=0; k<n; k++) cs[k]=k;
    iC( sample(ridx,cs,M2) );
    DblNumMat Mexact;
    iC( sample(rs,cs,Mexact) );
    oRSF Mexactfile("Mexact");
    Mexactfile.put("n1",N);
    Mexactfile.put("n2",N);
    Mexactfile << Mexact;
    DblNumMat Mlr,tmpM;
    iC(dgemm(1.0,mid,M2,0.0,tmpM);
    iC(dgemm(1.0,M1,tmpM,0.0,Mlr);
    oRSF Mlrfile("Mlr");
    Mlrfile.put("n1",N);
    Mlrfile.put("n2",N);
    Mlrfile << Mlr;

    int SIZE=5;
    double stmp[] = {0,1,2,3,4,5};
    DblNumMat s(SIZE,1); for(int k=0; k<SIZE; k++) s._data[k]=stmp[k];    
    DblNumMat ktmp(1,N); for(int k=0; k<N; k++) ktmp._data[k]=ks[k];
    DblNumMat B(SIZE,N);
    iC(dgemm(2*pi*dx,s,ktmp,0.0,B);
    for(int k=0; k<B._m*B._n; k++) B._data[k]=cos(B._data[k]);
    DblNumMat IB;    iC( pinv(B, 1e-16, IB) );
    DblNumMat coef(n,5);
    iC(dgemm(1.0,M2,IB,0.0,coef);

    DblNumMat G(N,5), tmpG(m,5);
    iC(dgemm(1.0,mid,coef,0.0,tmpG);
    iC(dgemm(1.0,M1,tmpG,0.0,G);

    out.put("n1",N);
    out.put("n2",SIZE);
    out << G;

    return 0;
}





