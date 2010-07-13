// Lowrank matrix decomposition

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

int sample(vector<int>& rs, vector<int>& cs, DblNumMat& res)
{
    int nr = rs.size();
    int nc = cs.size();
    res.resize(nr,nc);  
    setvalue(res,0.0);
    for(int a=0; a<nr; a++) {
	for(int b=0; b<nc; b++) {
	    res(a,b) = matrix(rs[a],cs[b]);
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

    int m,n;
    in.get("n1",m);
    in.get("n2",n);
    std::valarray<float> fdata(m*n);
    in >> fdata;
    
    std::valarray<double> data(m*n);
    for (int k=0; k < m*n; k++) {
	data[k] = fdata[k];
    }
    
    DblNumMat mat(m,n,false,&(data[0]));
    matrix = mat;
    
    vector<int> midx(m), nidx(n), lidx, ridx;
    DblNumMat mid;

    iC( lowrank(m,n,sample,(double) eps,npk,lidx,ridx,mid) );

    int m2=mid.m();
    int n2=mid.n();
    double *dmid = mid.data();

    std::valarray<float> fmid(m2*n2);
    for (int k=0; k < m2*n2; k++) {
	fmid[k] = dmid[k];
    }

    oRSF middle("mid");
    middle.put("n1",m2);
    middle.put("n2",n2);
    middle << fmid;

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

    DblNumMat tmp1(m,n2);
    iC( dgemm(1.0, lmat, mid, 0.0, tmp1) );
    iC( dgemm(1.0, tmp1, rmat, 0.0, mat) );
    for (int k=0; k < m*n; k++) 
	fdata[k] = data[k];
    out << fdata;

    return 0;
}





