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
    
    vector<int> cidx;
    vector<int> ridx;
    DblNumMat mid;

    iC( lowrank(m,n,sample,(double) eps,npk,cidx,ridx,mid) );

    int m2=mid.m();
    int n2=mid.n();
    double *dmid = mid.data();

    std::valarray<float> fmid(m2*n2);
    for (int k=0; k < m2*n2; k++) {
	fmid[k] = dmid[k];
    }

    out.put("n1",m2);
    out.put("n2",n2);
    out << fmid;

    for(int k=0; k<m2; k++)
	cerr<<cidx[k]<<" ";
    cerr<<endl;

    for(int k=0; k<n2; k++)
	cerr<<ridx[k]<<" ";
    cerr<<endl;

    // output left, right, and approximation
  
    return 0;
}





