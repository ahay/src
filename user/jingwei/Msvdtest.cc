//   Singular value decomposition to select rank-1 approximation
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

using namespace std;
using std::cerr;

//------------------------------------------------------------

int main(int argc, char** argv)
{   
    sf_init(argc,argv); // Initialize RSF

    iRSF par(0); // Get parameters
    int nz, nx, nkz, nkx, nzx, nkzx, m, n;
    par.get("nz",nz); 
    par.get("nx",nx);    
    par.get("nkz",nkz);
    par.get("nkx",nkx);
 
    nzx = nz*nx; 
    nkzx = nkz*nkx;

    iRSF matrix; // Get input
    matrix.get("n1",m);
    matrix.get("n2",n);

    if(m!=nzx) cerr<<"Error in size nzx"<<endl;
    if(n!=nkzx) cerr<<"Error in size nkzx"<<endl;

    // Read matrix
    std::valarray<sf_complex> mat(m*n);
    matrix >> mat;

    CpxNumMat M(m,n);
    for (int i=0; i<m; i++)
	for (int j=0; j<n; j++)
	    M(i,j) = cpx(crealf(mat[j*m+i]),cimagf(mat[j*m+i]));
 
    // svd for M = U*S*VT
    int k = min(m,n);  
    CpxNumMat U(m,k);
    FltNumVec S(k);
    CpxNumMat VT(k,n);
    {
	char jobu  = 'S';
        char jobvt = 'S';
        CpxNumMat MC(M);
	int lwork = 20*max(m,n);
	CpxNumVec work(lwork);
        FltNumVec rwork(lwork);
        int info;
        cgesvd_(&jobu, &jobvt, &m, &n, 
	        (MKL_Complex8*) MC.data(), &m, S.data(), 
	        (MKL_Complex8*) U.data(), &m, 
	        (MKL_Complex8*) VT.data(), &k, 
	        (MKL_Complex8*) work.data(), &lwork, 
	        rwork.data(), &info);    
	iA(info==0);
    }
    cerr<<"First five singular values are "<<S(0)<<" "<<S(1)<<" "<<S(2)<<" "<<S(3)<<" "<<S(4)<<endl;
    cerr<<"Last five singular values are "<<S(k-5)<<" "<<S(k-4)<<" "<<S(k-3)<<" "<<S(k-2)<<" "<<S(k-1)<<endl;
   
    // extract alpha and beta
    CpxNumVec aa, bb;
    aa.resize(m);
    bb.resize(n);
    for(int i=0; i<m; i++)  aa(i) = U(i,0)*sqrt(S(0));
    for(int i=0; i<n; i++)  bb(i) = sqrt(S(0))*VT(0,i);

    // check rank-1 error
    CpxNumMat Mapp(m,n), Merr(m,n);
    for(int i=0; i<m; i++)
	for(int j=0; j<n; j++) {
	    Mapp(i,j) = aa(i)*bb(j);
            Merr(i,j) = M(i,j) - Mapp(i,j);
	}
    cerr<<"Exact rank-1 rel err "<<sqrt(energy(Merr))/sqrt(energy(M))<<endl;
   
    
    // Write alpha and beta
    std::valarray<sf_complex> adata(m), bdata(n), appMdata(m*n);
    for (int i=0; i<m; i++)
        adata[i] = sf_cmplx(real(aa(i)),imag(aa(i)));
    for (int i=0; i<n; i++) 
        bdata[i] = sf_cmplx(real(bb(i)),imag(bb(i)));
    for (int i=0; i<m; i++) 
	for (int j=0; j<n; j++) 
	    appMdata[j*m+i] = sf_cmplx(real(Mapp(i,j)),imag(Mapp(i,j)));

    oRSF alpha, beta("beta"), appM("appM");
    alpha.type(SF_COMPLEX);
    alpha.put("n1",m);
    alpha.put("n2",1);
    alpha << adata;  

    beta.type(SF_COMPLEX);
    beta.put("n1",n);
    beta.put("n2",1);
    beta << bdata;
    
    appM.type(SF_COMPLEX);
    appM.put("n1",m);
    appM.put("n2",n);
    appM << appMdata;
    
    
    exit(0);
}
