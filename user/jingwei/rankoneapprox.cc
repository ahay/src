//   Get rank-1 approximation alpha, beta based on lowrank decomposition M1, M2
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

#include "rankoneapprox.hh"

//---------------------------------------
int rankoneapprox(const CpxNumMat& M1, const CpxNumMat& M2, CpxNumVec& alpha, CpxNumVec& beta, int npk)
{
    int m = M1.m();
    int n = M2.n();
    int k = M1.n();
    if(k==1) {
	alpha.resize(m);
	beta.resize(n);
	for(int i=0; i<m; i++)  alpha(i) = M1(i,0);
        for(int i=0; i<n; i++)  beta(i) = M2(0,i);
        cerr<<"Already a rank-1 approximation"<<endl;
    } else {

    // svd for M1 = U1*S1*VT1
    m = M1.m();
    n = M1.n();
    k = min(m,n);  
    CpxNumMat U1(m,k);
    FltNumVec S1(k);
    CpxNumMat VT1(k,n);
    {
	char jobu  = 'S';
        char jobvt = 'S';
        CpxNumMat MC(M1);
	int lwork = 20*max(m,n);
	CpxNumVec work(lwork);
        FltNumVec rwork(lwork);
        int info;
        cgesvd_(&jobu, &jobvt, &m, &n, 
	        (MKL_Complex8*) MC.data(), &m, S1.data(), 
	        (MKL_Complex8*) U1.data(), &m, 
	        (MKL_Complex8*) VT1.data(), &k, 
	        (MKL_Complex8*) work.data(), &lwork, 
	        rwork.data(), &info);    
	iA(info==0);
    }
    
    // svd for M2 = U2*S2*VT2
    m = M2.m();
    n = M2.n();
    k = min(m,n);  
    CpxNumMat U2(m,k);
    FltNumVec S2(k);
    CpxNumMat VT2(k,n);
    {
	char jobu  = 'S';
	char jobvt = 'S';
	CpxNumMat MC(M2);
	int lwork = 20*max(m,n);
	CpxNumVec work(lwork);
	FltNumVec rwork(lwork);
	int info;
	cgesvd_(&jobu, &jobvt, &m, &n, 
		(MKL_Complex8*) MC.data(), &m, S2.data(), 
		(MKL_Complex8*) U2.data(), &m, 
		(MKL_Complex8*) VT2.data(), &k, 
		(MKL_Complex8*) work.data(), &lwork, 
		rwork.data(), &info);    
	iA(info==0);
    }

    // VT1 = S1*VT1
    m = VT1.m();
    n = VT1.n();
    for(int i=0; i<m; i++)
	for(int j=0; j<n; j++)
	    VT1(i,j) *=S1(i); 

    // U2 = U2*S2
    m = U2.m();
    n = U2.n();
    for(int i=0; i<m; i++)
	for(int j=0; j<n; j++)
	    U2(i,j) *=S2(j);   
  
    // MM = VT1*U2
    m = VT1.m();
    n = U2.n();
    CpxNumMat MM(m,n);
    iC( zgemm(1.0, VT1, U2, 0.0, MM) );

    // svd for MM = UM*SM*VTM
    m = MM.m();
    n = MM.n();
    k = min(m,n);  
    CpxNumMat UM(m,k);
    FltNumVec SM(k);
    CpxNumMat VTM(k,n);
    {
	char jobu  = 'S';
	char jobvt = 'S';
	CpxNumMat MC(MM);
	int lwork = 20*max(m,n);
	CpxNumVec work(lwork);
	FltNumVec rwork(lwork);
	int info;
	cgesvd_(&jobu, &jobvt, &m, &n, 
		(MKL_Complex8*) MC.data(), &m, SM.data(), 
		(MKL_Complex8*) UM.data(), &m, 
	        (MKL_Complex8*) VTM.data(), &k, 
	        (MKL_Complex8*) work.data(), &lwork, 
		rwork.data(), &info);    
	iA(info==0);
    }

    // U = U1*UM    
    m = U1.m();
    n = UM.n();
    CpxNumMat U(m,n);
    iC( zgemm(1.0, U1, UM, 0.0, U) );

    // VT = VTM*VT2
    m = VTM.m();
    n = VT2.n();
    CpxNumMat VT(m,n);
    iC( zgemm(1.0, VTM, VT2, 0.0, VT) );

    // extract alpha and beta
    m = U.m();
    n = VT.n();
    alpha.resize(m);
    beta.resize(n); 
    for(int i=0; i<m; i++)  alpha(i) = U(i,0);
    for(int i=0; i<n; i++)  beta(i) = SM(0)*VT(0,i);
    cerr<<"Singular Values are "<<SM<<" "<<endl;

    }

    {	// check rank-1 error
        int m = M1.m();
        int n = M2.n();
        int k = M1.n();

        int nr = min(npk,m);
	vector<int> rs(nr);
	for(int i=0; i<nr; i++)      rs[i] = int( floor(drand48()*m) );

	int nc = min(npk,n);
	vector<int> cs(nc);
	for(int i=0; i<nc; i++)      cs[i] = int( floor(drand48()*n) );
	
	CpxNumMat sM1(nr,k);
        for(int i=0; i<nr; i++)
	    for(int j=0; j<k; j++)
		sM1(i,j) = M1(i,j);

	CpxNumMat sM2(k,nc);
        for(int i=0; i<k; i++)
	    for(int j=0; j<nc; j++)
		sM2(i,j) = M2(i,j);	

	CpxNumMat Mext(nr,nc), Mapp(nr,nc), Merr(nr,nc);
        iC( zgemm(1.0, sM1, sM2, 0.0, Mext) );

  	for(int i=0; i<nr; i++)
	    for(int j=0; j<nc; j++) {
		Mapp(i,j) = alpha(i)*beta(j);
                Merr(i,j) = Mext(i,j) - Mapp(i,j);
	    }
	cerr<<"rank-1 rel err "<<sqrt(energy(Merr))/sqrt(energy(Mext))<<endl;
    }
    
    return 0;
}
//------------------------------------------------------------
