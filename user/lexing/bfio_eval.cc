// Butterfly algorithm
//   Copyright (C) 2011 University of Texas at Austin
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

#include "bfio.hh"
#include "serialize.hh"

using std::istringstream;
using std::ifstream;
using std::ofstream;
using std::set;
using std::queue;
using std::cerr;

//---------------------------------------
int BFIO::eval(int N, const CpxNumMat& f, CpxNumMat& u)
{
    int N1 = f.m();
    int N2 = f.n();
    int M1 = u.m();
    int M2 = u.n();
    u.resize(M1,M2);
    setvalue(u,cpx(0,0));
    //--------
    Entry& entx1 = _e2dmap[_EPSx1];
    DblNumVec& gridx1 = entx1.grid();
    NumVec<CpxNumMat>& matsx1 = entx1.mats();
    // CpxNumMat& dir = ent.dir();
    //the transpose matrices
    NumVec<CpxNumMat> tmatsx1(2);
    iC( ztran(matsx1(0), tmatsx1(0)) );
    iC( ztran(matsx1(1), tmatsx1(1)) );
    // CpxNumMat tdir;
    // iC( ztran(dir, tdir) );
    ///
    Entry& entx2 = _e2dmap[_EPSx2]; 
    DblNumVec& gridx2 = entx2.grid();
    NumVec<CpxNumMat>& matsx2 = entx2.mats();
    NumVec<CpxNumMat> tmatsx2(2);
    iC( ztran(matsx2(0), tmatsx2(0)) );
    iC( ztran(matsx2(1), tmatsx2(1)) );
    //	
    Entry& entk1 = _e2dmap[_EPSk1]; 
    DblNumVec& gridk1 = entk1.grid();
    NumVec<CpxNumMat>& matsk1 = entk1.mats();
    NumVec<CpxNumMat> tmatsk1(2);
    iC( ztran(matsk1(0), tmatsk1(0)) );
    iC( ztran(matsk1(1), tmatsk1(1)) );
    //	
    Entry& entk2 = _e2dmap[_EPSk2]; 
    DblNumVec& gridk2 = entk2.grid();
    NumVec<CpxNumMat>& matsk2 = entk2.mats();
    NumVec<CpxNumMat> tmatsk2(2);
    iC( ztran(matsk2(0), tmatsk2(0)) );
    iC( ztran(matsk2(1), tmatsk2(1)) );
    //
    //	
    int NGx1 = _EPSx1;
    int NGx2 = _EPSx2;
    int NGk1 = _EPSk1;
    int NGk2 = _EPSk2;
    //
    int TL = int(round(log(double(N))/log(2)));
    int SL = TL;
    int EL = 0;
    int ML = int(floor((SL+EL)/2.0));
    int nz = pow2(EL);
    int zB = N/nz;
    //
    for(int z1=0; z1<nz; z1++)
	for(int z2=0; z2<nz; z2++) {
	    cerr<<z1<<" "<<z2<<endl;
	    int k1stt = z1*zB;      int k1end = (z1+1)*zB;
	    int k2stt = z2*zB;      int k2end = (z2+1)*zB;
	    NumMat< NumMat< CpxNumMat > > NOW;
	    //----------------------------------------------------------------------
	    for(int ell=SL; ell>=ML; ell--) {
		int nk = pow2(ell);	int nx = N/nk;
		int kB = N/nk;	int xB = N/nx;
		NumMat< NumMat< CpxNumMat > > PRE = NOW;
		NOW.resize(nk,nk);
		for(int k1=k1stt/kB; k1<k1end/kB; k1++)
		    for(int k2=k2stt/kB; k2<k2end/kB; k2++) {
			NOW(k1,k2).resize(nx,nx);
			for(int x1=0; x1<nx; x1++)
			    for(int x2=0; x2<nx; x2++) {
				NOW(k1,k2)(x1,x2).resize(NGk1,NGk2);		setvalue(NOW(k1,k2)(x1,x2),cpx(0,0));
			    }
		    }
		//
		vector<Point2> ko;
		for(int j=0; j<NGk2; j++)
		    for(int i=0; i<NGk1; i++)
			ko.push_back( Point2(gridk1(i)*kB-N/2, gridk2(j)*kB-N/2) );
		vector<Point2> co;
		for(int j=0; j<NGk2; j++)
		    for(int i=0; i<NGk1; i++)
			co.push_back( Point2(gridk1(i)*kB/2-N/2, gridk2(j)*kB/2-N/2) );
		//
		for(int k1=k1stt/kB; k1<k1end/kB; k1++)
		    for(int k2=k2stt/kB; k2<k2end/kB; k2++) {
			Point2 kc( (k1+0.5)*kB - N/2, (k2+0.5)*kB - N/2 );
			for(int x1=0; x1<nx; x1++)
			    for(int x2=0; x2<nx; x2++) {
				Point2 xc( (x1+0.5)*xB, (x2+0.5)*xB );
				//-------
				CpxNumMat all(NGk1,NGk2);		setvalue(all, cpx(0,0));
				if(ell==SL) {
				    //get
				    int kB1=N1/pow2(SL), kB2=N2/pow2(SL);
				    CpxNumMat ext(kB1, kB2);
				    for(int j=0; j<kB2; j++)
					for(int i=0; i<kB1; i++) {
					    ext(i,j) = f(i+k1*kB1, j+k2*kB2);
					}
				    //scale
				    vector<Point2> trg;		trg.push_back(xc);
				    vector<Point2> so;
				    for(int j=0; j<kB2; j++)
					for(int i=0; i<kB1; i++)
					    so.push_back( Point2(double(i*N)/double(N1)-N/2,double(j*N)/double(N2)-N/2) );
				    vector<Point2> src(so);
				    for(int g=0; g<(int)src.size(); g++)
					src[g] = src[g] + Point2(double(k1*kB1*N)/double(N1), double(k2*kB2*N)/double(N2));
				    CpxNumMat scl;		iC( kernel(N, trg, src, scl) );
				    CpxNumMat sclaux(kB1,kB2,false,scl.data());
				    for(int j=0; j<kB2; j++)
					for(int i=0; i<kB1; i++) {
					    ext(i,j) = ext(i,j) * sclaux(i,j);
					}
				    //transform
				    NumVec<CpxNumMat>& dirk1 = entk1.dirc();
				    CpxNumMat dirSLk1 = dirk1(int(round(log(double(kB1))/log(2))));
				    NumVec<CpxNumMat>& dirk2 = entk2.dirc();
				    CpxNumMat dirSLk2 = dirk2(int(round(log(double(kB2))/log(2))));
				    //CpxNumMat& dirSLk1 = entk1.dirc(kB1);
				    //CpxNumMat& dirSLk2 = entk2.dirc(kB2);
				    //the transpose matrices
				    CpxNumMat tdirSLk2;
				    iC( ztran(dirSLk2, tdirSLk2) );
				    CpxNumMat tmp(NGk1,kB2);		setvalue(tmp,cpx(0,0));
				    iC( zgemm(1, dirSLk1, ext, 0, tmp) );
				    iC( zgemm(1, tmp, tdirSLk2, 1, all) );
				} else {
				    int p1 = int(floor(x1/2));		int p2 = int(floor(x2/2));
				    for(int a1=0; a1<2; a1++)
					for(int a2=0; a2<2; a2++) {
					    int c1 = 2*k1+a1;		    int c2 = 2*k2+a2;
					    //get
					    CpxNumMat ext(PRE(c1,c2)(p1,p2));
					    //scale
					    vector<Point2> trg;		    trg.push_back(xc);
					    vector<Point2> src(co);
					    for(int g=0; g<(int)src.size(); g++)
						src[g] = src[g] + Point2(c1*kB/2, c2*kB/2);
					    CpxNumMat scl;		    iC( kernel(N, trg, src, scl) );
					    CpxNumMat sclaux(NGk1,NGk2,false,scl.data());
					    for(int j=0; j<NGk2; j++)
						for(int i=0; i<NGk1; i++)
						    ext(i,j) = ext(i,j) * sclaux(i,j);
					    //transform
					    CpxNumMat tmp(NGk1,NGk2);		setvalue(tmp,cpx(0,0));
					    iC( zgemm(1, matsk1(a1), ext, 0, tmp) );
					    iC( zgemm(1, tmp, tmatsk2(a2), 1, all) );
					}
				}
				//scale
				vector<Point2> trg;		trg.push_back(xc);
				vector<Point2> src(ko);
				for(int g=0; g<(int)src.size(); g++)
				    src[g] = src[g] + Point2(k1*kB,k2*kB);
				CpxNumMat scl;		iC( kernel(N, trg, src, scl) );
				CpxNumMat sclaux(NGk1,NGk2,false,scl.data());
				for(int j=0; j<NGk2; j++)
				    for(int i=0; i<NGk1; i++)
					all(i,j) = all(i,j) / sclaux(i,j);
				//put
				NOW(k1,k2)(x1,x2) = all;
				//cerr<<k1<<" "<<k2<<" "<<x1<<" "<<x2<<" "<<real(all(0,0))<<" "<<imag(all(0,0))<<endl;
			    }//x1 x2
		    }//k1k2
		PRE.resize(0,0);
	    }//ell
	    //----------------------------------------------------------------------
	    if(1) {
		int ell = ML;
		int nk = pow2(ell);	int nx = N/nk;
		int kB = N/nk;	int xB = N/nx;
		//
		vector<Point2> ko;
		for(int j=0; j<NGk2; j++)
		    for(int i=0; i<NGk1; i++)
			ko.push_back( Point2(gridk1(i)*kB-N/2, gridk2(j)*kB-N/2) );
		vector<Point2> xo;
		for(int j=0; j<NGx2; j++)
		    for(int i=0; i<NGx1; i++)
			xo.push_back( Point2(gridx1(i)*xB, gridx2(j)*xB) );
		//
		for(int k1=k1stt/kB; k1<k1end/kB; k1++)
		    for(int k2=k2stt/kB; k2<k2end/kB; k2++) {
			Point2 kc( (k1+0.5)*kB - N/2, (k2+0.5)*kB - N/2 );
			for(int x1=0; x1<nx; x1++)
			    for(int x2=0; x2<nx; x2++) {
				Point2 xc( (x1+0.5)*xB, (x2+0.5)*xB );
				//--------
				vector<Point2> src(ko);
				for(int g=0; g<(int)src.size(); g++)
				    src[g] = src[g] + Point2(k1*kB, k2*kB);
				vector<Point2> trg(xo);
				for(int g=0; g<(int)trg.size(); g++)
				    trg[g] = trg[g] + Point2(x1*xB, x2*xB);
				CpxNumMat evl(NGx1*NGx2,NGk1*NGk2);		iC( kernel(N, trg, src, evl) );
				CpxNumVec den(NGk1*NGk2,true,NOW(k1,k2)(x1,x2).data());
				CpxNumVec val(NGx1*NGx2,false,NOW(k1,k2)(x1,x2).data());
				iC( zgemv(1, evl, den, 0, val) );		//NOW(k1,k2)(x1,x2) = tmp;
				//cerr<<k1<<" "<<k2<<" "<<x1<<" "<<x2<<" "<<real(val(0))<<" "<<imag(val(0))<<endl;
			    }
		    }
	    }//ell
	    //----------------------------------------------------------------------
	    for(int ell=ML; ell>=EL; ell--) {
		int nk = pow2(ell);	int nx = N/nk;
		int kB = N/nk;	int xB = N/nx;
		NumMat< NumMat< CpxNumMat > > NXT;
		if(ell!=EL) {
		    NXT.resize(nk/2,nk/2);
		    for(int k1=k1stt/(2*kB); k1<k1end/(2*kB); k1++)
			for(int k2=k2stt/(2*kB); k2<k2end/(2*kB); k2++) {
			    NXT(k1,k2).resize(2*nx,2*nx);
			    for(int x1=0; x1<2*nx; x1++)
				for(int x2=0; x2<2*nx; x2++) {
				    NXT(k1,k2)(x1,x2).resize(NGx1,NGx2);		setvalue(NXT(k1,k2)(x1,x2),cpx(0,0));
				}
			}
		}
		//
		vector<Point2> xo;
		for(int j=0; j<NGx2; j++)
		    for(int i=0; i<NGx1; i++)
			xo.push_back( Point2(gridx1(i)*xB, gridx2(j)*xB) );
		vector<Point2> co;
		for(int j=0; j<NGx2; j++)
		    for(int i=0; i<NGx1; i++)
			co.push_back( Point2(gridx1(i)*xB/2, gridx2(j)*xB/2) );
		//
		for(int k1=k1stt/kB; k1<k1end/kB; k1++)
		    for(int k2=k2stt/kB; k2<k2end/kB; k2++) {
			Point2 kc( (k1+0.5)*kB - N/2, (k2+0.5)*kB - N/2 );
			for(int x1=0; x1<nx; x1++)
			    for(int x2=0; x2<nx; x2++) {
				Point2 xc( (x1+0.5)*xB, (x2+0.5)*xB );
				//-------
				//get
				CpxNumMat all(NOW(k1,k2)(x1,x2));
				//scale
				vector<Point2> src;		src.push_back(kc);
				vector<Point2> trg(xo);
				for(int g=0; g<(int)trg.size(); g++)
				    trg[g] = trg[g] + Point2(x1*xB,x2*xB);
				CpxNumMat scl;		iC( kernel(N, trg, src, scl) );
				CpxNumMat sclaux(NGx1,NGx2,false,scl.data());
				for(int j=0; j<NGx2; j++)
				    for(int i=0; i<NGx1; i++) {
					all(i,j) = all(i,j) / sclaux(i,j);
				    }
				//
				if(ell!=EL) {
				    int q1 = int(floor(k1/2));		  int q2 = int(floor(k2/2));
				    for(int a1=0; a1<2; a1++)
					for(int a2=0; a2<2; a2++) {
					    int c1 = 2*x1+a1;		      int c2 = 2*x2+a2;
					    //transform
					    CpxNumMat ext(NGx1,NGx2);		      setvalue(ext,cpx(0,0));
					    CpxNumMat tmp(NGx1,NGx2);		      setvalue(tmp,cpx(0,0));
					    iC( zgemm(1, tmatsx1(a1), all, 0, tmp) );
					    iC( zgemm(1, tmp, matsx2(a2), 1, ext) );
					    //scale
					    vector<Point2> src;		      src.push_back(kc);
					    vector<Point2> trg(co);
					    for(int g=0; g<(int)trg.size(); g++)
						trg[g] = trg[g] + Point2(c1*xB/2,c2*xB/2);
					    CpxNumMat scl;		      iC( kernel(N, trg, src, scl) );
					    CpxNumMat sclaux(NGx1,NGx2,false,scl.data());
					    for(int j=0; j<NGx2; j++)
						for(int i=0; i<NGx1; i++)
						    ext(i,j) = ext(i,j) * sclaux(i,j);
					    //put
					    for(int j=0; j<NGx2; j++)
						for(int i=0; i<NGx1; i++)
						    NXT(q1,q2)(c1,c2)(i,j) = NXT(q1,q2)(c1,c2)(i,j) + ext(i,j);
					}
				} else {
				    //transform
				    int xB1=M1/pow2(SL), xB2=M2/pow2(SL);
				    NumVec<CpxNumMat>& dirx1 = entx1.dirc();
				    CpxNumMat dirSLx1 = dirx1(int(round(log(double(xB1))/log(2))));
				    NumVec<CpxNumMat>& dirx2 = entx2.dirc();
				    CpxNumMat dirSLx2 = dirx2(int(round(log(double(xB2))/log(2))));			
				    //CpxNumMat& dirSLx1 = entx1.dirc(xB1);
				    //CpxNumMat& dirSLx2 = entx2.dirc(xB2);
				    //the transpose matrices
				    CpxNumMat tdirSLx1;
				    iC( ztran(dirSLx1, tdirSLx1) );
				    CpxNumMat ext(xB1,xB2);		  setvalue(ext,cpx(0,0));
				    CpxNumMat tmp(xB1,NGx2);		  setvalue(tmp,cpx(0,0));
				    iC( zgemm(1, tdirSLx1, all, 0, tmp) );
				    iC( zgemm(1, tmp, dirSLx2, 1, ext) );
				    //scale
				    vector<Point2> src;		  src.push_back(kc);
				    vector<Point2> to;
				    for(int j=0; j<xB2; j++)
					for(int i=0; i<xB1; i++)
					    to.push_back( Point2(double(i*N)/double(M1),double(j*N)/double(M2)) );
				    vector<Point2> trg(to);
				    for(int g=0; g<(int)trg.size(); g++)
					trg[g] = trg[g] + Point2(double(x1*xB1*N)/double(M1), double(x2*xB2*N)/double(M2));
				    CpxNumMat scl;		  iC( kernel(N, trg, src, scl) );
				    CpxNumMat sclaux(xB1,xB2,false,scl.data());
				    for(int j=0; j<xB2; j++)
					for(int i=0; i<xB1; i++)
					    ext(i,j) = ext(i,j) * sclaux(i,j);
				    //cerr<<k1<<" "<<k2<<" "<<x1<<" "<<x2<<" "<<real(ext(0,0))<<" "<<imag(ext(0,0))<<endl;
				    //put
				    for(int j=0; j<xB2; j++)
					for(int i=0; i<xB1; i++)
					    u(i+x1*xB1,j+x2*xB2) = u(i+x1*xB1,j+x2*xB2) + ext(i,j);
				}
			    }//x1x2
		    }//k1k2
		NOW = NXT;
		NXT.resize(0,0);
	    }//ell
	}//z1z2
    return 0;
}

