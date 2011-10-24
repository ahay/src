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
int BFIO::eval(const CpxNumMat& f, CpxNumMat& u)
{
  int N = f.m();
  u.resize(N,N);
  setvalue(u,cpx(0,0));
  //--------
  Entry& ent = _e2dmap[_EPS];
  DblNumVec& grid = ent.grid();
  NumVec<CpxNumMat>& mats = ent.mats();
  CpxNumMat& dir = ent.dir();
  //the transpose matrices
  NumVec<CpxNumMat> tmats(2);
  iC( ztran(mats(0), tmats(0)) );
  iC( ztran(mats(1), tmats(1)) );
  CpxNumMat tdir;
  iC( ztran(dir, tdir) );
  ///
  int NG = _EPS;
  //
  int TL = int(round(log(double(N))/log(2)));
  int SL = TL-3;
  int EL = 3;
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
		NOW(k1,k2)(x1,x2).resize(NG,NG);		setvalue(NOW(k1,k2)(x1,x2),cpx(0,0));
	      }
	  }
	//
	vector<Point2> so;
	for(int j=0; j<kB; j++)
	  for(int i=0; i<kB; i++)
	    so.push_back( Point2(i-N/2,j-N/2) );
	vector<Point2> ko;
	for(int j=0; j<NG; j++)
	  for(int i=0; i<NG; i++)
	    ko.push_back( Point2(grid(i)*kB-N/2, grid(j)*kB-N/2) );
	vector<Point2> co;
	for(int j=0; j<NG; j++)
	  for(int i=0; i<NG; i++)
	    co.push_back( Point2(grid(i)*kB/2-N/2, grid(j)*kB/2-N/2) );
	//
	for(int k1=k1stt/kB; k1<k1end/kB; k1++)
	  for(int k2=k2stt/kB; k2<k2end/kB; k2++) {
	    Point2 kc( (k1+0.5)*kB - N/2, (k2+0.5)*kB - N/2 );
	    for(int x1=0; x1<nx; x1++)
	      for(int x2=0; x2<nx; x2++) {
		Point2 xc( (x1+0.5)*xB, (x2+0.5)*xB );
		//-------
		CpxNumMat all(NG,NG);		setvalue(all, cpx(0,0));
		if(ell==SL) {
		  //get
		  CpxNumMat ext(kB, kB);
		  for(int j=0; j<kB; j++)
		    for(int i=0; i<kB; i++) {
		      ext(i,j) = f(i+k1*kB, j+k2*kB);
		    }
		  //scale
		  vector<Point2> trg;		trg.push_back(xc);
		  vector<Point2> src(so);
		  for(int g=0; g< (int)src.size(); g++)
		    src[g] = src[g] + Point2(k1*kB, k2*kB);
		  CpxNumMat scl;		iC( kernel(N, trg, src, scl) );
		  CpxNumMat sclaux(kB,kB,false,scl.data());
		  for(int j=0; j<kB; j++)
		    for(int i=0; i<kB; i++) {
		      ext(i,j) = ext(i,j) * sclaux(i,j);
		    }
		  //transform
		  CpxNumMat tmp(NG,kB);		setvalue(tmp,cpx(0,0));
		  iC( zgemm(1, dir, ext, 0, tmp) );
		  iC( zgemm(1, tmp, tdir, 1, all) );
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
		      CpxNumMat sclaux(NG,NG,false,scl.data());
		      for(int j=0; j<NG; j++)
			for(int i=0; i<NG; i++)
			  ext(i,j) = ext(i,j) * sclaux(i,j);
		      //transform
		      CpxNumMat tmp(NG,NG);		setvalue(tmp,cpx(0,0));
		      iC( zgemm(1, mats(a1), ext, 0, tmp) );
		      iC( zgemm(1, tmp, tmats(a2), 1, all) );
		    }
		}
		//scale
		vector<Point2> trg;		trg.push_back(xc);
		vector<Point2> src(ko);
		for(int g=0; g<(int)src.size(); g++)
		  src[g] = src[g] + Point2(k1*kB,k2*kB);
		CpxNumMat scl;		iC( kernel(N, trg, src, scl) );
		CpxNumMat sclaux(NG,NG,false,scl.data());
		for(int j=0; j<NG; j++)
		  for(int i=0; i<NG; i++)
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
	for(int j=0; j<NG; j++)
	  for(int i=0; i<NG; i++)
	    ko.push_back( Point2(grid(i)*kB-N/2, grid(j)*kB-N/2) );
	vector<Point2> xo;
	for(int j=0; j<NG; j++)
	  for(int i=0; i<NG; i++)
	    xo.push_back( Point2(grid(i)*xB, grid(j)*xB) );
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
		CpxNumMat evl(NG*NG,NG*NG);		iC( kernel(N, trg, src, evl) );
		CpxNumVec den(NG*NG,true,NOW(k1,k2)(x1,x2).data());
		CpxNumVec val(NG*NG,false,NOW(k1,k2)(x1,x2).data());
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
		  NXT(k1,k2)(x1,x2).resize(NG,NG);		setvalue(NXT(k1,k2)(x1,x2),cpx(0,0));
		}
	    }
	}
	//
	vector<Point2> to;
	for(int j=0; j<xB; j++)
	  for(int i=0; i<xB; i++)
	    to.push_back( Point2(i,j) );
	vector<Point2> xo;
	for(int j=0; j<NG; j++)
	  for(int i=0; i<NG; i++)
	    xo.push_back( Point2(grid(i)*xB, grid(j)*xB) );
	vector<Point2> co;
	for(int j=0; j<NG; j++)
	  for(int i=0; i<NG; i++)
	    co.push_back( Point2(grid(i)*xB/2, grid(j)*xB/2) );
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
		CpxNumMat sclaux(NG,NG,false,scl.data());
		for(int j=0; j<NG; j++)
		  for(int i=0; i<NG; i++) {
		    all(i,j) = all(i,j) / sclaux(i,j);
		  }
		//
		if(ell!=EL) {
		  int q1 = int(floor(k1/2));		  int q2 = int(floor(k2/2));
		  for(int a1=0; a1<2; a1++)
		    for(int a2=0; a2<2; a2++) {
		      int c1 = 2*x1+a1;		      int c2 = 2*x2+a2;
		      //transform
		      CpxNumMat ext(NG,NG);		      setvalue(ext,cpx(0,0));
		      CpxNumMat tmp(NG,NG);		      setvalue(tmp,cpx(0,0));
		      iC( zgemm(1, tmats(a1), all, 0, tmp) );
		      iC( zgemm(1, tmp, mats(a2), 1, ext) );
		      //scale
		      vector<Point2> src;		      src.push_back(kc);
		      vector<Point2> trg(co);
		      for(int g=0; g<(int)trg.size(); g++)
			trg[g] = trg[g] + Point2(c1*xB/2,c2*xB/2);
		      CpxNumMat scl;		      iC( kernel(N, trg, src, scl) );
		      CpxNumMat sclaux(NG,NG,false,scl.data());
		      for(int j=0; j<NG; j++)
			for(int i=0; i<NG; i++)
			  ext(i,j) = ext(i,j) * sclaux(i,j);
		      //put
		      for(int j=0; j<NG; j++)
			for(int i=0; i<NG; i++)
			  NXT(q1,q2)(c1,c2)(i,j) = NXT(q1,q2)(c1,c2)(i,j) + ext(i,j);
		    }
		} else {
		  //transform
		  CpxNumMat ext(xB,xB);		  setvalue(ext,cpx(0,0));
		  CpxNumMat tmp(xB,NG);		  setvalue(tmp,cpx(0,0));
		  iC( zgemm(1, tdir, all, 0, tmp) );
		  iC( zgemm(1, tmp, dir, 1, ext) );
		  //scale
		  vector<Point2> src;		  src.push_back(kc);
		  vector<Point2> trg(to);
		  for(int g=0; g<(int)trg.size(); g++)
		    trg[g] = trg[g] + Point2(x1*xB, x2*xB);
		  CpxNumMat scl;		  iC( kernel(N, trg, src, scl) );
		  CpxNumMat sclaux(xB,xB,false,scl.data());
		  for(int j=0; j<xB; j++)
		    for(int i=0; i<xB; i++)
		      ext(i,j) = ext(i,j) * sclaux(i,j);
		  //cerr<<k1<<" "<<k2<<" "<<x1<<" "<<x2<<" "<<real(ext(0,0))<<" "<<imag(ext(0,0))<<endl;
		  //put
		  for(int j=0; j<xB; j++)
		    for(int i=0; i<xB; i++)
		      u(i+x1*xB,j+x2*xB) = u(i+x1*xB,j+x2*xB) + ext(i,j);
		}
	      }//x1x2
	  }//k1k2
	NOW = NXT;
	NXT.resize(0,0);
      }//ell
    }//z1z2
  return 0;
}


