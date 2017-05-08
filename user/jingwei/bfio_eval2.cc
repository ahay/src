//   2to2 butterfly
//   BFIO::prep_aux
//   BFIO::eval2       if fi=1 or 3, switch at L/2 (even L), (L-1)/2 (odd L)
//                     if fi=2 or 4, switch at L/2 (even L), (L+1)/2 (odd L)
//
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
inline int BFIO::prep_aux(FltNumVec& grid, vector<float>& ts, CpxNumMat& tmp)
{
  int NG = grid.m();
  int NT = ts.size();
  tmp.resize(NG,NT);
  for(int b=0; b<NT; b++)
    for(int a=0; a<NG; a++) {
      float prd = 1.0;
      for(int c=0; c<NG; c++) {
	if(c!=a) {
	  prd *= (ts[b]-grid(c)) / (grid(a)-grid(c));
	}
      }
      tmp(a,b) = prd;
    }
  return 0;
}

//---------------------------------------
int BFIO::eval2(int N, const CpxNumMat& f, const FltNumVec& w, const FltNumVec& x, CpxNumMat& u, const FltNumVec& tau, const FltNumVec& p)
{
  //---------------
  clock_t time0 = clock();  
  //--------------

  int N1 = f.m();
  int N2 = f.n();
  int M1 = u.m();
  int M2 = u.n();
  u.resize(M1,M2);    setvalue(u,cpx(0,0));
  //--------
  Entry& entx1 = _e2dmap[_EPSx1];
  FltNumVec& gridx1 = entx1.grid();
  NumVec<CpxNumMat>& matsx1 = entx1.mats();
  //the transpose matrices
  NumVec<CpxNumMat> tmatsx1(2);
  iC( ztran(matsx1(0), tmatsx1(0)) );
  iC( ztran(matsx1(1), tmatsx1(1)) );
  //
  Entry& entx2 = _e2dmap[_EPSx2]; 
  FltNumVec& gridx2 = entx2.grid();
  NumVec<CpxNumMat>& matsx2 = entx2.mats();
  NumVec<CpxNumMat> tmatsx2(2);
  iC( ztran(matsx2(0), tmatsx2(0)) );
  iC( ztran(matsx2(1), tmatsx2(1)) );
  //	
  Entry& entk1 = _e2dmap[_EPSk1]; 
  FltNumVec& gridk1 = entk1.grid();
  NumVec<CpxNumMat>& matsk1 = entk1.mats();
  NumVec<CpxNumMat> tmatsk1(2);
  iC( ztran(matsk1(0), tmatsk1(0)) );
  iC( ztran(matsk1(1), tmatsk1(1)) );
  //	
  Entry& entk2 = _e2dmap[_EPSk2]; 
  FltNumVec& gridk2 = entk2.grid();
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
  cerr<<"EPSx1 "<<NGx1<<" EPSx2 "<<NGx2<<endl;
  cerr<<"EPSk1 "<<NGk1<<" EPSk2 "<<NGk2<<endl;
  //
  //
  int EL = _EL;
  int TL = int(round(log(float(N))/log(2)));
  int SL = TL-EL;
  int ML = 0;
  if ( (SL+EL)%2 == 0 ) {
    ML = (SL+EL)/2;
  } else {
    if (_fi==1 || _fi==3) {
      // if forward hyper Radon or x*k
      ML = floor((SL+EL)/2.0);
    } else if (_fi==2 || _fi==4) {
      // if adjoint hyper Radon or -x*k
      ML = floor((SL+EL)/2.0)+1;
    }
  }
  cerr<<"EL "<<EL<<" SL "<<SL<<endl;
  cerr<<"ML "<<ML<<endl;
  // grouping things for input
  NumMat<Point2> ks(N1,N2);
  for(int i=0; i<N1; i++)
    for(int j=0; j<N2; j++)
      ks(i,j) = Point2((w(i)-wmin)/(wmax-wmin), (x(j)-xmin)/(xmax-xmin));
  //
  int nk1 = pow2(SL);     int nk2 = nk1;
  float kB1 = 1.0/nk1;    float kB2 = 1.0/nk2;
  NumMat< vector<Index2> > grp(nk1,nk2);
  for(int j=0; j<N2; j++)
    for(int i=0; i<N1; i++) {
      Point2 tk = ks(i,j);
      int d1 = int(floor(tk(0)/kB1));    d1 = min(max(d1,0),nk1-1);
      int d2 = int(floor(tk(1)/kB2));    d2 = min(max(d2,0),nk2-1);
      grp(d1,d2).push_back( Index2(i,j) );
    }
  // number of entries in each box
  vector<IntNumMat> nbv(TL+1);
  for(int ell=SL; ell>=EL; ell--) {
    int nk1 = pow2(ell);    int nk2 = nk1;
    nbv[ell].resize(nk1,nk2);
    for(int k1=0; k1<nk1; k1++)
      for(int k2=0; k2<nk2; k2++) {
	if(ell==SL) {
	  nbv[ell](k1,k2) = grp(k1,k2).size();
	} else {
	  int ttl = 0;
	  for(int a1=0; a1<2; a1++)
	    for(int a2=0; a2<2; a2++) {
	      int c1 = 2*k1+a1;    int c2 = 2*k2+a2;
	      ttl = ttl + nbv[ell+1](c1,c2);
	    }
	  nbv[ell](k1,k2) = ttl;
	}
      }
  }
  // grouping things for output
  NumMat<Point2> xs(M1,M2);
  for(int i=0; i<M1; i++)
    for(int j=0; j<M2; j++)
      xs(i,j) = Point2((tau(i)-taumin)/(taumax-taumin), (p(j)-pmin)/(pmax-pmin));
  //
  int nx1 = pow2(SL);     int nx2 = nx1;
  float xB1 = 1.0/nx1;    float xB2 = 1.0/nx2;
  NumMat< vector<Index2> > grpx(nx1,nx2);
  for(int j=0; j<M2; j++)
    for(int i=0; i<M1; i++) {
      Point2 tx = xs(i,j);
      int d1 = int(floor(tx(0)/xB1));    d1 = min(max(d1,0),nx1-1);
      int d2 = int(floor(tx(1)/xB2));    d2 = min(max(d2,0),nx2-1);
      grpx(d1,d2).push_back( Index2(i,j) );
    }

  //----------------
  clock_t time1 = clock(); 
  float timediff = float(time1-time0)/CLOCKS_PER_SEC; 
  cerr<<"prepare data "<<timediff<<endl;
  //----------------


  //
  //
  //--------------------------------------------------------------------------------------
  int nz1 = pow2(EL);     int nz2 = nz1;
  float zB1 = 1.0/nz1;    float zB2 = 1.0/nz2;
  //
  for(int z1=0; z1<nz1; z1++)
    for(int z2=0; z2<nz2; z2++) {

      //-------------
      time0 = clock(); 
      //-------------

      //cerr<<" "<<z1<<" "<<z2<<endl;
      float k1stt = z1*zB1;    float k1end = (z1+1)*zB1;
      float k2stt = z2*zB2;    float k2end = (z2+1)*zB2;
      NumMat< NumMat< CpxNumMat > > NOW;
      //----------------------------------------------------------------------
      for(int ell=SL; ell>=ML; ell--) {
        int nk1 = pow2(ell);    int nk2 = nk1;
	float kB1 = 1.0/nk1;	float kB2 = 1.0/nk2;
	int nx1 = N/nk1;        int nx2 = N/nk2;
	float xB1 = 1.0/nx1;    float xB2 = 1.0/nx2;
	NumMat< NumMat< CpxNumMat > > PRE = NOW;
	NOW.resize(nk1,nk2);
	for(int k1=int(round(k1stt/kB1)); k1<int(round(k1end/kB1)); k1++)
	  for(int k2=int(round(k2stt/kB2)); k2<int(round(k2end/kB2)); k2++) {
	    NOW(k1,k2).resize(nx1,nx2);
	    for(int x1=0; x1<nx1; x1++)
	      for(int x2=0; x2<nx2; x2++) {
		NOW(k1,k2)(x1,x2).resize(NGk1,NGk2);		
                setvalue(NOW(k1,k2)(x1,x2),cpx(0,0));
	      }
	  }
	//
	vector<Point2> ko;
	for(int j=0; j<NGk2; j++)
	  for(int i=0; i<NGk1; i++)
	    ko.push_back( Point2(gridk1(i), gridk2(j)) );
	//
        for(int k1=int(round(k1stt/kB1)); k1<int(round(k1end/kB1)); k1++)
	  for(int k2=int(round(k2stt/kB2)); k2<int(round(k2end/kB2)); k2++) {
	    if(nbv[ell](k1,k2)==0)    continue;
	    //
	    Point2 kc( (k1+0.5)*kB1, (k2+0.5)*kB2 );
	    for(int x1=0; x1<nx1; x1++)
	      for(int x2=0; x2<nx2; x2++) {
		Point2 xc( (x1+0.5)*xB1, (x2+0.5)*xB2 );
		//-------
		CpxNumMat all(NGk1,NGk2);    setvalue(all, cpx(0,0));
                if(ell==SL) {
		  //get
		  vector<Index2>& gud = grp(k1,k2);
		  CpxNumVec ext(gud.size());
		  for(int g=0; g<int(gud.size()); g++)
		    ext(g) = f(gud[g](0),gud[g](1));
                  // Jingwei: good example for sample nonuniform points
		  //scale
		  vector<Point2> trg;    trg.push_back(xc);
		  vector<Point2> src(gud.size());
		  for(int g=0; g<int(gud.size()); g++)
		    src[g] = ks(gud[g](0),gud[g](1));
		  CpxNumMat scl;    iC( kernel2(N, trg, src, scl) );
		  for(int g=0; g<int(gud.size()); g++)
		    ext(g) = ext(g) * scl(0,g);
		  //transform
		  vector<float> s1(gud.size());
		  vector<float> s2(gud.size());
		  for(int g=0; g<int(gud.size()); g++) {
		    s1[g] = src[g](0)/kB1 - k1;
		    s2[g] = src[g](1)/kB2 - k2;
                    // Jingwei: pay attention to scaling
		  }
		  CpxNumMat tmp1(NGk1, gud.size());
		  iC( prep_aux(gridk1, s1, tmp1) );
		  CpxNumMat tmp2(NGk2, gud.size());
		  iC( prep_aux(gridk2, s2, tmp2) );
                  CpxNumMat ttmp2(gud.size(),NGk2);
		  iC( ztran(tmp2, ttmp2) );
                  //
                  CpxNumMat mid1(NGk1,gud.size());
                  {
                    for(int i=0; i<NGk1; i++)
		      for(int j=0; j<int(gud.size()); j++)
		        mid1(i,j) = tmp1(i,j) * ext(j);
                  }
		  iC( zgemm(1, mid1, ttmp2, 0, all) );
        	} else {
		  int q1 = int(floor(x1/2));    int q2 = int(floor(x2/2));
		  for(int a1=0; a1<2; a1++)
		    for(int a2=0; a2<2; a2++) {
		      int c1 = 2*k1+a1;    int c2 = 2*k2+a2;
                      if(nbv[ell+1](c1,c2)==0)    continue;
		      //
		      //get
		      CpxNumMat ext(PRE(c1,c2)(q1,q2));
		      //scale
		      vector<Point2> trg;    trg.push_back(xc);
		      vector<Point2> src(ko);
		      for(int g=0; g<int(src.size()); g++)
		        src[g] = ewmul(src[g] + Point2(c1,c2), Point2(kB1/2,kB2/2));
		      CpxNumMat scl;    iC( kernel2(N, trg, src, scl) );
		      CpxNumMat sclaux(NGk1,NGk2,false,scl.data());
		      for(int j=0; j<NGk2; j++)
		        for(int i=0; i<NGk1; i++)
			  ext(i,j) = ext(i,j) * sclaux(i,j);
		      //transform
		      CpxNumMat tmp(NGk1,NGk2);    setvalue(tmp,cpx(0,0));
		      iC( zgemm(1, matsk1(a1), ext, 0, tmp) );
		      iC( zgemm(1, tmp, tmatsk2(a2), 1, all) );
		    }
		}
		//scale
		vector<Point2> trg;    trg.push_back(xc);
		vector<Point2> src(ko);
		for(int g=0; g<int(src.size()); g++)
		  src[g] = ewmul(src[g] + Point2(k1,k2), Point2(kB1,kB2));
		CpxNumMat scl;    iC( kernel2(N, trg, src, scl) );
		CpxNumMat sclaux(NGk1,NGk2,false,scl.data());
		for(int j=0; j<NGk2; j++)
		  for(int i=0; i<NGk1; i++)
		    all(i,j) = all(i,j) / sclaux(i,j);
		//put
		NOW(k1,k2)(x1,x2) = all;
	      }//x1x2
	  }//k1k2
        PRE.resize(0,0);
      }//ell

      //-------------------
      time1 = clock(); 
      timediff = float(time1-time0)/CLOCKS_PER_SEC; 
      cerr<<"first half "<<timediff<<endl;
      //-------------------

      //----------------------------------------------------------------------

      //-------------------
      time0 = clock();
      //-------------------

      if(1) {
        int ell = ML;
	int nk1 = pow2(ell);    int nk2 = nk1;
	float kB1 = 1.0/nk1;	float kB2 = 1.0/nk2;
	int nx1 = N/nk1;	int nx2 = N/nk2;
	float xB1 = 1.0/nx1;	float xB2 = 1.0/nx2;
	//
	vector<Point2> ko;
	for(int j=0; j<NGk2; j++)
	  for(int i=0; i<NGk1; i++)
	    ko.push_back( Point2(gridk1(i), gridk2(j)) );
	vector<Point2> xo;
	for(int j=0; j<NGx2; j++)
	  for(int i=0; i<NGx1; i++)
	    xo.push_back( Point2(gridx1(i), gridx2(j)) );
	//
	for(int k1=int(round(k1stt/kB1)); k1<int(round(k1end/kB1)); k1++)
	  for(int k2=int(round(k2stt/kB2)); k2<int(round(k2end/kB2)); k2++) {
	    if(nbv[ell](k1,k2)==0)    continue;
            //
	    Point2 kc( (k1+0.5)*kB1, (k2+0.5)*kB2 );
	    for(int x1=0; x1<nx1; x1++)
	      for(int x2=0; x2<nx2; x2++) {
		Point2 xc( (x1+0.5)*xB1, (x2+0.5)*xB2 );
		//--------
		vector<Point2> src(ko);
		for(int g=0; g<int(src.size()); g++)
		  src[g] = ewmul(src[g] + Point2(k1,k2), Point2(kB1,kB2));
		vector<Point2> trg(xo);
		for(int g=0; g<int(trg.size()); g++)
		  trg[g] = ewmul(trg[g] + Point2(x1,x2), Point2(xB1,xB2));
		CpxNumMat evl(NGx1*NGx2,NGk1*NGk2);    iC( kernel2(N, trg, src, evl) );
		CpxNumVec den(NGk1*NGk2,true,NOW(k1,k2)(x1,x2).data());
		NOW(k1,k2)(x1,x2).resize(NGx1,NGx2);
		CpxNumVec val(NGx1*NGx2,false,NOW(k1,k2)(x1,x2).data());
		iC( zgemv(1, evl, den, 0, val) );
	      }
	  }    
      }//ell

      //---------------
      time1 = clock(); 
      timediff = float(time1-time0)/CLOCKS_PER_SEC; 
      cerr<<"middle step "<<timediff<<endl;
      //---------------
 
      //----------------------------------------------------------------------

      //-------------------
      time0 = clock();
      //-------------------

      for(int ell=ML; ell>=EL; ell--) {
        int nk1 = pow2(ell);    int nk2 = nk1;
	float kB1 = 1.0/nk1;	float kB2 = 1.0/nk2;
	int nx1 = N/nk1;        int nx2 = N/nk2;
	float xB1 = 1.0/nx1;	float xB2 = 1.0/nx2;
	NumMat< NumMat< CpxNumMat > > NXT;
	if(ell!=EL) {
	  NXT.resize(nk1/2,nk2/2);
	  for(int k1=int(round(k1stt/(2*kB1))); k1<int(round(k1end/(2*kB1))); k1++)
	    for(int k2=int(round(k2stt/(2*kB2))); k2<int(round(k2end/(2*kB2))); k2++) {
	      NXT(k1,k2).resize(2*nx1,2*nx2);
	      for(int x1=0; x1<2*nx1; x1++)
		for(int x2=0; x2<2*nx2; x2++) {
		  NXT(k1,k2)(x1,x2).resize(NGx1,NGx2);	    
                  setvalue(NXT(k1,k2)(x1,x2),cpx(0,0));
		}
	    }
	}
	//
	vector<Point2> xo;
	for(int j=0; j<NGx2; j++)
	  for(int i=0; i<NGx1; i++)
	    xo.push_back( Point2(gridx1(i), gridx2(j)) );
	//
	for(int k1=int(round(k1stt/kB1)); k1<int(round(k1end/kB1)); k1++)
	  for(int k2=int(round(k2stt/kB2)); k2<int(round(k2end/kB2)); k2++) {
	    if(nbv[ell](k1,k2)==0)    continue;
	    //
	    Point2 kc( (k1+0.5)*kB1, (k2+0.5)*kB2 );
	    for(int x1=0; x1<nx1; x1++)
	      for(int x2=0; x2<nx2; x2++) {
	        Point2 xc( (x1+0.5)*xB1, (x2+0.5)*xB2 );
		//-------
		//get
		CpxNumMat all(NOW(k1,k2)(x1,x2));
		//scale
		vector<Point2> src;    src.push_back(kc);
		vector<Point2> trg(xo);
		for(int g=0; g<int(trg.size()); g++)
		  trg[g] = ewmul(trg[g] + Point2(x1,x2), Point2(xB1,xB2));
		CpxNumMat scl;    iC( kernel2(N, trg, src, scl) );
		CpxNumMat sclaux(NGx1,NGx2,false,scl.data());
		for(int j=0; j<NGx2; j++)
		  for(int i=0; i<NGx1; i++) {
		    all(i,j) = all(i,j) / sclaux(i,j);
		  }
		//
		if(ell!=EL) {
		  int q1 = int(floor(k1/2));    int q2 = int(floor(k2/2));
		  for(int a1=0; a1<2; a1++)
		    for(int a2=0; a2<2; a2++) {
		      int c1 = 2*x1+a1;    int c2 = 2*x2+a2;
		      //transform
		      CpxNumMat ext(NGx1,NGx2);    setvalue(ext,cpx(0,0));
		      CpxNumMat tmp(NGx1,NGx2);	   setvalue(tmp,cpx(0,0));
		      iC( zgemm(1, tmatsx1(a1), all, 0, tmp) );
		      iC( zgemm(1, tmp, matsx2(a2), 1, ext) );
		      //scale
		      vector<Point2> src;    src.push_back(kc);
		      vector<Point2> trg(xo);
		      for(int g=0; g<int(trg.size()); g++)
		        trg[g] = ewmul(trg[g] + Point2(c1,c2), Point2(xB1/2,xB2/2));
		      CpxNumMat scl;    iC( kernel2(N, trg, src, scl) );
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
                  vector<Index2>& gudx = grpx(x1,x2);
                  vector<Point2> trg(gudx.size());
		  for(int g=0; g<int(gudx.size()); g++)
		    trg[g] = xs(gudx[g](0),gudx[g](1));
                  vector<float> s1(gudx.size());
		  vector<float> s2(gudx.size());
		  for(int g=0; g<int(gudx.size()); g++) {
		    s1[g] = trg[g](0)/xB1 - x1;
		    s2[g] = trg[g](1)/xB2 - x2;
		  }	   
                  CpxNumMat tmp1(NGx1, gudx.size());
		  iC( prep_aux(gridx1, s1, tmp1) );
                  CpxNumMat ttmp1(gudx.size(), NGx1);
                  iC( ztran(tmp1, ttmp1) );
		  CpxNumMat tmp2(NGx2, gudx.size());
		  iC( prep_aux(gridx2, s2, tmp2) );
                  //
                  CpxNumMat mid1(gudx.size(), NGx2);     
                  iC( zgemm(1, ttmp1, all, 0, mid1) );
                  //
                  CpxNumVec ext(gudx.size());    setvalue(ext, cpx(0,0));
                  {
                    for(int i=0; i<int(gudx.size()); i++)
		      for(int j=0; j<NGx2; j++)
		        ext(i) = ext(i) + mid1(i,j) * tmp2(j,i);
		  }
                  //scale
		  vector<Point2> src;    src.push_back(kc);
                  CpxNumMat scl;    iC( kernel2(N, trg, src, scl) );
                  for(int g=0; g<int(gudx.size()); g++)
		    ext(g) = ext(g) * scl(g,0);
	          //put
		  for(int g=0; g<int(gudx.size()); g++)
		    u(gudx[g](0),gudx[g](1)) = u(gudx[g](0),gudx[g](1)) + ext(g);  
		} 
	      }//x1x2
	  }//k1k2 
	NOW = NXT;
	NXT.resize(0,0);
      }//ell

      //----------------
      time1 = clock(); 
      timediff = float(time1-time0)/CLOCKS_PER_SEC; 
      cerr<<"second half "<<timediff<<endl;
      //-----------------

    }//z1z2
  return 0;
}


