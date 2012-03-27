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
int BFIO::eval(int N, const CpxNumMat& f, const FltNumVec& w, const FltNumVec& z, CpxNumMat& u, const FltNumVec& t, const FltNumVec& p)
{
  int N1 = f.m();
  int N2 = f.n();
  int M1 = u.m();
  int M2 = u.n();
  //cerr<<"N1 "<<N1<<" N2 "<<N2<<endl;
  //cerr<<"M1 "<<M1<<" M2 "<<M2<<endl;
  u.resize(M1,M2);
  setvalue(u,cpx(0,0));
  //int Nw = w.m();
  //int Nz = z.m();
  //cerr<<"Nw "<<Nw<<" Nz "<<Nz<<endl;
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
  Entry& entp1 = _e2dmap[_EPSp1]; 
  FltNumVec& gridp1 = entp1.grid();
  NumVec<CpxNumMat>& matsp1 = entp1.mats();
  NumVec<CpxNumMat> tmatsp1(2);
  iC( ztran(matsp1(0), tmatsp1(0)) );
  iC( ztran(matsp1(1), tmatsp1(1)) );
  //	
  Entry& entp2 = _e2dmap[_EPSp2]; 
  FltNumVec& gridp2 = entp2.grid();
  NumVec<CpxNumMat>& matsp2 = entp2.mats();
  NumVec<CpxNumMat> tmatsp2(2);
  iC( ztran(matsp2(0), tmatsp2(0)) );
  iC( ztran(matsp2(1), tmatsp2(1)) );
  //
  //	
  int NGx1 = _EPSx1;
  int NGx2 = _EPSx2;
  int NGp1 = _EPSp1;
  int NGp2 = _EPSp2;
  cerr<<"EPSx1 "<<NGx1<<" EPSx2 "<<NGx2<<endl;
  cerr<<"EPSp1 "<<NGp1<<" EPSp2 "<<NGp2<<endl;
  //
  int MUL = 1;
  //
  int EL = _EL;
  int TL = int(round(log(float(N))/log(2)));
  int SL = TL-EL;
  int ML = int(floor((SL+EL)/2.0));
  cerr<<"EL "<<EL<<" SL "<<SL<<endl;
  // grouping things for input
  NumMat<Point2> ps(N1,N2);
  for(int i=0; i<N1; i++)
    for(int j=0; j<N2; j++)
      ps(i,j) = Point2((w(i)-wmin)/(wmax-wmin), (z(j)-zmin)/(zmax-zmin));
  //
  int np1 = pow2(SL);  int np2 = np1*MUL;
  float pB1 = 1.0/np1;  float pB2 = 1.0/np2;
  NumMat< vector<Index2> > grp(np1,np2);
  for(int j=0; j<N2; j++)
    for(int i=0; i<N1; i++) {
      Point2 tp = ps(i,j);
      int d1 = int(floor(tp(0)/pB1));      d1 = min(max(d1,0),np1-1);
      int d2 = int(floor(tp(1)/pB2));      d2 = min(max(d2,0),np2-1);
      grp(d1,d2).push_back( Index2(i,j) );
    }
  // number of entries in each box
  vector<IntNumMat> nbv(TL+1);
  for(int ell=SL; ell>=EL; ell--) {
    int np1 = pow2(ell);	int np2 = np1*MUL;
    nbv[ell].resize(np1,np2);
    for(int p1=0; p1<np1; p1++)
      for(int p2=0; p2<np2; p2++) {
	if(ell==SL) {
	  nbv[ell](p1,p2) = grp(p1,p2).size();
	} else {
	  int ttl = 0;
	  for(int a1=0; a1<2; a1++)
	    for(int a2=0; a2<2; a2++) {
	      int c1 = 2*p1+a1;		    int c2 = 2*p2+a2;
	      ttl = ttl + nbv[ell+1](c1,c2);
	    }
	  nbv[ell](p1,p2) = ttl;
	}
      }
  }
  // grouping things for output
  NumMat<Point2> xs(M1,M2);
  for(int i=0; i<M1; i++)
    for(int j=0; j<M2; j++)
      xs(i,j) = Point2((t(i)-tmin)/(tmax-tmin), (p(j)-pmin)/(pmax-pmin));
  //
  np1 = pow2(EL);
  int nx1 = N/np1;	int nx2 = nx1;
  float xB1 = 1.0/nx1;	float xB2 = 1.0/nx2;
  NumMat< vector<Index2> > grpx(nx1,nx2);
  for(int j=0; j<M2; j++)
    for(int i=0; i<M1; i++) {
      Point2 tx = xs(i,j);
      int d1 = int(floor(tx(0)/xB1));      d1 = min(max(d1,0),nx1-1);
      int d2 = int(floor(tx(1)/xB2));      d2 = min(max(d2,0),nx2-1);
      grpx(d1,d2).push_back( Index2(i,j) );
    }
  //
  int nz1 = pow2(EL);  int nz2 = nz1*MUL;
  float zB1 = 1.0/nz1;  float zB2 = 1.0/nz2;
  //
  for(int z1=0; z1<nz1; z1++)
    for(int z2=0; z2<nz2; z2++) {
      cerr<<" "<<z1<<" "<<z2<<endl;
      float p1stt = z1*zB1;      float p1end = (z1+1)*zB1;
      float p2stt = z2*zB2;      float p2end = (z2+1)*zB2;
      NumMat< NumMat< CpxNumMat > > NOW;
      //----------------------------------------------------------------------
      for(int ell=SL; ell>=ML; ell--) {
        int np1 = pow2(ell);	int np2 = np1*MUL;
	float pB1 = 1.0/np1;	float pB2 = 1.0/np2;
	int nx1 = N/np1;	int nx2 = nx1;
	float xB1 = 1.0/nx1;	float xB2 = 1.0/nx2;
	NumMat< NumMat< CpxNumMat > > PRE = NOW;
	NOW.resize(np1,np2);
	for(int p1=int(round(p1stt/pB1)); p1<int(round(p1end/pB1)); p1++)
	  for(int p2=int(round(p2stt/pB2)); p2<int(round(p2end/pB2)); p2++) {
	    NOW(p1,p2).resize(nx1,nx2);
	    for(int x1=0; x1<nx1; x1++)
	      for(int x2=0; x2<nx2; x2++) {
		NOW(p1,p2)(x1,x2).resize(NGp1,NGp2);		setvalue(NOW(p1,p2)(x1,x2),cpx(0,0));
	      }
	  }
	//
	vector<Point2> po;
	for(int j=0; j<NGp2; j++)
	  for(int i=0; i<NGp1; i++)
	    po.push_back( Point2(gridp1(i), gridp2(j)) );
	//
        for(int p1=int(round(p1stt/pB1)); p1<int(round(p1end/pB1)); p1++)
	  for(int p2=int(round(p2stt/pB2)); p2<int(round(p2end/pB2)); p2++) {
	    if(nbv[ell](p1,p2)==0)	      continue;
	    //
	    Point2 pc( (p1+0.5)*pB1, (p2+0.5)*pB2 );
	    for(int x1=0; x1<nx1; x1++)
	      for(int x2=0; x2<nx2; x2++) {
		Point2 xc( (x1+0.5)*xB1, (x2+0.5)*xB2 );
		//-------
		CpxNumMat all(NGp1,NGp2);		setvalue(all, cpx(0,0));
                if(ell==SL) {
		  //get
		  vector<Index2>& gud = grp(p1,p2);
		  CpxNumVec ext(gud.size());
		  for(int g=0; g<int(gud.size()); g++)
		    ext(g) = f(gud[g](0),gud[g](1));
                  // Jingwei: good example for sample nonuniform points
		  //scale
		  vector<Point2> trg;		trg.push_back(xc);
		  vector<Point2> src(gud.size());
		  for(int g=0; g<int(gud.size()); g++)
		    src[g] = ps(gud[g](0), gud[g](1));
		  CpxNumMat scl;		iC( kernel(N, trg, src, scl) );
		  for(int g=0; g<int(gud.size()); g++)
		    ext(g) = ext(g) * scl(0,g);
		  //transform
		  vector<float> s1(gud.size());
		  vector<float> s2(gud.size());
		  for(int g=0; g<int(gud.size()); g++) {
		    s1[g] = src[g](0)/pB1 - p1;
		    s2[g] = src[g](1)/pB2 - p2;
                    // Jingwei: pay attention to scaling
		  }
		  CpxNumMat tmp1(NGp1, gud.size());
		  iC( prep_aux(gridp1, s1, tmp1) );
		  CpxNumMat tmp2(NGp2, gud.size());
		  iC( prep_aux(gridp2, s2, tmp2) );
		  for(int j=0; j<int(gud.size()); j++)
		    for(int i=0; i<NGp1; i++)
		      tmp1(i,j) = tmp1(i,j) * ext(j);
		  CpxNumMat tmp3(gud.size(),NGp2);
		  iC( ztran(tmp2, tmp3) );
		  iC( zgemm(1, tmp1, tmp3, 0, all) );
        	} else {
		  int q1 = int(floor(x1/2));		int q2 = int(floor(x2/2));
		  for(int a1=0; a1<2; a1++)
		    for(int a2=0; a2<2; a2++) {
		      int c1 = 2*p1+a1;		    int c2 = 2*p2+a2;
                      if(nbv[ell+1](c1,c2)==0)	      continue;
		      //
		      //get
		      CpxNumMat ext(PRE(c1,c2)(q1,q2));
		      //scale
		      vector<Point2> trg;		    trg.push_back(xc);
		      vector<Point2> src(po);
		      for(int g=0; g<int(src.size()); g++)
		        src[g] = ewmul(src[g] + Point2(c1,c2), Point2(pB1/2,pB2/2));
		      CpxNumMat scl;		    iC( kernel(N, trg, src, scl) );
		      CpxNumMat sclaux(NGp1,NGp2,false,scl.data());
		      for(int j=0; j<NGp2; j++)
		        for(int i=0; i<NGp1; i++)
			  ext(i,j) = ext(i,j) * sclaux(i,j);
		      //transform
		      CpxNumMat tmp(NGp1,NGp2);		setvalue(tmp,cpx(0,0));
		      iC( zgemm(1, matsp1(a1), ext, 0, tmp) );
		      iC( zgemm(1, tmp, tmatsp2(a2), 1, all) );
		    }
		}
		//scale
		vector<Point2> trg;		trg.push_back(xc);
		vector<Point2> src(po);
		for(int g=0; g<int(src.size()); g++)
		  src[g] = ewmul(src[g] + Point2(p1,p2), Point2(pB1,pB2));
		CpxNumMat scl;		iC( kernel(N, trg, src, scl) );
		CpxNumMat sclaux(NGp1,NGp2,false,scl.data());
		for(int j=0; j<NGp2; j++)
		  for(int i=0; i<NGp1; i++)
		    all(i,j) = all(i,j) / sclaux(i,j);
		//put
		NOW(p1,p2)(x1,x2) = all;
	      }//x1x2
	  }//p1p2
        PRE.resize(0,0);
      }//ell
      //----------------------------------------------------------------------
      if(1) {
        int ell = ML;
	int np1 = pow2(ell);	int np2 = np1*MUL;
	float pB1 = 1.0/np1;	float pB2 = 1.0/np2;
	int nx1 = N/np1;	int nx2 = nx1;
	float xB1 = 1.0/nx1;	float xB2 = 1.0/nx2;
	//
	vector<Point2> po;
	for(int j=0; j<NGp2; j++)
	  for(int i=0; i<NGp1; i++)
	    po.push_back( Point2(gridp1(i), gridp2(j)) );
	vector<Point2> xo;
	for(int j=0; j<NGx2; j++)
	  for(int i=0; i<NGx1; i++)
	    xo.push_back( Point2(gridx1(i), gridx2(j)) );
	//
	for(int p1=int(round(p1stt/pB1)); p1<int(round(p1end/pB1)); p1++)
	  for(int p2=int(round(p2stt/pB2)); p2<int(round(p2end/pB2)); p2++) {
	    if(nbv[ell](p1,p2)==0)	      continue;
            //
	    Point2 pc( (p1+0.5)*pB1, (p2+0.5)*pB2 );
	    for(int x1=0; x1<nx1; x1++)
	      for(int x2=0; x2<nx2; x2++) {
		Point2 xc( (x1+0.5)*xB1, (x2+0.5)*xB2 );
		//--------
		vector<Point2> src(po);
		for(int g=0; g<int(src.size()); g++)
		  src[g] = ewmul(src[g] + Point2(p1,p2), Point2(pB1,pB2));
		vector<Point2> trg(xo);
		for(int g=0; g<int(trg.size()); g++)
		  trg[g] = ewmul(trg[g] + Point2(x1,x2), Point2(xB1,xB2));
		CpxNumMat evl(NGx1*NGx2,NGp1*NGp2);		iC( kernel(N, trg, src, evl) );
		CpxNumVec den(NGp1*NGp2,true,NOW(p1,p2)(x1,x2).data());
		NOW(p1,p2)(x1,x2).resize(NGx1,NGx2);
		CpxNumVec val(NGx1*NGx2,false,NOW(p1,p2)(x1,x2).data());
		iC( zgemv(1, evl, den, 0, val) );
	      }
	  }    
      }//ell
      //----------------------------------------------------------------------
      for(int ell=ML; ell>=EL; ell--) {
        int np1 = pow2(ell);	int np2 = np1*MUL;
	float pB1 = 1.0/np1;	float pB2 = 1.0/np2;
	int nx1 = N/np1;	int nx2 = nx1;
	float xB1 = 1.0/nx1;	float xB2 = 1.0/nx2;
	NumMat< NumMat< CpxNumMat > > NXT;
	if(ell!=EL) {
	  NXT.resize(np1/2,np2/2);
	  for(int p1=int(round(p1stt/(2*pB1))); p1<int(round(p1end/(2*pB1))); p1++)
	    for(int p2=int(round(p2stt/(2*pB2))); p2<int(round(p2end/(2*pB2))); p2++) {
	      NXT(p1,p2).resize(2*nx1,2*nx2);
	      for(int x1=0; x1<2*nx1; x1++)
		for(int x2=0; x2<2*nx2; x2++) {
		  NXT(p1,p2)(x1,x2).resize(NGx1,NGx2);		setvalue(NXT(p1,p2)(x1,x2),cpx(0,0));
		}
	    }
	}
	//
	vector<Point2> xo;
	for(int j=0; j<NGx2; j++)
	  for(int i=0; i<NGx1; i++)
	    xo.push_back( Point2(gridx1(i), gridx2(j)) );
	//
	for(int p1=int(round(p1stt/pB1)); p1<int(round(p1end/pB1)); p1++)
	  for(int p2=int(round(p2stt/pB2)); p2<int(round(p2end/pB2)); p2++) {
	    if(nbv[ell](p1,p2)==0)	      continue;
	    //
	    Point2 pc( (p1+0.5)*pB1, (p2+0.5)*pB2 );
	    for(int x1=0; x1<nx1; x1++)
	      for(int x2=0; x2<nx2; x2++) {
	        Point2 xc( (x1+0.5)*xB1, (x2+0.5)*xB2 );
		//-------
		//get
		CpxNumMat all(NOW(p1,p2)(x1,x2));
		//scale
		vector<Point2> src;		src.push_back(pc);
		vector<Point2> trg(xo);
		for(int g=0; g<int(trg.size()); g++)
		  trg[g] = ewmul(trg[g] + Point2(x1,x2), Point2(xB1,xB2));
		CpxNumMat scl;		iC( kernel(N, trg, src, scl) );
		CpxNumMat sclaux(NGx1,NGx2,false,scl.data());
		for(int j=0; j<NGx2; j++)
		  for(int i=0; i<NGx1; i++) {
		    all(i,j) = all(i,j) / sclaux(i,j);
		  }
		//
		if(ell!=EL) {
		  int q1 = int(floor(p1/2));		  int q2 = int(floor(p2/2));
		  for(int a1=0; a1<2; a1++)
		    for(int a2=0; a2<2; a2++) {
		      int c1 = 2*x1+a1;		      int c2 = 2*x2+a2;
		      //transform
		      CpxNumMat ext(NGx1,NGx2);		      setvalue(ext,cpx(0,0));
		      CpxNumMat tmp(NGx1,NGx2);		      setvalue(tmp,cpx(0,0));
		      iC( zgemm(1, tmatsx1(a1), all, 0, tmp) );
		      iC( zgemm(1, tmp, matsx2(a2), 1, ext) );
		      //scale
		      vector<Point2> src;		      src.push_back(pc);
		      vector<Point2> trg(xo);
		      for(int g=0; g<int(trg.size()); g++)
		        trg[g] = ewmul(trg[g] + Point2(c1,c2), Point2(xB1/2,xB2/2));
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
                  vector<Index2>& gudx = grpx(x1,x2);
                  vector<Point2> trg(gudx.size());
		  for(int g=0; g<int(gudx.size()); g++)
		    trg[g] = xs(gudx[g](0), gudx[g](1));
                  vector<float> s1(gudx.size());
		  vector<float> s2(gudx.size());
		  for(int g=0; g<int(gudx.size()); g++) {
		    s1[g] = trg[g](0)/xB1 - x1;
		    s2[g] = trg[g](1)/xB2 - x2;
		  }
                  CpxNumVec ext(gudx.size());    setvalue(ext, cpx(0));	
		  CpxNumMat tmp(gudx.size(), NGx2);       
                  CpxNumMat tmp1(NGx1, gudx.size());
		  iC( prep_aux(gridx1, s1, tmp1) );
		  CpxNumMat tmp2(NGx2, gudx.size());
		  iC( prep_aux(gridx2, s2, tmp2) );
                  CpxNumMat tmp3(gudx.size(), NGx1);
		  iC( ztran(tmp1, tmp3) );
                  iC( zgemm(1, tmp3, all, 0, tmp) );
                  for(int g=0; g<int(gudx.size()); g++)
		    for(int j=0; j<NGx2; j++)
		      ext(g) = ext(g) + tmp(g,j)*tmp2(j,g);
                  //scale
		  vector<Point2> src;		  src.push_back(pc);
                  CpxNumMat scl;		iC( kernel(N, trg, src, scl) );
                  for(int g=0; g<int(gudx.size()); g++)
		    ext(g) = ext(g) * scl(g,0);
	          //put
		  for(int g=0; g<int(gudx.size()); g++)
		    u(gudx[g](0),gudx[g](1)) = u(gudx[g](0),gudx[g](1)) + ext(g);  
		} 
	      }//x1x2
	  }//p1p2
	NOW = NXT;
	NXT.resize(0,0);
      }//ell
    }//z1z2
  return 0;
}


