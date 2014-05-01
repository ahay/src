//   3to3 butterfly
//   BFIO::prep_aux
//   BFIO::eval_addaux
//   BFIO::eval3       if fi=1 or 2, switch at L/2 (even L), (L-1)/2 (odd L)
//                     if fi=3 or 4, switch at L/2 (even L), (L+1)/2 (odd L)
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
inline int BFIO::eval_addaux(const CpxNumTns& ext, CpxNumTns& all, CpxNumMat& m1, CpxNumMat& m2, CpxNumMat& m3)
{
  int N1 = ext.m();
  int N2 = ext.n();
  int N3 = ext.p();
  int NGk1 = all.m();
  int NGk2 = all.n();
  int NGk3 = all.p();
  //
  CpxNumTns mid1(NGk1,N2,N3);
  {
    CpxNumMat t1(N1,N2*N3,false,ext.data());
    CpxNumMat t2(NGk1,N2*N3,false,mid1.data());
    iC( zgemm(1.0, m1, t1, 0.0, t2) );
  }
  CpxNumTns mid2(N2,N3,NGk1);
  iC( shiftleft(mid1, mid2) );
  //
  CpxNumTns mid3(NGk2,N3,NGk1);
  {
    CpxNumMat t1(N2,N3*NGk1,false,mid2.data());
    CpxNumMat t2(NGk2,N3*NGk1,false,mid3.data());
    iC( zgemm(1.0, m2, t1, 0.0, t2) );
  }
  CpxNumTns mid4(N3,NGk1,NGk2);
  iC( shiftleft(mid3, mid4) );
  //
  CpxNumTns mid5(NGk3,NGk1,NGk2);
  {
    CpxNumMat t1(N3,NGk1*NGk2,false,mid4.data());
    CpxNumMat t2(NGk3,NGk1*NGk2,false,mid5.data());
    iC( zgemm(1.0, m3, t1, 0.0, t2) );
  }
  CpxNumTns mid6(NGk1,NGk2,NGk3);
  iC( shiftleft(mid5,mid6) );
  //
  for(int i=0; i<NGk1; i++)
    for(int j=0; j<NGk2; j++)
      for(int k=0; k<NGk3; k++)
	all(i,j,k) += mid6(i,j,k);
  return 0;
}

//---------------------------------------
int BFIO::eval3(int N, const CpxNumTns& f, const FltNumVec& w, const FltNumVec& x, const FltNumVec& y, CpxNumTns& u, const FltNumVec& tau, const FltNumVec& p, const FltNumVec& q)
{
  //---------------
  clock_t time0 = clock();  
  //--------------
  
  int N1 = f.m();
  int N2 = f.n();
  int N3 = f.p();
  int M1 = u.m();
  int M2 = u.n();
  int M3 = u.p();
  u.resize(M1,M2,M3);    setvalue(u,cpx(0,0));
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
  Entry& entx3 = _e2dmap[_EPSx3]; 
  FltNumVec& gridx3 = entx3.grid();
  NumVec<CpxNumMat>& matsx3 = entx3.mats();
  NumVec<CpxNumMat> tmatsx3(2);
  iC( ztran(matsx3(0), tmatsx3(0)) );
  iC( ztran(matsx3(1), tmatsx3(1)) );
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
  Entry& entk3 = _e2dmap[_EPSk3]; 
  FltNumVec& gridk3 = entk3.grid();
  NumVec<CpxNumMat>& matsk3 = entk3.mats();
  NumVec<CpxNumMat> tmatsk3(2);
  iC( ztran(matsk3(0), tmatsk3(0)) );
  iC( ztran(matsk3(1), tmatsk3(1)) );
  //
  //	
  int NGx1 = _EPSx1;
  int NGx2 = _EPSx2;
  int NGx3 = _EPSx3;
  int NGk1 = _EPSk1;
  int NGk2 = _EPSk2;
  int NGk3 = _EPSk3;
  cerr<<"EPSx1 "<<NGx1<<" EPSx2 "<<NGx2<<" EPSx3 "<<NGx3<<endl;
  cerr<<"EPSk1 "<<NGk1<<" EPSk2 "<<NGk2<<" EPSk3 "<<NGk3<<endl;
  //
  //
  int EL = _EL;
  int TL = int(round(log(float(N))/log(2)));
  int SL = TL-EL;
  int ML = 0;
  if ( (SL+EL)%2 == 0 ) {
    ML = (SL+EL)/2;
  } else {
    if (_fi==1 || _fi==2) {
      // if forward reflection/diffraction
      ML = floor((SL+EL)/2.0);
    } else if (_fi==3 || _fi==4) {
      // if adjoint reflection/diffraction
      ML = floor((SL+EL)/2.0)+1;
    }
  }
  cerr<<"EL "<<EL<<" SL "<<SL<<endl;
  cerr<<"ML "<<ML<<endl;
  // grouping things for input
  NumTns<Point3> ks(N1,N2,N3);
  for(int i=0; i<N1; i++)
    for(int j=0; j<N2; j++)
      for(int k=0; k<N3; k++)
	ks(i,j,k) = Point3((w(i)-wmin)/(wmax-wmin), (x(j)-xmin)/(xmax-xmin), (y(k)-ymin)/(ymax-ymin));
  //
  int nk1 = pow2(SL);     int nk2 = nk1;          int nk3 = nk1;
  float kB1 = 1.0/nk1;    float kB2 = 1.0/nk2;    float kB3 = 1.0/nk3;
  NumTns< vector<Index3> > grp(nk1,nk2,nk3);
  for(int k=0; k<N3; k++)
    for(int j=0; j<N2; j++)
      for(int i=0; i<N1; i++) {
	Point3 tk = ks(i,j,k);
        int d1 = int(floor(tk(0)/kB1));    d1 = min(max(d1,0),nk1-1);
        int d2 = int(floor(tk(1)/kB2));    d2 = min(max(d2,0),nk2-1);
        int d3 = int(floor(tk(2)/kB3));    d3 = min(max(d3,0),nk3-1);
        grp(d1,d2,d3).push_back( Index3(i,j,k) );
      }
  // number of entries in each box
  vector<IntNumTns> nbv(TL+1);
  for(int ell=SL; ell>=EL; ell--) {
    int nk1 = pow2(ell);    int nk2 = nk1;    int nk3 = nk1;
    nbv[ell].resize(nk1,nk2,nk3);
    for(int k1=0; k1<nk1; k1++)
      for(int k2=0; k2<nk2; k2++) 
        for(int k3=0; k3<nk3; k3++) {
	  if(ell==SL) {
	    nbv[ell](k1,k2,k3) = grp(k1,k2,k3).size();
	  } else {
	    int ttl = 0;
	    for(int a1=0; a1<2; a1++)
	      for(int a2=0; a2<2; a2++) 
                for(int a3=0; a3<2; a3++) {
	          int c1 = 2*k1+a1;    int c2 = 2*k2+a2;    int c3 = 2*k3+a3;
	          ttl = ttl + nbv[ell+1](c1,c2,c3);
	        }
	    nbv[ell](k1,k2,k3) = ttl;
	  }
        }
  }
  // grouping things for output
  NumTns<Point3> xs(M1,M2,M3);
  for(int i=0; i<M1; i++)
    for(int j=0; j<M2; j++)
      for(int k=0; k<M3; k++)
	xs(i,j,k) = Point3((tau(i)-taumin)/(taumax-taumin), (p(j)-pmin)/(pmax-pmin), (q(k)-qmin)/(qmax-qmin));
  //
  int nx1 = pow2(SL);     int nx2 = nx1;          int nx3 = nx1;
  float xB1 = 1.0/nx1;    float xB2 = 1.0/nx2;    float xB3 = 1.0/nx3;
  NumTns< vector<Index3> > grpx(nx1,nx2,nx3);
  for(int k=0; k<M3; k++)
    for(int j=0; j<M2; j++)
      for(int i=0; i<M1; i++) {
	Point3 tx = xs(i,j,k);
        int d1 = int(floor(tx(0)/xB1));    d1 = min(max(d1,0),nx1-1);
        int d2 = int(floor(tx(1)/xB2));    d2 = min(max(d2,0),nx2-1);
        int d3 = int(floor(tx(2)/xB3));    d3 = min(max(d3,0),nx3-1);
        grpx(d1,d2,d3).push_back( Index3(i,j,k) );
      }
  
  //----------------
  clock_t time1 = clock(); 
  float timediff = float(time1-time0)/CLOCKS_PER_SEC; 
  cerr<<"prepare data "<<timediff<<endl;
  //----------------

  //
  //
  //--------------------------------------------------------------------------------------
  int nz1 = pow2(EL);     int nz2 = nz1;          int nz3 = nz1;
  float zB1 = 1.0/nz1;    float zB2 = 1.0/nz2;    float zB3 = 1.0/nz3;
  //
  for(int z1=0; z1<nz1; z1++)
    for(int z2=0; z2<nz2; z2++) 
      for(int z3=0; z3<nz3; z3++) {

        //-------------
        time0 = clock(); 
        //-------------

        //cerr<<" "<<z1<<" "<<z2<<" "<<z3<<endl;
        float k1stt = z1*zB1;    float k1end = (z1+1)*zB1;
        float k2stt = z2*zB2;    float k2end = (z2+1)*zB2;
        float k3stt = z3*zB3;    float k3end = (z3+1)*zB3;
        NumTns< NumTns< CpxNumTns > > NOW;
        //----------------------------------------------------------------------
        for(int ell=SL; ell>=ML; ell--) {
          int nk1 = pow2(ell);    int nk2 = nk1;          int nk3 = nk1;
	  float kB1 = 1.0/nk1;    float kB2 = 1.0/nk2;    float kB3 = 1.0/nk3;
	  int nx1 = N/nk1;        int nx2 = N/nk2;        int nx3 = N/nk3;
	  float xB1 = 1.0/nx1;    float xB2 = 1.0/nx2;    float xB3 = 1.0/nx3;
	  NumTns< NumTns< CpxNumTns > > PRE = NOW;
	  NOW.resize(nk1,nk2,nk3);
	  for(int k1=int(round(k1stt/kB1)); k1<int(round(k1end/kB1)); k1++)
	    for(int k2=int(round(k2stt/kB2)); k2<int(round(k2end/kB2)); k2++) 
              for(int k3=int(round(k3stt/kB3)); k3<int(round(k3end/kB3)); k3++) {
	        NOW(k1,k2,k3).resize(nx1,nx2,nx3);
	        for(int x1=0; x1<nx1; x1++)
	          for(int x2=0; x2<nx2; x2++)
                    for(int x3=0; x3<nx3; x3++) {
		      NOW(k1,k2,k3)(x1,x2,x3).resize(NGk1,NGk2,NGk3);		
                      setvalue(NOW(k1,k2,k3)(x1,x2,x3),cpx(0,0));
	            }
	      }
	  //
	  vector<Point3> ko;
          for(int k=0; k<NGk3; k++)
	    for(int j=0; j<NGk2; j++)
	      for(int i=0; i<NGk1; i++)
	        ko.push_back( Point3(gridk1(i), gridk2(j), gridk3(k)) );
	  //
          for(int k1=int(round(k1stt/kB1)); k1<int(round(k1end/kB1)); k1++)
	    for(int k2=int(round(k2stt/kB2)); k2<int(round(k2end/kB2)); k2++) 
              for(int k3=int(round(k3stt/kB3)); k3<int(round(k3end/kB3)); k3++) {
	        if(nbv[ell](k1,k2,k3)==0)    continue;
	        //
	        Point3 kc( (k1+0.5)*kB1, (k2+0.5)*kB2, (k3+0.5)*kB3 );
	        for(int x1=0; x1<nx1; x1++)
	          for(int x2=0; x2<nx2; x2++) 
                    for(int x3=0; x3<nx3; x3++) {
		      Point3 xc( (x1+0.5)*xB1, (x2+0.5)*xB2, (x3+0.5)*xB3 );
		      //-------
		      CpxNumTns all(NGk1,NGk2,NGk3);    setvalue(all, cpx(0,0));
                      if(ell==SL) {
		        //get
		        vector<Index3>& gud = grp(k1,k2,k3);
		        CpxNumVec ext(gud.size());
		        for(int g=0; g<int(gud.size()); g++)
			  ext(g) = f(gud[g](0),gud[g](1),gud[g](2));
                        // Jingwei: good example for sample nonuniform points
		        //scale
		        vector<Point3> trg;    trg.push_back(xc);
		        vector<Point3> src(gud.size());
		        for(int g=0; g<int(gud.size()); g++)
			  src[g] = ks(gud[g](0),gud[g](1),gud[g](2));
		        CpxNumMat scl;    iC( kernel3(N, trg, src, scl) );
		        for(int g=0; g<int(gud.size()); g++)
		          ext(g) = ext(g) * scl(0,g);
		        //transform
		        vector<float> s1(gud.size());
		        vector<float> s2(gud.size());
                        vector<float> s3(gud.size());
		        for(int g=0; g<int(gud.size()); g++) {
		          s1[g] = src[g](0)/kB1 - k1;
		          s2[g] = src[g](1)/kB2 - k2;
                          s3[g] = src[g](2)/kB3 - k3;
                          // Jingwei: pay attention to scaling
		        }
		        CpxNumMat tmp1(NGk1, gud.size());
		        iC( prep_aux(gridk1, s1, tmp1) );
		        CpxNumMat tmp2(NGk2, gud.size());
		        iC( prep_aux(gridk2, s2, tmp2) );
                        CpxNumMat tmp3(NGk3, gud.size());
		        iC( prep_aux(gridk3, s3, tmp3) );
                        //
                        CpxNumMat mid1(NGk1,gud.size());
                        { 
                          for(int i=0; i<NGk1; i++)
                            for(int j=0; j<int(gud.size()); j++)   
		              mid1(i,j) = tmp1(i,j) * ext(j);
		        }
                        CpxNumMat mid2(gud.size(),NGk1);
		        iC( ztran(mid1, mid2) );
                        //
                        CpxNumTns mid3(NGk2,gud.size(),NGk1);
                        {
                          for(int i=0; i<NGk2; i++)
                            for(int j=0; j<int(gud.size()); j++)
                              for(int k=0; k<NGk1; k++)
			        mid3(i,j,k) = tmp2(i,j) * mid2(j,k);
                        }
                        CpxNumTns mid4(gud.size(),NGk1,NGk2);
                        iC( shiftleft(mid3,mid4) );
                        //
                        CpxNumTns mid5(NGk3,NGk1,NGk2);	
                        {	
                          CpxNumMat t1(gud.size(),NGk1*NGk2,false,mid4.data());
                          CpxNumMat t2(NGk3,NGk1*NGk2,false,mid5.data());
                          iC( zgemm(1.0, tmp3, t1, 0.0, t2) );       
                        }    
                        iC( shiftleft(mid5,all) );
        	      } else {
		        int q1 = int(floor(x1/2));    int q2 = int(floor(x2/2));    int q3 = int(floor(x3/2));
		        for(int a1=0; a1<2; a1++)
		          for(int a2=0; a2<2; a2++) 
                            for(int a3=0; a3<2; a3++) {
			      int c1 = 2*k1+a1;    int c2 = 2*k2+a2;    int c3 = 2*k3+a3;
                              if(nbv[ell+1](c1,c2,c3)==0)    continue;
		              //
		              //get
		              CpxNumTns ext(PRE(c1,c2,c3)(q1,q2,q3));
		              //scale
		              vector<Point3> trg;    trg.push_back(xc);
		              vector<Point3> src(ko);
		              for(int g=0; g<int(src.size()); g++)
			        src[g] = ewmul(src[g] + Point3(c1,c2,c3), Point3(kB1/2,kB2/2,kB3/2));
		              CpxNumMat scl;    iC( kernel3(N, trg, src, scl) );
		              CpxNumTns sclaux(NGk1,NGk2,NGk3,false,scl.data());
                              for(int k=0; k<NGk3; k++)
		                for(int j=0; j<NGk2; j++)
		                  for(int i=0; i<NGk1; i++)
				    ext(i,j,k) = ext(i,j,k) * sclaux(i,j,k);
		              //transform
                              iC( eval_addaux(ext, all, matsk1(a1), matsk2(a2), matsk3(a3)) );
		            }
		      }
		      //scale
		      vector<Point3> trg;    trg.push_back(xc);
		      vector<Point3> src(ko);
		      for(int g=0; g<int(src.size()); g++)
		        src[g] = ewmul(src[g] + Point3(k1,k2,k3), Point3(kB1,kB2,kB3));
		      CpxNumMat scl;    iC( kernel3(N, trg, src, scl) );
		      CpxNumTns sclaux(NGk1,NGk2,NGk3,false,scl.data());
                      for(int k=0; k<NGk3; k++)
		        for(int j=0; j<NGk2; j++)
		          for(int i=0; i<NGk1; i++)
			    all(i,j,k) = all(i,j,k) / sclaux(i,j,k);
		      //put
		      NOW(k1,k2,k3)(x1,x2,x3) = all;
	            }//x1x2x3
	      }//k1k2k3
	  PRE.resize(0,0,0);
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
	  int nk1 = pow2(ell);    int nk2 = nk1;          int nk3 = nk1;
	  float kB1 = 1.0/nk1;	  float kB2 = 1.0/nk2;    float kB3 = 1.0/nk3;
	  int nx1 = N/nk1;	  int nx2 = N/nk2;        int nx3 = N/nk3;
	  float xB1 = 1.0/nx1;	  float xB2 = 1.0/nx2;    float xB3 = 1.0/nx3;
	  //
	  vector<Point3> ko;
          for(int k=0; k<NGk3; k++)
	    for(int j=0; j<NGk2; j++)
	      for(int i=0; i<NGk1; i++)
		ko.push_back( Point3(gridk1(i), gridk2(j), gridk3(k)) );
	  vector<Point3> xo;
          for(int k=0; k<NGx3; k++)
	    for(int j=0; j<NGx2; j++)
	      for(int i=0; i<NGx1; i++)
		xo.push_back( Point3(gridx1(i), gridx2(j), gridx3(k)) );
	  //
	  for(int k1=int(round(k1stt/kB1)); k1<int(round(k1end/kB1)); k1++)
	    for(int k2=int(round(k2stt/kB2)); k2<int(round(k2end/kB2)); k2++) 
              for(int k3=int(round(k3stt/kB3)); k3<int(round(k3end/kB3)); k3++) {
		if(nbv[ell](k1,k2,k3)==0)    continue;
                //
		Point3 kc( (k1+0.5)*kB1, (k2+0.5)*kB2, (k3+0.5)*kB3 );
	        for(int x1=0; x1<nx1; x1++)
	          for(int x2=0; x2<nx2; x2++) 
                    for(int x3=0; x3<nx3; x3++) {
		      Point3 xc( (x1+0.5)*xB1, (x2+0.5)*xB2, (x3+0.5)*xB3 );
		      //--------
		      vector<Point3> src(ko);
		      for(int g=0; g<int(src.size()); g++)
			src[g] = ewmul(src[g] + Point3(k1,k2,k3), Point3(kB1,kB2,kB3));
		      vector<Point3> trg(xo);
		      for(int g=0; g<int(trg.size()); g++)
			trg[g] = ewmul(trg[g] + Point3(x1,x2,x3), Point3(xB1,xB2,xB3));
		      CpxNumMat evl(NGx1*NGx2*NGx3,NGk1*NGk2*NGk3);    iC( kernel3(N, trg, src, evl) );
		      CpxNumVec den(NGk1*NGk2*NGk3,true,NOW(k1,k2,k3)(x1,x2,x3).data());
		      NOW(k1,k2,k3)(x1,x2,x3).resize(NGx1,NGx2,NGx3);
		      CpxNumVec val(NGx1*NGx2*NGx3,false,NOW(k1,k2,k3)(x1,x2,x3).data());
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
          int nk1 = pow2(ell);    int nk2 = nk1;          int nk3 = nk1;
	  float kB1 = 1.0/nk1;	  float kB2 = 1.0/nk2;    float kB3 = 1.0/nk3;
	  int nx1 = N/nk1;        int nx2 = N/nk2;        int nx3 = N/nk3;
	  float xB1 = 1.0/nx1;	  float xB2 = 1.0/nx2;    float xB3 = 1.0/nx3;
	  NumTns< NumTns< CpxNumTns > > NXT;
	  if(ell!=EL) {
	    NXT.resize(nk1/2,nk2/2,nk3/2);
	    for(int k1=int(round(k1stt/(2*kB1))); k1<int(round(k1end/(2*kB1))); k1++)
	      for(int k2=int(round(k2stt/(2*kB2))); k2<int(round(k2end/(2*kB2))); k2++)
                for(int k3=int(round(k3stt/(2*kB3))); k3<int(round(k3end/(2*kB3))); k3++) {
		  NXT(k1,k2,k3).resize(2*nx1,2*nx2,2*nx3);
	          for(int x1=0; x1<2*nx1; x1++)
		    for(int x2=0; x2<2*nx2; x2++) 
                      for(int x3=0; x3<2*nx3; x3++) {
			NXT(k1,k2,k3)(x1,x2,x3).resize(NGx1,NGx2,NGx3);	    
			setvalue(NXT(k1,k2,k3)(x1,x2,x3),cpx(0,0));
		    }
	        }
	  }
	  //
	  vector<Point3> xo;
          for(int k=0; k<NGx3; k++)
	    for(int j=0; j<NGx2; j++)
	      for(int i=0; i<NGx1; i++)
		xo.push_back( Point3(gridx1(i), gridx2(j), gridx3(k)) );
	  //
	  for(int k1=int(round(k1stt/kB1)); k1<int(round(k1end/kB1)); k1++)
	    for(int k2=int(round(k2stt/kB2)); k2<int(round(k2end/kB2)); k2++) 
              for(int k3=int(round(k3stt/kB3)); k3<int(round(k3end/kB3)); k3++) {
		if(nbv[ell](k1,k2,k3)==0)    continue;
	        //
	        Point3 kc( (k1+0.5)*kB1, (k2+0.5)*kB2, (k3+0.5)*kB3 );
	        for(int x1=0; x1<nx1; x1++)
	          for(int x2=0; x2<nx2; x2++) 
                    for(int x3=0; x3<nx3; x3++) {
		      Point3 xc( (x1+0.5)*xB1, (x2+0.5)*xB2, (x3+0.5)*xB3 );
		      //-------
		      //get
	  	      CpxNumTns all(NOW(k1,k2,k3)(x1,x2,x3));
		      //scale
		      vector<Point3> src;    src.push_back(kc);
		      vector<Point3> trg(xo);
		      for(int g=0; g<int(trg.size()); g++)
			trg[g] = ewmul(trg[g] + Point3(x1,x2,x3), Point3(xB1,xB2,xB3));
		      CpxNumMat scl;    iC( kernel3(N, trg, src, scl) );
		      CpxNumTns sclaux(NGx1,NGx2,NGx3,false,scl.data());
                      for(int k=0; k<NGx3; k++)
	                for(int j=0; j<NGx2; j++)
		          for(int i=0; i<NGx1; i++) {
			    all(i,j,k) = all(i,j,k) / sclaux(i,j,k);
		          }
		      //
		      if(ell!=EL) {
		        int q1 = int(floor(k1/2));    int q2 = int(floor(k2/2));    int q3 = int(floor(k3/2));
		        for(int a1=0; a1<2; a1++)
		          for(int a2=0; a2<2; a2++) 
                            for(int a3=0; a3<2; a3++) {
		              int c1 = 2*x1+a1;    int c2 = 2*x2+a2;    int c3 = 2*x3+a3;
		              //transform
		              CpxNumTns ext(NGx1,NGx2,NGx3);    setvalue(ext,cpx(0,0));
                              iC( eval_addaux(all, ext, tmatsx1(a1), tmatsx2(a2), tmatsx3(a3)) );
		              //scale
		              vector<Point3> src;    src.push_back(kc);
		              vector<Point3> trg(xo);
		              for(int g=0; g<int(trg.size()); g++)
				trg[g] = ewmul(trg[g] + Point3(c1,c2,c3), Point3(xB1/2,xB2/2,xB3/2));
		              CpxNumMat scl;    iC( kernel3(N, trg, src, scl) );
		              CpxNumTns sclaux(NGx1,NGx2,NGx3,false,scl.data());
                              for(int k=0; k<NGx3; k++)
		                for(int j=0; j<NGx2; j++)
			          for(int i=0; i<NGx1; i++)
				    ext(i,j,k) = ext(i,j,k) * sclaux(i,j,k);
		              //put
                              for(int k=0; k<NGx3; k++)
		                for(int j=0; j<NGx2; j++)
			          for(int i=0; i<NGx1; i++)
				    NXT(q1,q2,q3)(c1,c2,c3)(i,j,k) = NXT(q1,q2,q3)(c1,c2,c3)(i,j,k) + ext(i,j,k);
		            }
		      } else {
                        //transform
			vector<Index3>& gudx = grpx(x1,x2,x3);
                        vector<Point3> trg(gudx.size());
		        for(int g=0; g<int(gudx.size()); g++)
			  trg[g] = xs(gudx[g](0),gudx[g](1),gudx[g](2));
                        vector<float> s1(gudx.size());
		        vector<float> s2(gudx.size());
                        vector<float> s3(gudx.size());
		        for(int g=0; g<int(gudx.size()); g++) {
		          s1[g] = trg[g](0)/xB1 - x1;
		          s2[g] = trg[g](1)/xB2 - x2;
                          s3[g] = trg[g](2)/xB3 - x3;
		        }   
                        CpxNumMat tmp1(NGx1, gudx.size());
		        iC( prep_aux(gridx1, s1, tmp1) );
                        CpxNumMat ttmp1(gudx.size(), NGx1);
                        iC( ztran(tmp1, ttmp1) );
		        CpxNumMat tmp2(NGx2, gudx.size());
		        iC( prep_aux(gridx2, s2, tmp2) );
                        CpxNumMat tmp3(NGx3, gudx.size());
		        iC( prep_aux(gridx3, s3, tmp3) );
                        //
                        CpxNumTns mid1(gudx.size(),NGx2,NGx3);	
                        {	
                          CpxNumMat t1(NGx1,NGx2*NGx3,false,all.data());
                          CpxNumMat t2(gudx.size(),NGx2*NGx3,false,mid1.data());
                          iC( zgemm(1.0, ttmp1, t1, 0.0, t2) );       
                        }    
                        CpxNumTns mid2(NGx2,NGx3,gudx.size());
                        iC( shiftleft(mid1,mid2) );
                        //
                        CpxNumMat mid3(gudx.size(),NGx3);    setvalue(mid3, cpx(0,0));
                        {
                          for(int i=0; i<int(gudx.size()); i++)
                            for(int j=0; j<NGx2; j++)
                              for(int k=0; k<NGx3; k++)
				mid3(i,k) = mid3(i,k) + tmp2(j,i) * mid2(j,k,i);
			}
                        //
                        CpxNumVec ext(gudx.size());    setvalue(ext, cpx(0,0));	
                        {
                          for(int i=0; i<int(gudx.size()); i++)
                            for(int j=0; j<NGx3; j++)
			      ext(i) = ext(i) + mid3(i,j) * tmp3(j,i);
			}
                        //scale
		        vector<Point3> src;    src.push_back(kc);
                        CpxNumMat scl;    iC( kernel3(N, trg, src, scl) );
                        for(int g=0; g<int(gudx.size()); g++)
		          ext(g) = ext(g) * scl(g,0);
	                //put
		        for(int g=0; g<int(gudx.size()); g++)
			  u(gudx[g](0),gudx[g](1),gudx[g](2)) = u(gudx[g](0),gudx[g](1),gudx[g](2)) + ext(g);
		      } 
	            }//x1x2x3
	      }//k1k2k3 
	  NOW = NXT;
	  NXT.resize(0,0,0);
        }//ell

        //----------------
        time1 = clock(); 
        timediff = float(time1-time0)/CLOCKS_PER_SEC; 
        cerr<<"second half "<<timediff<<endl;
        //-----------------
        
      }//z1z2z3
  return 0;
}


