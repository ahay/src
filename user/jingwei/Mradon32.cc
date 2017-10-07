//   azimuthally isotropic 3to2 Radon transform (using 2to2 butterfly)
//   complex f(w,x1,x2) --> complex u(tau,p)
//   BFIO::setup32    
//   BFIO::kernel2    fi=1 hyper Radon
//   BFIO::check2
//   BFIO::eval2
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

using namespace std;
using std::cerr;

int main(int argc, char** argv)
{
  srand48(time(NULL));
  clock_t ck0, ck1;
  time_t tt0, tt1;
    
  sf_init(argc,argv); // Initialize RSF
    
  // Get input
  iRSF input;

  int nw, nx1, nx2;
  input.get("n1",nw);
  input.get("n2",nx1);
  input.get("n3",nx2);

  float w0, dw;
  input.get("o1",w0);
  input.get("d1",dw);  

  float x10, dx1;
  input.get("o2",x10);
  input.get("d2",dx1);

  float x20, dx2;
  input.get("o3",x20);
  input.get("d3",dx2);

  
  std::valarray<sf_complex> fdata(nw*nx1*nx2);
  //fdata.resize(nw*nx1*nx2);

  input >> fdata;

  CpxNumMat f(nw,nx1*nx2);  setvalue(f,cpx(0,0));
  for (int i=0; i<nw; i++)
    for (int j=0; j<nx1; j++)
      for (int k=0; k<nx2; k++)
	f(i,nx1*k+j) = cpx(crealf(fdata[nw*nx1*k+nw*j+i]),cimagf(fdata[nw*nx1*k+nw*j+i]));
  FltNumVec w(nw);
  for (int i=0; i<nw; i++)  w(i) = w0+i*dw;
  FltNumVec x(nx1*nx2);   
  for (int j=0; j<nx1; j++)
    for (int k=0; k<nx2; k++) {
      float xx1 = x10+j*dx1;
      float xx2 = x20+k*dx2; 
      x(nx1*k+j) = sqrt(xx1*xx1+xx2*xx2);
    }

  // Set output
  iRSF par(0);
  int ntau, np;
  par.get("ntau",ntau); 
  par.get("np",np); 

  float tau0, dtau;
  par.get("tau0",tau0);
  par.get("dtau",dtau);  

  float p0, dp;
  par.get("p0",p0);
  par.get("dp",dp);
 
  CpxNumMat u(ntau,np);  setvalue(u,cpx(0,0));
  FltNumVec tau(ntau);   
  for (int i=0; i<ntau; i++)  tau(i) = tau0+i*dtau;
  FltNumVec p(np);
  for (int j=0; j<np; j++)  p(j) = p0+j*dp;

  oRSF output;
  output.put("n1",ntau);
  output.put("n2",np);
  output.put("n3",1);

  output.put("o1",tau0);
  output.put("d1",dtau);

  output.put("o2",p0);
  output.put("d2",dp);

  //output.type(SF_FLOAT);
  // this has be there if the input and output types are different

  // BFIO setup
  BFIO bfio("bfio_");
  iC( bfio.setup32(par,input) );

  int N;
  par.get("N",N); // number of partitions
  cerr<<"N "<<N<<endl;


  float time_eval;

  if(N<=256) {
    ck0 = clock();
    iC( bfio.eval2(N,f,w,x,u,tau,p) );
    ck1 = clock();    
    time_eval = float(ck1-ck0)/CLOCKS_PER_SEC;
  } else {
    tt0 = time(0);
    iC( bfio.eval2(N,f,w,x,u,tau,p) );
    tt1 = time(0);    
    time_eval = difftime(tt1,tt0);
  }
  //
  float relerr = 0;
  int NC = 128;
  ck0 = clock();
  iC( bfio.check2(N,f,w,x,u,tau,p,NC,relerr) );
  ck1 = clock();
  float time_chck = float(ck1-ck0)/CLOCKS_PER_SEC*float(ntau)*float(np)/float(NC);
  //
  cerr<<"Ta "<<time_eval<<endl;
  cerr<<"Td "<<time_chck<<endl;
  cerr<<"Rt "<<time_chck/time_eval<<endl;
  cerr<<"Ea "<<relerr<<endl;
  //
  
  std::valarray<sf_complex> udata(ntau*np);
  //std::valarray<float> udata(ntau*np);
  //udata.resize(ntau*np);
    
  for (int i=0; i<ntau; i++)
    for (int j=0; j<np; j++)
      udata[ntau*j+i]=sf_cmplx(real(u(i,j)),imag(u(i,j)));
      //udata[ntau*j+i]=real(u(i,j));

  output << udata;

  exit(0);
}
