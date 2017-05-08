//   apex shifted 2to3 Radon transform (using 2to2 butterfly)
//   complex f(w,x) --> complex u(tau,p,x)
//   BFIO::setup2 
//   BFIO::apkernel2   fi=1 apex shifted hyper Radon
//   BFIO::apcheck2
//   BFIO::apeval2
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
    
  sf_init(argc,argv); // Initialize RSF
    
  // Get input
  iRSF input;

  int nw, nx;
  input.get("n1",nw);
  input.get("n2",nx);

  float w0, dw;
  input.get("o1",w0);
  input.get("d1",dw);  

  float x0, dx;
  input.get("o2",x0);
  input.get("d2",dx);


  std::valarray<sf_complex> fdata(nw*nx);
  //fdata.resize(nw*nx);

  input >> fdata;

  CpxNumMat f(nw,nx);  setvalue(f,cpx(0,0));
  for (int i=0; i<nw; i++)
    for (int j=0; j<nx; j++)
      f(i,j) = cpx(crealf(fdata[nw*j+i]),cimagf(fdata[nw*j+i]));
  FltNumVec w(nw);
  for (int i=0; i<nw; i++)  w(i) = w0+i*dw;
  FltNumVec x(nx);   
  for (int j=0; j<nx; j++)  x(j) = x0+j*dx;
 

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


  CpxNumTns u(ntau,np,nx);  setvalue(u,cpx(0,0));
  FltNumVec tau(ntau);   
  for (int i=0; i<ntau; i++)  tau(i) = tau0+i*dtau;
  FltNumVec p(np);
  for (int j=0; j<np; j++)  p(j) = p0+j*dp;

  oRSF output;
  output.put("n1",ntau);
  output.put("n2",np);
  output.put("n3",nx);

  output.put("o1",tau0);
  output.put("d1",dtau);

  output.put("o2",p0);
  output.put("d2",dp);

  output.put("o3",x0);
  output.put("d3",dx);
 
  //output.type(SF_FLOAT);
  // this has be there if the input and output types are different

  // BFIO setup
  BFIO bfio("bfio_");
  iC( bfio.setup2(par,input) );

  int N;
  par.get("N",N); // number of partitions
  cerr<<"N "<<N<<endl;

  float time_eval;
  CpxNumMat utmp(ntau,np);   

  {
    int k=0;
    float xx=x(k);
    setvalue(utmp,cpx(0,0));

    ck0 = clock();
    iC( bfio.apeval2(N,f,w,x,utmp,tau,p,xx) );
    ck1 = clock();    
    time_eval = float(ck1-ck0)/CLOCKS_PER_SEC;
    //
    float relerr = 0;
    int NC = 128;
    ck0 = clock();
    iC( bfio.apcheck2(N,f,w,x,utmp,tau,p,xx,NC,relerr) );
    ck1 = clock();
    float time_chck = float(ck1-ck0)/CLOCKS_PER_SEC*float(ntau)*float(np)/float(NC);
    //
    cerr<<"Ta "<<time_eval<<endl;
    cerr<<"Td "<<time_chck<<endl;
    cerr<<"Rt "<<time_chck/time_eval<<endl;
    cerr<<"Ea "<<relerr<<endl;
    //
    for (int i=0; i<ntau; i++)
      for (int j=0; j<np; j++)
	u(i,j,k)=utmp(i,j);
  }

  for (int k=1; k<nx; k++) {
    float xx=x(k);
    setvalue(utmp,cpx(0,0));
    iC( bfio.apeval2(N,f,w,x,utmp,tau,p,xx) );
    cerr<<"k = "<<k<<endl;
    //
   
    for (int i=0; i<ntau; i++)
      for (int j=0; j<np; j++)
	u(i,j,k)=utmp(i,j);
  }

  std::valarray<sf_complex> udata(ntau*np*nx);
  //std::valarray<float> udata(ntau*np*nx);
  //udata.resize(ntau*np*nx);
    
  for (int i=0; i<ntau; i++)
    for (int j=0; j<np; j++)
      for (int k=0; k<nx; k++)
        udata[ntau*np*k+ntau*j+i]=sf_cmplx(real(u(i,j,k)),imag(u(i,j,k)));
        //udata[ntau*np*k+ntau*j+i]=real(u(i,j,k));

  output << udata;

  exit(0);
}
