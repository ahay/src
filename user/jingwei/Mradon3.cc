//   3to3 Radon transform (using 3to3 butterfly)
//   complex f(w,x,y) --> complex u(tau,p,q)
//   BFIO::setup3
//   BFIO::kernel3    fi=0 linear Radon
//                    fi=1 reflection Radon             
//                    fi=2 diffraction Radon
//                    fi=3 adjoint of reflection Radon  
//                    fi=4 adjoint of diffraction Radon
//   BFIO::check3
//   BFIO::eval3
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

  int nw, nx, ny;
  input.get("n1",nw);
  input.get("n2",nx);
  input.get("n3",ny);

  float w0, dw;
  input.get("o1",w0);
  input.get("d1",dw);  

  float x0, dx;
  input.get("o2",x0);
  input.get("d2",dx);

  float y0, dy;
  input.get("o3",y0);
  input.get("d3",dy);


  std::valarray<sf_complex> fdata(nw*nx*ny);
  //fdata.resize(nw*nx*ny);

  input >> fdata;

  CpxNumTns f(nw,nx,ny);  setvalue(f,cpx(0,0));
  for (int i=0; i<nw; i++)
    for (int j=0; j<nx; j++)
      for (int k=0; k<ny; k++)
        f(i,j,k) = cpx(crealf(fdata[nw*nx*k+nw*j+i]),cimagf(fdata[nw*nx*k+nw*j+i]));
  FltNumVec w(nw);
  for (int i=0; i<nw; i++)  w(i) = w0+i*dw;
  FltNumVec x(nx);   
  for (int j=0; j<nx; j++)  x(j) = x0+j*dx;
  FltNumVec y(ny);   
  for (int k=0; k<ny; k++)  y(k) = y0+k*dy;


  // Set output
  iRSF par(0);
  int ntau, np, nq;
  par.get("ntau",ntau); 
  par.get("np",np); 
  par.get("nq",nq); 

  float tau0, dtau;
  par.get("tau0",tau0);
  par.get("dtau",dtau);  

  float p0, dp;
  par.get("p0",p0);
  par.get("dp",dp);

  float q0, dq;
  par.get("q0",q0);
  par.get("dq",dq);
 
  CpxNumTns u(ntau,np,nq);  setvalue(u,cpx(0,0));
  FltNumVec tau(ntau);   
  for (int i=0; i<ntau; i++)  tau(i) = tau0+i*dtau;
  FltNumVec p(np);
  for (int j=0; j<np; j++)  p(j) = p0+j*dp;
  FltNumVec q(nq);
  for (int k=0; k<nq; k++)  q(k) = q0+k*dq;

  oRSF output;
  output.put("n1",ntau);
  output.put("n2",np);
  output.put("n3",nq);

  output.put("o1",tau0);
  output.put("d1",dtau);

  output.put("o2",p0);
  output.put("d2",dp);

  output.put("o3",q0);
  output.put("d3",dq);

  //output.type(SF_FLOAT);
  // this has be there if the input and output types are different

  // BFIO setup
  BFIO bfio("bfio_");
  iC( bfio.setup3(par,input) );

  int N;
  par.get("N",N); // number of partitions
  cerr<<"N "<<N<<endl;


  float time_eval;

  if(N<=256) {
    ck0 = clock();
    iC( bfio.eval3(N,f,w,x,y,u,tau,p,q) );
    ck1 = clock();    
    time_eval = float(ck1-ck0)/CLOCKS_PER_SEC;
  } else {
    tt0 = time(0);
    iC( bfio.eval3(N,f,w,x,y,u,tau,p,q) );
    tt1 = time(0);    
    time_eval = difftime(tt1,tt0);
  }
  //
  float relerr = 0;
  int NC = 128;
  ck0 = clock();
  iC( bfio.check3(N,f,w,x,y,u,tau,p,q,NC,relerr) );
  ck1 = clock();
  float time_chck = float(ck1-ck0)/CLOCKS_PER_SEC*float(ntau)*float(np*nq)/float(NC);
  //
  cerr<<"Ta "<<time_eval<<endl;
  cerr<<"Td "<<time_chck<<endl;
  cerr<<"Rt "<<time_chck/time_eval<<endl;
  cerr<<"Ea "<<relerr<<endl;
  //
  
  std::valarray<sf_complex> udata(ntau*np*nq);
  //std::valarray<float> udata(ntau*np*nq);
  //udata.resize(ntau*np*nq);
    
  for (int i=0; i<ntau; i++)
    for (int j=0; j<np; j++)
      for (int k=0; k<nq; k++)
	udata[ntau*np*k+ntau*j+i]=sf_cmplx(real(u(i,j,k)),imag(u(i,j,k)));
        //udata[ntau*np*k+ntau*j+i]=real(u(i,j,k));

  output << udata;

  exit(0);
}
