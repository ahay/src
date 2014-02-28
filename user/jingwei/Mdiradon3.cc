//   direct 3to3 reflection/diffraction Radon transform
//   real f(t,x,y) --> real u(tau,p,q)
//   BFIO::dikernel3    fi=1 direct reflection Radon      
//                      fi=2 direct diffraction Radon
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

  int nt, nx, ny;
  input.get("n1",nt);
  input.get("n2",nx);
  input.get("n3",ny);

  float t0, dt;
  input.get("o1",t0);
  input.get("d1",dt);  

  float x0, dx;
  input.get("o2",x0);
  input.get("d2",dx);

  float y0, dy;
  input.get("o3",y0);
  input.get("d3",dy);

  std::valarray<float> fdata(nt*nx*ny);

  input >> fdata;

  // Set output
  iRSF par(0);
  int fi;
  par.get("fi",fi);
  cerr<<"fi "<<fi<<endl;

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

  std::valarray<float> udata(ntau*np*nq);

  cerr<<"tmin "<<t0<<" tmax "<<t0+nt*dt<<endl;
  cerr<<"xmin "<<x0<<" xmax "<<x0+nx*dx<<endl;
  cerr<<"ymin "<<y0<<" ymax "<<y0+ny*dy<<endl;
  cerr<<"taumin "<<tau0<<" taumax "<<tau0+ntau*dtau<<endl;
  cerr<<"pmin "<<p0<<" pmax "<<p0+np*dp<<endl;
  cerr<<"qmin "<<q0<<" qmax "<<q0+nq*dq<<endl;

  // BFIO setup
  BFIO bfio("bfio_");

  float time_eval1, time_eval2;
  float t, x, y, tau, p, q; 
  int l;

  ck0 = clock(); 
  tt0 = time(0);
  for (int i=0; i<ntau; i++)
    for (int j=0; j<np; j++)
      for (int k=0; k<nq; k++) {
        udata[ntau*np*k+ntau*j+i] = 0.0;
	for (int m=0; m<nx; m++)
	  for (int n=0; n<ny; n++) {
            tau = tau0 + i*dtau;
            p = p0 + j*dp;
            q = q0 + k*dq;
            x = x0 + m*dx;
            y = y0 + n*dy;
            iC( bfio.dikernel3(fi,tau,p,q,x,y,t) );
            l = int(round((t-t0)/dt));
            if (l>=0 && l<nt) {
              udata[ntau*np*k+ntau*j+i] += fdata[nt*nx*n+nt*m+l];    
            }
	  }
      }
  ck1 = clock();    
  tt1 = time(0);  
  time_eval1 = float(ck1-ck0)/CLOCKS_PER_SEC;
  time_eval2 = difftime(tt1,tt0);
  cerr<<"Teval1 "<<time_eval1<<endl;
  cerr<<"Teval2 "<<time_eval2<<endl;
  
  output << udata;

  exit(0);
}
