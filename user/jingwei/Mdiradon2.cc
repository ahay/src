//   direct 2to2 hyper Radon transform (single integral, nearest point interpolation)
//   real f(t,x) --> real u(tau,p)
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

  int nt, nx;
  input.get("n1",nt);
  input.get("n2",nx);

  float t0, dt;
  input.get("o1",t0);
  input.get("d1",dt);  

  float x0, dx;
  input.get("o2",x0);
  input.get("d2",dx);

  std::valarray<float> fdata(nt*nx);

  input >> fdata;

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

  oRSF output;
  output.put("n1",ntau);
  output.put("n2",np);

  output.put("o1",tau0);
  output.put("d1",dtau);

  output.put("o2",p0);
  output.put("d2",dp);

  std::valarray<float> udata(ntau*np);

  cerr<<"tmin "<<t0<<" tmax "<<t0+nt*dt<<endl;
  cerr<<"xmin "<<x0<<" xmax "<<x0+nx*dx<<endl;
  cerr<<"taumin "<<tau0<<" taumax "<<tau0+ntau*dtau<<endl;
  cerr<<"pmin "<<p0<<" pmax "<<p0+np*dp<<endl;


  float time_eval1, time_eval2;
  float t, x, tau, p; 
  int m;

  ck0 = clock(); 
  tt0 = time(0);
  for (int i=0; i<ntau; i++)
    for (int j=0; j<np; j++) {
      udata[ntau*j+i] = 0.0;
      for (int n=0; n<nx; n++) {
        tau = tau0 + i*dtau;
        p = p0 + j*dp;
        x = x0 + n*dx;
        t = sqrt(tau*tau+p*p*x*x);
        m = int(round((t-t0)/dt));
        if (m>=0 && m<nt)
          udata[ntau*j+i] += fdata[nt*n+m];   
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
