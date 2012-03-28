// fft by hu perform on 1st axis of 1d data
// input complex; output real
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
  sf_init(argc,argv); // Initialize RSF
    
  // Get input
  iRSF input;

  int nw;
  input.get("n1",nw);

  float w0, dw;
  input.get("o1",w0);
  input.get("d1",dw);  

  std::valarray<sf_complex> fdata(nw);
  input >> fdata;
   
  CpxNumVec f(nw);  setvalue(f,cpx(0,0));
  for (int k=0; k<nw; k++)
    f(k) = cpx(crealf(fdata[k]),cimagf(fdata[k]));


  // Set output
  iRSF par(0);
  int inv;
  par.get("inv",inv);
  
  oRSF output;

  int nt;
  par.get("nt",nt); 

  float t0, dt;
  par.get("t0",t0);
  par.get("dt",dt);  

  output.put("n1",nt);
 
  output.put("o1",t0);
  output.put("d1",dt);

  output.type(SF_FLOAT);
  

  CpxNumVec u(nt);  setvalue(u,cpx(0,0));
    for (int j=0; j<nt; j++) {
      for (int k=0; k<nw; k++) {
	//float phs = 2*M_PI*(k*dw)*(j*dt);
        float phs = 2*M_PI*(w0+k*dw)*(t0+j*dt);
        cpx res(0,0);
        float cc=cos(phs);
        float ss=sin(phs);
        if (inv==0) {
          res = cpx(cc,-ss);
        } else {
          res = cpx(cc,ss);
        } 
        u(j)=u(j)+res*f(k);
      }
    }

  std::valarray<float> udata(nt);
  for (int j=0; j<nt; j++)
    udata[j]=real(u(j));

  output << udata;

  exit(0);
}
