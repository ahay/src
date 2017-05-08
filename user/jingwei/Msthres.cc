//   soft thresholding function
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
  sf_init(argc,argv); // Initialize RSF
    
  // Get input
  iRSF input;
  int nt, nx;
  input.get("n1",nt);
  input.get("n2",nx);

  std::valarray<float> fdata(nt*nx);

  input >> fdata;

  // Set output
  iRSF par(0);
  float miu;
  par.get("miu",miu);

  oRSF output;
  
  for (int i=0; i<nt*nx; i++) {
    if (fdata[i]<=-miu/2.0) {
      fdata[i]=fdata[i]+miu/2.0;
    } else if (fdata[i]>=miu/2.0) {
      fdata[i]=fdata[i]-miu/2.0;
    } else {
      fdata[i]=0;
    }
  }
  
  output << fdata;

  exit(0);
}
