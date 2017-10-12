//   adjoint test
//   check <f1,f2>?=<g1,g2>
//   complex f1(t,x,y), f2(t,x,y), g1(t,x,y), g2(t,x,y)
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
  iRSF input1;
  int nt1, nx1, ny1;
  input1.get("n1",nt1);
  input1.get("n2",nx1);
  input1.get("n3",ny1);

  iRSF input2("input2");
  int nt2, nx2, ny2;
  input2.get("n1",nt2);
  input2.get("n2",nx2);
  input2.get("n3",ny2);

  if (nt1!=nt2 || nx1!=nx2 || ny1!=ny2)
    cerr<<"wrong size of f1, f2"<<endl; 

  iRSF input3("input3");
  int nt3, nx3, ny3;
  input3.get("n1",nt3);
  input3.get("n2",nx3);
  input3.get("n3",ny3);

  iRSF input4("input4");
  int nt4, nx4, ny4;
  input4.get("n1",nt4);
  input4.get("n2",nx4);
  input4.get("n3",ny4);

  if (nt3!=nt4 || nx3!=nx4 || ny3!=ny4)
    cerr<<"wrong size of g1, g2"<<endl; 


  std::valarray<sf_complex> f1data(nt1*nx1*ny1);
  std::valarray<sf_complex> f2data(nt2*nx2*ny2);

  input1 >> f1data;
  input2 >> f2data;
 
  double sumreal1=0.0;
  double sumimag1=0.0;
  for (int i=0; i<nt1; i++)
    for (int j=0; j<nx1; j++) 
      for (int k=0; k<ny1; k++) {
        double a=crealf(f1data[nt1*nx1*k+nt1*j+i]);
        double b=cimagf(f1data[nt1*nx1*k+nt1*j+i]);
        double c=crealf(f2data[nt1*nx1*k+nt1*j+i]);
        double d=cimagf(f2data[nt1*nx1*k+nt1*j+i]);
        sumreal1=sumreal1+a*c+b*d;
        sumimag1=sumimag1+b*c-a*d;
      }
  double sumabs1=sqrt(sumreal1*sumreal1+sumimag1*sumimag1);
  cerr<<"sumreal1 "<<sumreal1<<" sumimag1 "<<sumimag1<<endl;    
  cerr<<"sumabs1 "<<sumabs1<<endl;  


  std::valarray<sf_complex> g1data(nt3*nx3*ny3);
  std::valarray<sf_complex> g2data(nt4*nx4*ny4);

  input3 >> g1data;
  input4 >> g2data;
 
  double sumreal2=0.0;
  double sumimag2=0.0;
  for (int i=0; i<nt3; i++)
    for (int j=0; j<nx3; j++) 
      for (int k=0; k<ny3; k++) {
        double a=crealf(g1data[nt3*nx3*k+nt3*j+i]);
        double b=cimagf(g1data[nt3*nx3*k+nt3*j+i]);
        double c=crealf(g2data[nt3*nx3*k+nt3*j+i]);
        double d=cimagf(g2data[nt3*nx3*k+nt3*j+i]);
        sumreal2=sumreal2+a*c+b*d;
        sumimag2=sumimag2+b*c-a*d;
      }
  double sumabs2=sqrt(sumreal2*sumreal2+sumimag2*sumimag2);
  cerr<<"sumreal2 "<<sumreal2<<" sumimag2 "<<sumimag2<<endl;  
  cerr<<"sumabs2 "<<sumabs2<<endl;  

  double diffreal=sumreal1-sumreal2;
  double diffimag=sumimag1-sumimag2;
  double relerr=sqrt(diffreal*diffreal+diffimag*diffimag)/sumabs1;
  cerr<<"relerr "<<relerr<<endl;  

  // Set output
  //oRSF output;
  //output << f1data;

  exit(0);
}
