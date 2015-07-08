/* 3D mutter (only linear). */
/*
  Copyright (C) 2015 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/
#include <rsf.h>
int main(int argc, char* argv[])
{
  int  i1, i2,i3, n1,n2,n3,i1_max;
  float o1,o2,o3,d1,d2,d3,x2,x3,v,*dd,*pp,t_max;
  sf_file inp, out, mask;

  sf_init(argc,argv);
  inp = sf_input("in");
  out = sf_output("out"); 
  
  if(NULL!=sf_getstring("mask")) mask=sf_output("mask");	
    
  if (!sf_getfloat("v",&v)) v=1.5;
  if (!sf_getfloat("x2",&x2)) x2=0.0;
  if (!sf_getfloat("x3",&x3)) x3=0.0;
  	
  if (!sf_histint(inp,"n1",&n1)) sf_error("No n1= in input");
  if (!sf_histint(inp,"n2",&n2)) sf_error("No n2= in input");
  if (!sf_histint(inp,"n3",&n3)) sf_error("No n3= in input");
  if (!sf_histfloat(inp,"d1",&d1)) sf_error("No d1= in input");
  if (!sf_histfloat(inp,"d2",&d2)) sf_error("No d2= in input");  
  if (!sf_histfloat(inp,"d3",&d3)) sf_error("No d3= in input");
  if (!sf_histfloat(inp,"o1",&o1)) sf_error("No o1= in input");
  if (!sf_histfloat(inp,"o2",&o2)) sf_error("No o2= in input");
  if (!sf_histfloat(inp,"o3",&o3)) sf_error("No o3= in input");    

       
  dd = sf_floatalloc(n1);
  pp = sf_floatalloc(n1);
	memset(pp,0,n1*sizeof(float));  
	  
  for (i3=0;i3 < n3;i3++) {
    for (i2=0; i2 < n2; i2++) {
      t_max=sqrtf((o3+i3*d3-x3)*(o3+i3*d3-x3)+(o2+i2*d2-x2)*(o2+i2*d2-x2))/v;
      i1_max=(t_max-o1)/d1+1;
      if(i1_max<0) i1_max=0; if(i1_max>n1) i1_max=n1; 
      if(NULL!=sf_getstring("mask"))	
	 {
      	memset(dd,0,n1*sizeof(float));    
      	for (i1=i1_max; i1 < n1; i1++) {
			dd[i1]=1.0;
      	}
      	sf_floatwrite(dd,n1,mask);
      }
/*      	memset(pp,0,n1*sizeof(float));       */         
      	sf_floatread(pp,n1,inp);
       	for (i1=0; i1 < i1_max; i1++) {
			pp[i1]=0.0;
      	}     	
      	sf_floatwrite(pp,n1,out);      
      
      
    }
  }
  
  exit(0);
}



