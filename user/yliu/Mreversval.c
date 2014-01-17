/* Reverse data value 
   Data' = Max(Data) - Data
*/
/*
  Copyright (C) 2014 Jilin University
  
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
#include <math.h>

#include <rsf.h>

int main(int argc,char *argv[])
{
    int i,j,n1,n2,n3;
    float o1,d1,o2,d2,max;
    float *orig,*output;
    sf_file in,out;
    
    sf_init (argc, argv); 
    in = sf_input("in");
    out = sf_output("out");
    
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    if (!sf_histfloat(in,"d2",&d2)) d2=1;
    if (!sf_histfloat(in,"o2",&o2)) o2=1;
    n3=sf_leftsize(in,2);
    
    orig = sf_floatalloc(n1*n2);
    output = sf_floatalloc(n1*n2);
    for(j = 0; j < n3; j++){
	
	sf_floatread(orig,n1*n2,in);
	max=0;
	for(i=0;i<n1*n2;i++){
	    if(orig[i]>max) max=orig[i];
	}
	
	for(i=0;i<n1*n2;i++){
	    output[i]=max-orig[i];
	}
	
	sf_floatwrite(output,n1*n2,out);
    }
    exit(0);
    
}
/* 	$Id$	 */
