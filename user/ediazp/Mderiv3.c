/* Second order derivative along axis
  
   int axis=[1] axis to differentiate
   
   int operator=[2] 1 backward, 2 centered, 3 forward

   
 */
/*
  Copyright (C) 2011 Colorado School of Mines

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

int main (int argc, char* argv[])
{
    int n1,n2,n3;
	int i3,i2,i1;
	int axis,n_axis;
	int operator;
    int arr_size, t2loop;
	
	float *uo,*u1,*u2,*der;
    float c0,c1,c2,dt;
    sf_file in=NULL, out=NULL;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

	
	
	//=====================================================	
	//Get parameters from command line:

	//along which axis you want to derivate?
	if(! sf_getint("axis",&axis)) axis=1; 
    //1=backward,2=centered,3=forward der 
	if(! sf_getint("operator",&operator)) operator=2; 
	
	//=====================================================

	
	
	
	
	
	
	
	//=====================================================
	//Get parameters from input file
    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    
	
	if(axis==3) {
		if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
		if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
		if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");

		if (!sf_histfloat(in,"d3",&dt)) sf_error("No d3= in input");

		arr_size=n1*n2;
		t2loop=sf_leftsize(in,3);
		n_axis=n3;
		
	}else if (axis==2) {
		if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
		if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
		
		if (!sf_histfloat(in,"d2",&dt)) sf_error("No d2= in input");
		
		arr_size=n1;
		t2loop=sf_leftsize(in,2);
		n_axis=n2;
	}else {
		if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
		
		if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
		
		arr_size=1;
		t2loop=sf_leftsize(in,1);
		n_axis=n1;
	}
	
	//=====================================================

	
	

	
	
	// Now define differential operator
	// der[i] = c0*trace[i-1] + c1*trace[i] +c2*trace[i+1]   
	
	if ( operator==1){

		c0=(1/dt)*-1.0; c1=-c0; c2=0.0;

	}else if (operator==2) {
	
		c0=-(1/dt); c1=0.0; c2=-c0;		
	
	}else if (operator==3) {
	
		c0=0.0; c1=-(1/dt); c2=-c1;		
	
	} else {
	    c0=c1=c2=0.0;
		sf_error("operator can only be 1,2 or 3");
	
	}

	
	
	
	
    uo  = sf_floatalloc(arr_size);
    u1  = sf_floatalloc(arr_size);
    u2  = sf_floatalloc(arr_size);
    der = sf_floatalloc(arr_size);



	
	
	
	
	
	

	
	
	for (i3=0; i3<t2loop; i3++) {
		
		for (i1=0; i1<arr_size; i1++){
			uo[i1]=0.0;
			u1[i1]=0.0;
			u2[i1]=0.0; 
			der[i1]=0.0;	
		}

		sf_floatread(u1,arr_size,in);
		sf_floatread(u2,arr_size,in);
		
		for (i1=0; i1<arr_size; i1++){
			der[i1]=c0*uo[i1]+c1*u1[i1]+c2*u2[i1];	
		}

		sf_floatwrite(der,arr_size,out);

		for (i2=1; i2 < n_axis-1; i2++) {

			for (i1=0; i1<arr_size; i1++){
				uo[i1]= u1[i1];
				u1[i1]= u2[i1];
			} 	
			sf_floatread(u2,arr_size,in);

			for (i1=0; i1<arr_size; i1++){ 
				der[i1]=c0*uo[i1] +c1*u1[i1]+c2*u2[i1];	
			}		
			sf_floatwrite(der,arr_size,out);
   
		}

		for (i1=0; i1<arr_size; i1++){
			uo[i1]= u1[i1];
			u1[i1]= u2[i1];		
			u2[i1]=0.0;
			
			der[i1]=c0*uo[i1]+c1*u1[i1]+c2*u2[i1];	
		}
	
		sf_floatwrite(der,arr_size,out);

    }
	
	exit(0);
}
