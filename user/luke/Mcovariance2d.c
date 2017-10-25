/* determine covariance from 2d data*/
/* input data must be mean zero*/
/*
  Copyright (C) 2017 University of Texas at Austin

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


// make a vector from data file
void makevec(float *vec, float** data, int i, int nx){
	int j;
	for( j=0; j < nx; j++){
		vec[j] = data[j][i];				
	}	
	return;	
}

void dotproduct(float *dot, float *vec1, float *vec2, int nx){
	int i;
	*dot = 0.0; //initialize
	for (i=0; i<nx ; i++){
		*dot += vec1[i]*vec2[i];
	}
	return;
}

int main(int argc, char *argv[])
{	//declare variables
	int i,j;
	int nt, nx;

	float dt, dx, t0, x0;
	float dot;

	float **data=NULL;
	float **cov=NULL;
	float *vec1=NULL;
	float *vec2=NULL;


	sf_file in=NULL, out=NULL;
	sf_init (argc,argv);
	in = sf_input("in");
	out = sf_output("out");

	//read data
	if (!sf_histint(in,"n1",&nt)) sf_error("No n1= in input");
	if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");


	if (!sf_histfloat(in,"d1",&dt)) sf_error("No d1= in input");
	if (!sf_histfloat(in,"d2",&dx)) sf_error("No d2= in input");


	if (!sf_histfloat(in,"o1",&t0)) sf_error("No o1= in input");
	if (!sf_histfloat(in,"o2",&x0)) sf_error("No o1= in input");

	// allocate arrays
	data = sf_floatalloc2(nt,nx);
	cov = sf_floatalloc2(nt,nt);
	vec1 = sf_floatalloc(nx);
	vec2 = sf_floatalloc(nx);

	// setup output
	sf_putint(out,"n1",nt);
	sf_putfloat(out,"d1",dt);
	sf_putfloat(out,"o1",t0);
	sf_putint(out,"n2",nt);
	sf_putfloat(out,"d2",dt);
	sf_putfloat(out,"o2",t0);



	// read data
	sf_floatread (data[0],nt*nx,in);
/*
	// zero out covariance 
	for (i=0; i<nt; i++){
		for (j=0; j<nt; j++){
			cov[i][j] = 0.;
		}
	}
*/

	// dot products for covariance		
	for (i=0; i<nt; i++){
		makevec(vec1,data,i,nx);
		for (j=0; j<i+1; j++){ // calculate up to diagonal
			makevec(vec2,data,j,nx);
			dotproduct(&dot,vec1,vec2,nx);
			cov[i][j] = dot;
			cov[j][i] = dot; // make use of symmetry
		}
	}

	// write data
	sf_floatwrite(cov[0],nt*nt,out);
	exit(0);
}
