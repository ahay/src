/* 1D Least-Sqaure Adaptive Matched-Filter for Multiple Suppression 
	x = argmin || d - M x ||^2 

	The Program uses internal (icaf.c) or transient convolution (tcaf.c)
*/
/*
  Copyright (C) 2010 Politecnico di Milano

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


//FIXME: bug with tcaf

#include <rsf.h>
#include <math.h>
#include "match1d.h"

#define MY_MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MY_MAX(X,Y) ((X) > (Y) ? (X) : (Y))


int main (int argc, char *argv[])
{
    char *method;
	bool verb,transient;
	int n1,n2,n3,i3;
	int order,w1,w2;
	float eps;
	float *data=NULL,*model=NULL,*output=NULL;
	sf_file in=NULL, out=NULL, multiple=NULL;
	
	sf_init (argc,argv);
	
	in=sf_input("in");
	out=sf_output("out");


	if (NULL != sf_getstring("multiple")) multiple=sf_input("multiple");
	else sf_error("No multiple file specified");

//FIXME: to be reactivated when tcaf issue solved
    //if (NULL == (method = sf_getstring("method")))
    	method="new";
    /* method to use (old,new) */

	if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
	if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");

	if (!sf_getint("w1",&w1)) w1=9;
	/* data window length along 1st dimentions (time/depth) */
	if (!sf_getint("order",&order)) order=w1-2;
	/* matchied-filter order */
	if (order>w1-2)
		sf_error("order=%d must be less equal than %d (w1-2)",order,w1-2); 
	
	if (!sf_getint("w2",&w2)) w2=3;
	/* data window length along 1st dimentions (time/depth) */

	if (!sf_getfloat("eps",&eps)) eps=0.01;
	/* dumping parameter  */

	if (!sf_getbool("verb",&verb)) verb=false;

//FIXME: to be reactivated when tcaf issue solved
//	if (!sf_getbool("transient",&transient))
		transient=false;
	/* transient convolution [y/n] */


	data= sf_floatalloc(n1*n2);
	model= sf_floatalloc(n1*n2);
	output= sf_floatalloc(n1*n2);
	// adaptative subtraction
	n3 = sf_leftsize(in,2);	
	


	for (i3=0;i3<n3;i3++) { /*gahters loop */
	    
		sf_warning("Gather %d/%d",i3+1,n3);
		sf_floatread(data,n1*n2,in);
		sf_floatread(model,n1*n2,multiple);
		
		if (transient)
			adaptsub_tcaf(data,model, output, n1, n2, order, w1, w2,eps,method,verb);
		else
			adaptsub_icaf(data,model, output, n1, n2, order, w1, w2,eps,method,verb);

		sf_floatwrite(output,n1*n2,out);
		
	}
	free(data);
	free(model);
	free(output);

    exit (0);
}
