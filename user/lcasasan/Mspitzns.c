/* Missing data interpolation in 2-D using F-X Prediction Error Filters
   based on Seismic trace interpolation in the F-X domain
   S.Spitz Geophysics 56, 785(1991). 

   Uses 2D Patching. 
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

#include <rsf.h>
#include <rsfgee.h>

#include "patching1.h"

int main(int argc, char* argv[])
{
    bool verb, norm;


    int n[2], w[2], k[2], a[2], l[2], w_out[2], n_out[2];   
    int  n1, n2, n12, i3, n3;
	int  n12_out, w12_out ;
    float d2; 			   /* data parameters */    
    int ntraces, order, L; /* input parameters */
 	float *wall, *data, *windwt;
    
    
    sf_file in, out;// mask;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    
    if (SF_FLOAT !=sf_gettype(in)) sf_error("Need float type");    

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
 
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input");


    n3 = sf_leftsize(in,2); /*number of gathers in the line*/
	if (n3==0) n3=1;    

	n[0]=n1;
	n[1]=n2;	
		
	//TENT PARAMETERS
	a[0] = 1; a[1] = 1; l[0] = 1;	l[1] = 1;
	
	if (!sf_getint("w1",&w[0])) w[0] = n1; /*lenght of patch along the first dimension */
    if (!sf_getint("w2",&w[1])) w[1] = n2; /*lenght of patch along the second dimension */

    if (!sf_getint("k1",&k[0])) k[0] = 1;  /*number of patches along the first dimension */
    if (!sf_getint("k2",&k[1])) k[1] = 1;  /*number of patches along the second dimension */

    if (!sf_getint("order",&order)) order=3;
    /* linear PEF order*/
    if (!sf_getint("ntraces",&ntraces)) ntraces=1;
    /* number of traces to be interpolated */
    
    if (!sf_getbool("verb",&verb)) verb = false;
    /* verbosity flag */

	if (!sf_getbool("norm",&norm)) norm = true;
    /* output normalization flag */


    n12 = n[0]*n[1];
    
	L = ntraces+1;
	
	n_out[0] = n[0];
	n_out[1] = L * (n2-1) + 1;

	n12_out= n_out[0]*n_out[1];
		
	w_out[0] = w[0];
	w_out[1] = L * (w[1]-1) + 1;
	
	w12_out= w_out[0]*w_out[1];

	
	/*sf_warning("n12 =%d",n12);
	sf_warning("n12_out =%d",n12_out);
	sf_warning("w12 =%d",w12);    
	sf_warning("w12_out =%d",w12_out);  
	sf_warning("patches =%d x %d",k[0],k[1]); */
	
	wall = sf_floatalloc(n12);
    data = sf_floatalloc(n12_out);
    windwt = sf_floatalloc(w12_out);


	
	sf_putint  (out,"n2",n_out[1]);
	sf_putfloat  (out,"d2",d2/(ntraces+1));	
  	tent (2, w_out, l, a, windwt);
	for (i3=0;i3<n3;i3++) { 	/*GATHERS LOOP [3rd dimension]*/
	    sf_warning("Gather %d/%d",i3+1,n3);
		
		sf_floatread(wall,n12,in);
 
		patching1(wall, data, 2, k, n, w, n_out, w_out, windwt, order, ntraces ,verb, norm);
	
    	sf_floatwrite(data,n12_out,out);
	} /* END GATHERS LOOP [3rd dimension]*/


	free(wall);
	free(data);
	free(windwt);

    exit(0);
}



