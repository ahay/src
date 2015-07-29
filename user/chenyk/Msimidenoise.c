/* Random noise attenuation using local similarity */
/*
  Copyright (C) 2013 University of Texas at Austin
   
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
#include "fxdeconutil.h"


int main(int argc, char *argv[])
{
	int i1,i2,i3,n1,n11,n2,n22,n3,n33;
	float *tracein,*ref,thr;
	sf_file in=NULL, simi=NULL, out=NULL;
    	sf_init(argc,argv);
    	in = sf_input("in");
    	out = sf_output("out");	
	simi=sf_input("similarity");

	if(!sf_histint(in,"n1",&n1)) sf_error("No n1 in input");
	if(!sf_histint(in,"n2",&n2)) sf_error("No n2 in input");
	n3=sf_leftsize(in,2);	
	if(!sf_histint(simi,"n1",&n11)) sf_error("No n1 in input");
	if(!sf_histint(simi,"n2",&n22)) sf_error("No n2 in input");
	n33=sf_leftsize(simi,2);	
	if(n1!=n11 || n2!=n22 || n3!=n33) sf_error("Dimension dismatch");

	tracein=sf_floatalloc(n1);
	ref=sf_floatalloc(n1);

    	if (!sf_getfloat("thr",&thr)) sf_error("Need thr=");
    	/* thresholding level */ 

	for(i3=0;i3<n3;i3++)
	{
	   for(i2=0;i2<n2;i2++)
		{

		sf_floatread(tracein,n1,in);
		sf_floatread(ref,n1,simi);
		for(i1=0;i1<n1;i1++)
			{
				if(ref[i1] < thr)
				tracein[i1] = 0;
				else
				tracein[i1]=tracein[i1];
			}
		sf_floatwrite(tracein,n1,out);
		}
	}

	exit(0);
}







