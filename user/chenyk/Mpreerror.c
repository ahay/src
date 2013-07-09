/* Prediction error */
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

int main(int argc, char* argv[])
{
    int n1,n2,n3,n11,n22,n33,i1,i2,i3,ifrela;
    float *trace1,*trace2, sum, dif,error;
    sf_file in, pre, out;

    sf_init(argc,argv);
    in = sf_input("in");
    pre = sf_input("predict");
    out = sf_output("out");

    if(!sf_histint(in,"n1",&n1)) sf_error("No n1 in input");
    if(!sf_histint(in,"n2",&n2)) sf_error("No n2 in input");
    n3=sf_leftsize(in,2); 

    if(!sf_histint(pre,"n1",&n11)) sf_error("No n1 in predict");
    if(!sf_histint(pre,"n2",&n22)) sf_error("No n2 in predict");
    n33=sf_leftsize(in,2);    

    if(n1!=n11 || n2!=n22 || n3!=n33) sf_error("Dimension mismatch between two input and predict");

    if(!sf_getint("type",&ifrela)) ifrela=1;
    /* if compute relative error, 1: yes, 0: no, default is yes. */
    
    trace1=sf_floatalloc(n1*n2);
    trace2=sf_floatalloc(n1*n2);
    sf_putint(out,"n1",1);
    sf_putint(out,"n2",n3);
    sf_putint(out,"n3",1);		

    for (i3=0; i3 < n3; i3++)
    {
	sf_floatread(trace1,n1*n2,in);
	sf_floatread(trace2,n1*n2,pre);
	sum=0;dif=0;
	for(i2=0;i2<n2;i2++)
		for(i1=0;i1<n1;i1++)
		{
			sum+=trace1[i2*n1+i1]*trace1[i2*n1+i1]; dif+=(trace1[i2*n1+i1]-trace2[i2*n1+i1])*(trace1[i2*n1+i1]-trace2[i2*n1+i1]);
		}
	if(ifrela==1)error=dif;
	else error=dif/sum;
        sf_floatwrite(&error,1,out);
    }



    exit(0);
}
