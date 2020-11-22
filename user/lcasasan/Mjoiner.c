/* Join two selected points along the first dimension */
/*
  Copyright (C) 2010 University of Texas at Austin

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


void sigmoid (int w, float* coeff);
void interp_joint (float* data, const float* coeff, int m, int k, int datalength);


int main (int argc, char* argv[])
{		

    float o1,d1;	
    int n1,	n2, n3,nw, nindex;
    int i2,i3;	
    float *column=NULL,*window=NULL;
    int  *index=NULL;
    sf_file in=NULL, index_FILE=NULL, out=NULL;
	
    sf_init (argc,argv);
	
    in=sf_input("in");
    index_FILE=sf_input("index");	
    out=sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
	
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
	
    if (!sf_getint("nw",&nw)) sf_error("Need nw=");     /* length of joining window */

	
    if (!sf_histint(index_FILE,"n1",&nindex)) sf_error("No n1= in index file");
	
    if (!(n2==nindex)) sf_error("Error index file lenght nindex=%d differs from input n2=%d",nindex,n2);
	
    column = sf_floatalloc(n1);
	
    index  = sf_intalloc(n2);
    window = sf_floatalloc(nw);
    sf_intread(index,n2,index_FILE);	
	
    sigmoid(nw,window);

    /*for (int i=0;i<nw;i++)
      sf_warning("window[%d]=%f",i,window[i]); */
	
    /* reading the number of gahters in data*/
    n3 = sf_leftsize(in,2);	

    if (n3==0) n3=1;
	
    /* sf_warning("I'm here nw=%d",nw); */
    for (i3=0;i3<n3;i3++) { /*gahters loop */	
	sf_warning("Gather %d/%d",i3+1,n3);
	for (i2=0;i2<n2;i2++) {
	    sf_floatread(column,n1,in);
		
	    interp_joint (column, window, index[i2], (nw+1)/2+1, n1);

	    sf_floatwrite(column,n1,out);
	    
	}
    } /* END gahters loop */
    sf_close();
    exit (0);
}




void sigmoid (int w, float* coeff){
    int i;
    for (i = 0; i < w; i++)
    {
	coeff[i] = cos( (SF_PI*i) / (2*(w-1)));
	coeff[i] = coeff[i]* coeff[i];
    }	
}

/* m = punto di ingresso nel dato
   k = 2 * w */
void interp_joint (float* data, const float* coeff, int m, int k, int datalength) {

    int i, indStart, indEnd;
    float a1, a2;
	
    if (k > m) indStart = 0;
    else indStart = m-k;
    if (m+k >= datalength) indEnd = datalength-1;
    else indEnd = m+k;

    a1 = data[indStart];
    /* sf_warning("Son qua,m=%d a1[%d]=%f",m,indStart,a1); */
    a2 = data[indEnd];
    /* sf_warning("Son qua,m=%d a2[%d]=%f",m,indEnd,a2); */
	
    for (i = indStart+1; i < indEnd; i++)
	data[i] = a1 * coeff[i-(m-k)] + a2 * coeff[(m+k)-i];

    /*for (int i=0;i<datalength;i++)
      sf_warning("data[%d]=%f",i,data[i]); */
}
