/* Create a data mask using multiple muting curve from MRKE */
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


void sigmoid (int w, float* coeff);
void interp_joint (float* data, const float* coeff, int m, int k, int datalength);

int main (int argc, char* argv[])
{
	bool verb,start,shift,smooth;
	float o1,o2,d1,d2;	
	int n1,n2,n3;
	int i1,i2,i3;
	int index,nw,nws=0,sign=1;	
	float *trace=NULL,*window=NULL;
	
	sf_file in=NULL, out=NULL, mask=NULL;
	sf_init (argc,argv);
	
	in=sf_input("in");


	mask=sf_input("mask");	
	if (SF_INT != sf_gettype(mask)) sf_error("Need integer mask file");

	out=sf_output("out");
	sf_settype (out, SF_FLOAT);
	
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");
	
	if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&o2)) sf_error("No o2= in input");

	

    if (!sf_getbool("start",&start)) start=false; /*mask from starting sample to index value in mask */

    if (!sf_getbool("smooth",&smooth)) smooth=false; /*smoothed mask [raised cosine] */
	
    if (!sf_getint("nw",&nw)) nw=0; /*smoothing window length must be odd*/
	if (nw<0)
		sign=-1; 
	nw=nw*sign;   

	if (!nw%2)
	nw=nw+1;
    
    if (!sf_getbool("shift",&shift)) shift=false; /*shift */
    
    if (shift)
	nws=(nw-1)/2;	
    
    if (!sf_getbool("verb",&verb)) verb=false;
    
    if (smooth) {
	window = sf_floatalloc(abs(nw));
	sigmoid(nw,window);
    }
 
	/* reading the number of gahters in data*/
    n3 = sf_leftsize(in,2);	


	for (i3=0;i3<n3;i3++) { /*gahters loop */
	    sf_warning("Gather %d/%d",i3+1,n3);
		for (i2=0;i2<n2;i2++) {
			trace = sf_floatalloc(n1);
			
			sf_intread(&index,1,mask);
			
			index= index - ( sign * nws);
			//index= index - nws;

			for (i1=0;i1<n1;i1++) {
				if (start) {
					if (i1<  index)
							trace[i1]=1.0;
					else   
							trace[i1]=0.0;
				}
				else {/*start=false*/ 
					if (i1>=  index)
							trace[i1]=1.0;
					else   
							trace[i1]=0.0;
				}
			}
	    	//sf_warning("Son qua, max=%f index=%d",max[i2],index[i2]);

			if (smooth) {
				if (start)
					interp_joint (trace, window, ( index) + 1, (nw+1)/2, n1);
				else
					interp_joint (trace, window, ( index) - 1, (nw+1)/2, n1);
			}

			sf_floatwrite(trace,n1,out);
			free(trace);
		}




	} /* END gahters loop */
	sf_close();
}

void sigmoid (int w, float* coeff){
	int i;
	for (i = 0; i < w; i++)
        {
                coeff[i] = cos( (SF_PI*i) / (2*(w-1)));
                coeff[i] = coeff[i]* coeff[i];
        }	
}

void interp_joint (float* data, const float* coeff, int m, int k, int datalength) {

        int i, indStart, indEnd;
        float a1, a2;
	
	if (k > m) indStart = 0;
        else indStart = m-k;
        if (m+k >= datalength) indEnd = datalength-1;
        else indEnd = m+k;

        a1 = data[indStart];
		//sf_warning("Son qua,m=%d a1[%d]=%f",m,indStart,a1);
        a2 = data[indEnd];
		//sf_warning("Son qua,m=%d a2[%d]=%f",m,indEnd,a2);
	
	for (i = indStart+1; i < indEnd; i++)
                data[i] = a1 * coeff[i-(m-k)] + a2 * coeff[(m+k)-i];

	//for (int i=0;i<datalength;i++)
	//	sf_warning("data[%d]=%f",i,data[i]);
}
