/* Integer header attributes. */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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
#include <stdio.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int i1, i2, n1, n2, *max, *min, *inp, imax, imin;
    float *mean;
    char pad[] = "                    ", out[21];
    sf_file head;
    
    sf_init (argc,argv);
    head = sf_input("in");

    if (SF_INT != sf_gettype(head)) sf_error("Need int input");
    if (!sf_histint(head,"n1",&n1)) sf_error("No n1= in input");
    n2 = sf_leftsize(head,1);

    printf("******************************************* \n");
    printf("%d headers, %d traces\n",n1,n2);
 
    inp = sf_intalloc(n1);
    max = sf_intalloc(n1);
    min = sf_intalloc(n1);
    mean = sf_floatalloc(n1);

    sf_intread(inp,n1,head);
    imax = imin = 0;
    for (i1=0; i1 < n1; i1++) {
	min[i1] = max[i1] = inp[i1];
	mean[i1] = inp[i1];
    }


    for (i2=1; i2 < n2; i2++) {
	sf_intread(inp,n1,head);
	for (i1=0; i1 < n1; i1++) {
	    if (min[i1] > inp[i1]) {
		min[i1] = inp[i1];
		imin = i2;
	    }
	    if (max[i1] < inp[i1]) {
		max[i1] = inp[i1];
		imax = i2;
	    }
	    mean[i1] += inp[i1];
	}
    }

    for (i1=0; i1 < n1; i1++) {
	snprintf(out,21,"key[%d]=\"%s\"",i1,
		 i1 < SF_NKEYS? sf_segykeyword(i1): "?");
	printf("%s%s",out,pad+strlen(out));
	snprintf(out,21,"min[%d]=%d",imin,min[i1]);
	printf("%s%s",out,pad+strlen(out));
	snprintf(out,21,"max[%d]=%d",imax,max[i1]);
	printf("%s%s",out,pad+strlen(out));
	printf("mean=%g\n",mean[i1]/n2);
    }
   
    printf("******************************************* \n");

    exit(0);
}








