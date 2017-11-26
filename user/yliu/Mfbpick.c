/* First break picking from instantaneous traveltime attribute. */
/*
  Copyright (C) 2017 Jilin University
  
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
#include <stdio.h>
#include <math.h>
 
int main (int argc, char* argv[]) 
{
    int n1, n2, n3;
    int i1, i2, i3;
    float d1, d2, o1, o2, travelt;
    char *type;
    float *trace,*travelp=NULL;
    sf_file in, out;


    sf_init (argc, argv);
    in=sf_input("in");
    out=sf_output("out");
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"d2",&d2)) d2=1.;
    if (!sf_histfloat(in,"d1",&o1)) sf_error("No o1= in input");
    if (!sf_histfloat(in,"d2",&o2)) o2=1.;

    n3 = sf_leftsize(in,2);

    if (NULL == (type=sf_getstring("type"))) type="traveltime";
    /* [traveltime,position] type, the default is traveltime  */
    
    trace=sf_floatalloc(n1);

    switch (type[0]) {
	case 't':
	    sf_unshiftdim(in, out, 1);
	    break;
	case 'p':
	    travelp=sf_floatalloc(n1);
	    break;
	default:
	    sf_error("Unknown calculation type=%c",type);
	    break;
    }
  
    for(i3=0; i3 < n3; i3++) {
        for(i2=0; i2 < n2; i2++) {
	    sf_floatread(trace,n1,in);

	    switch (type[0]) {
		case 't':
		    for(i1=0; i1 < n1; i1++) {
			if(trace[i1]==1) {  
			    travelt=i1*d1+o1;
			    break;
			}
		    }
		    sf_floatwrite(&travelt,1,out);
		    break;
		case 'p':
		    for(i1=0; i1 < n1; i1++) {
			travelp[i1]=0.;
		    }
		    for(i1=0; i1 < n1; i1++) {
			if(trace[i1]==1) {  
			    travelp[i1]=1.;
			    break;
			}
		    }
		    sf_floatwrite(travelp,n1,out);
		    break;
		default:
		    sf_error("Unknown calculation type=%c",type);
		    break;
	    }
	}
    }
    exit(0);
}

/* 	$Id$	 */
