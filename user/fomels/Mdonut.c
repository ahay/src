/* Donut filter */
/*
  Copyright (C) 2022 University of Texas at Austin
   
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

int main(int argc, char *argv[])
{
    bool adj;
    int i1, n1, ia, na;
    float *trace, *a, *holed;
    sf_file hole, data, filt;

    sf_init(argc,argv);

    if (!sf_getbool("adj",&adj)) adj=false;

    data = sf_input("data");

    if (adj) {
	hole = sf_input("in");
	if (!sf_histint(hole,"n1",&n1)) sf_error("No n1= in input");

	if (!sf_getint("na",&na)) sf_error("Need na=");	
	
	filt = sf_output("out");
	sf_putint(filt,"n1",na);
    } else {
	filt = sf_input("in");
	if (!sf_histint(filt,"n1",&na)) sf_error("No n1= in input");

	if (!sf_histint(data,"n1",&n1)) sf_error("No n1= in data");
	
	hole = sf_output("out");
	sf_putint(hole,"n1",n1);
    }
      
    trace = sf_floatalloc(n1);
    holed = sf_floatalloc(n1);
    a = sf_floatalloc(na);

    sf_floatread(trace,n1,data);

    if (adj) {
	sf_floatread(holed,n1,hole);
	for (ia=0; ia < na; ia++) {
	    a[ia] = 0.0f;
	}
    } else {
	sf_floatread(a,na,filt);	    
	for (i1=0; i1 < n1; i1++) {
	    holed[i1] = 0.0f;
	}
    }
    
    for (i1=na; i1 < n1-na; i1++) {
	for (ia=0; ia < na; ia++) {
	    if (adj) {
		a[ia] += holed[i1]*(trace[i1-ia]+trace[i1+ia]); 
	    } else {
		holed[i1] += a[ia]*(trace[i1-ia]+trace[i1+ia]); 
	    }
	}
    }

    if (adj) {
	sf_floatwrite(a,na,filt);
    } else {
	sf_floatwrite(holed,n1,hole);
    }

    exit(0);
}
