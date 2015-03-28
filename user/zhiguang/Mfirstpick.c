/* Pick first arrival by recognizing the first nonzero value at each trace */
/*
 Copyright (C) 2014 University of Texas at Austin
 
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
	int is,i,j,n1,n2,n3;
	float threshold, d1, d2, o1, o2, o3, d3;
	float **data, *time;
	sf_file in,out;

	sf_init(argc,argv);

	in=sf_input("in");
	out=sf_output("out");

	if(!sf_getfloat("threshold", &threshold)) sf_error("Need thresholding.");
	if(!sf_histint(in, "n1", &n1)) sf_error("No n1= in input.");
	if(!sf_histint(in, "n2", &n2)) sf_error("No n2= in input.");
	if(!sf_histint(in, "n3", &n3)) sf_error("No n3= in input.");
	if(!sf_histfloat(in, "d1", &d1)) sf_error("No d1= in input.");
	if(!sf_histfloat(in, "d2", &d2)) sf_error("No d2= in input.");
	if(!sf_histfloat(in, "d3", &d3)) sf_error("No d3= in input.");
	if(!sf_histfloat(in, "o1", &o1)) sf_error("No o1= in input.");
	if(!sf_histfloat(in, "o2", &o2)) sf_error("No o2= in input.");
	if(!sf_histfloat(in, "o3", &o3)) sf_error("No o3= in input.");

	data=sf_floatalloc2(n1, n2);
	time=sf_floatalloc(n2);

	sf_putint(out, "n1", n2);
	sf_putint(out, "n2", n3);
	sf_putint(out, "n3", 1);
	sf_putfloat(out, "d1", d2);
	sf_putfloat(out, "d2", d3);
	sf_putfloat(out, "o1", o2);
	sf_putfloat(out, "o2", o3);
	sf_putstring(out, "label1", "Offset");
	sf_putstring(out, "unit1", "m");
	sf_putstring(out, "label2", "Shots");
	sf_putstring(out, "unit2", "m");

	sf_warning("threshold is %g", threshold);

	for(is=0; is<n3; is++){
		sf_floatread(data[0], n1*n2, in);
		for(i=0; i<n2; i++){
			for(j=0; j<n1; j++){
				if(data[i][j]>threshold){
					time[i]=o1+j*d1;
					break;
				}
			}
		}
		sf_floatwrite(time, n2, out);
	}

	exit(0);
}
