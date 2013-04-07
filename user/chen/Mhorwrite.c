/* save rsf data into horizon ASCII format, eq to dd form=ascii */

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WA:RRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/


#include <rsf.h>


int main(int argc, char*argv[])
{
	sf_file in;
	int n1, n2, n3;
	int i1, i2, i3;
	float o3, d3, o2, d2, *u, *p;
	char* head, *tail, *format;
	sf_init(argc, argv);

	in  = sf_input("in");
	
	head = sf_getstring("head");
	/* header */
	tail = sf_getstring("tail");
	/* tail */
	if((format = sf_getstring("fromat"))==NULL) format="%g ";
	/* format */

	if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

	if (!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
	if (!sf_histint(in, "n2", &n2)) sf_error("No n2= in input");
	if (!sf_histint(in, "n3", &n3)) n3=1;
	if (!sf_histfloat(in, "o2", &o2)) o2=1;
	if (!sf_histfloat(in, "o3", &o3)) o3=1;
	if (!sf_histfloat(in, "d2", &d2)) d2=1;
	if (!sf_histfloat(in, "d3", &d3)) d3=1;

	u = sf_floatalloc(n1*n2*n3);
	p = u;

	sf_floatread(u, n1*n2*n3, in);

	fprintf(stdout, "%s\n", head);
	for(i3=0; i3<n3; i3++)
	for(i2=0; i2<n2; i2++)
	{
		fprintf(stdout, format, o3+d3*i3);
		fprintf(stdout, format, o2+d2*i2);
		for(i1=0; i1<n1; i1++)
			fprintf(stdout, format, *p++);
		fprintf(stdout, "\n");
	}
	fprintf(stdout, "%s\n", tail);

	free(u);
	return 0;
}



