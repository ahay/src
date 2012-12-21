/* generate simultaneous sources grid from delay file */

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
	sf_file in, out ;
	int i1, i2, n1, n2, *v;
	float o1, d1, **u;
	char *l1, *u1;
	sf_axis ax;

	sf_init(argc, argv);

	in  = sf_input("in");  /* delay file (int) */
	out = sf_output("out");

	if(!sf_histint(in, "n1", &n2)) sf_error("n1 needed");
	sf_shiftdim(in, out, 1);


	if(!sf_getint("n1", &n1)) n1=1000; /* samples */
	if(!sf_getfloat("o1", &o1)) o1=0.0; /* sampling interval */
	if(!sf_getfloat("d1", &d1)) d1=0.004; /* original  */
	if((l1=sf_getstring("l1")) == NULL) l1="Time"; /* label "Time" */
	if((u1=sf_getstring("u1")) == NULL) u1="s"; /* unit "s" */

	ax = sf_maxa(n1, o1, d1);
	sf_setlabel(ax, l1);
	sf_setunit(ax, u1);
	sf_oaxa(out, ax, 1);
	sf_putint(out, "n2", n2);
	sf_settype(out, SF_FLOAT);

	v = sf_intalloc(n2);
	u = sf_floatalloc2(n1, n2);
	sf_intread(v, n2, in);

	for(i2=0; i2<n2; i2++) 
	for(i1=0; i1<n1; i1++) 
	if(i1==v[i2])	u[i2][i1] = 1;
	else u[i2][i1] = 0;

	sf_floatwrite(u[0], n1*n2, out);

	free(v);
	free(u[0]);
	free(u);

	return 0;

}


