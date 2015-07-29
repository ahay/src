/* Compare 2 data set */

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
#include "comp.h"

int main(int argc,char*argv[])
{
	sf_file in,out,ref;
	int dim, dimr, n[SF_MAX_DIM], nr[SF_MAX_DIM];
	int n1, n2, id, i2,  mode;
	float *u,*ur,*p;

	sf_init(argc,argv);

	in  = sf_input("in");
	out = sf_output("out");

	ref = sf_input("ref");

	dimr = sf_filedims(ref,nr);
	dim  = sf_filedims(in,n);
	if(dim < dimr) sf_error("dim < dim_ref");

	for(id=0, n1=1; id<dimr; id++)
	{
		if(n[id] != nr[id]) 
			sf_error("Dimension %d not match!", id);
		n1 *= nr[id];
	}
	sf_unshiftdim(in,out,1);
	for(id=1; id<dimr; id++)	sf_unshiftdim(out,out,1);

	for(id=dimr, n2=1; id<dim; id++)	n2 *= n[id];

	if (!sf_getint("mode",&mode)) mode=0;
	/* compare method: 
	0 - normalized xcorrelation; 
	1 - mean square error */

	u  = sf_floatalloc(n1);
	ur = sf_floatalloc(n1);
	p  = sf_floatalloc(n2);

	sf_floatread(ur,n1,ref);

	for(i2=0; i2<n2; i2++)
	{
		sf_floatread(u,n1,in);
		switch(mode)
		{
		case 1:
			p[i2]=comp_mae(u, ur, n1);
			break;
		case 2:
			p[i2]=comp_mse(u, ur, n1);
			break;
		case 3:
			p[i2]=comp_cor(u, ur, n1);
			break;
		default:
			p[i2]=comp_ncc(u, ur, n1);
		}
	}
	
	sf_floatwrite(p, n2, out);

	exit(0);
}

