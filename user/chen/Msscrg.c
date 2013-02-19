/* Extract common reciever gathers from simultaneous data */

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
	sf_file in, delay, out ;
	int n1, nt, nr, ns;
	int ir, is;
	int *pdelay, nn[SF_MAX_DIM], ndim;
	float *pi, **po;
	/* char buf[4]; */
	sf_axis ax;
	int jt;

	sf_init(argc, argv);

	in  = sf_input("in");
	delay  = sf_input("delay");
	out = sf_output("out");

	ndim = sf_filedims(in, nn);
	sf_shiftdim(in, out, 2);

	if(!sf_histint(delay, "n1", &ns)) sf_error("n1 needed in delay");
	ax = sf_iaxa(delay, 1);
	sf_oaxa(out, ax, 2);

	nr = 1;
	for(ir=1; ir<ndim; ir++) nr *= nn[ir];

	pdelay = sf_intalloc(ns);
	sf_intread(pdelay, ns, delay);

	if(!sf_getint("jt", &jt)) jt = 1;
	/* subsampling [nps] in observation */
	for(is=0, n1=0; is<ns; is++)
	{
		pdelay[is] /= jt;
		if(pdelay[is]>n1) n1 = pdelay[is];
	}
	n1 = nn[0] - n1;
	if(!sf_getint("nt", &nt)) nt = n1;
	n1 = nn[0];

	pi = sf_floatalloc(n1);
	po = sf_floatalloc2(nt, nr);

	sf_putint(out, "n1", nt);

	for(ir=0; ir<nr; ir++)
	{
		sf_floatread(pi, n1, in);
		for(is=0; is<ns; is++)
		memcpy(po[is], pi+pdelay[is], nt*sizeof(float));
		sf_floatwrite(po[0], nt*ns, out);
	}

	free(po[0]);
	free(po);
	free(pi);
	free(pdelay);

	return 0;

}


