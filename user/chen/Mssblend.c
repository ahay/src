/* blend reciever gathers (T-S-R) to generate simultaneous data */

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
	int ir, is, jt, inx;
	int *pdelay, nn[SF_MAX_DIM], ndim;
	float **pi, *po;
#ifndef NO_BLAS
	float alpha;
#endif

	sf_init(argc, argv);

	in  = sf_input("in");
	delay  = sf_input("delay");
	out = sf_output("out");

	ndim = sf_filedims(in, nn);
	if(ndim<3) sf_error("input dimensions < 3");

	if(!sf_histint(delay, "n1", &ns)) sf_error("n1 needed in delay");
	n1 = nn[0];
	if(nn[1] != ns) sf_error("delay and input has different shots");

	nr = 1;
	for(ir=2; ir<ndim; ir++) nr *= nn[ir];

	if(!sf_getint("jt", &jt)) jt=1;
	/* subsampling nps */

	pdelay = sf_intalloc(ns);
	sf_intread(pdelay, ns, delay);
	for(is=0, nt=0; is<ns; is++)
	{
		pdelay[is] /= jt;
		if(pdelay[is]>nt) nt = pdelay[is];
	}
	nt += n1;

	pi = sf_floatalloc2(n1, ns);
	po = sf_floatalloc(nt);
	sf_putint(out, "n1", nt);
	sf_unshiftdim(in,out,2);

#ifndef NO_BLAS
	inx =1;
	alpha = 1.0; 
#endif
	
	for(ir=0; ir<nr; ir++)
	{
		sf_floatread(pi[0], ns*n1, in);
		memset(po, 0, sizeof(float)*nt);
		for(is=0; is<ns; is++)
		{
#ifndef NO_BLAS
		saxpy_(&n1, &alpha, pi[is], &inx, po+pdelay[is], &inx);
#else
		for(inx=0; inx<n1; inx++)
			po[pdelay[is]+inx] += pi[is][inx];
#endif
		}
		sf_floatwrite(po, nt, out);
	}

	free(pi[0]);
	free(pi);
	free(po);
	free(pdelay);

	return 0;

}


