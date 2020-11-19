/* Edge preserving (smoothing) filter by window selection */

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
#ifdef _OPENMP
#include <omp.h>
#endif

#include "epsmedian.h"
#include "epsmean.h"
#include "epspoly.h"


int main(int argc, char* argv[])
{
    int i2, n1, n2;
	int order, nfw;
    sf_file in, out;
	float *u1;
	char *filter;
	void *h=NULL;

    sf_init(argc, argv);

    in  = sf_input("in");	
    out = sf_output("out");	

    if (!sf_getint("nfw",&nfw)) nfw=5;
    /* window size */
    if ((filter=sf_getstring("filter"))==NULL) filter="mean";
    /* filter: mean,median,poly,fir */
    if (!sf_getint("order",&order)) order=nfw-1;
    /* filter order (<= nfw, only for polynomial and fir filters) */

	if(!sf_histint(in, "n1", &n1)) sf_error("n1 needed in input");
	n2 = sf_leftsize(in, 1);

	u1 = sf_floatalloc(n1);

	if(strcmp(filter, "mean")==0) h = epsmean_init(n1, nfw);
	else if(strcmp(filter, "median")==0) h = epsmedian_init(n1, nfw);
	else if(strcmp(filter, "poly")==0) h = epspoly_init(n1, nfw, order);
	else sf_error("filter %s not support", filter);

#ifdef _OPENMP
#pragma omp parallel for  ordered       \
    schedule(dynamic,5)          \
    private(i2)                  
#endif
    for(i2=0; i2<n2; i2++)
    {
		sf_floatread(u1, n1, in);
		if(strcmp(filter, "mean")==0) epsmean(h, u1, 1);
		else if(strcmp(filter, "median")==0) epsmedian(h, u1, 1);	
		else if(strcmp(filter, "poly")==0) epspoly(h, u1, 1);	
		sf_floatwrite(u1, n1, out);
    }


    free(u1);
    return 0;
}

