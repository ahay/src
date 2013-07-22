/* dip estimation by line interpolating pwd */

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
#include "ldip.h"


int main(int argc, char*argv[])
{
    sf_file in, out, p0;
    int m, n, n1, n2, n3, rect[2], niter, liter;
    int i3;
    bool verb;
    float **wav, **dip, eta, **idip;
    char *interp;

    sf_init(argc, argv);

    in  = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float type");

    if (!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
    if (!sf_histint(in, "n2", &n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(in, 2);

    if(!sf_getint("m", &m)) m=1;
    /* b[-m, ... ,n] */
    if(!sf_getint("n", &n)) n=1;
    /* b[-m, ... ,n] */
    if ((interp=sf_getstring("interp"))==NULL) interp="maxflat";
    /* interpolation method: maxflat lagrange bspline */

    if (!sf_getint("rect1",&rect[0])) rect[0]=0;
    /* dip smoothness on 1st axis */
    if (!sf_getint("rect2",&rect[1])) rect[1]=0;
    /* dip smoothness on 2nd axis */
	
    if (!sf_getint("niter",&niter)) niter=5;
    /* number of iterations */
    if (!sf_getint("liter",&liter)) liter=20;
    /* number of linear iterations */
    if (!sf_getfloat("eta", &eta)) eta = 1.0;
    /* steps for iteration */

    if(sf_getstring("idip")!=NULL)
    {
	idip = sf_floatalloc2(n1, n2);
	p0 = sf_input("idip");
	sf_histint(p0, "n1", &i3);
	if(i3!=n1) sf_error("p0.n1=%d and in.n1=%d",i3,n1);
	sf_histint(p0, "n2", &i3);
	if(i3!=n2) sf_error("p0.n2=%d and in.n2=%d",i3,n2);
	sf_floatread(idip[0], n1*n2, p0);
    } else	{
	idip = NULL;
	p0=NULL;
    }
	
    /* starting dip */
    if (!sf_getbool("verb", &verb)) verb = false;
    /* verbosity flag */

    wav = sf_floatalloc2(n1, n2);
    dip = sf_floatalloc2(n1, n2);

	
    /* initialize dip estimation */
    ldip_init(interp, m, n, n1, n2, rect, liter, verb);


    for(i3=0; i3<n3; i3++)
    {
	sf_floatread(wav[0], n1*n2, in);
	if(p0!=NULL) memcpy(dip[0], idip[0], n1*n2*sizeof(float));
	else memset(dip[0], 0, n1*n2*sizeof(float));
	ldip(wav, dip, niter, eta);
	sf_floatwrite(dip[0], n1*n2, out);
    }

    if(p0!=NULL)
    {
	free(idip[0]);
	free(idip);
    }
    ldip_close();
    free(dip[0]);
    free(wav[0]);
    free(dip);
    free(wav);
    return 0;
}



