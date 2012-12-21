/* Generate signal series */

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
#include "wavelet.h"


int main(int argc, char*argv[])
{

	sf_file  out;
	int i1, n1;
	float o1, d1, *u1, par[4];
	char *waveform;
	sf_init(argc, argv);

	out = sf_output("out");

	if ((waveform=sf_getstring("waveform"))==NULL) waveform="ricker";
	/* waveform: ricker,sinc,harmonic,randn,rand */
	if (!sf_getfloats("para", par, 4)) {par[0] = 25.0; par[1] = 0.0;}
/* parameters of waveform\n
  ricker(freq)
  sinc(freq)
  harmonic(freq,phase)
  randn(seed)
  rand(seed)
*/
	if (!sf_getint("n", &n1)) n1 = 100;
	/* length */
	if (!sf_getfloat("o", &o1)) o1 = 0.0;
	/* original */
	if (!sf_getfloat("d", &d1)) d1 = 0.004;
	/* interval */

	sf_putint(out, "n1", n1);
	sf_putfloat(out, "o1", o1);
	sf_putfloat(out, "d1", d1);

	u1 = sf_floatalloc(n1);

	if(strcmp(waveform, "ricker")==0)
	{
		for(i1=0; i1<n1; i1++)	u1[i1] = o1+d1*i1;
		sf_wvlt_rck(n1, u1, par);
	}else if(strcmp(waveform, "harmonic")==0)
	{
		for(i1=0; i1<n1; i1++)	u1[i1] = o1+d1*i1;
		sf_wvlt_harmonic(n1, u1, par);
	}else if(strcmp(waveform, "sinc")==0)
	{
		for(i1=0; i1<n1; i1++)	u1[i1] = o1+d1*i1;
		sf_wvlt_sinc(n1, u1, par);
	}else if(strcmp(waveform, "randn")==0)
	{
		init_genrand(par[0]);
		sf_randn(n1, u1);
	}else if(strcmp(waveform, "rand")==0)
	{
		init_genrand(par[0]);
		sf_random(n1, u1);
	}else sf_error("waveform %s not implemented", waveform);

	sf_floatwrite(u1, n1, out);

	free(u1);
	return 0;

}


