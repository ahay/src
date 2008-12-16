/* Aliasing test. 
*/
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include <math.h>

#include <rsf.h>

int main(int argc, char* argv[])
{
    int nt,nx,it,ix,ix0;
    float *wave, *data, cycles, slow;
    sf_file out;
    
    sf_init(argc,argv);
    out = sf_output("out");
    sf_setformat(out,"native_float");

    if (!sf_getint("n1",&nt)) nt=600;
    if (!sf_getint("n2",&nx)) nx=24;
    /* dimensions */

    if (!sf_getfloat("cycles",&cycles)) cycles=10.;
    /* wave frequency */

    if (!sf_getint("ix0",&ix0)) ix0=0; 
    /* central trace */
    /* try ix0=2 */

    if (!sf_getfloat("slow",&slow)) slow=0.1;
    /* slowness */

    sf_putint(out,"n1",nt);
    sf_putint(out,"n2",nx);
    sf_putint(out,"d1",1.);
    sf_putint(out,"d2",1.);
    sf_putfloat(out,"o1",0.);
    sf_putfloat(out,"o2",0.);

    data = sf_floatalloc(nt);
    wave = sf_floatalloc(nt);

    for (it=0; it < nt; it++) {
	data[it]=0.;
	wave[it] = sinf(2.*SF_PI*it*cycles/nt) * 
	    expf(- 3.*(it+1.)/nt);
    }

    for (ix=0; ix < nx; ix++) {
	it = nt*hypotf(1.,slow*(ix-ix0))/3.;
	if (it > nt) it=nt;
	if (it > 0)  sf_floatwrite(data,it,out);
	if (nt > it) sf_floatwrite(wave,nt-it,out);
    }
 
    exit(0);
}
