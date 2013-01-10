/* 2D divn by stationary LS */

/*
  Copyright (C) 2013 Zhonghuan Chen, UT Austin, Tsinghua University
  
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
#include "runtime.h"
#include "tls2.h"

int main(int argc, char* argv[])
{
    bool tls, angle;
    int i1, i3, n1, n2, n3, n12, rect[3];
    float norm, *u1, *u2;
    sf_file inp, out, den;
	void *h=NULL;

    sf_init(argc,argv);
    inp = sf_input("in");
    den = sf_input("den");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(inp)) sf_error("Need float input");
    
    if (!sf_histint(inp,"n1", &n1)) sf_error("No n1= in input");
    if (!sf_histint(inp,"n2", &n2)) sf_error("No n2= in input");
    if (!sf_histint(inp,"n3", &n3)) n3=1;
    n12 = n1*n2;

    if (!sf_getint("rect1",&rect[0])) rect[0]=0;
    if (!sf_getint("rect2",&rect[1])) rect[1]=0;
    if (!sf_getint("rect3",&rect[2])) rect[2]=0;
    /* smoothing radius */

    if (!sf_getbool("tls",&tls)) tls=false;
    /* total least squares */
    if (!sf_getbool("angle",&angle)) angle=false;
    /* angle or slope */

    u1 = sf_floatalloc(n12);
    u2 = sf_floatalloc(n12);

	runtime_init(n12*sizeof(float));
	h = tls2_init(n1, n2, rect, tls, angle);
    for (i3=0; i3 < n3+rect[2]; i3++)
	{
		if(i3<n3)
		{ 
			sf_floatread(u1, n12, inp);
			sf_floatread(u2, n12, den);
		}else{
			for (i1=0; i1 < n12; i1++) 
			{ u1[i1] = 0.0; u2[i1] = 0.0; }
		}
		tls2(h, u1, u2);

		if(i3>=rect[2])
		{
			sf_floatwrite(u1, n12, out);
			norm = runtime(1);
			sf_warning("%d of %d, %f MB/sec;", i3-rect[2], n3, norm);
		}
	}
	tls2_close(h);

	free(u1);
	free(u2);
    exit(0);
}


