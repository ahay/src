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
    bool tls;
    int i1, i3, n1, n2, n3, n12, rect[3];
    float norm, *u1, *u2, *u3;
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

    if (!sf_getint("rect1",&rect[0])) rect[0]=1;
    if (!sf_getint("rect2",&rect[1])) rect[1]=1;
    if (!sf_getint("rect3",&rect[2])) rect[2]=1;
    /* smoothing radius */

    if (!sf_getbool("tls",&tls)) tls=true;
    /* total least squares */

    u1 = sf_floatalloc(n12);
    u2 = sf_floatalloc(n12);
    u3 = sf_floatalloc(n12*2);

	runtime_init(n12*sizeof(float));
	h = tls2_init(n1, n2, rect, false);
    for (i3=0; i3 < n3+rect[2]; i3++)
	{ 
		sf_floatread(u1, n12, inp);
		sf_floatread(u2, n12, den);
			/* smooth division */
		for (i1=0; i1 < n12; i1++) 
		{
			u3[i1*2] = u1[i1];
			u3[i1*2+1] = u1[i1];
		}
		tls2(h, u3);

		if(i3>=rect[2])
		{
			sf_floatwrite(u3, n12, out);
			norm = runtime(1);
			sf_warning("%d of %d, %f MB/sec;", i3-rect[2], n3, norm);
		}
	}
	tls2_close(h);

	free(u1);
	free(u2);
	free(u3);
    exit(0);
}


