/* Max/Min curvatures by azimuth curvature cube */

/*
  Copyright (C) 2013 Zhonghuan Chen
  
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
#include "dim3axis.h"


int main(int argc, char* argv[])
{
    int n[3], n12, i1, i2, i3, dim, axis;
    float **u1, *u2;
	char * mode;

	int m1, m2;

    sf_file in, out;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output ("out"); 	/*  */

	if ((mode=sf_getstring("mode"))==NULL) mode="max";
	/* max/min/mean/gauss/mpo/mne/dip/strike */
	if (!sf_getint("axis",&axis)) axis=2;
	/* azimuth dimension */
	dim = dim3axis(in, axis, n);
	if(dim < axis) sf_error("dim < axis");
	n12 = n[0]*n[1];
	sf_warning("n1=%d n2=%d n3=%d\n", n[0], n[1], n[2]);

	sf_unshiftdim(in,out,axis);

    u1 = sf_floatalloc2(n[0], n[1]);
    u2 = sf_floatalloc(n[0]);

    for(i3=0; i3<n[2]; i3++)
    {
	sf_floatread(u1[0], n12, in);
	for(i1=0; i1<n[0]; i1++)
	{
		if(strcmp(mode, "max")==0)
		{
			m1 = 0;
			for(i2=0; i2<n[1]; i2++)
				if(fabs(u1[i2][i1])>fabs(u1[m1][i1])) m1 = i2;
			u2[i1] = u1[m1][i1];
		}else if(strcmp(mode, "min")==0)
		{
			m1 = 0;
			for(i2=0; i2<n[1]; i2++)
				if(fabs(u1[i2][i1])<fabs(u1[m1][i1])) m1 = i2;
			u2[i1] = u1[m1][i1];
		}else if(strcmp(mode, "mpo")==0)
		{
			m1 = 0;
			for(i2=0; i2<n[1]; i2++)
				if(u1[i2][i1]>u1[m1][i1]) m1 = i2;
			u2[i1] = u1[m1][i1];
		}else if(strcmp(mode, "mne")==0)
		{
			m1 = 0;
			for(i2=0; i2<n[1]; i2++)
				if(u1[i2][i1]<u1[m1][i1]) m1 = i2;
			u2[i1] = u1[m1][i1];
		}else if(strcmp(mode, "mean")==0)
		{
			u2[i1] = 0.0;
			for(i2=0; i2<n[1]; i2++)
				u2[i1] += u1[i2][i1];
			u2[i1] /= n[1];
		}else if(strcmp(mode, "gauss")==0)
		{
			m1 = 0;	m2 = 0;
			for(i2=0; i2<n[1]; i2++)
			{
				if(u1[i2][i1]<u1[m1][i1]) m1 = i2;
				if(u1[i2][i1]>u1[m2][i1]) m2 = i2;
			}
			u2[i1] = u1[m1][i1]*u1[m2][i1];
		}
	}
	sf_floatwrite(u2, n[0], out);
    }

    free(u1[0]);
    free(u1);
    free(u2);

    return 0;
}



