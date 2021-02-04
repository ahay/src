/* 2D data resampling with ability to extrapolate. */

/*
  Copyright (C) 2012 Zhonghuan Chen, UT Austin, Tsinghua University
  
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

float resample(float **u, int n1, int n2, float x1, float x2)
{
	int y1, y2;
	float a1, b1, a2, b2;

	y1 = x1;
	y2 = x2;
	a1 = x1-y1;
	b1 = 1.0 - a1;
	a2 = x2-y2;
	b2 = 1.0 - a2;
	return (b2*(b1*u[y2][y1]+a1*u[y2][y1+1]) +
			a2*(b1*u[y2+1][y1]+a1*u[y2+1][y1+1]));
}

int main(int argc, char*argv[])
{
	sf_file in, out;
	int n[SF_MAX_DIM], dim, nr;
	float d[2], o[2];

	int n1, n2;
	float d1, d2, o1, o2;

	int j, ir, i1, i2;
	float **ibuf, **obuf, x1, x2;

	sf_init(argc, argv);

	in = sf_input ("in");
	out = sf_output ("out");

	dim = sf_filedims(in,n);
	for(nr=1, j=2; j < dim; j++) 	nr *= n[j];

	if (!sf_histint(in,"n1", n)) sf_error("n1");
	if (!sf_histfloat(in,"o1", o)) sf_error("o1");
	if (!sf_histfloat(in,"d1", d)) sf_error("d1");
	if (!sf_histint(in,"n2", n+1)) sf_error("n2");
	if (!sf_histfloat(in,"o2", o+1)) sf_error("o2");
	if (!sf_histfloat(in,"d2", d+1)) sf_error("d2");
	
	if (!sf_getfloat("o1",&o1)) o1=o[0];
	/* first sample sample on 1st axis */
	if (!sf_getfloat("d1",&d1)) d1=d[0];
	/* sample interval on 1st axis */
	if (!sf_getfloat("o2",&o2)) o2=o[1];
	/* first sample on 2nd axis */
	if (!sf_getfloat("d2",&d2)) d2=d[1];
	/* sample interval on 2nd axis */

	
    if (!sf_getint("n1",&n1)) n1 = (d[0]*(n[0]-1)+o[0]-o1)/d1;
	/* number of samples on first axis */
    if (!sf_getint("n2",&n2)) n2 = (d[1]*(n[1]-1)+o[1]-o2)/d2;
	/* number of samples on second axis */

	sf_putint(out,"n1", n1);
	sf_putfloat(out,"o1", o1);
	sf_putfloat(out,"d1", d1);
	sf_putint(out,"n2", n2);
	sf_putfloat(out,"o2", o2);
	sf_putfloat(out,"d2", d2);

	ibuf = sf_floatalloc2(n[0], n[1]);
	obuf = sf_floatalloc2(n1, n2);

	for(ir=0;ir<nr;ir++)
	{
		sf_floatread(ibuf[0], n[0]*n[1], in);
		for(i2=0; i2<n2; i2++)
		{
			x2 = (i2*d2 + o2 - o[1])/d[1];
			if (x2 < 0.0) x2 = 0.0;
			if (x2 > (float)n[1]-1.01) x2 = (float)n[1]-1.01;
			for(i1=0; i1<n1; i1++)
			{
				x1 = (i1*d1 + o1 - o[0])/d[0];
				if (x1 < 0.0) x1 = 0.0;
				if (x1 > (float)n[0]-1.01) x1=(float)n[0]-1.01;
				obuf[i2][i1] = resample(ibuf, n1, n2, x1, x2);
			}
		}
		sf_floatwrite(obuf[0], n1*n2, out);
	}	

	return 0;
}

