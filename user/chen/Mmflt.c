/* 3D Recursive median filter */

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


#include<rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "mflt.h"

int main(int argc, char*argv[])
{
	sf_file in, out;
	int n1, n2, n3, n4, i1, i2, i3, i4, n123;
	int rect[3];
	float *u1;
	void *h1=NULL, *h2=NULL, *h3=NULL;

	sf_init(argc, argv);

	in  = sf_input("in");
	out = sf_output("out");

	if (!sf_getint("rect1",&rect[0])) rect[0]=1;
	/* filter length on 1st axis */
	if (!sf_getint("rect2",&rect[1])) rect[1]=0;
	/* filter length on 2nd axis */
	if (!sf_getint("rect3",&rect[2])) rect[2]=0;
	/* filter length on 3nd axis */

	n4 = sf_leftsize(in, 3);
	if (!sf_histint(in, "n1", &n1)) sf_error("No n1= in input");
	if (!sf_histint(in, "n2", &n2)) 
	{
		n2=1; n3=1; n4=1;
		rect[1]=0;	rect[2]=0;
		sf_warning("No n2= in input");
	}else if (!sf_histint(in, "n3", &n3))
	{
		n3=1; n4=1;
		rect[2]=0;
		sf_warning("No n3= in input");
	}



	n123 = n1*n2*n3;
	u1 = sf_floatalloc(n123);

	if(rect[0]>0) h1 = mflt_init(n1, rect[0]);
	if(rect[1]>0) h2 = mflt_init(n2, rect[1]);
	if(rect[2]>0) h3 = mflt_init(n3, rect[2]);

#ifdef _OPENMP
#pragma omp parallel for  ordered       \
    schedule(dynamic,8)          \
    private(i1, i2, i3, i4)                  
#endif
	for(i4=0; i4<n4; i4++)
	{
		sf_floatread(u1, n123, in);
		if(rect[0]>0) 
			for(i2=0; i2<n2; i2++) 
			for(i3=0; i3<n3; i3++) 
				mflt(h1, u1+i3*n1*n2+i2*n1, 1);
		if(rect[1]>0)
			for(i1=0; i1<n1; i1++) 
			for(i3=0; i3<n3; i3++) 
				mflt(h2, u1+i3*n1*n2+i1, n1);
		if(rect[2]>0)
			for(i1=0; i1<n1; i1++) 
			for(i2=0; i2<n2; i2++) 
				mflt(h3, u1+i2*n1+i1, n1*n2);
		sf_floatwrite(u1, n123, out);
	}

	if(rect[0]>0) mflt_close(h1);
	if(rect[1]>0) mflt_close(h2);
	if(rect[2]>0) mflt_close(h3);
	
	free(u1);
}


