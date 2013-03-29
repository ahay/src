

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
#include "fcoh2.h"


int main(int argc, char*argv[])
{
	float **in, **p1;
	void *h;	
	h = fcoh2_init(5,3,-1,1,1);

	in = sf_floatalloc2(5,3);
	p1 = sf_floatalloc2(5,3);

	in[0][0] = 0; in[1][0] = 1; in[2][0] = 2;
	in[0][1] = 1; in[1][1] = 2; in[2][1] = 3;
	in[0][2] = 2; in[1][2] = 3; in[2][2] = 2;
	in[0][3] = 3; in[1][3] = 2; in[2][3] = 1;
	in[0][4] = 2; in[1][4] = 1; in[2][4] = 0;

	fcoh2_2d(h, in, p1, 1);


	return 0;
}


