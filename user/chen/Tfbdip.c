

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
#include "fbdip.h"


int main(int argc, char*argv[])
{
	int rect[2]={0,0};
	float ****in, **p1;

	in = sf_floatalloc4(1, 1, 3, 3);
	p1 = sf_floatalloc2(1,1);

/*  0.1   0.2   0.3
//  1.1   1.2   1.3
//  2.1   2.2   2.3 */
	in[0][0][0][0] = 1.2;  in[1][0][0][0] = 0.05; in[2][0][0][0] = 0.0;
	in[0][1][0][0] = 0.5;  in[1][1][0][0] = 0.0;  in[2][1][0][0] = 0.0;
	in[0][2][0][0] = 0.0;  in[1][2][0][0] = 0.0;  in[2][2][0][0] = 0.0;

	fbdip_init(1.0, 1, 1,
    	rect, 1, 0, false);

	fbdip(3, 3, in, p1, 1, 0.5);
/* large angle model
//  1	2	1.5
//  1	2	1.5
//  1	2	1.5 */
	
	in[0][0][0][0] = 1.75;  in[1][0][0][0] = 0.0;  in[2][0][0][0] = 0.0;
	in[0][1][0][0] = 0.125; in[1][1][0][0] = 0.0;  in[2][1][0][0] = 0.0;
	in[0][2][0][0] =-4.7/36;in[1][2][0][0] = 0.0;  in[2][2][0][0] = 0.0;

	fbdip(3, 3, in, p1, 1, 0.5);
	exit(0);
}


