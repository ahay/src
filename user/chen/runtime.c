/* runtime */

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
#include <time.h>

static	long nblk, szblk;
static time_t start, end;

void runtime_init(int size)
/*< initialize >*/
{
	szblk = size;
	nblk = 0;
	start = time(NULL);
}


float runtime(long db)
/*< data flow >*/
{
	float f;
	double d;
	nblk += db; 
	end = time(NULL);
	d = difftime(end, start);
	d = nblk / d;
	d *= szblk;
	f = d/1024/1024;
	return f;
}

