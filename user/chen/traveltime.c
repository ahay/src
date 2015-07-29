#include <rsf.h>

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



float sf_traveltime(int type, float t0, float x, float v)
/*< travel time:  type =	0 line, 1 parabola, x hyperbola >*/
{
	float deltat;
	deltat= x/v;
	switch(type)
	{
	case 0:
		return (t0+deltat);
		break;
	case 1:
		return (t0+deltat*deltat);
		break;
	default:
		return sqrtf(t0*t0+deltat*deltat);
	}
}
