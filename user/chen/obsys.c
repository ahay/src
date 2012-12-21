/* Observation system */

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


typedef struct tagVel
{
	sf_axis z, x, y;
	int nd, nxyz;
	float *p;
}Vel, *HVel;
/* velocity file */
/*^*/

HVel obs_vel(sf_file fp)
/*< velocity file >*/
{
	int nn[SF_MAX_DIM], nd;
	HVel p;

	p = sf_alloc(1,sizeof(Vel));

	nd = sf_filedims(fp, nn);
	p->nxyz = nn[0];
	p->z = sf_iaxa(fp, 1);
	p->nd = 1;
	if(nd > 1 && nn[1]>1)
	{ 
		p->x = sf_iaxa(fp, 2);
		p->nxyz *= nn[1];
		p->nd = 2;
	}
	if(nd > 2 && nn[2]>1) 
	{
		p->y = sf_iaxa(fp, 3);	
		p->nxyz *= nn[2];
		p->nd = 3;
	}

	p->p = sf_floatalloc(p->nxyz);
	sf_floatread(p->p, p->nxyz, fp);
	return p;
}


typedef struct tagCoord
{
	sf_axis a1, a2;
	int *p;
}Coord, *HCoord;
/* coordinate file */
/*^*/


HCoord obs_coord(sf_file fp, int *m, int nm)
/*< read coordinate file >*/
{
	int n, i1, i2, pos, n1, n2;
	HCoord p;
	int **t;

	p = sf_alloc(1,sizeof(Coord));
	p->a1 = sf_iaxa(fp, 1);
	p->a2 = sf_iaxa(fp, 2);
	n1 = sf_n(p->a1);
	n2 = sf_n(p->a2);

	p->p = sf_intalloc(n2);
	t = sf_intalloc2(n1, n2);

	sf_intread(t[0], n1*n2, fp);
	
	n = (n2<nm)? n1: nm;
	for(i2=0; i2<n2; i2++)
	{
		pos=0;
		for(i1=n-1; i1>=0; i1--)
		{
			if(t[i2][i1]<0 || t[i2][i1] >= m[i1]) 
				sf_error("Coord (%d, %d) incorrect", i1, i2);
			pos *= m[i1];
			pos += t[i2][i1];
		}
		p->p[i2] = pos;
	}
	return p;
}



