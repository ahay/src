/* Upwind update sequence */
/*
  Copyright (C) 2010 University of Texas at Austin
  
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

static float *t0;
static int *visit, n1, n2;

static int fermat(const void *a, const void *b)
/* comparison for traveltime sorting from small to large */
{
    float ta, tb;

    ta = t0[*(int *)a];
    tb = t0[*(int *)b];

    if (ta >  tb) return 1;
    if (ta == tb) return 0;
    return -1;
}


void update_init(int m1, int m2 /* dimensions */,
		 float *t       /* [m1,m2] traveltime */)
/*< initialize >*/
{
    int i, n12;

    n1 = m1;
    n2 = m2;
    n12 = n1*n2;
    t0 = t;

    /* sort from small to large traveltime */
    visit = sf_intalloc(n12);
    for (i = 0; i < n12; i++) {
	visit[i] = i;
    }
    qsort(visit, n12, sizeof(int), fermat);
}

void update_close(void)
/*< free allocated storage >*/
{
    free(visit);
}

unsigned char get_update(int i, bool *up1, bool *up2, int *j)
/*< next update step >*/
{
    float t1;
    int i1, i2, a1, a2, b1, b2, c1, c2;
    unsigned char update;

    *j = visit[i];
    t1 = t0[*j];

    i1 = *j%n1;
    i2 = *j/n1;
    
    update = 0;
    
    a1 = *j-1;
    b1 = *j+1;
    *up1 = (bool) (i1 && (i1 == n1-1 || 1 != fermat(&a1,&b1)));
    c1 = *up1? a1:b1;
    if (t1 > t0[c1]) update |= 1;
    
    a2 = *j-n1;
    b2 = *j+n1;
    *up2 = (bool) (i2 && (i2 == n2-1 || 1 != fermat(&a2,&b2)));
    c2 = *up2? a2:b2;
    if (t1 > t0[c2]) update |= 2;

    return update;
}
