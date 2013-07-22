/* Handling of caustics. */
/*
  Copyright (C) 1993 The Board of Trustees of Stanford University
  
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

#include "norgen.h"

static int crosschk (struct point pt0, struct point pt1, struct point pp0, struct point pp1);

void caustics_up (struct heptagon *cube, int start, int nr)
/*< caustics up >*/
{
    int ii, jj;
    int sp, sp2;
    struct point pt[4], pt0, pt2;

    for(jj=start;jj<nr;jj++) {
	ii = (jj+nr)%nr;

        if(cube[ii].cf==END) continue;
        if(cube[ii].cf==CAUSTIC) continue;
        pt[0] = cube[ii].x0; pt[2] = cube[ii].x1;

        sp = ii+1;
	while(cube[(sp+nr)%nr].cf==CAUSTIC) sp++;
	sp = (sp+nr) % nr;	
	if(cube[sp].cf==END) continue;
        pt[1] = cube[sp].x0; pt[3] = cube[sp].x1;

        if(crosschk(pt[0], pt[2], pt[1], pt[3])) {
            cube[ii].cf = CAUSTIC;
	    cube[sp].cf = CAUSTIC;
	    caustics_up (cube, jj-1, nr);
	    break;
        }

        else if(crosschk(pt[0], pt[1], pt[2], pt[3])) {
	    cube[ii].cf = CAUSTIC;
	    cube[sp].cf = CAUSTIC;
            caustics_up (cube, jj-1, nr);
	    break;
        }

        sp2 = sp + 1;
	while(cube[(sp2+nr)%nr].cf==CAUSTIC) sp2++;
	sp2 = (sp2+nr) % nr;	
	if(cube[sp2].cf==END) continue;
        pt0 = cube[sp2].x0;  pt2 = cube[sp2].x1;

        if(crosschk(pt0, pt[1], pt[0], pt[2])) {
	    cube[ii].cf = CAUSTIC;
            cube[sp].cf = CAUSTIC;
	    cube[sp2].cf = CAUSTIC;		
            caustics_up (cube, jj-1, nr);
	    break;
        }

        else if(crosschk(pt0, pt[1], pt[2], pt[3])) {
            cube[sp].cf = CAUSTIC;
	    cube[sp2].cf = CAUSTIC;		
            caustics_up (cube, jj-1, nr);
	    break;
        }

	else if(crosschk(pt2, pt[3], pt[0], pt[2])) {
            cube[ii].cf = CAUSTIC;
            cube[sp].cf = CAUSTIC;
	    cube[sp2].cf = CAUSTIC;		
            caustics_up (cube, jj-1, nr);
	    break;
        }

	else if(crosschk(pt2, pt[3], pt[0], pt[1])) {
            cube[sp].cf = CAUSTIC;
	    cube[sp2].cf = CAUSTIC;		
            caustics_up (cube, jj-1, nr);
	    break;
        }

        else if(crosschk(pt2, pt0, pt[2], pt[0])) {
	    cube[ii].cf = CAUSTIC;
	    cube[sp].cf = CAUSTIC;
	    cube[sp2].cf = CAUSTIC;		
            caustics_up (cube, jj-1, nr);
	    break;
        }

        else if(crosschk(pt2, pt0, pt[0], pt[1])) {
            cube[sp].cf = CAUSTIC;
	    cube[sp2].cf = CAUSTIC;		
            caustics_up (cube, jj-1, nr);
	    break;
        }
        else if(crosschk(pt2, pt0, pt[2], pt[3])) {
            cube[sp].cf = CAUSTIC;
	    cube[sp2].cf = CAUSTIC;		
            caustics_up (cube, jj-1, nr);
	    break;
        }

    }
    return;
}

/*----------------------------------------------------------------------------*/

void caustics_down (struct heptagon *cube, int start, int nr)
/*< caustics down >*/
{
    int ii, jj;
    int sp, sp2;
    struct point pt[4], pt0, pt2;

    for(jj=start;jj>=0;jj--) {
	ii=(jj+nr)%nr;	

        if(cube[ii].cf==CAUSTIC) continue;
        pt[0] = cube[ii].x0; pt[2] = cube[ii].x1;

        sp = ii-1;
        while(cube[(sp+nr)%nr].cf==CAUSTIC) sp--;
	sp = (sp+nr)%nr;
	if(cube[sp].cf==END) continue;
        pt[1] = cube[sp].x0; pt[3] = cube[sp].x1;

        if(crosschk(pt[0], pt[2], pt[1], pt[3])) {
            cube[ii].cf = (cube[ii].cf==END)? END:CAUSTIC;
            cube[sp].cf = CAUSTIC;
            caustics_down (cube, jj+1, nr);
	    break;
        }
        else if(crosschk(pt[0], pt[1], pt[2], pt[3])) {
            cube[ii].cf = (cube[ii].cf==END)? END:CAUSTIC;
            cube[sp].cf = CAUSTIC;
            caustics_down (cube, jj+1, nr);
	    break;
        }

        sp2 = sp - 1;
        while(cube[(sp2+nr)%nr].cf==CAUSTIC) sp2--;
	sp2 = (sp2+nr)%nr;
	if(cube[sp2].cf==END) continue;
        pt0 = cube[sp2].x0;  pt2 = cube[sp2].x1;

        if(crosschk(pt0, pt[1], pt[0], pt[2])) {
            cube[ii].cf = (cube[ii].cf==END)? END:CAUSTIC;
            cube[sp].cf = CAUSTIC;
	    cube[sp2].cf = CAUSTIC;		
            caustics_down (cube, jj+1, nr);
	    break;
        }

        else if(crosschk(pt0, pt[1], pt[2], pt[3])) {
            cube[sp].cf = CAUSTIC;
	    cube[sp2].cf = CAUSTIC;		
            caustics_down (cube, jj+1, nr);
	    break;
        }

        else if(crosschk(pt2, pt[3], pt[0], pt[2])) {
            cube[ii].cf = (cube[ii].cf==END)? END:CAUSTIC;
            cube[sp].cf = CAUSTIC;
	    cube[sp2].cf = CAUSTIC;		
            caustics_down (cube, jj+1, nr);
	    break;
        }

        else if(crosschk(pt2, pt[3], pt[0], pt[1])) {
            cube[sp].cf = CAUSTIC;
	    cube[sp2].cf = CAUSTIC;		
            caustics_down (cube, jj+1, nr);
	    break;
        }

        if(crosschk(pt2, pt0, pt[2], pt[0])) {
	    cube[ii].cf = (cube[ii].cf==END)? END:CAUSTIC;
	    cube[sp].cf = CAUSTIC;
	    cube[sp2].cf = CAUSTIC;
            caustics_down (cube, jj+1, nr);
	    break;
        }

        else if(crosschk(pt2, pt0, pt[0], pt[1])) {
            cube[sp].cf = CAUSTIC;
	    cube[sp2].cf = CAUSTIC;
            caustics_down (cube, jj+1, nr);
	    break;
        }
        else if(crosschk(pt2, pt0, pt[2], pt[3])) {
            cube[sp].cf = CAUSTIC;
	    cube[sp2].cf = CAUSTIC;
            caustics_down (cube, jj+1, nr);
	    break;
        }

    }
    return;
}

/*----------------------------------------------------------------------------*/

static int crosschk (struct point pt0, struct point pt1, struct point pp0, struct point pp1)
{
    float mt=0.0f, bt=0.0f, mp=0.0f, bp=0.0f;
    float X, Z;
    int inft, infp;

    inft = infp = 0;

    /* Obtain line eqn. that defines both segments */
    if(pt0.x == pt1.x) inft = 1;
    else {
        mt = slope (pt0, pt1);
        bt = ordinate (pt0, mt);
    }

    if(pp0.x == pp1.x) infp = 1;

    else {
        mp = slope (pp0, pp1);
        bp = ordinate (pp0, mp);
    }

    /* Obtain intersection point of the two lines */

    if(inft && infp) return 0;
    if(inft) {
        X = pt0.x;
        Z = mp * X + bp;
    } else if(infp) {
        X = pp0.x;
        Z = mt * X + bt;
    } else {
        if (mt==mp) return 0;
        X = (bt - bp) / (mp - mt);
        Z = mt * X + bt;
    }

    /* Check if it belongs to both segments */
    if(!bel(X, pt0.x, pt1.x)) return 0;
    if(!bel(Z, pt0.z, pt1.z)) return 0;
    if(!bel(X, pp0.x, pp1.x)) return 0;
    if(!bel(Z, pp0.z, pp1.z)) return 0;

    /* At this point we know that both segments intersect */
    return 1;
}
