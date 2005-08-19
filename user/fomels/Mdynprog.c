/* Automatic picking by dynamical programming */
/*
  Copyright (C) 2004 University of Texas at Austin
  
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

#include "dynprog.h"

int main(int argc, char* argv[])
{
    int n1, n2, n3, i3, *pick, gate;
    float **scan, an;
    sf_file scn, pik;

    sf_init(argc,argv);
    scn = sf_input("in");
    pik = sf_output("out");

    if (SF_FLOAT != sf_gettype(scn)) sf_error("Need float input");
    sf_settype(pik,SF_INT);

    if (!sf_histint(scn,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(scn,"n2",&n2)) sf_error("No n2= in input");
    n3 = sf_leftsize(scn,2);

    sf_putint(pik,"n1",n2);
    sf_putint(pik,"n2",1);

    if (!sf_getfloat("an",&an)) an=1.; 
    /* axes anisotropy */
    if (!sf_getint("gate",&gate)) gate=n1; 
    /* picking gate */

    scan = sf_floatalloc2(n1,n2);
    pick = sf_intalloc(n2);

    dynprog_init(n2,n1,an);

    for (i3=0; i3 < n3; i3++) {
	sf_floatread(scan[0],n1*n2,scn);
	dynprog(gate,scan,pick);
	sf_intwrite(pick,n2,pik);
    }

    exit(0);
}
