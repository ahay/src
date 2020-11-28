/* Generate receiver file for Helmholtz solver. */
/*
  Copyright (C) 2012 University of Texas at Austin
  
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

int main(int argc, char* argv[])
{
    sf_file out;
    int n1,n2,i,j;
    int recz, recx0, recdx;
    float d1,d2;
    float **f;

    sf_init(argc, argv);
    out = sf_output("out");

    if (!sf_getint("n1",&n1)) n1=1;
    if (!sf_getint("n2",&n2)) n2=1;
    if (!sf_getfloat("d1",&d1)) d1=0.1;
    if (!sf_getfloat("d2",&d2)) d2=0.1;
    
    if (!sf_getint("recz",&recz)) sf_error("No srcz=.");
    if (!sf_getint("recx0",&recx0)) sf_error("No srcx0=.");
    if (!sf_getint("recdx",&recdx)) sf_error("No srcdx=.");

    f=sf_floatalloc2(n1,n2);

    sf_putint(out,"n1",n1);
    sf_putint(out,"n2",n2);
    sf_putfloat(out,"d1",d1);
    sf_putfloat(out,"d2",d2);

    for (j=0; j<n2; j++) { 
    for (i=0; i<n1; i++) {
        f[j][i]=-1.0;
        if ( i == recz ) {
        if ( j > recx0 && j < n2-recx0 ) {
        if ( j % recdx == 0 ) { 
            f[j][i]=1.0;
        }
        }
        }
    } /* for i */
    }  /* for j */

    sf_floatwrite(f[0],n1*n2,out);

    exit(0);
}
