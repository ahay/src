/* 2-D Kirchhoff pre-stack time migration and demigration with antialiasing. */ 

/*
  Copyright (C) 2012 China University of Petrolum (East China)
  
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

#include "tkirmig.h"

int main(int argc, char* argv[])
{
    int n1, n2, n3, n123;
    float o1, o2, o3, d1, d2, d3, dip;
    bool adj, verb, half;
    float *data, *modl, **vrms, **mask;
    sf_file in, out, vel;

    sf_init(argc,argv);
    in = sf_input("in");
    vel = sf_input("vel");
    out = sf_output("out");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histfloat(in,"d1",&d1)) sf_error("No d1= in input");
    if (!sf_histfloat(in,"o1",&o1)) sf_error("No o1= in input");

    if (!sf_histint(in,"n2",&n2)) sf_error("No n2= in input");
    if (!sf_histfloat(in,"d2",&d2)) sf_error("No d2= in input");
    if (!sf_histfloat(in,"o2",&o2)) sf_error("No o2= in input");

    if (!sf_histint(in,"n3",&n3)) sf_error("No n3= in input");
    if (!sf_histfloat(in,"d3",&d3)) sf_error("No d3= in input");
    if (!sf_histfloat(in,"o3",&o3)) sf_error("No o3= in input");

    if (!sf_getbool("adj",&adj)) adj=true;
    /* yes: migration, no: modeling */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */
    if (!sf_getfloat("dip",&dip)) dip=45.;
    /* the dip range of migration aperture */
    if (!sf_getbool("half",&half)) half=false;
    /* half offset flag */

    n123 = n1*n2*n3;
    vrms = sf_floatalloc2(n1,n2);
    mask = sf_floatalloc2(n2,n3);
    data = sf_floatalloc(n123);
    modl = sf_floatalloc(n123);
    /* read velocity file */
    sf_floatread(vrms[0],n1*n2,vel);

//    sf_triangle1_init(10,n1);

    for (int i3=0; i3 < n3; i3++) {
        for (int i2=0; i2 < n2; i2++) {
            mask[i3][i2]=1.;
        }
    }

    tkirmig_init(n1,d1,o1,n2,d2,o2,n3,d3,o3,dip,vrms,mask,half,verb);

    sf_warning("n123=%d,n3=%d,n2=%d,n1=%d",n123,n3,n2,n1);

    if (adj) {
       sf_floatread(data,n123,in);
       tkirmig_lop(true,false,n123,n123,modl,data);
       sf_floatwrite(modl,n123,out);
       
    } else {
      sf_floatread(modl,n123,in);
      tkirmig_lop(false,false,n123,n123,modl,data);
      sf_floatwrite(data,n123,out);
    }

    sf_close();

    exit(0);
}
