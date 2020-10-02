/* TF Weights Preconditioner for Complex input as linear oper. 

December 2013 program of the month:
http://ahay.org/blog/2013/12/01/program-of-the-month-sfcausint/
*/
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
#include <math.h>
#include <rsf.h>
#include "ctf2dprec.h"
#include "cfft2w.h"

int main(int argc, char* argv[])
{
    int nz, nx, nz2, nx2, nk, nzx, n_left, i;
    bool adj;
    sf_complex **pp, **qq;
    sf_file in=NULL, out=NULL, w=NULL, wf=NULL;

    sf_init(argc,argv);

    in = sf_input("in");
    out = sf_output("out");
    w = sf_input("w");
    wf = sf_input("wf");

    if (!sf_histint(in,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&nx)) sf_error("No n2= in input");

    nzx = nz*nx;
    n_left = sf_leftsize(in,2);
    
    nk = cfft2_init(1,nz,nx,&nz2,&nx2);


    ctf2dprec_init(nz, nx, nk, nz2, nx2, w, wf);


    pp = sf_complexalloc2(nz,nx);
    qq = sf_complexalloc2(nz,nx);

    sf_complexread(pp[0],nzx,in);

    if (!sf_getbool("adj",&adj)) adj=false;

    for (i=0; i < n_left; i++) {
      if (adj) {
          ctf2dprec_lop(true,false,nzx,nzx,qq,pp);
      } else {
          ctf2dprec_lop(false,false,nzx,nzx,pp,qq);
      }
    }

    sf_complexwrite(qq,nzx,out);


    exit(0);
}
