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
    int nk_rfft, nx_rfft;
    bool adj;
    float *ww, *ff;
    sf_complex *pp, *qq;
    sf_file src, out, w, wf;

    sf_init(argc,argv);

    src = sf_input("in");
    out = sf_output("out");
    w = sf_input("w");
    wf = sf_input("wf");

    if (!sf_histint(src,"n1",&nz)) sf_error("No n1= in input");
    if (!sf_histint(src,"n2",&nx)) sf_error("No n2= in input");
    /* dim nk from frequency weight - derived from real fft2 */
    if (!sf_histint(wf,"n1",&nk_rfft)) sf_error("No n1= in wf");
    if (!sf_histint(wf,"n2",&nx_rfft)) sf_error("No n2= in wf");


    nzx = nz*nx;
    n_left = sf_leftsize(src,2);
    
    nk = cfft2_init(1,nz,nx,&nz2,&nx2);
    sf_warning("nz is %d nx is %d", nz,nx);
    sf_warning("(Padded Dom) nz2 is %d nx2 is %d", nz2,nx2);
    sf_warning("(FFT Dom) nk is %d nx2 is %d nk(comb) is %d", nk/nx2, nx2,nk);

    pp = sf_complexalloc(nzx);
    qq = sf_complexalloc(nzx);

    ww = sf_floatalloc(nz*nx);
    sf_floatread(ww,nzx,w);
    sf_fileclose(w);

    ff = sf_floatalloc(nk);
    sf_floatread(ff,nk_rfft,wf);
    sf_fileclose(wf);

    ctf2dprec_init(nz, nx, nk, nz2, nx2, ww, ff);


    sf_complexread(pp,nzx,src);



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