/* Inverse velocity spectrum with interpolation by modeling from inversion result (C version) */
/*
  Copyright (C) 2020 University of Texas at Austin
   
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
    int i, id, dim, n[SF_MAX_DIM], nd;
    kiss_fft_cpx *din;
    float *dout;
    sf_file fin, fout;

    sf_init(argc,argv);
    fin = sf_input("in");
    fout = sf_output("out");

    if (SF_COMPLEX != sf_gettype(fin)) sf_error("Need complex input");

    dim = sf_filedims (fin,n);
    nd = 1;

    for (i=0; i < dim; i++) nd *= n[i];

    din = (kiss_fft_cpx*) sf_complexalloc(nd);
    dout = sf_floatalloc(nd);

    sf_floatread((float*)din,nd*2,fin);

    for (id=0; id < nd; id++) {
	dout[id] = hypotf(din[id].r,din[id].i);
    }

    sf_settype(fout, SF_FLOAT);
    sf_setform(fout, SF_NATIVE);

    sf_floatwrite(dout,nd,fout);

    exit(0);
}
