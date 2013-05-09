/* Direct division for complex data. */
/*
  Copyright (C) 2013 University of Texas at Austin
   
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
    sf_complex *num, *den, *rat;
    sf_file fnum, fden, frat;

    sf_init(argc,argv);
    fnum = sf_input("in");
    fden = sf_input("den");
    frat = sf_output("out");

    if (SF_COMPLEX != sf_gettype(fnum) ||
	SF_COMPLEX != sf_gettype(fden)) sf_error("Need complex input");

    dim = sf_filedims (fnum,n);
    nd = 1;

    for (i=0; i < dim; i++) nd *= n[i];

    num = sf_complexalloc(nd);
    den = sf_complexalloc(nd);
    rat = sf_complexalloc(nd);

    sf_complexread(num,nd,fnum);
    sf_complexread(den,nd,fden);

    for (id=0; id < nd; id++) {
#ifdef SF_HAS_COMPLEX_H
	rat[id] = num[id]/den[id]; 
#else
	rat[id] = sf_cdiv(num[id],den[id]);
#endif
    }

    sf_complexwrite(rat,nd,frat);

    exit(0);
}
