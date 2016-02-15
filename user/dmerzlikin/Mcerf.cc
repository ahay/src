// Complex error function.

//   Copyright (C) 2010 University of Texas at Austin
//  
//   This program is free software; you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation; either version 2 of the License, or
//   (at your option) any later version.
//  
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//  
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#include <rsf.hh>

#include "Faddeeva3.hh"

int main(int argc, char** argv)
{   
    sf_init(argc,argv); // Initialize RSF

    iRSF inp;
    oRSF out;

    if (SF_COMPLEX != inp.type()) sf_error("Need complex input");

    int nbuf = out.bufsiz()/sizeof(sf_complex);
    std::valarray<sf_complex> cbuf(nbuf);
    std::complex<double> c;

    for (int nsiz = inp.size(); nsiz > 0; nsiz -= nbuf) {
	if (nbuf > nsiz) {
	    nbuf = nsiz;
	    cbuf.resize(nbuf);
	}

	inp >> cbuf;

	for (int i=0; i < nbuf; i++) {
	    c = std::complex<double>(crealf(cbuf[i]),cimagf(cbuf[i]));
	    c = Faddeeva::erf(c);
	    cbuf[i] = sf_cmplx(std::real(c),std::imag(c));
	}

	out << cbuf;
    }

    exit(0);
}
