/* Compute a histogram of data values. */
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

int main(int argc, char* argv[]) 
{
    int i, n, i1, n1, nbuf, *hist;
    float o1, d1, *fbuf;
    char key[5];
    sf_file in=NULL, out=NULL;

    sf_init(argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");
    n = sf_filesize(in);

    nbuf = BUFSIZ/sizeof(float);
    fbuf = sf_floatalloc(nbuf);

    if (!sf_getint("n1",&n1)) sf_error("Need n1=");
    /* number of histogram samples */
    if (!sf_getfloat("o1",&o1)) sf_error("Need o1=");
    /* histogram origin */
    if (!sf_getfloat("d1",&d1)) sf_error("Need d1=");
    /* histogram sampling */

    sf_settype(out,SF_INT);
    sf_putint(out,"n1",n1);
    sf_putfloat(out,"o1",o1);
    sf_putfloat(out,"d1",d1);
    for (i=1; i < SF_MAX_DIM; i++) {
	sprintf(key,"n%d",i+1);

	if (!sf_histint(in,key,&i1)) break;
	if (i1 > 1) sf_putint(out,key,1);
    }

    hist = sf_intalloc(n1);
    for (i1=0; i1 < n1; i1++) {
	hist[i1]=0;
    }

    for (; n > 0; n -= nbuf) {
	if (nbuf > n) nbuf = n;

	sf_floatread(fbuf,nbuf,in);
	
	for (i=0; i < nbuf; i++) {
	    i1 = (int) floorf((fbuf[i]-o1)/d1);
	    if (i1 >= 0 && i1 < n1) hist[i1]++;
	}
    }

    sf_intwrite(hist,n1,out);
    if (in != NULL) sf_fileclose(in);
    exit(0);
}
