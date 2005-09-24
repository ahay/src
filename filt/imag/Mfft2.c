/* Fast Fourier Transform on axis 1 and/or axis 2 */
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

#include "fft1.h"
#include "fft2.h"

#define X2K(a,b,p) b.n=a.n+p; \
                   b.d=2.0*SF_PI/(b.n*a.d); \
                   b.o=(1==b.n)?0:-SF_PI/a.d;

#define LOOP(a) for(i2=0;i2<a2.n;i2++){ \
                for(i1=0;i1<a1.n;i1++){ {a} }}

int main (int argc, char *argv[])
{
    axa a1,a2,a3;
    axa b1,b2;
    int i1,i2,i3;
    int s1,s2; /* transformation sign */
    int p1,p2;

    char *mode;
    bool  flag;
    float complex **dx, **dk;
   
    sf_file Fi,Fo;

    sf_init(argc,argv);

    if(! sf_getint("s1",&s1)) s1=0;
    if(! sf_getint("s2",&s2)) s2=0;
    if(! sf_getint("p1",&p1)) p1=0;
    if(! sf_getint("p2",&p2)) p2=0;

    Fi = sf_input ( "in"); if (SF_COMPLEX !=sf_gettype(Fi)) sf_error("Need complex data");
    Fo = sf_output("out"); sf_settype(Fo,SF_COMPLEX);

    iaxa(Fi,&a1,1);
    iaxa(Fi,&a2,2);
    iaxa(Fi,&a3,3); oaxa(Fo,&a3,3);

    if(s1!=0) {  
	X2K(a1,b1,p1);
	oaxa(Fo,&b1,1);
    }
    if(s2!=0) {
	X2K(a2,b2,p2);
	oaxa(Fo,&b2,2);
    }

    dx=sf_complexalloc2(a1.n,a2.n);

    mode = sf_charalloc(2);
    if(s1 != 0 && s2 != 0) mode = "b"; /* both */
    if(s1 != 0 && s2 == 0) mode = "1";
    if(s1 == 0 && s2 != 0) mode = "2";

    switch(mode[0]) {
	case 'b':
	    sf_warning("FFT on 1 and 2");
	    if(s1 != s2) sf_error("%s: s1 != s2",__FILE__);
	    flag = (s1>0);

	    fft2_init(b1.n,b2.n);
	    sft2_init(a1.o,a1.d,a2.o,a2.d);

	    dk=sf_complexalloc2(b1.n,b2.n);

	    oaxa(Fo,&b1,1);
	    oaxa(Fo,&b2,2);
	    break;
	case '1':
	    sf_warning("FFT on 1");
	    flag = (s1>0);

	    fft1_init(b1.n,a2.n,1);
	    sft1_init(a1.o,a1.d);

	    dk=sf_complexalloc2(b1.n,a2.n);

	    oaxa(Fo,&b1,1);
	    oaxa(Fo,&a2,2);
	    break;
	case '2':
	    sf_warning("FFT on 2");
	    flag = (s2>0);

	    fft1_init(a1.n,b2.n,2);
	    sft1_init(a2.o,a2.d);

	    dk=sf_complexalloc2(a1.n,b2.n);

	    oaxa(Fo,&a1,1);
	    oaxa(Fo,&b2,2);
	    break;
	default:
	    sf_warning("no FFT");
	    dk=sf_complexalloc2(a1.n,a2.n);
	    flag = false;

	    oaxa(Fo,&a1,1);
	    oaxa(Fo,&a2,2);
	    break;
    }

    for( i3=0; i3<a3.n; i3++) {
	sf_warning("i3=%d of %d",i3+1,a3.n);

	sf_complexread(dx[0],a1.n*a2.n,Fi);
	LOOP( dk[i2][i1] = dx[i2][i1]; );

	switch(mode[0]) {
	    case 'b':
		cnt2(     dk);
		fft2(flag,dk);
		sft2(     dk);
		sf_complexwrite(dk[0],b1.n*b2.n,Fo);
		break;
	    case '1':
		cnt1a1(      dk);
		fft1a1(false,dk);
		sft1a1(      dk);
		sf_complexwrite(dk[0],b1.n*a2.n,Fo);
		break;
	    case '2':
		cnt1a2(      dk);
		fft1a2(false,dk);
		sft1a2(      dk);
		sf_complexwrite(dk[0],a1.n*b2.n,Fo);
		break;
	    default:
		sf_complexwrite(dk[0],a1.n*a2.n,Fo);
		break;
	}
    }

    switch(mode[0]) {
	case 'b':
	    fft2_close();
	    sft2_close();
	    break;
	case '1':
	case '2':
	    fft1_close();
	    sft1_close();
	    break;
    }

    exit (0);
}

