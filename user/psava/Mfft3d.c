/* FFT in 3D */
/*
  Copyright (C) 2008 Colorado School of Mines
  
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

#include "fft3.h"

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
    bool verb;
    int  ompchunk; 
    bool inv;
    bool cnt;
    int  axis;

    sf_file Fi=NULL; /* input   */
    sf_file Fo=NULL; /* output  */

    sf_complex ***cc=NULL;

    /* cube axes */
    sf_axis a1,a2,a3;
    sf_axis b1=NULL,b2=NULL,b3=NULL;

    float scale1,scale2,scale3;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;  /* OpenMP data chunk size */
    if(! sf_getbool("verb",&verb)) verb=false;         /* verbosity flag */
    if(! sf_getbool("inv",&inv)) inv=false;            /* forward/inverse */
    if(! sf_getbool("cnt",&cnt)) cnt=false;            /* centering */
    if(! sf_getint("axis",&axis)) axis=0;              /* FFT axis or axes */

    Fi = sf_input ("in" );
    Fo = sf_output("out");

    /* input axes */
    a1 = sf_iaxa(Fi,1);
    a2 = sf_iaxa(Fi,2);
    a3 = sf_iaxa(Fi,3);
    if(verb) {
	sf_raxa(a1);
	sf_raxa(a2);
	sf_raxa(a3);
    }
    b3=a3;
    b2=a2;
    b1=a1;

    if(cnt) {
	switch(axis) {
	    case 3:
		b3= sf_maxa(sf_n(a3),-1./(2.*sf_d(a3)),1./(sf_n(a3)*sf_d(a3)));
		break;		
	    case 2:
		b2= sf_maxa(sf_n(a2),-1./(2.*sf_d(a2)),1./(sf_n(a2)*sf_d(a2)));
		break;		
	    case 1:
		b1= sf_maxa(sf_n(a1),-1./(2.*sf_d(a1)),1./(sf_n(a1)*sf_d(a1)));
		break;		
	    case 0:
	    default:
		b1= sf_maxa(sf_n(a1),-1./(2.*sf_d(a1)),1./(sf_n(a1)*sf_d(a1)));
		b2= sf_maxa(sf_n(a2),-1./(2.*sf_d(a2)),1./(sf_n(a2)*sf_d(a2)));
		b3= sf_maxa(sf_n(a3),-1./(2.*sf_d(a3)),1./(sf_n(a3)*sf_d(a3)));
		break;
	}
    } else {
	switch(axis) {
	    case 3:
		b3= sf_maxa(sf_n(a3),                0,1./(sf_n(a3)*sf_d(a3)));
		break;		
	    case 2:
		b2= sf_maxa(sf_n(a2),                0,1./(sf_n(a2)*sf_d(a2)));
		break;		
	    case 1:
		b1= sf_maxa(sf_n(a1),                0,1./(sf_n(a1)*sf_d(a1)));
		break;		
	    case 0:
	    default:
		b1= sf_maxa(sf_n(a1),                0,1./(sf_n(a1)*sf_d(a1)));
		b2= sf_maxa(sf_n(a2),                0,1./(sf_n(a2)*sf_d(a2)));
		b3= sf_maxa(sf_n(a3),                0,1./(sf_n(a3)*sf_d(a3)));
		break;
	}
    }
    sf_oaxa(Fo,b1,1);
    sf_oaxa(Fo,b2,2);
    sf_oaxa(Fo,b3,3);

    /*------------------------------------------------------------*/
    /* allocate arrays */
    cc = sf_complexalloc3(sf_n(a1),sf_n(a2),sf_n(a3)); /*    data */

    /*------------------------------------------------------------*/
    sf_complexread(cc[0][0],sf_n(a1)*sf_n(a2)*sf_n(a3),Fi); /* read */
    /*------------------------------------------------------------*/

    switch(axis) {
	case 3:
	    sf_warning("FFT on axis 3");
	    scale3=fft3a3_init(sf_n(a1),sf_n(a2),sf_n(a3));
	    if(cnt && !inv) cnt3a3(cc);
	    fft3a3(inv,(kiss_fft_cpx***) cc,scale3);
	    if(cnt &&  inv) cnt3a3(cc);
	    break;

	case 2:
	    sf_warning("FFT on axis 2");
	    scale2=fft3a2_init(sf_n(a1),sf_n(a2),sf_n(a3));
	    if(cnt && !inv) cnt3a2(cc);
	    fft3a2(inv,(kiss_fft_cpx***) cc,scale2);
	    if(cnt &&  inv) cnt3a2(cc);
	    break;

	case 1:
	    sf_warning("FFT on axis 1");
	    scale1=fft3a1_init(sf_n(a1),sf_n(a2),sf_n(a3));
	    if(cnt && !inv) cnt3a1(cc);
	    fft3a1(inv,(kiss_fft_cpx***) cc,scale1);
	    if(cnt &&  inv) cnt3a1(cc);
	    break;

	case 0:
	default:
	    sf_warning("FFT on axes 123");
	    fft3_init(sf_n(a1),sf_n(a2),sf_n(a3));
	    if(cnt && !inv) cnt3(cc);
	    fft3(inv,(kiss_fft_cpx***) cc);
	    if(cnt &&  inv) cnt3(cc);
	    break;
    }


    /*------------------------------------------------------------*/
    sf_complexwrite(cc[0][0],sf_n(a1)*sf_n(a2)*sf_n(a3),Fo); /* write */
    /*------------------------------------------------------------*/

    /* close FFT */
    fft3_close();

    exit (0);
}
