/* 3D FFT with centering and Hermitian scaling  */
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

    sf_fft3d ft1=NULL,ft2=NULL,ft3=NULL; /* FT structures */

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
	    case 23:
	    case 32:
		b2= sf_maxa(sf_n(a2),-1./(2.*sf_d(a2)),1./(sf_n(a2)*sf_d(a2)));
		b3= sf_maxa(sf_n(a3),-1./(2.*sf_d(a3)),1./(sf_n(a3)*sf_d(a3)));
		break;
	    case 13:
	    case 31:
		b1= sf_maxa(sf_n(a1),-1./(2.*sf_d(a1)),1./(sf_n(a1)*sf_d(a1)));
		b3= sf_maxa(sf_n(a3),-1./(2.*sf_d(a3)),1./(sf_n(a3)*sf_d(a3)));
		break;
	    case 12:
	    case 21:
		b1= sf_maxa(sf_n(a1),-1./(2.*sf_d(a1)),1./(sf_n(a1)*sf_d(a1)));
		b2= sf_maxa(sf_n(a2),-1./(2.*sf_d(a2)),1./(sf_n(a2)*sf_d(a2)));
		break;
	    case 3:
		b3= sf_maxa(sf_n(a3),-1./(2.*sf_d(a3)),1./(sf_n(a3)*sf_d(a3)));
		break;		
	    case 2:
		b2= sf_maxa(sf_n(a2),-1./(2.*sf_d(a2)),1./(sf_n(a2)*sf_d(a2)));
		break;		
	    case 1:
		b1= sf_maxa(sf_n(a1),-1./(2.*sf_d(a1)),1./(sf_n(a1)*sf_d(a1)));
		break;		
	    case 123:
	    case 321:
	    case 0:
	    default:
		b1= sf_maxa(sf_n(a1),-1./(2.*sf_d(a1)),1./(sf_n(a1)*sf_d(a1)));
		b2= sf_maxa(sf_n(a2),-1./(2.*sf_d(a2)),1./(sf_n(a2)*sf_d(a2)));
		b3= sf_maxa(sf_n(a3),-1./(2.*sf_d(a3)),1./(sf_n(a3)*sf_d(a3)));
		break;
	}
    } else {
	switch(axis) {
	    case 23:
	    case 32:
		b2= sf_maxa(sf_n(a2),                0,1./(sf_n(a2)*sf_d(a2)));
		b3= sf_maxa(sf_n(a3),                0,1./(sf_n(a3)*sf_d(a3)));
		break;
	    case 13:
	    case 31:
		b1= sf_maxa(sf_n(a1),                0,1./(sf_n(a1)*sf_d(a1)));
		b3= sf_maxa(sf_n(a3),                0,1./(sf_n(a3)*sf_d(a3)));
		break;
	    case 12:
	    case 21:
		b1= sf_maxa(sf_n(a1),                0,1./(sf_n(a1)*sf_d(a1)));
		b2= sf_maxa(sf_n(a2),                0,1./(sf_n(a2)*sf_d(a2)));
		break;
	    case 3:
		b3= sf_maxa(sf_n(a3),                0,1./(sf_n(a3)*sf_d(a3)));
		break;		
	    case 2:
		b2= sf_maxa(sf_n(a2),                0,1./(sf_n(a2)*sf_d(a2)));
		break;		
	    case 1:
		b1= sf_maxa(sf_n(a1),                0,1./(sf_n(a1)*sf_d(a1)));
		break;		
	    case 123:
	    case 321:
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
    cc = sf_complexalloc3(sf_n(a1),sf_n(a2),sf_n(a3));

    /*------------------------------------------------------------*/
    /* read */
    sf_complexread(cc[0][0],sf_n(a1)*sf_n(a2)*sf_n(a3),Fi);
    /*------------------------------------------------------------*/

    switch(axis) {
	case 23:
	case 32:
	    sf_warning("FFT on axes 2 and 3");

	    ft2=sf_fft3a2_init(sf_n(a1),sf_n(a2),sf_n(a3));
	    if(cnt && !inv)           sf_cnt3a2(cc,ft2);
	    sf_fft3a2(    inv,(kiss_fft_cpx***) cc,ft2);
	    if(cnt &&  inv)           sf_cnt3a2(cc,ft2);

	    ft3=sf_fft3a3_init(sf_n(a1),sf_n(a2),sf_n(a3));
	    if(cnt && !inv)           sf_cnt3a3(cc,ft3);
	    sf_fft3a3(    inv,(kiss_fft_cpx***) cc,ft3);
	    if(cnt &&  inv)           sf_cnt3a3(cc,ft3);
	    
	    break;
	case 13:
	case 31:
	    sf_warning("FFT on axes 1 and 3");
	    
	    ft1=sf_fft3a1_init(sf_n(a1),sf_n(a2),sf_n(a3));
	    if(cnt && !inv)           sf_cnt3a1(cc,ft1);
	    sf_fft3a1(    inv,(kiss_fft_cpx***) cc,ft1);
	    if(cnt &&  inv)           sf_cnt3a1(cc,ft1);

	    ft3=sf_fft3a3_init(sf_n(a1),sf_n(a2),sf_n(a3));
	    if(cnt && !inv)           sf_cnt3a3(cc,ft3);
	    sf_fft3a3(    inv,(kiss_fft_cpx***) cc,ft3);
	    if(cnt &&  inv)           sf_cnt3a3(cc,ft3);
	    
	    break;
	case 12:
	case 21:
	    sf_warning("FFT on axes 1 and 2");
	    
	    ft1=sf_fft3a1_init(sf_n(a1),sf_n(a2),sf_n(a3));
	    if(cnt && !inv)           sf_cnt3a1(cc,ft1);
	    sf_fft3a1(    inv,(kiss_fft_cpx***) cc,ft1);
	    if(cnt &&  inv)           sf_cnt3a1(cc,ft1);

	    ft2=sf_fft3a2_init(sf_n(a1),sf_n(a2),sf_n(a3));
	    if(cnt && !inv)           sf_cnt3a2(cc,ft2);
	    sf_fft3a2(    inv,(kiss_fft_cpx***) cc,ft2);
	    if(cnt &&  inv)           sf_cnt3a2(cc,ft2);
	    
	    break;
	case 3:
	    sf_warning("FFT on axis 3");
	    
	    ft3=sf_fft3a3_init(sf_n(a1),sf_n(a2),sf_n(a3));
	    if(cnt && !inv)           sf_cnt3a3(cc,ft3);
	    sf_fft3a3(    inv,(kiss_fft_cpx***) cc,ft3);
	    if(cnt &&  inv)           sf_cnt3a3(cc,ft3);

	    break;
	case 2:
	    sf_warning("FFT on axis 2");
	    
	    ft2=sf_fft3a2_init(sf_n(a1),sf_n(a2),sf_n(a3));
	    if(cnt && !inv)           sf_cnt3a2(cc,ft2);
	    sf_fft3a2(    inv,(kiss_fft_cpx***) cc,ft2);
	    if(cnt &&  inv)           sf_cnt3a2(cc,ft2);
	    
	    break;
	case 1:
	    sf_warning("FFT on axis 1");

	    ft1=sf_fft3a1_init(sf_n(a1),sf_n(a2),sf_n(a3));
	    if(cnt && !inv)           sf_cnt3a1(cc,ft1);
	    sf_fft3a1(    inv,(kiss_fft_cpx***) cc,ft1);
	    if(cnt &&  inv)           sf_cnt3a1(cc,ft1);

	    break;
	case 123:
	case 321:
	case 0:
	default:
	    sf_warning("FFT on axes 1,2 and 3");

	    ft1=sf_fft3a1_init(sf_n(a1),sf_n(a2),sf_n(a3));
	    if(cnt && !inv)           sf_cnt3a1(cc,ft1);
	    sf_fft3a1(    inv,(kiss_fft_cpx***) cc,ft1);
	    if(cnt &&  inv)           sf_cnt3a1(cc,ft1);

	    ft2=sf_fft3a2_init(sf_n(a1),sf_n(a2),sf_n(a3));
	    if(cnt && !inv)           sf_cnt3a2(cc,ft2);
	    sf_fft3a2(    inv,(kiss_fft_cpx***) cc,ft2);
	    if(cnt &&  inv)           sf_cnt3a2(cc,ft2);

	    ft3=sf_fft3a3_init(sf_n(a1),sf_n(a2),sf_n(a3));
	    if(cnt && !inv)           sf_cnt3a3(cc,ft3);
	    sf_fft3a3(    inv,(kiss_fft_cpx***) cc,ft3);
	    if(cnt &&  inv)           sf_cnt3a3(cc,ft3);

	    break;
    }

    /*------------------------------------------------------------*/
    /* write */
    sf_complexwrite(cc[0][0],sf_n(a1)*sf_n(a2)*sf_n(a3),Fo); 
    /*------------------------------------------------------------*/

    /* close FFT */
    switch(axis) {
	case 23:
	case 32:
	    sf_fft3a2_close(ft2);
	    sf_fft3a3_close(ft3);
	    break;
	case 13:
	case 31:
	    sf_fft3a1_close(ft1);
	    sf_fft3a3_close(ft3);
	    break;
	case 12:
	case 21:
	    sf_fft3a1_close(ft1);
	    sf_fft3a2_close(ft2);
	    break;
	case 3:
	    sf_fft3a3_close(ft3);
	    break;		
	case 2:
	    sf_fft3a2_close(ft2);
	    break;		
	case 1:
	    sf_fft3a1_close(ft1);
	    break;	
	case 123:
	case 321:	
	case 0:
	default:
	    sf_fft3a1_close(ft1);
	    sf_fft3a2_close(ft2);
	    sf_fft3a3_close(ft3);
	    break;
    }
    /*------------------------------------------------------------*/


    exit (0);
}
