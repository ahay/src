/* 2D convolution with arbitrary filter */
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
    bool centered;

    sf_file Fx=NULL; /* input  */
    sf_file Fy=NULL; /* output */
    sf_file Ff=NULL; /* filter - assumes centered filter */

    float  **x=NULL; /* data in  */
    float  **y=NULL; /* data out */
    float  **f=NULL; /* filter   */
    
    /* cube axes */
    sf_axis a1,a2;
    sf_axis f1,f2;
    
    int i1,i2;
    int j1,j2;
    int k1,k2;
    int m1,m2;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;      /* OpenMP data chunk size */
    if(! sf_getbool("verb",&verb)) verb=false;             /* verbosity flag */
    if(! sf_getbool("centered",&centered)) centered=false; /* centered filter */

    Fx = sf_input ("in" );
    Fy = sf_output("out");
    Ff = sf_input ("flt");
    
    /* input axes */
    a1 = sf_iaxa(Fx,1);
    a2 = sf_iaxa(Fx,2);
    f1 = sf_iaxa(Ff,1);
    f2 = sf_iaxa(Ff,2);

    if(verb) { 
	sf_raxa(a1); 
	sf_raxa(a2); 
	sf_raxa(f1); 
	sf_raxa(f2); 
    }

    /*------------------------------------------------------------*/
    /* allocate arrays */
    x=sf_floatalloc2(sf_n(a1),sf_n(a2)); /* data in  */
    y=sf_floatalloc2(sf_n(a1),sf_n(a2)); /* data out */ 
    f=sf_floatalloc2(sf_n(f1),sf_n(f2)); /* filter */

    /*------------------------------------------------------------*/
    sf_floatread(x[0],sf_n(a1)*sf_n(a2),Fx); /* read data         */
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    sf_floatread(f[0],sf_n(f1)*sf_n(f2),Ff); /* read filter       */
    /*------------------------------------------------------------*/

    /* zero output */
    for(i2=0; i2<sf_n(a2); i2++) {
	for(i1=0; i1<sf_n(a1); i1++) {
	    y[i2][i1]=0;
	}
    }
    
    /*------------------------------------------------------------*/
    if(centered) {
	m1= (sf_n(f1)+1)/2-1; m1=sf_n(f1)>1?m1:0;
	m2= (sf_n(f2)+1)/2-1; m2=sf_n(f2)>1?m2:0;
    }
    if(verb) sf_warning("m1=%d m2=%d",m1,m2);

    /*------------------------------------------------------------*/
    for(    k2=0; k2<sf_n(f2); k2++) {
	for(k1=0; k1<sf_n(f1); k1++) {
	    
	    for(    j2=0; j2<sf_n(a2); j2++) {
		i2=j2+k2-m2;		
		if( (i2>=0 && i2<sf_n(a2))) {
		    
		    for(j1=0; j1<sf_n(a1); j1++) {
			i1=j1+k1-m1;
			if( (i1>=0 && i1<sf_n(a1)) ) {
			    
			    y[j2][j1] += x[i2][i1] * f[k2][k1];
			} /* if i1 */
		    } /* j1 loop*/
		} /* if i2 */
	    } /* j2 loop*/
	    
	} /* k1 loop*/
    } /* k2 loop*/

    /*------------------------------------------------------------*/
    sf_floatwrite(y[0],sf_n(a1)*sf_n(a2),Fy); /* write data       */
    /*------------------------------------------------------------*/

    exit (0);
}
