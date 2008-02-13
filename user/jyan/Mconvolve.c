/* 3D convolution with arbitrary filter */
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

    sf_file Fx=NULL; /* input  */
    sf_file Fy=NULL; /* output */
    sf_file Ff=NULL; /* filter - assumes centered filter */

    float ***x=NULL; /* data in  */
    float ***y=NULL; /* data out */
    float ***f=NULL; /* filter   */
    
    /* cube axes */
    sf_axis a1,a2,a3;
    sf_axis f1,f2,f3;
    
    int i1,i2,i3;
    int j1,j2,j3;
    int k1,k2,k3;
    int m1,m2,m3;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getint("ompchunk",&ompchunk)) ompchunk=1;      /* OpenMP data chunk size */
    if(! sf_getbool("verb",&verb)) verb=false;             /* verbosity flag */

    Fx = sf_input ("in" );
    Fy = sf_output("out");
    Ff = sf_input ("flt");
    
    /* input axes */
    a1 = sf_iaxa(Fx,1);
    a2 = sf_iaxa(Fx,2);
    a3 = sf_iaxa(Fx,3);

    f1 = sf_iaxa(Ff,1);
    f2 = sf_iaxa(Ff,2);
    f3 = sf_iaxa(Ff,3);

    if(verb) { 
	sf_raxa(a1); 
	sf_raxa(a2); 
	sf_raxa(a3); 

	sf_raxa(f1); 
	sf_raxa(f2); 
	sf_raxa(f3); 
    }

    /*------------------------------------------------------------*/
    /* allocate arrays */
    x=sf_floatalloc3(sf_n(a1),sf_n(a2),sf_n(a3)); /* data in  */
    y=sf_floatalloc3(sf_n(a1),sf_n(a2),sf_n(a3)); /* data out */ 
    f=sf_floatalloc3(sf_n(f1),sf_n(f2),sf_n(f3)); /* filter */

    /*------------------------------------------------------------*/
    /* read data  */
    sf_floatread(x[0][0],sf_n(a1)*sf_n(a2)*sf_n(a3),Fx);
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* read filter*/
    sf_floatread(f[0][0],sf_n(f1)*sf_n(f2)*sf_n(f3),Ff);
    /*------------------------------------------------------------*/

    /* zero output */
    for        (i3=0; i3<sf_n(a3); i3++) {
	for    (i2=0; i2<sf_n(a2); i2++) {
	    for(i1=0; i1<sf_n(a1); i1++) {
		y[i3][i2][i1]=0;
	    }
	}
    }
    
    /*------------------------------------------------------------*/
	/* center filter */
    m1=(sf_n(f1)+1)/2-1; m1=sf_n(f1)>1?m1:0;
    m2=(sf_n(f2)+1)/2-1; m2=sf_n(f2)>1?m2:0;
    m3=(sf_n(f3)+1)/2-1; m3=sf_n(f3)>1?m3:0;
    if(verb) sf_warning("m1=%d m2=%d m3=%d",m1,m2,m3);

    /*------------------------------------------------------------*/
    for(        k3=0; k3<sf_n(f3); k3++) {
	for(    k2=0; k2<sf_n(f2); k2++) {
	    for(k1=0; k1<sf_n(f1); k1++) {
		
		for(j3=0; j3<sf_n(a3); j3++) {
		    i3=SF_MIN(SF_MAX(j3+k3-m3,0),sf_n(a3)-1);
		    
		    for(j2=0; j2<sf_n(a2); j2++) {
			i2=SF_MIN(SF_MAX(j2+k2-m2,0),sf_n(a2)-1);
			
			for(j1=0; j1<sf_n(a1); j1++) {
			    i1=SF_MIN(SF_MAX(j1+k1-m1,0),sf_n(a1)-1);
			    
			    y[j3][j2][j1] += x[i3][i2][i1] * f[k3][k2][k1];
			    
			} /* j1 loop */				
		    } /* j2 loop */
		} /* j3 loop */
		
	    } /* k1 loop */
	} /* k2 loop */
    } /* k3 loop */

    /*------------------------------------------------------------*/
    /* write data */
    sf_floatwrite(y[0][0],sf_n(a1)*sf_n(a2)*sf_n(a3),Fy);
    /*------------------------------------------------------------*/

    exit (0);
}
