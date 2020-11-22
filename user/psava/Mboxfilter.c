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

#define INBOUND(imin,imax,i) ((i>=imin && i<imax)?true:false)

int main(int argc, char* argv[])
{
    bool verb,stat;
    int ompnth=1,ompith = 0;

    sf_file Fx=NULL; /* input  */
    sf_file Fy=NULL; /* output */
    sf_file Ff=NULL; /* filter - assumes centered filter */

    float ***x=NULL; /* data in  */
    float ***y=NULL; /* data out */
    float ***fs=NULL; /*     stationary filter */
    float****fn=NULL; /* non-stationary filter */
    
    /* cube axes */
    sf_axis a1,a2,a3;
    sf_axis f1,f2,f3;
    
    int i1,i2,i3;
    int j1,j2,j3;
    int k1,k2,k3;
    int m1,m2,m3;
    
    off_t start;
    float f;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);
    
    /* OMP parameters */
#ifdef _OPENMP
    ompnth=omp_init();
#endif

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("stat",&stat)) stat=true;  /* stationary flag */

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
    x =sf_floatalloc3(sf_n(a1),sf_n(a2),sf_n(a3)); /* data in  */
    y =sf_floatalloc3(sf_n(a1),sf_n(a2),sf_n(a3)); /* data out */ 

    if(stat) 
	fs=sf_floatalloc3(       sf_n(f1),sf_n(f2),sf_n(f3)); /* filter */
    else
	fn=sf_floatalloc4(ompnth,sf_n(f1),sf_n(f2),sf_n(f3)); /* filter */

    /*------------------------------------------------------------*/
    /* read data  */
    sf_floatread(x[0][0],sf_n(a1)*sf_n(a2)*sf_n(a3),Fx);
    /*------------------------------------------------------------*/

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

    /*------------------------------------------------------------*/
    m1=(sf_n(f1)+1)/2-1; m1=sf_n(f1)>1?m1:0;
    m2=(sf_n(f2)+1)/2-1; m2=sf_n(f2)>1?m2:0;
    m3=(sf_n(f3)+1)/2-1; m3=sf_n(f3)>1?m3:0;
    if(verb) sf_warning("m3=%d m2=%d m1=%d",m3,m2,m1);

    /*------------------------------------------------------------*/
    if(stat) {
	/* read filter*/
	sf_floatread(fs[0][0],sf_n(f1)*sf_n(f2)*sf_n(f3),Ff);

	for(k3=0; k3<sf_n(f3); k3++) {   
	    for(k2=0; k2<sf_n(f2); k2++) {		
		for(k1=0; k1<sf_n(f1); k1++) {

		    f = fs[k3][k2][k1];

/*#ifdef _OPENMP*/
/*#pragma omp parallel					\*/
/*    private(j3,j2,j1,i3,i2,i1)		\*/
/*    shared (a3,a2,a1,k3,k2,k1,m3,m2,m1,y,x,f)*/
/*#endif*/
		    for        (j3=0; j3<sf_n(a3); j3++) { i3=j3-k3+m3;
			for    (j2=0; j2<sf_n(a2); j2++) { i2=j2-k2+m2;	    
			    for(j1=0; j1<sf_n(a1); j1++) { i1=j1-k1+m1;
		    
				if( INBOUND(0,sf_n(a3),i3) &&
				    INBOUND(0,sf_n(a2),i2) &&
				    INBOUND(0,sf_n(a1),i1) ) {
				    y[j3][j2][j1] += x[i3][i2][i1] * f;
				}
			    }
			}
		    }
		}
	    }
	}

    } else {
	
	start = sf_tell(Ff);
	
	for        (j3=0; j3<sf_n(a3); j3++) {

#ifdef _OPENMP
#pragma omp parallel	       \
    private(j2,j1,k3,k2,k1,i3,i2,i1)			\
    shared (j3,   f3,f2,f1,a3,a2,a1,m3,m2,m1,y,x,fn)
#endif
	    for    (j2=0; j2<sf_n(a2); j2++) {
#ifdef _OPENMP	    
		ompith = omp_get_thread_num();
#endif
		for(j1=0; j1<sf_n(a1); j1++) {
		    
#ifdef _OPENMP	    
#pragma omp critical
#endif
		    {
			/* read filter*/
			sf_seek(Ff,start+(j3*sf_n(a1)*sf_n(a2)+
					  j2*sf_n(a1)         +
					  j1)*sizeof(float),SEEK_SET);
			sf_floatread(fn[ompith][0][0],sf_n(f1)*sf_n(f2)*sf_n(f3),Ff);
		    }

		    for(k3=0; k3<sf_n(f3); k3++) {
			i3=j3-k3+m3;
			if( INBOUND(0,sf_n(a3),i3)) {	
			    
			    for(k2=0; k2<sf_n(f2); k2++) {
				i2=j2-k2+m2;
				if( INBOUND(0,sf_n(a2),i2)) {			
				    
				    for(k1=0; k1<sf_n(f1); k1++) {
					i1=j1-k1+m1;
					if( INBOUND(0,sf_n(a1),i1)) {   

					    y[j3][j2][j1] += x[i3][i2][i1] * fn[ompith][k3][k2][k1];
  
					}
				    } /* k1 loop */	    
				}
			    } /* k2 loop */	    
		    	}
		    } /* k3 loop */
		    
		} /* j1 loop */
	    } /* j2 loop */	
	} /* j3 loop */	

    }

    /*------------------------------------------------------------*/
    /* write data */
    sf_floatwrite(y[0][0],sf_n(a1)*sf_n(a2)*sf_n(a3),Fy);
    /*------------------------------------------------------------*/


    exit (0);
}
