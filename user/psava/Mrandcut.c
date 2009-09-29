/* cut a random dataset from a 3D cube */

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

#include <math.h>
#include <rsf.h>

#define BOUND(i,n) (i<0) ? 0 : ( (i>n-1) ? n-1 : i );

int main(int argc, char* argv[])
{
    bool verb;
    int  axis;

    sf_file Fr,Fo,Fi;    /* I/O files */
    sf_axis a1,a2,a3,ar; /* cube axes */
    int     i1,i2,i3,ir,jr;
    int     n1,n2,n3,nr;
    
    int     *rr=NULL;
    float ***ui=NULL;
    float ***uo=NULL;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getint ("axis",&axis)) axis=2;     /* stack axis */

    Fi = sf_input ("in" );
    Fo = sf_output("out");
    Fr = sf_input ("rr" );

    /* read axes */
    a1=sf_iaxa(Fi,1); if(verb) sf_raxa(a1); n1=sf_n(a1);
    a2=sf_iaxa(Fi,2); if(verb) sf_raxa(a2); n2=sf_n(a2);
    a3=sf_iaxa(Fi,3); if(verb) sf_raxa(a3); n3=sf_n(a3);
    ar=sf_iaxa(Fr,1); if(verb) sf_raxa(ar); nr=sf_n(ar);

    /*------------------------------------------------------------*/
    ui=sf_floatalloc3(n1,n2,n3);
    sf_floatread(ui[0][0],n1*n2*n3,Fi);

    rr=sf_intalloc(nr);
    sf_intread(rr,nr,Fr);
    /*------------------------------------------------------------*/
    switch(axis) {
	case 3:
	    if(verb) sf_warning("cut on axis 3");
	    sf_oaxa(Fo,ar,3);
	    uo=sf_floatalloc3(n1,n2,nr);

	    for        (ir=0; ir<nr; ir++) { jr=BOUND( rr[ir] ,n3);
		for    (i2=0; i2<n2; i2++) {
		    for(i1=0; i1<n1; i1++) {
			uo[ir][i2][i1] = ui[jr][i2][i1];			
		    }
		}
	    }
	    sf_floatwrite(uo[0][0],n1*n2*nr,Fo);
    
	    break;
	case 2:
	    if(verb) sf_warning("cut on axis 2");
	    sf_oaxa(Fo,ar,2);
	    uo=sf_floatalloc3(n1,nr,n3);

	    for        (i3=0; i3<n3; i3++) {
		for    (ir=0; ir<nr; ir++) { jr=BOUND( rr[ir] ,n2);
		    for(i1=0; i1<n1; i1++) {
			uo[i3][ir][i1] = ui[i3][jr][i1];			
		    }
		}
	    }
	    sf_floatwrite(uo[0][0],n1*nr*n3,Fo);

	    break;
	case 1:
	default:
	    if(verb) sf_warning("cut on axis 1");
	    sf_oaxa(Fo,ar,1);
	    uo=sf_floatalloc3(nr,n2,n3);
    
	    for        (i3=0; i3<n3; i3++) {
		for    (i2=0; i2<nr; i2++) {
		    for(ir=0; ir<nr; ir++) { jr=BOUND( rr[ir] ,n1);
			uo[i3][i2][ir] = ui[i3][i2][jr];			
		    }
		}
	    }
	    sf_floatwrite(uo[0][0],nr*n2*n3,Fo);

	    break;
    }

    free(**ui); free(*ui); free(ui);
    free(**uo); free(*uo); free(uo);
    /*------------------------------------------------------------*/


    exit (0);
}


