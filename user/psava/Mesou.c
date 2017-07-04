/* source for quasistatic electric potential */
/*
  Copyright (C) 2012 Colorado School of Mines
  
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
    bool verb,zbnd;
    int ompnth=1;

    /* I/O files */
    sf_file Fpre=NULL; /* pressure */
    sf_file Feso=NULL; /* electric source */
    sf_file Fske=NULL; /* Skempton coef */
    sf_file Fqke=NULL; /* Q k/etaw */

    /* cube axes */
    sf_axis az,ax,at;
    int     nz,nx,nt;
    int     iz,ix,it;
    float   dz,dx,idx,idz;

    /* I/O arrays */
    float **pre=NULL; /* pressure */
    float **eso=NULL; /* electric source */
    float **ske=NULL; /* Skempton coef */
    float **qke=NULL; /* Q k/etaw */

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    /* init OMP */
    ompnth=omp_init();
    /*------------------------------------------------------------*/

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("zbnd",&zbnd)) zbnd=false; /*  boundary flag */

    /*------------------------------------------------------------*/
    /* I/O files */
    Fpre = sf_input ("in" ); /* pressure */
    Feso = sf_output("out"); /* electric potential */
    Fske = sf_input ("ske"); /* Skempton coef */
    Fqke = sf_input ("qke"); /* Q k/etaw */

    /*------------------------------------------------------------*/
    /* axes */
    az = sf_iaxa(Fpre,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* depth */
    ax = sf_iaxa(Fpre,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* space */
    at = sf_iaxa(Fpre,3); sf_setlabel(at,"t"); if(verb) sf_raxa(at); /* time  */

    nx = sf_n(ax); dx = sf_d(ax); dx *= dx; idx=0.25/dx; /* 1/(2dx)^2 */
    nz = sf_n(az); dz = sf_d(az); dz *= dz; idz=0.25/dz; /* 1/(2dz)^2 */
    nt = sf_n(at);
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* Skempton coef */
    if(verb) fprintf(stderr,"Skempton");

    ske = sf_floatalloc2(nz,nx); 
    sf_floatread(ske[0],nz*nx,Fske);

    if(verb) fprintf(stderr," OK\n");
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* permeability */
    if(verb) fprintf(stderr,"Qk/e");
    
    qke = sf_floatalloc2(nz,nx); 
    sf_floatread(qke[0],nz*nx,Fqke); 

    if(verb) fprintf(stderr," OK\n");
    /*------------------------------------------------------------*/

    
    /*------------------------------------------------------------*/
    /* initialize electric source */
    eso = sf_floatalloc2(nz,nx);
    for    (ix=0; ix<nx; ix++) {
	for(iz=0; iz<nz; iz++) {
	    eso[ix][iz]=0;
	}
    }

    /* allocate pressure */
    pre = sf_floatalloc2(nz,nx); 

    for(it=0;it<nt;it++) {

	/*------------------------------------------------------------*/
	/* P: pressure */
	sf_floatread(pre[0],nz*nx,Fpre); 

	/* B P */
	for    (ix=0; ix<nx; ix++) {
	    for(iz=0; iz<nz; iz++) {
		pre[ix][iz]*=ske[ix][iz];
	    }
	}
	/*------------------------------------------------------------*/
	
	/*------------------------------------------------------------*/
	/* electric source */
	for    (ix=1; ix<nx-1; ix++) {
	    for(iz=1; iz<nz-1; iz++) {
		eso[ix][iz]
		    = idx* (qke[ix+1][iz  ]-qke[ix-1][iz  ])*(pre[ix+1][iz  ]-pre[ix-1][iz  ])
		    + idz* (qke[ix  ][iz+1]-qke[ix  ][iz-1])*(pre[ix  ][iz+1]-pre[ix  ][iz-1])
		    + idx* qke[ix][iz] * (pre[ix-1][iz  ]-2*pre[ix][iz]+pre[ix+1][iz  ]) 
		    + idz* qke[ix][iz] * (pre[ix  ][iz-1]-2*pre[ix][iz]+pre[ix  ][iz+1]);
	    }
	}
	/*------------------------------------------------------------*/
	
	/*------------------------------------------------------------*/
	/* write electric source */
	sf_floatwrite(eso[0],nz*nx,Feso);
	/*------------------------------------------------------------*/
    }

    /*------------------------------------------------------------*/
    /* deallocate arrays */
    free(*pre); free(pre);
    free(*eso); free(eso);
    free(*qke); free(qke);
    free(*ske); free(ske);

    exit (0);
}

