/* compute quasi-static electric potential */
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
    bool verb,cbnd,csrf;
    int meth;
    int ompnth=1, ompith=0;
    int nit;

    /* I/O files */
    sf_file Feso=NULL; /* electric source */
    sf_file Fepo=NULL; /* electric field  */
    sf_file Fcon=NULL; /* conductivity    */

    /* cube axes */
    sf_axis az,ax,as;
    int     nz,nx,ns;
    int     iz,ix,is,it;
    float   idz,idx;

    /* I/O arrays */
    float***els=NULL; /* electric source */
    float***epo=NULL; /* electric field  */
    float***tmp=NULL; /* temporary array */
    float **jnk=NULL;

    float **con=NULL; /* conductivity    */
    float   cnx_, cnz_, cox_, coz_, cxz_;
    float **cxp,**cxm,**czp,**czm,**cxz;
    float *epo1d,*els1d,*tmp1d;
    float *cxp1d,*cxm1d,*czp1d,*czm1d;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    /* init OMP */
    ompnth=omp_init();
    /*------------------------------------------------------------*/

    if(! sf_getbool("verb",&verb)) verb=false; /*      verbosity flag */
    if(! sf_getbool("cbnd",&cbnd)) cbnd=true;  /* conductive boundary */
    if(! sf_getbool("csrf",&csrf)) csrf=false; /* conductive  surface */
    if(! sf_getint ("meth",&meth)) meth=0;     /*         method flag */
    if(! sf_getint ( "nit",&nit))   nit=100000;/*   Jacobi iterations */ 

    /*------------------------------------------------------------*/
    /* I/O files */
    Feso = sf_input ("in" ); /* electric source */
    Fepo = sf_output("out"); /* electric field  */
    Fcon = sf_input ("con"); /* conductivity    */

    /*------------------------------------------------------------*/
    /* axes */
    az = sf_iaxa(Feso,1); sf_setlabel(az,"z"); if(verb) sf_raxa(az); /* z */
    ax = sf_iaxa(Feso,2); sf_setlabel(ax,"x"); if(verb) sf_raxa(ax); /* x */
    as = sf_iaxa(Feso,3); sf_setlabel(as,"s"); if(verb) sf_raxa(as); /* sources */

    nx = sf_n(ax); idx = 1./sf_d(ax); idx*=idx; /* 1/dx^2 */
    nz = sf_n(az); idz = 1./sf_d(az); idz*=idz; /* 1/dz^2 */
    ns = sf_n(as); /* number of sources */
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* prepare conductivity */
    con = sf_floatalloc2(nz,nx); 
    cxp = sf_floatalloc2(nz,nx); 
    cxm = sf_floatalloc2(nz,nx); 
    czp = sf_floatalloc2(nz,nx); 
    czm = sf_floatalloc2(nz,nx); 
    cxz = sf_floatalloc2(nz,nx); 

    /* read conductivity */
    sf_floatread(con[0],nz*nx,Fcon);

    /* compute conductivity gradient */
    /* use centered FD stencil */
    for    (ix=1; ix<nx-1; ix++) {
	for(iz=1; iz<nz-1; iz++) {	    
	    cnx_ = con[ix][iz]*idx;/* con(ix,iz)/dx^2 */
	    cnz_ = con[ix][iz]*idz;/* con(ix,iz)/dz^2 */	    
	    cox_ = (con[ix+1][iz  ]-con[ix-1][iz  ])*(0.25*idx);/* (con(ix+1,iz)-con(ix-1,iz))/(4dx^2) */
	    coz_ = (con[ix  ][iz+1]-con[ix  ][iz-1])*(0.25*idz);/* (con(ix,iz+1)-con(ix,iz-1))/(4dz^2) */
	    cxz_ = 2*(cnx_+cnz_);

	    cxm[ix][iz] = (cnx_-cox_)/cxz_;
	    cxp[ix][iz] = (cnx_+cox_)/cxz_;
	    czm[ix][iz] = (cnz_-coz_)/cxz_;
	    czp[ix][iz] = (cnz_+coz_)/cxz_;
	    cxz[ix][iz] = cxz_;
	}
    }

    /*------------------------------------------------------------*/ 
    if(verb) sf_warning("allocate arrays");
    tmp = sf_floatalloc3(nz,nx,ompnth); 
    els = sf_floatalloc3(nz,nx,ompnth);
    epo = sf_floatalloc3(nz,nx,ompnth); 
    if(verb) sf_warning("OK");
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    if(verb) sf_warning("write output file");
    for    (ix=0; ix<nx; ix++) {
	for(iz=0; iz<nz; iz++) {
	    tmp[0][ix][iz] = 0;
	}
    }
    for(is=0;is<ns;is++) {
	sf_seek(Fepo,is*nz*nx*sizeof(float),SEEK_SET);
	sf_floatwrite(tmp[0][0],nz*nx,Fepo);
    }
    if(verb) sf_warning("OK");
    /*------------------------------------------------------------*/

    if(verb) fprintf(stderr,"\n");
#ifdef _OPENMP
#pragma omp parallel							\
    for									\
    schedule(dynamic)							\
    private(ompith, it,is,ix,iz,tmp1d,els1d,epo1d,czm1d,czp1d,cxm1d,cxp1d,jnk) \
    shared( ompnth,nit,ns,nx,nz,tmp,  els,  epo,  czm,  czp,  cxm,  cxp,  cxz,csrf,cbnd)
#endif
    /* loop over sources */
    for(is=0;is<ns;is++) {
	ompith=omp_get_thread_num();

	/* initialize arrays */
	for    (ix=0; ix<nx; ix++) {
	    for(iz=0; iz<nz; iz++) {
		epo[ompith][ix][iz] = 0;
	    }
	}

	/* read electric source */
#ifdef _OPENMP
#pragma omp critical
#endif 
	{
	    sf_seek(Feso,is*nz*nx*sizeof(float),SEEK_SET);
	    sf_floatread(els[ompith][0],nz*nx,Feso);
	}	

	for    (ix=1; ix<nx-1; ix++) {
	    for(iz=1; iz<nz-1; iz++) {
		els[ompith][ix][iz] /= cxz[ix][iz];
	    }
	}

	// start Jacobi method
	it=0;
	while(it<nit) {
	    if(verb & (it%1000==0)) 
		fprintf(stderr,"is=%03d (ith=%02d) | %5.1f%% \n",is,ompith,(100.0*it/nit));

	    for(ix=1;ix<nx-1;ix++) {
		tmp1d = tmp[ompith][ix];
		els1d = els[ompith][ix];
		
		czm1d = czm[ix];
		czp1d = czp[ix];
		cxm1d = cxm[ix];
		cxp1d = cxp[ix];
		    
		epo1d = epo[ompith][ix];
		for(iz=1;iz<nz-1;iz++) {
		    tmp1d[iz] 
			= czm1d[iz] * epo1d[iz-1]
			+ czp1d[iz] * epo1d[iz+1] 
			+ els1d[iz];
		}
		
		epo1d = epo[ompith][ix-1];
		for(iz=1; iz<nz-1; iz++) {
		    tmp1d[iz] += cxm1d[iz] * epo1d[iz];
		}
		    
		epo1d = epo[ompith][ix+1];
		for(iz=1; iz<nz-1; iz++) {
		    tmp1d[iz] += cxp1d[iz] * epo1d[iz];
		}    
	    }
	    
	    // circulate arrays
	    jnk=epo[ompith];
	    epo[ompith]=tmp[ompith];
	    tmp[ompith]=jnk;
	    
	    // conductive surface boundary
	    if(csrf) {
		for(ix=0; ix<nx; ix++) {
		    epo[ompith][ix][0]=epo[ompith][ix][1];
		}	
	    }
	    
	    // conductive boundaries
	    if(cbnd) {
		for(ix=0; ix<nx; ix++) {
		    epo[ompith][ix][nz-1]=epo[ompith][ix][nz-2];
		}	
		for(iz=0; iz<nz; iz++) {
		    epo[ompith][   0][iz]=epo[ompith][   1][iz];
		    epo[ompith][nx-1][iz]=epo[ompith][nx-2][iz];
		}
	    }
	    	 
	    it++;
	}
	// end Jacobi method 
	
	/* write electric potential */
#ifdef _OPENMP
#pragma omp critical
#endif 
	{
	    sf_seek(Fepo,is*nz*nx*sizeof(float),SEEK_SET);
	    sf_floatwrite(epo[ompith][0],nz*nx,Fepo);
	}

    } // end source iteration
    if(verb) fprintf(stderr,"\n");

    /*------------------------------------------------------------*/
    /* deallocate arrays */
    free(**tmp); free(*tmp); free(tmp);
    free(**els); free(*els); free(els);
    free(**epo); free(*epo); free(epo);

    free(*con); free(con);
    free(*cxp); free(cxp);
    free(*cxm); free(cxm);
    free(*czp); free(czp);
    free(*czm); free(czm);
    free(*cxz); free(cxz);

    exit (0);
}

