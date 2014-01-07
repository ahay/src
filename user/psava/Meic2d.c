/* Extended IC 3D */

/*
  Copyright (C) 2013 Colorado School of Mines
  
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

/* 
 * inputs: wavefields organized as z-x-t(w)
 * output:      image organized as hz-hx-ht-cc
 */

#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>
#endif

int main(int argc, char* argv[])
{
    bool verb,isreversed;

    sf_file Fs,Fr,Fi,Fc;        /* I/O files */
    sf_axis az,ax,at,ac,aa;     /* cube axes */
    int     nz,nx,nt, nhx,  nhz, nht,nc;
    int           it, ihx,  ihz, iht,ic;
    int               nhx2,nhz2,nht2;
    off_t iseek;

    float ***us=NULL,***ur=NULL,****ii=NULL; 

    pt2d *cc=NULL;
    bool *ccin=NULL;
    float cxmin,czmin;
    float cxmax,czmax;
    int  icx, icz;
    int  mcx, mcz, mct;
    int  pcx, pcz, pct;

    int **mcxall, **pcxall;
    int **mczall, **pczall;
    int  *mctall,  *pctall;
    int lht;

    float scale;

    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    /* OMP parameters */
#ifdef _OPENMP
    omp_init();
#endif

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(! sf_getbool("isreversed",&isreversed)) isreversed=false; /* received wavefield */
    
    Fs = sf_input ("in" ); /*   source wavefield */
    Fr = sf_input ("ur" ); /* receiver wavefield */
    Fc = sf_input ("cc" ); /* CIP coordinates    */
    Fi = sf_output("out"); /* image */

    /*------------------------------------------------------------*/
    /* read axes */
    az=sf_iaxa(Fs,1); nz = sf_n(az);
    ax=sf_iaxa(Fs,2); nx = sf_n(ax);
    at=sf_iaxa(Fs,3); nt = sf_n(at);

    /* CIP coordinates */
    ac = sf_iaxa(Fc,2); sf_setlabel(ac,"cc"); sf_setunit(ac,"");
    nc = sf_n(ac); 

    if(! sf_getint("nhz",&nhz)) nhz=0; nhz2=2*nhz+1; /* z lags */
    if(! sf_getint("nhx",&nhx)) nhx=0; nhx2=2*nhx+1; /* x lags */
    if(! sf_getint("nht",&nht)) nht=0; nht2=2*nht+1; /* t lags */

    lht=2*nht;

    if(verb) {
	sf_warning("nhx=%3d nhz=%3d nht=%3d",nhx2,nhz2,nht2);

	sf_raxa(az);
	sf_raxa(ax);
	sf_raxa(at);
	sf_raxa(ac);
    }

    /* set output axes */
    aa=sf_maxa(nhz2,-nhz*sf_d(az),sf_d(az));
    sf_setlabel(aa,"hz"); sf_setunit(aa,"");
    if(verb) sf_raxa(aa);
    sf_oaxa(Fi,aa,1);
    
    aa=sf_maxa(nhx2,-nhx*sf_d(ax),sf_d(ax)); 
    sf_setlabel(aa,"hx"); sf_setunit(aa,"");
    if(verb) sf_raxa(aa);
    sf_oaxa(Fi,aa,2);

    aa=sf_maxa(nht2,-nht*sf_d(at),sf_d(at));
    sf_setlabel(aa,"ht"); sf_setunit(aa,"");
    if(verb) sf_raxa(aa);
    sf_oaxa(Fi,aa,3);

    sf_oaxa(Fi,ac,4);

    /*------------------------------------------------------------*/
    /* allocate work arrays */
    us=sf_floatalloc3(nz,nx,nht2);
    ur=sf_floatalloc3(nz,nx,nht2);
    ii=sf_floatalloc4(nhz2,nhx2,nht2,nc);
    /* zero output */
    for(ic=0; ic<nc; ic++) {
	for        (iht=0; iht<nht2; iht++) {
	    for    (ihx=0; ihx<nhx2; ihx++) {
		for(ihz=0; ihz<nhz2; ihz++) {
		    ii[ic][iht][ihx][ihz] = 0;
		}
	    }
	}
    }

    /*------------------------------------------------------------*/
    /* CIP coordinates */
    cc= (pt2d*) sf_alloc(nc,sizeof(*cc));
    pt2dread1(Fc,cc,nc,2);

    mcxall=sf_intalloc2(nhx2,nc);
    pcxall=sf_intalloc2(nhx2,nc);
    mczall=sf_intalloc2(nhz2,nc);
    pczall=sf_intalloc2(nhz2,nc);

    ccin=sf_boolalloc(nc);

    cxmin = sf_o(ax) +             nhx *sf_d(ax);
    cxmax = sf_o(ax) + (sf_n(ax)-1-nhx)*sf_d(ax);
    czmin = sf_o(az) +             nhz *sf_d(az);
    czmax = sf_o(az) + (sf_n(az)-1-nhz)*sf_d(az);
    if(verb) {
	sf_warning("cxmin=%f,cxmax=%f",cxmin,cxmax);
	sf_warning("czmin=%f,czmax=%f",czmin,czmax);
    }

    for(ic=0; ic<nc; ic++) {
	ccin[ic]=(cc[ic].x>=cxmin && cc[ic].x<=cxmax &&
		  cc[ic].z>=czmin && cc[ic].z<=czmax)?true:false;
	
	if(ccin[ic]) {

	    icx = 0.5+(cc[ic].x-sf_o(ax))/sf_d(ax);
	    for(ihx=-nhx; ihx<nhx+1; ihx++) {
		mcxall[ic][nhx+ihx] = icx-ihx;
		pcxall[ic][nhx+ihx] = icx+ihx;
	    }

	    icz = 0.5+(cc[ic].z-sf_o(az))/sf_d(az);
	    for(ihz=-nhz; ihz<nhz+1; ihz++) {
		mczall[ic][nhz+ihz] = icz-ihz;
		pczall[ic][nhz+ihz] = icz+ihz;
	    }

	}
    }
       
    mctall=sf_intalloc(nht2);
    pctall=sf_intalloc(nht2);
    for (iht=0; iht<nht2; iht++) { 
	mctall[iht]=iht;
	pctall[iht]=2*nht-iht;
    }
    
    /*------------------------------------------------------------*/
    if(isreversed) { /* receiver wavefield is reversed */

	/* read wavefield @ [0...2nht-1]*/
	for(iht=0;iht<2*nht;iht++) {
	    sf_floatread(us[iht][0],nz*nx,Fs);
	    sf_floatread(ur[iht][0],nz*nx,Fr);
	}

	if(verb) fprintf(stderr,"nt\n");
	for(it=nht;it<nt-nht;it++) {
	    if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b%04d",it);
	    
	    /* read wavefield @ [2nht]*/
	    sf_floatread(us[ lht ][0],nz*nx,Fs);
	    sf_floatread(ur[ lht ][0],nz*nx,Fr);
	    
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)				\
    private(ic,								\
	    ihx,ihz,iht,						\
	    mcx,mcz,mct,						\
	    pcx,pcz,pct) 						\
    shared (nc,ii,us,ur,						\
	    nhx2,  nhz2,  nht2,						\
	    mcxall,mczall,mctall,					\
	    pcxall,pczall,pctall,ccin)
#endif
	    for(ic=0; ic<nc; ic++) {
		if(ccin[ic]) {
		    
		    for        (iht=0; iht<nht2; iht++) { mct=mctall    [iht]; pct=pctall    [iht];
			for    (ihx=0; ihx<nhx2; ihx++) { mcx=mcxall[ic][ihx]; pcx=pcxall[ic][ihx];
			    for(ihz=0; ihz<nhz2; ihz++) { mcz=mczall[ic][ihz]; pcz=pczall[ic][ihz];
				
				ii[ic][iht][ihx][ihz] += us[mct][mcx][mcz]*ur[pct][pcx][pcz];
				
			    } /* ihz */
			}     /* ihx */
		    }         /* iht */
		    
		}
	    } /* ic */
	    
	    /* update wavefield index (cycle) */
	    for(iht=0;iht<nht2;iht++) {
		mctall[iht] = (mctall[iht]+1) % nht2;
		pctall[iht] = (pctall[iht]+1) % nht2;
	    }
	    lht = (lht+1) % nht2;
	    
	} /* it */
	if(verb) fprintf(stderr,"\n");

    } else { /* receiver wavefield is NOT reversed */

	/* read wavefield @ [0...2nht-1]*/
	for(iht=0;iht<2*nht;iht++) {
	    sf_floatread(us[iht][0],nz*nx,Fs);
	    iseek = (off_t)(nt-1-iht)*nz*nx*sizeof(float);
	    sf_seek(Fr,iseek,SEEK_SET);
	    sf_floatread(ur[iht][0],nz*nx,Fr);
	}

	if(verb) fprintf(stderr,"nt\n");
	for(it=nht;it<nt-nht;it++) {
	    if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b%04d",nt-nht-1-it);
	   
	    /* read wavefield @ [2nht]*/
	    sf_floatread(us[ lht ][0],nz*nx,Fs);
	    iseek=(off_t)(nt-nht-1-it)*nz*nx*sizeof(float);
	    sf_seek(Fr,iseek,SEEK_SET);
	    sf_floatread(ur[ lht ][0],nz*nx,Fr);
	    
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)				\
    private(ic,								\
	    ihx,ihz,iht,						\
	    mcx,mcz,mct,						\
	    pcx,pcz,pct) 						\
    shared (nc,ii,us,ur,						\
	    nhx2,  nhz2,  nht2,						\
	    mcxall,mczall,mctall,					\
	    pcxall,pczall,pctall,ccin)
#endif
	    for(ic=0; ic<nc; ic++) {
		if(ccin[ic]) {
		    
		    for        (iht=0; iht<nht2; iht++) { mct=mctall    [iht]; pct=pctall    [iht];
			for    (ihx=0; ihx<nhx2; ihx++) { mcx=mcxall[ic][ihx]; pcx=pcxall[ic][ihx];
			    for(ihz=0; ihz<nhz2; ihz++) { mcz=mczall[ic][ihz]; pcz=pczall[ic][ihz];
				
				ii[ic][iht][ihx][ihz] += us[mct][mcx][mcz]*ur[pct][pcx][pcz];
				
			    } /* ihz */
			}     /* ihx */
		    }         /* iht */
		    
		}
	    } /* ic */
	    
	    /* update wavefield index (cycle) */
	    for(iht=0;iht<nht2;iht++) {
		mctall[iht] = (mctall[iht]+1) % nht2;
		pctall[iht] = (pctall[iht]+1) % nht2;
	    }
	    lht = (lht+1) % nht2;

	} /* it */
	if(verb) fprintf(stderr,"\n");


	
    } /* end "is reversed" */
    /*------------------------------------------------------------*/

    /*------------------------------------------------------------*/
    /* scale image */
    scale = 1./nt;
    for(ic=0; ic<nc; ic++) {
	for        (iht=0; iht<nht2; iht++) {
	    for    (ihx=0; ihx<nhx2; ihx++) {
		for(ihz=0; ihz<nhz2; ihz++) {
		    ii[ic][iht][ihx][ihz] *= scale;
		}
	    }
	}
    }

    /*------------------------------------------------------------*/
    /* write image */
    sf_floatwrite(ii[0][0][0],nc*(nhx2)*(nhz2)*(nht2),Fi);
    
    /*------------------------------------------------------------*/
    /* deallocate arrays */
    free(***ii); free(**ii); free(*ii); free(ii);
    free(*us); free(us);
    free(*ur); free(ur);
 
    free(cc);
    free(ccin);
    
    free(*mcxall); free(mcxall);
    free(*pcxall); free(pcxall);
    free(*mczall); free(mczall);
    free(*pczall); free(pczall);
    /*------------------------------------------------------------*/   

    exit (0);
}
