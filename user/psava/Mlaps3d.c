/* OpenMP lagged-products in the time-domain */

/*
  Copyright (C) 2011 Colorado School of Mines
  
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
 * inputs:      wavefields organized as z-x-y-t(w)
 * output: lagged products organized as hz-hx-hy-ht @ irregular points
 */

#include <rsf.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define rCOR(a,b) (a*b)

int main(int argc, char* argv[])
{
    bool verb,buf;

    sf_file Fs,Fr,Fi,Fc;           /* I/O files */
    sf_axis ax,ay,az,at,ac,aa;     /* cube axes */
    int     nx,ny,nz,nt,nhx, nhy, nhz, nht,nc;
    int              it,ihx, ihy, ihz, iht,ic;
    int                 nhx2,nhy2,nhz2,nht2;
    
    float *****imag=NULL;
    float  ****r_us=NULL;
    float  ****r_ur=NULL; 

    pt3d *cc=NULL;
    bool *ccin=NULL;
    float cxmin,cymin,czmin;
    float cxmax,cymax,czmax;
    int  icx, icy, icz;
    int  mcx, mcy, mcz, mct;
    int  pcx, pcy, pcz, pct;

    int **mcxall, **pcxall;
    int **mcyall, **pcyall;
    int **mczall, **pczall;
    int  *mctall,  *pctall;
    int lht;


    /*------------------------------------------------------------*/
    /* init RSF */
    sf_init(argc,argv);

    /* OMP parameters */
#ifdef _OPENMP
    int ompnth=1;
    ompnth=omp_init();
    if(!ompnth)
        abort();
#endif

    if(! sf_getbool("verb",&verb)) verb=false;         /* verbosity flag */
    if(! sf_getbool("buf", &buf)) buf=false;
    
    Fs = sf_input ("in" ); /*   source wavefield */
    Fr = sf_input ("ur" ); /* receiver wavefield */
    Fc = sf_input ("cc" ); /* CIP coordinates    */
    Fi = sf_output("out"); /* image */

    /* read axes */
    az=sf_iaxa(Fs,1); sf_setlabel(az,"z"); nz = sf_n(az);
    ax=sf_iaxa(Fs,2); sf_setlabel(ax,"x"); nx = sf_n(ax);
    ay=sf_iaxa(Fs,3); sf_setlabel(ay,"y"); ny = sf_n(ay);
    at=sf_iaxa(Fs,4); sf_setlabel(at,"t"); nt = sf_n(at);

    /* CIP coordinates */
    ac = sf_iaxa(Fc,2); sf_setlabel(ac,"cc"); sf_setunit(ac,"");
    nc = sf_n(ac);

    if(! sf_getint("nhx",&nhx)) nhx=0; /* number of lags on the x axis */
    if(! sf_getint("nhy",&nhy)) nhy=0; /* number of lags on the y axis */
    if(! sf_getint("nhz",&nhz)) nhz=0; /* number of lags on the z axis */
    if(! sf_getint("nht",&nht)) nht=0; /* number of lags on the t axis */

    nhx2=2*nhx+1;
    nhy2=2*nhy+1;
    nhz2=2*nhz+1;
    nht2=2*nht+1;

    if(verb) {
	sf_warning("nhx=%3d nhy=%3d nhz=%3d nht=%3d",nhx2,nhy2,nhz2,nht2);

	sf_raxa(ax);
	sf_raxa(ay);
	sf_raxa(az);
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

    aa=sf_maxa(nhy2,-nhy*sf_d(ay),sf_d(ay)); 
    sf_setlabel(aa,"hy"); sf_setunit(aa,"");
    if(verb) sf_raxa(aa);
    sf_oaxa(Fi,aa,3);

    aa=sf_maxa(nht2,-nht*sf_d(at),sf_d(at));
    sf_setlabel(aa,"ht"); sf_setunit(aa,"");
    if(verb) sf_raxa(aa);
    sf_oaxa(Fi,aa,4);

    sf_oaxa(Fi,ac,5);

    /* work arrays */
    r_us=sf_floatalloc4(nz,nx,ny,nht2);
    r_ur=sf_floatalloc4(nz,nx,ny,nht2);
    imag=sf_floatalloc5(nhz2,nhx2,nhy2,nht2,nc);

    /*------------------------------------------------------------*/
    /* CIP coordinates */
    cc= (pt3d*) sf_alloc(nc,sizeof(*cc));
    pt3dread1(Fc,cc,nc,3);

    mcxall=sf_intalloc2(nhx2,nc);
    pcxall=sf_intalloc2(nhx2,nc);

    mcyall=sf_intalloc2(nhy2,nc);
    pcyall=sf_intalloc2(nhy2,nc);

    mczall=sf_intalloc2(nhz2,nc);
    pczall=sf_intalloc2(nhz2,nc);

    ccin=sf_boolalloc(nc);

    cxmin = sf_o(ax) +             nhx *sf_d(ax);
    cxmax = sf_o(ax) + (sf_n(ax)-1-nhx)*sf_d(ax);

    cymin = sf_o(ay) +             nhy *sf_d(ay);
    cymax = sf_o(ay) + (sf_n(ay)-1-nhy)*sf_d(ay);

    czmin = sf_o(az) +             nhz *sf_d(az);
    czmax = sf_o(az) + (sf_n(az)-1-nhz)*sf_d(az);

    if(verb) {
	sf_warning("cxmin=%f,cxmax=%f",cxmin,cxmax);
	sf_warning("cymin=%f,cymax=%f",cymin,cymax);
	sf_warning("czmin=%f,czmax=%f",czmin,czmax);
    }
	
    for(ic=0; ic<nc; ic++) {
	ccin[ic]=(cc[ic].x>=cxmin && cc[ic].x<=cxmax &&
		  cc[ic].y>=cymin && cc[ic].y<=cymax &&
		  cc[ic].z>=czmin && cc[ic].z<=czmax)?true:false;

	if(ccin[ic]) {

	    icx = 0.5+(cc[ic].x-sf_o(ax))/sf_d(ax);
	    for(ihx=-nhx; ihx<nhx+1; ihx++) {
		mcxall[ic][nhx+ihx] = icx-ihx;
		pcxall[ic][nhx+ihx] = icx+ihx;
	    }

	    icy = 0.5+(cc[ic].y-sf_o(ay))/sf_d(ay);
	    for(ihy=-nhy; ihy<nhy+1; ihy++) {
		mcyall[ic][nhy+ihy] = icy-ihy;
		pcyall[ic][nhy+ihy] = icy+ihy;
	    }

	    icz = 0.5+(cc[ic].z-sf_o(az))/sf_d(az);
	    for(ihz=-nhz; ihz<nhz+1; ihz++) {
		mczall[ic][nhz+ihz] = icz-ihz;
		pczall[ic][nhz+ihz] = icz+ihz;
	    }

/*	    sf_warning("ic=%d icz=%d mczall=%d pczall=%d",*/
/*		       ic,icz,mczall[ic][0],pczall[ic][0]);*/

	}
    }
    
    mctall=sf_intalloc(nht2);
    pctall=sf_intalloc(nht2);
    for (iht=0; iht<nht2; iht++) { 
	mctall[iht]=iht;
	pctall[iht]=2*nht-iht;
    }
    

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)	\
    private(ic,ihy, ihx, ihz, iht)		\
    shared (nc,nhy2,nhx2,nhz2,nht2,imag)
#endif	

    /* zero image */
    for(ic=0; ic<nc; ic++) {
	for            (iht=0; iht<nht2; iht++) {
	    for        (ihy=0; ihy<nhy2; ihy++) {
		for    (ihx=0; ihx<nhx2; ihx++) {
		    for(ihz=0; ihz<nhz2; ihz++) {
			imag[ic][iht][ihy][ihx][ihz] = 0;
		    }
		}
	    }
	}
    }
    
    /* read wavefield @ [0...2nht-1]*/
    sf_floatread(r_us[0][0][0],nx*ny*nz*(2*nht),Fs);
    sf_floatread(r_ur[0][0][0],nx*ny*nz*(2*nht),Fr);
    
    lht=2*nht;

    if(verb) fprintf(stderr,"nt\n");
    for(it=nht;it<nt-nht;it++) { /* loop over time */
	if(verb) fprintf(stderr,"\b\b\b\b\b\b\b\b\b\b%04d",nt-nht-1-it);
	
	/* read wavefield @ [2nht]*/
	sf_floatread(r_us[ lht ][0][0],nx*ny*nz,Fs);
	sf_floatread(r_ur[ lht ][0][0],nx*ny*nz,Fr);	

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)				\
    private(ic,								\
	    ihx,ihy,ihz,iht,						\
	    mcx,mcy,mcz,mct,						\
	    pcx,pcy,pcz,pct) 						\
    shared (nc,imag,r_us,r_ur,						\
	    nhx2,  nhy2,  nhz2,  nht2,					\
	    mcxall,mcyall,mczall,mctall,				\
	    pcxall,pcyall,pczall,pctall,ccin)
#endif
	for(ic=0; ic<nc; ic++) {
	    if(ccin[ic]) {
		
		for            (iht=0; iht<nht2; iht++) { mct=mctall    [iht]; pct=pctall    [iht];
		    for        (ihy=0; ihy<nhy2; ihy++) { mcy=mcyall[ic][ihy]; pcy=pcyall[ic][ihy];
			for    (ihx=0; ihx<nhx2; ihx++) { mcx=mcxall[ic][ihx]; pcx=pcxall[ic][ihx];
			    for(ihz=0; ihz<nhz2; ihz++) { mcz=mczall[ic][ihz]; pcz=pczall[ic][ihz];
				
				imag[ic][iht][ihy][ihx][ihz] += rCOR(r_us[mct][mcy][mcx][mcz],
								     r_ur[pct][pcy][pcx][pcz]);

/*				sf_warning("pct=%d pcy=%d pcx=%d pcz=%d",pct,pcy,pcx,pcz);*/
/*				imag[ic][iht][ihy][ihx][ihz] += r_ur[0][0][0][0];*/

				
			    } // ihz
			}     // ihx
		    }         // ihy
		}             // iht
		
	    }
	} // ic
	
	/* update wavefield index (cycle) */
	for(iht=0;iht<nht2;iht++) {
	    mctall[iht] = (mctall[iht]+1) % nht2;
	    pctall[iht] = (pctall[iht]+1) % nht2;
	}
	lht = (lht+1) % nht2;
	
    } // it
    if(verb) fprintf(stderr,"\n");
    
    /* write image */
    sf_floatwrite(imag[0][0][0][0],nc*(nhx2)*(nhy2)*(nhz2)*(nht2),Fi);
    
    /*------------------------------------------------------------*/
    free(****imag); free(***imag); free(**imag); free(*imag); free(imag);
    
    free(**r_us); free(*r_us); free(r_us);
    free(**r_ur); free(*r_ur); free(r_ur);
    
    free(cc);
    free(ccin);
    
    free(*mcxall); free(mcxall);
    free(*pcxall); free(pcxall);
    free(*mcyall); free(mcyall);
    free(*pcyall); free(pcyall);
    free(*mczall); free(mczall);
    free(*pczall); free(pczall);

    /*------------------------------------------------------------*/
    exit (0);
}
