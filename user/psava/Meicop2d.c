/* Extended IC 2D */
/* 
 * ----------------------------------------
 * FORWARD OP: opr * wfl = img
 *     inputs:   z- x- t     opr
 *               z- x- t     wfl
 *              (x,z)   -cc  cip
 *     output:  hz-hx-ht-cc  img
 * ----------------------------------------
 * ADJOINT OP: opr * img = wfl
 *     inputs:   z- x- t      opr
 *              hz-hx-ht-cc   img
 *              (x,z)   -cc   cip
 *     output:   z- x- t      wfl
 * ----------------------------------------
 */

#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define CICLOOP(a)				\
    for    (ix=0;ix<nx;ix++){                   \
        for(iz=0;iz<nz;iz++){                   \
            {a}                                 \
        }                                       \
    }

#define EICLOOP(a) \
    for        (iht=0; iht<sf_n(aht); iht++) { mct=mctall    [iht]; pct=pctall    [iht]; \
        for    (ihx=0; ihx<sf_n(ahx); ihx++) { mcx=mcxall[ic][ihx]; pcx=pcxall[ic][ihx]; \
            for(ihz=0; ihz<sf_n(ahz); ihz++) { mcz=mczall[ic][ihz]; pcz=pczall[ic][ihz]; \
                {a}                                                     \
            }                                                           \
        }                                                               \
    }

/*------------------------------------------------------------*/
void applyGaussian(float ****img,
		   sf_axis ac,
		   sf_axis aht,sf_axis ahx,sf_axis ahz,
		   float   gst,float   gsx,float   gsz)
{
    int     ic,iht,ihx,ihz;
    float       gt, gx, gz;
    int        nht,nhx,nhz;

    nht=(sf_n(aht)-1)/2;
    nhx=(sf_n(ahx)-1)/2;
    nhz=(sf_n(ahz)-1)/2;

    for(ic=0;ic<sf_n(ac);ic++) {
	for        (iht=0;iht<sf_n(aht);iht++){gt=(iht-nht)*sf_d(aht); gt*=gt;
	    for    (ihx=0;ihx<sf_n(ahx);ihx++){gx=(ihx-nhx)*sf_d(ahx); gx*=gx;
		for(ihz=0;ihz<sf_n(ahz);ihz++){gz=(ihz-nhz)*sf_d(ahz); gz*=gz;
		    img[ic][iht][ihx][ihz]*=exp(-gt*gst-gx*gsx-gz*gsz);
		}
	    }
	}
    }
}
/*------------------------------------------------------------*/
void applyScaling(float ****img,
		  sf_axis ac,
		  sf_axis aht,sf_axis ahx,sf_axis ahz,
		  float scale)
{
    int     ic,iht,ihx,ihz;

    for(ic=0;ic<sf_n(ac);ic++) {
	for        (iht=0;iht<sf_n(aht);iht++){
	    for    (ihx=0;ihx<sf_n(ahx);ihx++){
		for(ihz=0;ihz<sf_n(ahz);ihz++){
		    img[ic][iht][ihx][ihz]*=scale;
		}
	    }
	}
    }
}

/*------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    bool verb;     /* verbosity flag */
    bool pos; /* direction of spraying */
    bool adj;      /* adjoint operator flag */
    bool wflcausal, oprcausal; /* causal wfl?, opr? */

    sf_file   Fopr,       Fwfl,        Fimg,     Fcip; /* I/O files */
    float   ***opr=NULL,***wfl=NULL,****img=NULL; 
    int     itO,itW;
    
    sf_axis az,ax,at,ac,aa; /* wfld axes */
    int     nz,nx,nt,nc;
    int     iz,ix,it,ic;

    sf_axis ahx, ahz, aht; /* EIC axes */
    int     nhx, nhz, nht;
    int     ihx, ihz, iht;
    float   dhx, dhz, dht;

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
    int lht,fht; /* last buffer index */

    float scale; /* time summation scaling */
    int nslice;  /* wavefield slice size */

    bool gaus;         /* gaussian taper */
    float gsx,gsz,gst; /* std dev */

    /*------------------------------------------------------------*/
    sf_init(argc,argv);
#ifdef _OPENMP
    omp_init();
#endif

    if(! sf_getbool(    "verb",&verb    ))        verb=false; /* verbosity flag */
    if(! sf_getbool(    "positive",&pos ))        pos=true; /* if positive sprays opr to positive shits, else, sprays to negative shifts */
    if(! sf_getbool(     "adj",&adj     ))         adj=false; /* adjoint flag */
    if(! sf_getbool("wflcausal",&wflcausal)) wflcausal=false; /* causal wfl? */
    if(! sf_getbool("oprcausal",&oprcausal)) oprcausal=false; /* causal opr? */

    /*------------------------------------------------------------*/
    Fopr = sf_input ("opr" ); /* operator */
    az=sf_iaxa(Fopr,1); if(verb) sf_raxa(az); nz = sf_n(az);
    ax=sf_iaxa(Fopr,2); if(verb) sf_raxa(ax); nx = sf_n(ax);
    at=sf_iaxa(Fopr,3); if(verb) sf_raxa(at); nt = sf_n(at);

    scale = 1./nt;                /* time summation scaling */
    nslice = nz*nx*sizeof(float); /* wavefield slice */

    Fcip = sf_input ("cip" ); /* CIP coordinates    */
    ac = sf_iaxa(Fcip,2); 
    sf_setlabel(ac,"c"); sf_setunit(ac,""); if(verb) sf_raxa(ac); nc = sf_n(ac); 
    
    /*------------------------------------------------------------*/
    /* setup output */
    if(adj) {
	Fimg = sf_input ("in");  /*  read img */
	ahz=sf_iaxa(Fimg,1); nhz=(sf_n(ahz)-1)/2; if(verb) sf_raxa(ahz); 
	ahx=sf_iaxa(Fimg,2); nhx=(sf_n(ahx)-1)/2; if(verb) sf_raxa(ahx);
	aht=sf_iaxa(Fimg,3); nht=(sf_n(aht)-1)/2; if(verb) sf_raxa(aht); 

	aa=sf_maxa(1,0,1); sf_setlabel(aa,""); sf_setunit(aa,""); 

	/* set output axes */
	Fwfl = sf_output("out"); /* write wfl */
	sf_oaxa(Fwfl,az,1);
	sf_oaxa(Fwfl,ax,2);
	sf_oaxa(Fwfl,at,3);
	sf_oaxa(Fwfl,aa,4);

    } else {
	Fwfl = sf_input ( "in"); /*  read wfl */
	
	if(! sf_getint("nhz",&nhz)) nhz=0; /* z lags */
	dhz=2*sf_d(az);
	ahz=sf_maxa(2*nhz+1,-nhz*dhz,dhz); sf_setlabel(ahz,"hz"); sf_setunit(ahz,""); 
	if(verb) sf_raxa(ahz);
	
	if(! sf_getint("nhx",&nhx)) nhx=0; /* x lags */
	dhx=2*sf_d(ax);
	ahx=sf_maxa(2*nhx+1,-nhx*dhx,dhx); sf_setlabel(ahx,"hx"); sf_setunit(ahx,""); 
	if(verb) sf_raxa(ahx);

	if(! sf_getint("nht",&nht)) nht=0; /* t lags */
	dht=2*sf_d(at);
	aht=sf_maxa(2*nht+1,-nht*dht,dht); sf_setlabel(aht,"ht"); sf_setunit(aht,""); 
	if(verb) sf_raxa(aht);

	Fimg = sf_output("out"); /* write img */
	sf_oaxa(Fimg,ahz,1);
	sf_oaxa(Fimg,ahx,2);
	sf_oaxa(Fimg,aht,3);
	sf_oaxa(Fimg, ac,4);
    }

    /*------------------------------------------------------------*/
    if(! sf_getbool("gaus",&gaus)) gaus=false; /* Gaussian taper */
    if(gaus) {
	if(! sf_getfloat("gsx",&gsx)) gsx=0.25*sf_n(ahx)*sf_d(ahx); gsx=(nhx==0)?1:1./(2*gsx*gsx);
        if(! sf_getfloat("gsz",&gsz)) gsz=0.25*sf_n(ahz)*sf_d(ahz); gsz=(nhz==0)?1:1./(2*gsz*gsz);
        if(! sf_getfloat("gst",&gst)) gst=0.25*sf_n(aht)*sf_d(aht); gst=(nht==0)?1:1./(2*gst*gst);
    }
    
    /*------------------------------------------------------------*/
    /* allocate arrays */
    opr=sf_floatalloc3(nz,nx,sf_n(aht));
    wfl=sf_floatalloc3(nz,nx,sf_n(aht));
    img=sf_floatalloc4(sf_n(ahz),sf_n(ahx),sf_n(aht),sf_n(ac));

    /*------------------------------------------------------------*/
    /* CIP coordinates */
    cc= (pt2d*) sf_alloc(nc,sizeof(*cc));
    pt2dread1(Fcip,cc,nc,2);

    mcxall=sf_intalloc2(sf_n(ahx),sf_n(ac));
    pcxall=sf_intalloc2(sf_n(ahx),sf_n(ac));
    mczall=sf_intalloc2(sf_n(ahz),sf_n(ac));
    pczall=sf_intalloc2(sf_n(ahz),sf_n(ac));
    ccin=sf_boolalloc(sf_n(ac));

    cxmin = sf_o(ax) +             nhx *sf_d(ax);
    cxmax = sf_o(ax) + (sf_n(ax)-1-nhx)*sf_d(ax);
    czmin = sf_o(az) +             nhz *sf_d(az);
    czmax = sf_o(az) + (sf_n(az)-1-nhz)*sf_d(az);

    for(ic=0;ic<nc;ic++) {
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
       
    mctall=sf_intalloc(sf_n(aht));
    pctall=sf_intalloc(sf_n(aht));
    for (iht=0; iht<sf_n(aht); iht++) { 
	mctall[iht]=            iht;
	pctall[iht]=sf_n(aht)-1-iht;
    }
    
    if(adj) { /* ADJIONT OPERATOR */

	for(iht=0;iht<sf_n(aht);iht++)
	    CICLOOP( wfl[iht][ix][iz]=0; );                       /* zero wfl */
	for(it=0;it<nt;it++) sf_floatwrite(wfl[0][0],nz*nx,Fwfl); /* reserve wfl */ 
	sf_seek(Fwfl,0,SEEK_SET);                                 /* seek back */

	sf_floatread(img[0][0][0],sf_n(ac)*sf_n(ahx)*sf_n(ahz)*sf_n(aht),Fimg); /* read img */
	;        applyScaling (img,ac,aht,ahx,ahz,scale);                       /* scaling  */
	if(gaus) applyGaussian(img,ac,aht,ahx,ahz,gst,gsx,gsz);                 /* Gaussian */

	lht=0; itO=-999999; itW=-999999;
	for(it=-nht;it<nt+nht;it++) { if(verb) fprintf(stderr,"\b\b\b\b\b\b%04d",it);
	    fht=(lht+1) % sf_n(aht);

	    if(it<nt-nht) {
		itO = it + nht;
		if( !oprcausal ) sf_seek(Fopr,(off_t)(nt-2-itO)*nslice,SEEK_SET);
		else             sf_seek(Fopr,(off_t)      itO *nslice,SEEK_SET);
		sf_floatread(opr[ lht ][0],nz*nx,Fopr);
	    }

	    if(it>=0+nht && 
	       it<nt-nht) { 

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)				\
    private(ic,     ihx,ihz,iht,mcx,   mcz,   mct,   pcx,   pcz,   pct)	\
    shared (nc,ccin,ahx,ahz,aht,mcxall,mczall,mctall,pcxall,pczall,pctall)
#endif
		for(ic=0;ic<nc;ic++){ if(ccin[ic]) { /* sum over c only! */
      if(pos){
			  EICLOOP( wfl    [mct][mcx][mcz] +=
				   opr    [pct][pcx][pcz] *
				  img[ic][iht][ihx][ihz]; );
		  }else{
			  EICLOOP( wfl    [pct][pcx][pcz] +=
				   opr    [mct][mcx][mcz] *
				  img[ic][iht][ihx][ihz]; );
      }
     }
		}
	    }
            for(iht=0;iht<sf_n(aht);iht++) {
                mctall[iht] = (mctall[iht]+1) % sf_n(aht); /* cycle iht index */
                pctall[iht] = (pctall[iht]+1) % sf_n(aht);
            }

	    if(it>=0+nht) {
		itW = it - nht;
		if( !wflcausal ) sf_seek(Fwfl,(off_t)(nt-1-itW)*nslice,SEEK_SET);
		else             sf_seek(Fwfl,(off_t)      itW *nslice,SEEK_SET);
		sf_floatwrite(wfl[ fht ][0],nz*nx,Fwfl);/* write wfl */
		CICLOOP( wfl[ fht ][ix][iz]=0; );       /* reset wfl */
	    }

	    lht=(lht+1) % sf_n(aht);
	}
	if(verb) fprintf(stderr,"\n");

    } else { /* FORWARD OPERATOR */

	for(ic=0;ic<nc;ic++) EICLOOP( img[ic][iht][ihx][ihz] = 0; );            /* zero img */
	sf_floatwrite(img[0][0][0],sf_n(ac)*sf_n(ahx)*sf_n(ahz)*sf_n(aht),Fimg);/* reserve img */
	sf_seek(Fimg,0,SEEK_SET);                                               /* seek back */

	lht=0; itO=-999999; itW=-999999;
	for(it=-nht;it<nt+nht;it++) { if(verb) fprintf(stderr,"\b\b\b\b\b\b%04d",it);
	    fht=(lht+1) % sf_n(aht);

	    if(it<nt-nht) {
		itO = it + nht;
		if( !oprcausal ) sf_seek(Fopr,(off_t)(nt-1-itO)*nslice,SEEK_SET);
		else             sf_seek(Fopr,(off_t)      itO *nslice,SEEK_SET);
		sf_floatread(opr[ lht ][0],nz*nx,Fopr);

		itW = it + nht;
		if( !wflcausal ) sf_seek(Fwfl,(off_t)(nt-1-itW)*nslice,SEEK_SET);
		else             sf_seek(Fwfl,(off_t)      itW *nslice,SEEK_SET);
		sf_floatread(wfl[ lht ][0],nz*nx,Fwfl);
	    }

	    if(it>=0+nht && 
	       it<nt-nht) { 

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)				\
    private(ic,     ihx,ihz,iht,mcx,   mcz,   mct,   pcx,   pcz,   pct)	\
    shared (nc,ccin,ahx,ahz,aht,mcxall,mczall,mctall,pcxall,pczall,pctall)
#endif
		for(ic=0;ic<nc;ic++){ if(ccin[ic]) { /* sum over c and t! */
      if(pos){
		   	EICLOOP( img[ic][iht][ihx][ihz] +=
				   opr    [pct][pcx][pcz] *
			  	 wfl    [mct][mcx][mcz]; );
      }else{
		  	EICLOOP( img[ic][iht][ihx][ihz] +=
		  		 opr    [mct][mcx][mcz] *
		  		 wfl    [pct][pcx][pcz]; );
      }
    }
		}
	    }
            for(iht=0;iht<sf_n(aht);iht++) {
                mctall[iht] = (mctall[iht]+1) % sf_n(aht); /* cycle iht index */
                pctall[iht] = (pctall[iht]+1) % sf_n(aht);
            }

	    lht=(lht+1) % sf_n(aht);
	}
	if(verb) fprintf(stderr,"\n");

	if(gaus) applyGaussian(img,ac,aht,ahx,ahz,gst,gsx,gsz);                  /* Gaussian  */
	;        applyScaling (img,ac,aht,ahx,ahz,scale);                        /* scaling   */
	sf_floatwrite(img[0][0][0],sf_n(ac)*sf_n(ahx)*sf_n(ahz)*sf_n(aht),Fimg); /* write img */

    } /* end if adj */

    /*------------------------------------------------------------*/
    /* deallocate arrays */
    free(***img); free(**img); free(*img); free(img);
    ;             free(**opr); free(*opr); free(opr);
    ;             free(**wfl); free(*wfl); free(wfl);
 
    free(cc);
    free(ccin);
    free(*mcxall); free(mcxall);
    free(*pcxall); free(pcxall);
    free(*mczall); free(mczall);
    free(*pczall); free(pczall);
    /*------------------------------------------------------------*/   

    exit (0);
}


