/* Extended IC 2D */
/* 
 * ----------------------------------------
 * FORWARD OP: opr * wfl = img
 *     inputs:   z- x- y- t     opr
 *               z- x- y- t     wfl
 *              (x,y,z)    -cc  cip
 *     output:  hz-hx-hy-ht-cc  img
 * ----------------------------------------
 * ADJOINT OP: opr * img = wfl
 *     inputs:   z- x- y- t      opr
 *              hz-hx-hy-ht-cc   img
 *              (x,y,z)    -cc   cip
 *     output:   z- x- y- t      wfl
 * ----------------------------------------
 */

#include <rsf.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#define CICLOOP(a)				    \
    for    (iy=0;iy<ny;iy++){			    \
        for    (ix=0;ix<nx;ix++){                   \
	    for(iz=0;iz<nz;iz++){                   \
		{a}                                 \
	    }                                       \
	}					    \
    }

#define EICLOOP(a)							\
    for            (iht=0; iht<sf_n(aht); iht++) { mct=mctall    [iht]; pct=pctall    [iht]; \
	for        (ihy=0; ihy<sf_n(ahy); ihy++) { mcy=mcyall[ic][ihy]; pcy=pcyall[ic][ihy]; \
	    for    (ihx=0; ihx<sf_n(ahx); ihx++) { mcx=mcxall[ic][ihx]; pcx=pcxall[ic][ihx]; \
		for(ihz=0; ihz<sf_n(ahz); ihz++) { mcz=mczall[ic][ihz]; pcz=pczall[ic][ihz]; \
		    {a}							\
		}							\
	    }								\
	}								\
    }

/*------------------------------------------------------------*/
void applyGaussian(float *****img,
		   sf_axis ac,
		   sf_axis aht,sf_axis ahx,sf_axis ahy,sf_axis ahz,
		   float   gst,float   gsx,float   gsy,float   gsz)
{
    int     ic,iht,ihx,ihy,ihz;
    float       gt, gx, gy, gz;
    int        nht,nhx,nhy,nhz;

    nht=(sf_n(aht)-1)/2;
    nhx=(sf_n(ahx)-1)/2;
    nhy=(sf_n(ahy)-1)/2;
    nhz=(sf_n(ahz)-1)/2;

    for(ic=0;ic<sf_n(ac);ic++) {
	for            (iht=0;iht<sf_n(aht);iht++){gt=(iht-nht)*sf_d(aht); gt*=gt;
	    for        (ihy=0;ihy<sf_n(ahy);ihy++){gy=(ihy-nhy)*sf_d(ahy); gy*=gy;
		for    (ihx=0;ihx<sf_n(ahx);ihx++){gx=(ihx-nhx)*sf_d(ahx); gx*=gx;
		    for(ihz=0;ihz<sf_n(ahz);ihz++){gz=(ihz-nhz)*sf_d(ahz); gz*=gz;
			img[ic][iht][ihy][ihx][ihz]*=exp(-gt*gst-gx*gsx-gy*gsy-gz*gsz);
		    }
		}
	    }
	}
    }
}
/*------------------------------------------------------------*/
void applyScaling(float *****img,
		  sf_axis ac,
		  sf_axis aht,sf_axis ahx,sf_axis ahy,sf_axis ahz,
		  float scale)
{
    int     ic,iht,ihx,ihy,ihz;

    for(ic=0;ic<sf_n(ac);ic++) {
	for            (iht=0;iht<sf_n(aht);iht++){
	    for        (ihy=0;ihy<sf_n(ahy);ihy++){
		for    (ihx=0;ihx<sf_n(ahx);ihx++){
		    for(ihz=0;ihz<sf_n(ahz);ihz++){
			img[ic][iht][ihy][ihx][ihz]*=scale;
		    }
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

    sf_file    Fopr,        Fwfl,         Fimg,     Fcip; /* I/O files */
    float   ****opr=NULL,****wfl=NULL,*****img=NULL; 
    int     itO,itW;
    
    sf_axis az,ax,ay,at,ac,aa; /* wfld axes */
    int     nz,nx,ny,nt,nc;
    int     iz,ix,iy,it,ic;

    sf_axis ahx, ahy, ahz, aht; /* EIC axes */
    int     nhx, nhy, nhz, nht;
    int     ihx, ihy, ihz, iht;
    float   dhx, dhy, dhz, dht;

    pt3d *cc=NULL;
    bool *ccin=NULL;
    float cxmin,czmin,cymin;
    float cxmax,czmax,cymax;
    int  icx, icz, icy;
    int  mcx, mcz, mcy, mct;
    int  pcx, pcz, pcy, pct;
    int **mcxall, **pcxall;
    int **mcyall, **pcyall;
    int **mczall, **pczall;
    int  *mctall,  *pctall;
    int lht,fht; /* last buffer index */

    float scale; /* time summation scaling */
    int nslice;  /* wavefield slice size */

    bool gaus;         /* gaussian taper */
    float gsx,gsy,gsz,gst; /* std dev */

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
    ay=sf_iaxa(Fopr,3); if(verb) sf_raxa(ay); ny = sf_n(ay);
    at=sf_iaxa(Fopr,4); if(verb) sf_raxa(at); nt = sf_n(at);

    scale = 1./nt;                /* time summation scaling */
    nslice = nz*nx*ny*sizeof(float); /* wavefield slice */

    Fcip = sf_input ("cip" ); /* CIP coordinates    */
    ac = sf_iaxa(Fcip,2); 
    sf_setlabel(ac,"c"); sf_setunit(ac,""); if(verb) sf_raxa(ac); nc = sf_n(ac); 
    
    /*------------------------------------------------------------*/
    /* setup output */
    if(adj) {
	Fimg = sf_input ("in");  /*  read img */
	ahz=sf_iaxa(Fimg,1); nhz=(sf_n(ahz)-1)/2; if(verb) sf_raxa(ahz); 
	ahx=sf_iaxa(Fimg,2); nhx=(sf_n(ahx)-1)/2; if(verb) sf_raxa(ahx);
	ahy=sf_iaxa(Fimg,3); nhy=(sf_n(ahy)-1)/2; if(verb) sf_raxa(ahy);
	aht=sf_iaxa(Fimg,4); nht=(sf_n(aht)-1)/2; if(verb) sf_raxa(aht); 

	aa=sf_maxa(1,0,1); sf_setlabel(aa,""); sf_setunit(aa,""); 

	/* set output axes */
	Fwfl = sf_output("out"); /* write wfl */
	sf_oaxa(Fwfl,az,1);
	sf_oaxa(Fwfl,ax,2);
	sf_oaxa(Fwfl,ay,3);
	sf_oaxa(Fwfl,at,4);
	sf_oaxa(Fwfl,aa,5);

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

	if(! sf_getint("nhy",&nhy)) nhy=0; /* y lags */
	dhy=2*sf_d(ay);
	ahy=sf_maxa(2*nhy+1,-nhy*dhy,dhy); sf_setlabel(ahy,"hy"); sf_setunit(ahy,""); 
	if(verb) sf_raxa(ahy);

	if(! sf_getint("nht",&nht)) nht=0; /* t lags */
	dht=2*sf_d(at);
	aht=sf_maxa(2*nht+1,-nht*dht,dht); sf_setlabel(aht,"ht"); sf_setunit(aht,""); 
	if(verb) sf_raxa(aht);

	Fimg = sf_output("out"); /* write img */
	sf_oaxa(Fimg,ahz,1);
	sf_oaxa(Fimg,ahx,2);
	sf_oaxa(Fimg,ahy,3);
	sf_oaxa(Fimg,aht,4);
	sf_oaxa(Fimg, ac,5);
    }

    /*------------------------------------------------------------*/
    if(! sf_getbool("gaus",&gaus)) gaus=false; /* Gaussian taper */
    if(gaus) {
	if(! sf_getfloat("gsx",&gsx)) gsx=0.25*sf_n(ahx)*sf_d(ahx); gsx=(nhx==0)?1:1./(2*gsx*gsx);
	if(! sf_getfloat("gsy",&gsy)) gsy=0.25*sf_n(ahy)*sf_d(ahy); gsy=(nhy==0)?1:1./(2*gsy*gsy);
        if(! sf_getfloat("gsz",&gsz)) gsz=0.25*sf_n(ahz)*sf_d(ahz); gsz=(nhz==0)?1:1./(2*gsz*gsz);
        if(! sf_getfloat("gst",&gst)) gst=0.25*sf_n(aht)*sf_d(aht); gst=(nht==0)?1:1./(2*gst*gst);
    }
    
    /*------------------------------------------------------------*/
    /* allocate arrays */
    opr=sf_floatalloc4(nz,nx,ny,sf_n(aht));
    wfl=sf_floatalloc4(nz,nx,ny,sf_n(aht));
    img=sf_floatalloc5(sf_n(ahz),sf_n(ahx),sf_n(ahy),sf_n(aht),sf_n(ac));

    /*------------------------------------------------------------*/
    /* CIP coordinates */
    cc= (pt3d*) sf_alloc(nc,sizeof(*cc));
    pt3dread1(Fcip,cc,nc,3);

    mcxall=sf_intalloc2(sf_n(ahx),sf_n(ac));
    pcxall=sf_intalloc2(sf_n(ahx),sf_n(ac));
    mcyall=sf_intalloc2(sf_n(ahy),sf_n(ac));
    pcyall=sf_intalloc2(sf_n(ahy),sf_n(ac));
    mczall=sf_intalloc2(sf_n(ahz),sf_n(ac));
    pczall=sf_intalloc2(sf_n(ahz),sf_n(ac));
    ccin=sf_boolalloc(sf_n(ac));

    cxmin = sf_o(ax) +             nhx *sf_d(ax);
    cxmax = sf_o(ax) + (sf_n(ax)-1-nhx)*sf_d(ax);
    cymin = sf_o(ay) +             nhy *sf_d(ay);
    cymax = sf_o(ay) + (sf_n(ay)-1-nhy)*sf_d(ay);
    czmin = sf_o(az) +             nhz *sf_d(az);
    czmax = sf_o(az) + (sf_n(az)-1-nhz)*sf_d(az);

    for(ic=0;ic<nc;ic++) {
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
	    CICLOOP( wfl[iht][iy][ix][iz]=0; );                         /* zero wfl */
	for(it=0;it<nt;it++) sf_floatwrite(wfl[0][0][0],nz*nx*ny,Fwfl); /* reserve wfl */ 
	sf_seek(Fwfl,0,SEEK_SET);                                       /* seek back */

	sf_floatread(img[0][0][0][0],sf_n(ac)*sf_n(ahy)*sf_n(ahx)*sf_n(ahz)*sf_n(aht),Fimg); /* read img */
	;        applyScaling (img,ac,aht,ahx,ahy,ahz,scale);                       /* scaling  */
	if(gaus) applyGaussian(img,ac,aht,ahx,ahy,ahz,gst,gsx,gsy,gsz);             /* Gaussian */

	lht=0; itO=-999999; itW=-999999;
	for(it=-nht;it<nt+nht;it++) { if(verb) fprintf(stderr,"\b\b\b\b\b\b%04d",it);
	    fht=(lht+1) % sf_n(aht);

	    if(it<nt-nht) {
		itO = it + nht;
		if( !oprcausal ) sf_seek(Fopr,(off_t)(nt-1-itO)*nslice,SEEK_SET);
		else             sf_seek(Fopr,(off_t)      itO *nslice,SEEK_SET);
		sf_floatread(opr[ lht ][0][0],nz*nx*ny,Fopr);
	    }

	    if(it>=0+nht && 
	       it<nt-nht) { 

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)				\
    private(ic,     ihx,ihy,ihz,iht,mcx,   mcy,   mcz,   mct,   pcx,   pcy,   pcz,   pct) \
    shared (nc,ccin,ahx,ahy,ahz,aht,mcxall,mcyall,mczall,mctall,pcxall,pcyall,pczall,pctall)
#endif
		for(ic=0;ic<nc;ic++){ if(ccin[ic]) { /* sum over c only! */
      if(pos){

			EICLOOP( wfl    [mct][mcy][mcx][mcz] +=
				 opr    [pct][pcy][pcx][pcz] *
				 img[ic][iht][ihy][ihx][ihz]; );
      }else{
			EICLOOP( wfl    [pct][pcy][pcx][pcz] +=
				 opr    [mct][mcy][mcx][mcz] *
				 img[ic][iht][ihy][ihx][ihz]; );


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
		sf_floatwrite(wfl[ fht ][0][0],nz*nx*ny,Fwfl);/* write wfl */
		CICLOOP( wfl[ fht ][iy][ix][iz]=0; );       /* reset wfl */
	    }

	    lht=(lht+1) % sf_n(aht);
	}
	if(verb) fprintf(stderr,"\n");

    } else { /* FORWARD OPERATOR */

	for(ic=0;ic<nc;ic++) EICLOOP( img[ic][iht][ihy][ihx][ihz] = 0; );            /* zero img */
	sf_floatwrite(img[0][0][0][0],sf_n(ac)*sf_n(ahx)*sf_n(ahy)*sf_n(ahz)*sf_n(aht),Fimg);/* reserve img */
	sf_seek(Fimg,0,SEEK_SET);                                               /* seek back */

	lht=0; itO=-999999; itW=-999999;
	for(it=-nht;it<nt+nht;it++) { if(verb) fprintf(stderr,"\b\b\b\b\b\b%04d",it);
	    fht=(lht+1) % sf_n(aht);

	    if(it<nt-nht) {
		itO = it + nht;
		if( !oprcausal ) sf_seek(Fopr,(off_t)(nt-1-itO)*nslice,SEEK_SET);
		else             sf_seek(Fopr,(off_t)      itO *nslice,SEEK_SET);
		sf_floatread(opr[ lht ][0][0],nz*nx*ny,Fopr);

		itW = it + nht;
		if( !wflcausal ) sf_seek(Fwfl,(off_t)(nt-1-itW)*nslice,SEEK_SET);
		else             sf_seek(Fwfl,(off_t)      itW *nslice,SEEK_SET);
		sf_floatread(wfl[ lht ][0][0],nz*nx*ny,Fwfl);
	    }

	    if(it>=0+nht && 
	       it<nt-nht) { 

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)				\
    private(ic,     ihx,ihy,ihz,iht,mcx,   mcy,   mcz,   mct,   pcx,   pcy,   pcz,   pct) \
    shared (nc,ccin,ahx,ahy,ahz,aht,mcxall,mcyall,mczall,mctall,pcxall,pcyall,pczall,pctall)
#endif
		for(ic=0;ic<nc;ic++){ if(ccin[ic]) { /* sum over c and t! */
      if(pos){
			  EICLOOP( img[ic][iht][ihy][ihx][ihz] +=
				   opr    [pct][pcy][pcx][pcz] *
			  	 wfl    [mct][mcy][mcx][mcz]; );
      }else{
  			EICLOOP( img[ic][iht][ihy][ihx][ihz] +=
	  			 opr    [mct][mcy][mcx][mcz] *
		  		 wfl    [pct][pcy][pcx][pcz]; );


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

	if(gaus) applyGaussian(img,ac,aht,ahx,ahy,ahz,gst,gsx,gsy,gsz);                       /* Gaussian  */
	;        applyScaling (img,ac,aht,ahx,ahy,ahz,scale);                                 /* scaling   */
	sf_floatwrite(img[0][0][0][0],sf_n(ac)*sf_n(ahx)*sf_n(ahy)*sf_n(ahz)*sf_n(aht),Fimg); /* write img */

    } /* end if adj */

    /*------------------------------------------------------------*/
    /* deallocate arrays */
    free(****img); free(***img); free(**img); free(*img); free(img);
    ;              free(***opr); free(**opr); free(*opr); free(opr);
    ;              free(***wfl); free(**wfl); free(*wfl); free(wfl);
 
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


