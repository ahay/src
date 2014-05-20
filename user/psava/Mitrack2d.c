/* Datuming by 2D Green functions in constant media */

#include <rsf.h>
#include <math.h>

#ifdef _OPENMP
#include <omp.h>
#endif

static float mymem=0;
void memreport(float addmem)
{
    mymem += addmem;
    sf_warning("MEMORY=%5d Mb",(int)mymem);
}

/*------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    bool  verb; /* verbosity flag */
    float velo; /* medium velocity */
    bool  fast;

    sf_file  Fnn,     Fff; /* coordinates */
    pt2d    *nn=NULL,*ff=NULL;
    sf_axis  an,      af;
    int      jn,      jf;
    
    sf_file     Fwin,       Fwou; /* wavefield */
    sf_complex **win=NULL, **wou=NULL, tmp;
    sf_axis  aw;
    int      jw;

    float slow; /* slowness */
    float dist, **tim=NULL,**amp=NULL;; /* time & amplitude */
    float   d,t,a;
    sf_complex iomega;

    int ompnth=1,ompith=0;

    float *gau,gxc,gxs; /* Gaussian,center and sigma */

    /*------------------------------------------------------------*/
    sf_init(argc,argv);
#ifdef _OPENMP
    ompnth=omp_init();
#endif

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(!sf_getfloat("velo",&velo)) velo=1.0;   /* medium velocity */
    slow=1./velo; /* slowness */
    if(! sf_getbool("fast",&fast)) fast=true;  /* fast execution */

    /*------------------------------------------------------------*/
    /* coordinates */
    Fnn = sf_input ("cnn"); /* "near" coordinates  */
    Fff = sf_input ("cff"); /* "far" array */

    an = sf_iaxa(Fnn,2); sf_setlabel(an,"n"); if(verb) sf_raxa(an);
    af = sf_iaxa(Fff,2); sf_setlabel(af,"f"); if(verb) sf_raxa(af);

    nn = (pt2d*) sf_alloc(sf_n(an),sizeof(*nn)); 
    ff = (pt2d*) sf_alloc(sf_n(af),sizeof(*ff)); 

    memreport(sf_n(an)/1024./1024.*sizeof(*nn));
    memreport(sf_n(af)/1024./1024.*sizeof(*ff));

    pt2dread1(Fnn,nn,sf_n(an),2); /* read (x,z) */
    pt2dread1(Fff,ff,sf_n(af),2); /* read (x,z) */
      
    /*------------------------------------------------------------*/
    /* wavefield */
    Fwin = sf_input ("in" ); 
    Fwou = sf_output("out"); 
    aw = sf_iaxa(Fwin,2); sf_setlabel(aw,"w"); if(verb) sf_raxa(aw);  /* freq axis */

    win = sf_complexalloc2(sf_n(an),ompnth);
    memreport(sf_n(an)/1024.*ompnth/1024.*sizeof(sf_complex));    

    wou = sf_complexalloc2(sf_n(an),ompnth);
    memreport(sf_n(an)/1024.*ompnth/1024.*sizeof(sf_complex));

    /*------------------------------------------------------------*/
    if(fast) {

	/* precompute time */
	tim = sf_floatalloc2(sf_n(an),sf_n(af));
	memreport(sf_n(af)/1024.*sf_n(an)/1024*sizeof(float));
	
	/* precompute amplitude */
	amp = sf_floatalloc2(sf_n(an),sf_n(af));
	memreport(sf_n(af)/1024.*sf_n(an)/1024*sizeof(float));
	
	for(jf=0; jf<sf_n(af); jf++) {
	    for(jn=0; jn<sf_n(an); jn++) {
		dist = sqrt( pow((ff[jf].x-nn[jn].x),2) + 
			     pow((ff[jf].z-nn[jn].z),2) ); 
		amp[jf][jn] = dist==0?1.0:1.0/(2*SF_PI*dist);  
		tim[jf][jn] = dist*slow;                   
	    }
	}
    }
    
    /*------------------------------------------------------------*/
    if(!sf_getfloat("gxc",&gxc)) gxc=0.;   /* Gaussian center x */
    if(!sf_getfloat("gxs",&gxs)) gxs=0.;   /* Gaussian stdev  x */
    gau = sf_floatalloc(sf_n(af));
    memreport(sf_n(af)/1024.*1./1024.*sizeof(float));   

    if(gxs==0.)
	for(jf=0; jf<sf_n(af); jf++)
	    gau[jf] = 1.0;
    else {
	gxs = 1./ (2.*pow(gxs,2));
	for(jf=0; jf<sf_n(af); jf++)
	    gau[jf] = exp( - pow(ff[jf].x-gxc,2) * gxs);
    }

    /*------------------------------------------------------------*/
    /* reserve output binary */
    for(jn=0; jn<sf_n(an); jn++)
	wou[0][jn] = sf_cmplx(0.0,0.0);
    for(jw=0; jw<sf_n(aw); jw++)
	sf_complexwrite(wou[0],sf_n(an),Fwou);
    sf_seek(Fwou,0,SEEK_SET);

    /*------------------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)				\
    private(ompith,jw,jn,jf,iomega,tmp,d,t,a)				\
    shared (       aw,an,af,win,wou,tim,amp,gau)
#endif
    for(jw=0; jw<sf_n(aw); jw++) {
#ifdef _OPENMP
        ompith = omp_get_thread_num();
#endif
#ifdef _OPENMP
#pragma omp critical
#endif
	{
	    if(verb) sf_warning ("(ith=%2d) ... <jw=%4d of %4d>",ompith,jw,sf_n(aw)-1);
	    sf_seek(Fwin,(off_t)(jw)*sf_n(an)*sizeof(sf_complex),SEEK_SET);
	    sf_complexread (win[ompith],sf_n(an),Fwin); /* read win */
	}
	for(jn=0; jn<sf_n(an); jn++)        /* init wou */
	    wou[ompith][jn] = sf_cmplx(0.0,0.0);

	iomega = sf_cmplx(0.0, 2.*SF_PI* (sf_o(aw) + jw*sf_d(aw)));

	if(fast) {

	    for(jf=0; jf<sf_n(af); jf++) {
		/* stack */
		tmp = sf_cmplx(0.0,0.0);
		for(jn=0; jn<sf_n(an); jn++)
		    tmp += win[ompith][jn] * amp[jf][jn]*cexpf(-iomega*tim[jf][jn]);
		
		tmp *= gau[jf];  /* Gaussian taper */
		
		/* spray */
		for(jn=0; jn<sf_n(an); jn++)
		    wou[ompith][jn] += tmp * amp[jf][jn]*cexpf(+iomega*tim[jf][jn]);
	    }
	} else {
	    	    for(jf=0; jf<sf_n(af); jf++) {

		/* stack */
		tmp = sf_cmplx(0.0,0.0);
		for(jn=0; jn<sf_n(an); jn++) {
		    d = sqrt( pow((ff[jf].x-nn[jn].x),2) + 
			      pow((ff[jf].z-nn[jn].z),2) ); 
		    a = d==0?1.0:1.0/(2*SF_PI*d);  
		    t = d*slow; 
		    tmp += win[ompith][jn] * a*cexpf(-iomega*t);
		}

		tmp *= gau[jf]; /* Gaussian taper */
		
		/* spray */
		for(jn=0; jn<sf_n(an); jn++) {
		    d = sqrt( pow((ff[jf].x-nn[jn].x),2) + 
			      pow((ff[jf].z-nn[jn].z),2) ); 
		    a = d==0?1.0:1.0/(2*SF_PI*d);  
		    t = d*slow; 
		    wou[ompith][jn] += tmp * a*cexpf(+iomega*t);
		}
	    }
	}

#ifdef _OPENMP
#pragma omp critical
#endif
	{
	    sf_seek(Fwou,(off_t)(jw)*sf_n(an)*sizeof(sf_complex),SEEK_SET);
	    sf_complexwrite(wou[ompith],sf_n(an),Fwou); /* write wou */
	}
    }

    /*------------------------------------------------------------*/
    /* deallocate arrays */

    free(gau);             memreport(sf_n(af)/1024./1024.*sizeof(float));    
    if(fast) {
	free(*amp); free(amp); memreport(-sf_n(af)/1024.*sf_n(an)/1024.*sizeof(float));
	free(*tim); free(tim); memreport(-sf_n(af)/1024.*sf_n(an)/1024.*sizeof(float));
    }
    free(*wou); free(wou); memreport(-  ompnth/1024.*sf_n(an)/1024.*sizeof(sf_complex));
    free(*win); free(win); memreport(-  ompnth/1024.*sf_n(an)/1024.*sizeof(sf_complex));
    free(ff);              memreport(-sf_n(af)/1024./1024.*sizeof(*ff));
    free(nn);              memreport(-sf_n(an)/1024./1024.*sizeof(*nn));

    exit (0);
}
