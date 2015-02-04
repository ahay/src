/* Datuming by 3D Green functions in constant media */

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
    bool  fast; /* precompute TAA */
    bool nin; /*  input on the near array */
    bool nou; /* output on the near array */    
	
    sf_file  Fnn,     Fff; /* coordinates */
    pt3d    *nn=NULL,*ff=NULL;
    sf_axis  an,      af;
    int      jn,      jf;
    
    sf_file     Fwin,       Fwou; /* wavefield */
    sf_complex **win=NULL, **wnn=NULL, **wff=NULL;
    sf_axis  aw;
    int      jw;

    float slow; /* slowness */
    float dist, **tim=NULL,**amp=NULL,**ang=NULL; /* time & amplitude & angle */
    float   d,t,a;
    sf_complex iomega;

    float ox,oy,oz;
    pt3d oo;
    vc3d vecON, vecOF;
    float angFON,angMAX;
    float gauANG, g;

    int ompnth=1,ompith=0;

    size_t eseek;

    /*------------------------------------------------------------*/
    sf_init(argc,argv);
#ifdef _OPENMP
    ompnth=omp_init();
#endif

    if(! sf_getbool("verb",&verb)) verb=false; /* verbosity flag */
    if(!sf_getfloat("velo",&velo)) velo=1.0;   /* medium velocity */
    slow=1./velo; /* slowness */
    if(! sf_getbool("fast",&fast)) fast=true;  /* fast execution */

    if(! sf_getbool("nin",&nin)) nin=true;
    if(! sf_getbool("nou",&nou)) nou=true;
    
    /*------------------------------------------------------------*/
    if(! sf_getfloat("ox",&ox)) ox=0.0; oo.x=ox;
    if(! sf_getfloat("oy",&oy)) oy=0.0; oo.y=oy;
    if(! sf_getfloat("oz",&oz)) oz=0.0; oo.z=oz;
    if(! sf_getfloat("angMAX",&angMAX)) angMAX=90.0;
    if(! sf_getfloat("gauANG",&gauANG)) gauANG=0.3*angMAX;
    gauANG = 1./ (2.*pow(gauANG,2));

    /*------------------------------------------------------------*/
    /* coordinates */
    Fnn = sf_input ("cnn"); /* "near" coordinates  */
    Fff = sf_input ("cff"); /* "far" array */

    an = sf_iaxa(Fnn,2); sf_setlabel(an,"n"); if(verb) sf_raxa(an);
    af = sf_iaxa(Fff,2); sf_setlabel(af,"f"); if(verb) sf_raxa(af);

    nn = (pt3d*) sf_alloc(sf_n(an),sizeof(*nn)); 
    ff = (pt3d*) sf_alloc(sf_n(af),sizeof(*ff)); 

    memreport(sf_n(an)/1024./1024.*sizeof(*nn));
    memreport(sf_n(af)/1024./1024.*sizeof(*ff));

    pt3dread1(Fnn,nn,sf_n(an),3); /* read (x,y,z) */
    pt3dread1(Fff,ff,sf_n(af),3); /* read (x,y,z) */
      
    /*------------------------------------------------------------*/
    /* wavefield */
    Fwin = sf_input ("in" ); 
    Fwou = sf_output("out"); 
    aw = sf_iaxa(Fwin,2); sf_setlabel(aw,"w"); if(verb) sf_raxa(aw);  /* freq axis */

    if( nou ) sf_oaxa(Fwou,an,1);
    else      sf_oaxa(Fwou,af,1);
    
    sf_warning("allocate memory");

    if( nin ) {
	win = sf_complexalloc2(sf_n(an),ompnth);
	memreport(sf_n(an)/1024.*ompnth/1024.*sizeof(sf_complex));    
    }
    
    wff = sf_complexalloc2(sf_n(af),ompnth);
    memreport(sf_n(af)/1024.*ompnth/1024.*sizeof(sf_complex));
    
    if( nou ) {
	wnn = sf_complexalloc2(sf_n(an),ompnth);
	memreport(sf_n(an)/1024.*ompnth/1024.*sizeof(sf_complex));
    }
    sf_warning("OK");

    /*------------------------------------------------------------*/
    if(fast) {
	sf_warning("precompute T & A & A");

	/* precompute time */
	tim = sf_floatalloc2(sf_n(af),sf_n(an));
	memreport(sf_n(af)/1024.*sf_n(an)/1024*sizeof(float));

	/* precompute amplitude */
	amp = sf_floatalloc2(sf_n(af),sf_n(an));
	memreport(sf_n(af)/1024.*sf_n(an)/1024*sizeof(float));

	/* precompute angle */
	ang = sf_floatalloc2(sf_n(af),sf_n(an));
	memreport(sf_n(af)/1024.*sf_n(an)/1024*sizeof(float));
	
	for    (jn=0; jn<sf_n(an); jn++) { vecON = vec3d(&oo, &nn[jn]);
	    for(jf=0; jf<sf_n(af); jf++) { vecOF = vec3d(&oo, &ff[jf]);
		ang[jn][jf] = ang3d(&vecON, &vecOF);
		
		dist = sqrt( pow((ff[jf].x-nn[jn].x),2) + 
			     pow((ff[jf].y-nn[jn].y),2) + 
			     pow((ff[jf].z-nn[jn].z),2) );
		
		amp[jn][jf] = dist==0?1.0:1.0/(4*SF_PI*dist);  
		amp[jn][jf] *= exp( - pow(ang[jn][jf],2) * gauANG );

		tim[jn][jf] = dist*slow;
	    }
	}

	sf_warning("OK");
    }

    /*------------------------------------------------------------*/
    sf_warning("reserve output");

    if( nou ) {
	for(jn=0; jn<sf_n(an); jn++)
	    wnn[0][jn] = sf_cmplx(0.0,0.0);
	for(jw=0; jw<sf_n(aw); jw++)
	    sf_complexwrite(wnn[0],sf_n(an),Fwou);
    } else {
	for(jf=0; jf<sf_n(af); jf++)
	    wff[0][jf] = sf_cmplx(0.0,0.0);
	for(jw=0; jw<sf_n(aw); jw++)
	    sf_complexwrite(wff[0],sf_n(af),Fwou);
    }
    sf_seek(Fwou,0,SEEK_SET);
    sf_warning("OK");

    /*------------------------------------------------------------*/
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic)				\
    private(ompith,jw,jn,jf,iomega,d,t,a,g,vecON,vecOF,angFON)		\
    shared (       aw,an,af,win,wnn,wff,tim,amp,ang,oo,nn,ff,nin,nou)
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
	    if( nin ) {
		eseek =  jw*sf_n(an); /* seek elements */
		sf_seek(Fwin,(off_t)(eseek*sizeof(sf_complex)),SEEK_SET);
		sf_complexread (win[ompith],sf_n(an),Fwin);
	    } else {
		eseek =  jw*sf_n(af); /* seek elements */
		sf_seek(Fwin,(off_t)(eseek*sizeof(sf_complex)),SEEK_SET);
		sf_complexread (wff[ompith],sf_n(af),Fwin); 
	    }
	}
	
	iomega = sf_cmplx(0.0, 2.*SF_PI* (sf_o(aw) + jw*sf_d(aw)));

	if(fast) {

	    if( nin ) {
		
		/* spray nn -> ff */
		for(jf=0; jf<sf_n(af); jf++) wff[ompith][jf] = sf_cmplx(0.0,0.0);
		for    (jn=0; jn<sf_n(an); jn++) {
		    for(jf=0; jf<sf_n(af); jf++) {
			if( ang[jn][jf] < angMAX )
			    wff[ompith][jf] += win[ompith][jn] * amp[jn][jf]*cexpf(-iomega*tim[jn][jf]);
		    }
		}

	    }

	    if( nou ) {

		/* stack nn <- ff */
		for(jn=0; jn<sf_n(an); jn++) wnn[ompith][jn] = sf_cmplx(0.0,0.0);
		for    (jn=0; jn<sf_n(an); jn++) {
		    for(jf=0; jf<sf_n(af); jf++) {
			if( ang[jn][jf] < angMAX )
			    wnn[ompith][jn] += wff[ompith][jf] * amp[jn][jf]*cexpf(+iomega*tim[jn][jf]);
		    }
		}
		
	    }

	} else {

	    if( nin ) {
		
		/* spray nn -> ff */
		for(jf=0; jf<sf_n(af); jf++)
		    wff[ompith][jf] = sf_cmplx(0.0,0.0);
		for(jn=0; jn<sf_n(an); jn++) {
		    
		    vecON = vec3d(&oo, &nn[jn]);         /* vector O-N */
		    for(jf=0; jf<sf_n(af); jf++) {
			vecOF = vec3d(&oo, &ff[jf]);     /* vector O-F */
			angFON = ang3d(&vecON, &vecOF);  /* angle F-O-N */
			
			if( angFON < angMAX ) {
			    g = exp( - pow(angFON,2) * gauANG ); /* gaussian scaling */
			    d = sqrt( pow((ff[jf].x-nn[jn].x),2) +
				      pow((ff[jf].y-nn[jn].y),2) +
				      pow((ff[jf].z-nn[jn].z),2) );
			    a = d==0?1.0:1.0/(2*SF_PI*d);     /* spherical divergence */
			    t = d*slow;                       /* traveltime */
			    wff[ompith][jf] += win[ompith][jn] * g * a*cexpf(-iomega*t);
			}
		    }
		}

	    }

	    if( nou ) {

		/* stack nn <- ff */
		for(jn=0; jn<sf_n(an); jn++)        
		    wnn[ompith][jn] = sf_cmplx(0.0,0.0);
		
		for(jn=0; jn<sf_n(an); jn++) {
		    
		    vecON = vec3d(&oo, &nn[jn]);          /* vector O-N */
		    for(jf=0; jf<sf_n(af); jf++) {
			vecOF = vec3d(&oo, &ff[jf]);      /* vector O-F */
			angFON = ang3d(&vecON, &vecOF);  /* angle F-O-N */
			
			if( angFON < angMAX ) {
			    g = exp( - pow(angFON,2) * gauANG ); /* gaussian scaling */
			    d = sqrt( pow((ff[jf].x-nn[jn].x),2) +
				      pow((ff[jf].y-nn[jn].y),2) +
				      pow((ff[jf].z-nn[jn].z),2) );
			    a = d==0?1.0:1.0/(2*SF_PI*d);        /* spherical divergence */
			    t = d*slow;                          /* traveltime */
			    wnn[ompith][jn] += wff[ompith][jf] * g * a*cexpf(+iomega*t);
			}
		    }
		}

	    }
	    
	}
	
#ifdef _OPENMP
#pragma omp critical
#endif
	{/* write wou */
	    if( nou ) {
		eseek = jw*sf_n(an);
		sf_seek(Fwou,(off_t)(eseek*sizeof(sf_complex)),SEEK_SET);
		sf_complexwrite(wnn[ompith],sf_n(an),Fwou); 
	    } else {
		eseek = jw*sf_n(af);
		sf_seek(Fwou,(off_t)(eseek*sizeof(sf_complex)),SEEK_SET);
		sf_complexwrite(wff[ompith],sf_n(af),Fwou);
	    }
	}
    }

    /*------------------------------------------------------------*/
    /* deallocate arrays */
    
    if(fast) {
	free(*ang); free(ang); memreport(-sf_n(af)/1024.*sf_n(an)/1024.*sizeof(float));
	free(*amp); free(amp); memreport(-sf_n(af)/1024.*sf_n(an)/1024.*sizeof(float));
	free(*tim); free(tim); memreport(-sf_n(af)/1024.*sf_n(an)/1024.*sizeof(float));
    }
    if( nou ) {
	free(*wnn); free(wnn);
	memreport(-ompnth/1024.*sf_n(an)/1024.*sizeof(sf_complex));
    }
    
    free(*wff); free(wff);
    memreport(-ompnth/1024.*sf_n(af)/1024.*sizeof(sf_complex));

    if( nin ) {
	free(*win); free(win);
	memreport(-ompnth/1024.*sf_n(an)/1024.*sizeof(sf_complex));
    }
    
    free(ff);              memreport(-sf_n(af)/1024./1024.*sizeof(*ff));
    free(nn);              memreport(-sf_n(an)/1024./1024.*sizeof(*nn));

    exit (0);
}
