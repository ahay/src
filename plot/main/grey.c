/* Generate raster plot.

Takes: < data.rsf > plot.vpl
*/

#include <math.h>
#include <float.h>

#include <rsf.h>
#include <rsfplot.h>

#define TSIZE 4096
#define BSIZE 256

int main(int argc, char* argv[])
{
    int n1, n2, n3, gainstep, panel, it, nreserve, i1, i2, i3, j, orient;
    float o1, o2, o3, d1, d2, d3, gpow, clip, pclip, phalf, bias=0.;
    float pbias, gain=0., x1, y1, x2, y2, **data, f, barmin, barmax, dat;
    bool transp, yreverse, xreverse, allpos, polarity, verb;
    bool scalebar, nomin=true, nomax=true, framenum;
    char *gainpanel, *color;
    unsigned char tbl[TSIZE+1], **buf, tmp, *barbuf[1];
    enum {GAIN_EACH=-3,GAIN_ALL=-2,NO_GAIN=-1};
    sf_file in;
    
    sf_init(argc,argv);
    in = sf_input("in");
    vp_init();

    if (SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

    if (!sf_histint(in,"n1",&n1)) sf_error("No n1= in input");
    if (!sf_histint(in,"n2",&n2)) n2=1;
    n3 = sf_leftsize(in,2);

    if (!sf_histfloat(in,"o1",&o1)) o1=0.;
    if (!sf_histfloat(in,"o2",&o2)) o2=0.;
    if (!sf_histfloat(in,"o3",&o3)) o3=0.;

    if (!sf_histfloat(in,"d1",&d1)) d1=1.;
    if (!sf_histfloat(in,"d2",&d2)) d2=1.;
    if (!sf_histfloat(in,"d3",&d3)) d3=1.;
    
    if (!sf_getbool("transp",&transp)) transp=true;
    /* if y, transpose the display axes */
    if (!sf_getbool("yreverse",&yreverse)) yreverse=true;
    /* if y, reverse the vertical axis */
    if (!sf_getbool("xreverse",&xreverse)) xreverse=false;
    /* if y, reverse the horizontal axis */

    /* I don't understand this */
    if (transp) {
	orient = 3;
    } else {	
	orient = (xreverse==yreverse)? 0:2;
    }

    panel = NO_GAIN; /* no need for gain */

    phalf=85.;
    if (!sf_getfloat("gpow",&gpow)) {
	gpow=1.;
	/* raise data to gpow power for display */
    } else if (gpow <= 0.) {
	gpow=0.;
	sf_getfloat("phalf",&phalf);
	/* percentage for estimating gpow */
	if (phalf <=0. || phalf > 100.)
	    sf_error("phalf=%g should be > 0 and <= 100",phalf);
	panel = 0;
    }

    pclip=99.;
    if (!sf_getfloat("clip",&clip)) {
	/* data clip */
	clip = 0.;
	sf_getfloat("pclip",&pclip);
	/* data clip percentile */
	if (pclip <=0. || pclip > 100.)
	    sf_error("pclip=%g should be > 0 and <= 100",pclip);
	panel = 0;
    } else if (clip <= 0.) {
	sf_warning("clip=%g <= 0",clip);
	clip = FLT_EPSILON;
    }

    if (0==panel) {
	if (!sf_getint("gainstep",&gainstep)) gainstep=0.5+n1/256.;
	/* subsampling for gpow and clip estimation */
	if (gainstep <= 0) gainstep=1;

	gainpanel = sf_getstring("gainpanel");
	/* gain reference: 'a' for all, 'e' for each, or number */
	if (NULL != gainpanel) {
	    switch (gainpanel[0]) {
		case 'a': 
		    panel=GAIN_ALL; 
		    break;
		case 'e': 
		    panel=GAIN_EACH;
		    break;
		default:
		    if (0 ==sscanf(gainpanel,"%d",&panel) || 
			panel < 1 || panel > n3) 
			sf_error("gainpanel= should be all, each, or a number"
				 " between 1 and %d",n3);
		    panel--;
		    break;
	    }
	    free (gainpanel);
	} 

	sf_unpipe(in,sf_filesize(in)*sizeof(float));
    } 

    if (!sf_getbool("allpos",&allpos)) allpos=false;
    /* if y, assume positive data */
    if (!sf_getfloat("bias",&pbias)) pbias=0.;
    /* subtract bias from data */
    if (!sf_getbool("polarity",&polarity)) polarity=false;
    /* if y, reverse polarity (white is high by default) */
    if (!sf_getbool("verb",&verb)) verb=false;
    /* verbosity flag */

    if (!sf_getbool ("wantscalebar",&scalebar)) scalebar = false;
    /* if y, draw scalebar */
    if (scalebar) {
	nomin = !sf_getfloat("minval",&barmin);
	/* minimum value for scalebar (default is the data minimum) */
	nomax = !sf_getfloat("maxval",&barmax);
	/* maximum value for scalebar (default is the data maximum) */
	barbuf[0] = (unsigned char*) sf_alloc(BSIZE,sizeof(unsigned char));
    }

    if (!sf_getbool("wantframenum",&framenum)) framenum = (n3 > 1);
    /* if y, display third axis position in the corner */

    x1 = o1-0.5*d1;
    x2 = o1+(n1-1)*d1+0.5*d1;
    y1 = o2-0.5*d2;
    y2 = o2+(n2-1)*d2+0.5*d2;

    vp_stdplot_init (x1, x2, y1, y2, transp, false, yreverse, false);
    vp_frame_init(in,"tlb");
    if (scalebar && !nomin && !nomax) vp_barframe_init (barmin,barmax);

    if (transp) {
	f=x1; x1=y1; y1=f;
	f=x2; x2=y2; y2=f;
    }

    if (yreverse) {
	f=y1; y1=y2; y2=f;
    }

    if (xreverse) {
	f=x1; x1=x2; x2=f;
    }

    data = sf_floatalloc2(n1,n2);
    buf = sf_ucharalloc2(n1,n2);

    if (GAIN_ALL==panel || panel >= 0) {
	gainpar (in,data,n1,n2,gainstep,
		 o1,pclip,phalf,&clip,&gpow,pbias,d1,n3,panel);
	if (verb) sf_warning("panel=%d clip=%g gpow=%g",panel,clip,gpow);
	sf_seek(in,0L,SEEK_SET); /* rewind */
    }

    /* initialize color table */
    if (NULL == (color = sf_getstring("color"))) color="I";
    /* color scheme */
    if (!sf_getint ("nreserve",&nreserve)) nreserve = 8;
    /* reserved colors */
    vp_rascoltab (nreserve, color);

    for (i3=0; i3 < n3; i3++) {	
	if (i3 > 0) vp_erase (); 	
	if (GAIN_EACH == panel) {
	    gainpar (in,data,n1,n2,gainstep,
		     o1,pclip,phalf,&clip,&gpow,pbias,d1,n3,0);
	    if (verb) sf_warning("clip=%g gpow=%g",clip,gpow);
	} else {
	    sf_read(data[0],sizeof(float),n1*n2,in);
	}
	
	if (1 == panel || 0==i3) { 
	    /* initialize the conversion table */
	    if(!allpos) { /* negative and positive values */
		for (it=1; it<=TSIZE/2; it++) {
		    tbl[TSIZE-it] = (gpow != 1.)?
			252*(pow(((TSIZE-2.0*it)/TSIZE),gpow)+1.)/2.+3.:
			252*(    ((TSIZE-2.0*it)/TSIZE)      +1.)/2.+3.;
		    tbl[it] = 255 - tbl[TSIZE-it] + 2.0;
		}
		bias = TSIZE/2.;
		gain = TSIZE/(2.*clip);
	    } else { /* all positive */
		for (it=1; it < TSIZE ; it++) {
		    tbl[it] = 256*((it-1.0)/TSIZE);
		}
		bias = 0.;
		gain = TSIZE/clip;		
	    }
	    tbl[0] = tbl[1];
	    tbl[TSIZE] = tbl[TSIZE-1];
	    if (polarity) { /* switch polarity */
		for (it=0; it<=TSIZE; it++) {
		    tbl[it]=255-tbl[it];
		}
	    }
	}
	    
	/* convert to bytes */
	for (i2=0; i2 < n2; i2++) {
	    for (i1=0; i1 < n1; i1++) {
		j = (data[i2][i1]-pbias)*gain + bias;
		if      (j < 0) j=0;
		else if (j > TSIZE) j=TSIZE;
		buf[i2][i1] = tbl[j];
	    }
	}

	if (yreverse) {
	    for (i2=0; i2 < n2; i2++) {
		for (i1=0; i1 < n1/2; i1++) {			
		    tmp = buf[i2][i1];
		    buf[i2][i1] = buf[i2][n1-1-i1];
		    buf[i2][n1-1-i1] = tmp;
		}
	    }
	} 

	if ((xreverse && transp) || (!xreverse && !transp)) {
	    for (i2=0; i2 < n2/2; i2++) {
		for (i1=0; i1 < n1; i1++) {
		    tmp = buf[i2][i1];
		    buf[i2][i1] = buf[n2-1-i2][i1];
		    buf[n2-1-i2][i1] = tmp;
		}
	    }
	}

	if (framenum) vp_framenum(o3+i3*d3);
	vp_frame(); 
	vp_uraster (buf, false, 256, n1, n2, 
		    x1, y1, x2, y2, orient);
	vp_simpleframe();
	
	if (scalebar) {
	    if (nomin) barmin = data[0][0];
	    if (nomax) barmax = data[0][0];
	    if (nomin || nomax) {
		for (i2=0; i2 < n2; i2++) {
		    for (i1=0; i1 < n1; i1++) {
			dat = data[i2][i1];
			if (nomin && barmin > dat) barmin = dat;
			if (nomax && barmax < dat) barmax = dat;
		    }
		}
		vp_barframe_init (barmin,barmax);
	    }

	    for (it=0; it < BSIZE; it++) {
		dat = (barmax*it + barmin*(BSIZE-1-it))/(BSIZE-1);
		j = (dat-pbias)*gain + bias;
		if      (j < 0) j=0;
		else if (j > TSIZE) j=TSIZE;
		barbuf[0][it] = tbl[j];
	    }

	    vp_barraster(BSIZE, barbuf);
	}

	vp_purge(); 
    } 

    exit (0);
}

/* 	$Id: grey.c,v 1.14 2003/10/18 18:22:04 fomels Exp $	 */
