/* Generate raster plot.
   
Takes:  > (plot.vpl | char.rsf)

Can input char values.
If called "byte", outputs char values.

If called "bar", outputs scalebar data.

Run "sfdoc stdplot" for more parameters.
*/
/*
  Copyright (C) 2004 University of Texas at Austin
  
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
#include <float.h>
#include <string.h>

#include <rsf.h>
#include <rsfplot.h>

#define TSIZE 4096

int main(int argc, char* argv[])
{
    int n1, n2, n3, gainstep, panel, it, nreserve, i1, i2, i3, j, orient;
    float o1, o2, o3, d1, d2, d3, gpow, clip, pclip, phalf, bias=0., minmax[2];
    float pbias, gain=0., x1, y1, x2, y2, **data=NULL, f, barmin, barmax, dat;
    bool transp, yreverse, xreverse, allpos, polarity, symcp, verb;
    bool eclip=false, egpow=false, barreverse, mean=false;
    bool scalebar, nomin=true, nomax=true, framenum, sfbyte, sfbar, charin;
    char *gainpanel, *color, *barfile;
    unsigned char tbl[TSIZE+1], **buf, tmp, *barbuf[1];
    enum {GAIN_EACH=-3,GAIN_ALL=-2,NO_GAIN=-1};
    off_t pos;
    sf_file in, out=NULL, bar=NULL;
    
    sf_init(argc,argv);
    in = sf_input("in");

    sfbyte = (bool) (NULL != strstr (sf_getprog(),"byte"));
    sfbar = (bool) (NULL != strstr (sf_getprog(),"bar"));

    if (sfbyte) {
	out = sf_output("out");
	sf_settype(out,SF_UCHAR);
    } else if (sfbar) {
	bar = sf_output("out");
	sf_settype(bar,SF_UCHAR);
    } else {
	vp_init();
    }

    charin = (bool) (SF_UCHAR == sf_gettype(in));

    if (charin && sfbyte) sf_error("Cannot input uchar to byte");

    if (!charin && SF_FLOAT != sf_gettype(in)) sf_error("Need float input");

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

    if (transp) {
	orient = 3;
    } else {	
	orient = (xreverse==yreverse)? 0:2;
    }

    if (!charin) {
	panel = NO_GAIN; /* no need for gain */
	
	phalf=85.;
	egpow = false;
	if (!sf_getfloat("gpow",&gpow)) {
	    gpow=1.;
	    /*( gpow=1 raise data to gpow power for display )*/
	} else if (gpow <= 0.) {
	    gpow=0.;
	    egpow = true;
	    sf_getfloat("phalf",&phalf);
	    /* percentage for estimating gpow */
	    if (phalf <=0. || phalf > 100.)
		sf_error("phalf=%g should be > 0 and <= 100",phalf);
	    panel = 0;
	}
	
	pclip=99.;
	eclip = (bool) (!sf_getfloat("clip",&clip));
	/* data clip */
	if (eclip) {	    
	    clip = 0.;
	    sf_getfloat("pclip",&pclip);
	    /* data clip percentile (default is 99) */
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
			    sf_error("gainpanel= should be all,"
				     " each, or a number"
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
	if (!sf_getbool("mean",&mean)) mean=false;
	/* if y, bias on the mean value */
	if (!sf_getfloat("bias",&pbias)) pbias=0.;
	/* value mapped to the center of the color table */
	if (!sf_getbool("polarity",&polarity)) polarity=false;
	/* if y, reverse polarity (white is high by default) */
	if (!sf_getbool("symcp",&symcp)) symcp=false;
	/* if y, assume symmetric color palette of 255 colors */
	if (!sf_getbool("verb",&verb)) verb=false;
	/* verbosity flag */
    } /* if !charin */

    barfile = sf_getstring("bar");
    /* file for scalebar data */

    if (sfbyte) {
	scalebar = (bool) (NULL != barfile);
	if (scalebar) sf_putstring(out,"bar",barfile);
    } else if (sfbar) {
	scalebar = true;
    } else {
	if (!sf_getbool ("wantscalebar",&scalebar) && 
	    !sf_getbool ("scalebar",&scalebar)) scalebar = false;
	/* if y, draw scalebar */	
    }
    if (scalebar) {
	nomin = (bool) (!sf_getfloat("minval",&barmin));
	/* minimum value for scalebar (default is the data minimum) */
	nomax = (bool) (!sf_getfloat("maxval",&barmax));
	/* maximum value for scalebar (default is the data maximum) */
	
	barbuf[0] = sf_ucharalloc(VP_BSIZE);

	if (!sf_getbool("barreverse",&barreverse)) barreverse=false;
	/* if y, go from small to large on the bar scale */

	if (sfbyte || sfbar) {
	    if (sfbyte) {
		bar = sf_output("bar");
		sf_settype(bar,SF_UCHAR);
	    }
	    sf_putint(bar,"n1",VP_BSIZE+2*sizeof(float));
	    sf_putint(bar,"n2",1);
	    sf_putint(bar,"n3",n3);

	    if (!nomin) sf_putfloat(bar,"minval",barmin);
	    if (!nomax) sf_putfloat(bar,"maxval",barmax);
	} else if (charin) {
	    if (NULL == barfile) {
		barfile=sf_histstring(in,"bar");
		if (NULL == barfile) sf_error("Need bar=");
	    }

	    bar = sf_input(barfile);
	    if (SF_UCHAR != sf_gettype(bar)) sf_error("Need uchar in bar");

	    if (nomin) nomin = (bool) (!sf_histfloat(bar,"minval",&barmin));
	    if (nomax) nomax = (bool) (!sf_histfloat(bar,"maxval",&barmax));
	}
    }

    if (!sf_getbool("wantframenum",&framenum)) framenum = (bool) (n3 > 1);
    /* if y, display third axis position in the corner */

    x1 = o1-0.5*d1;
    x2 = o1+(n1-1)*d1+0.5*d1;
    y1 = o2-0.5*d2;
    y2 = o2+(n2-1)*d2+0.5*d2;

    if (!sfbyte && !sfbar) {
	vp_stdplot_init (x1, x2, y1, y2, transp, false, yreverse, false);
	vp_frame_init(in,"tlb",false);
/*	if (scalebar && !nomin && !nomax) 
	vp_barframe_init (in,barmin,barmax); */
    }

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

    buf = sf_ucharalloc2(n1,n2);

    if (!charin) {
	data = sf_floatalloc2(n1,n2);

	if (GAIN_ALL==panel || panel >= 0) {
	    pos = sf_tell(in);
	    if (panel > 0) sf_seek(in,
				   pos+panel*n1*n2*sizeof(float),
				   SEEK_SET);
	    vp_gainpar (in,data,n1,n2,gainstep,
			pclip,phalf,&clip,&gpow,mean,&pbias,
			n3,panel,panel);
	    if (verb) sf_warning("panel=%d bias=%g clip=%g gpow=%g",
				 panel,pbias,clip,gpow);
	    if (sfbyte) sf_putfloat(out,"clip",clip);
	    sf_seek(in,pos,SEEK_SET); /* rewind */
	}
    }

    if (!sfbyte && !sfbar) {
	/* initialize color table */
	if (NULL == (color = sf_getstring("color"))) color="i";
	/* color scheme (default is i) */
	if (!sf_getint ("nreserve",&nreserve)) nreserve = 8;
	/* reserved colors */
	vp_rascoltab (nreserve, color);
    }

    for (i3=0; i3 < n3; i3++) {	
	if (!charin) {
	    if (GAIN_EACH == panel) {
		if (eclip) clip=0.;
		if (egpow) gpow=0.;
		vp_gainpar (in,data,n1,n2,gainstep,
			    pclip,phalf,&clip,&gpow,
			    mean,&pbias,n3,0,n3);
		if (verb) sf_warning("bias=%g clip=%g gpow=%g",pbias,clip,gpow);
	    } else {
		sf_floatread(data[0],n1*n2,in);
	    }
	    
	    if (1 == panel || GAIN_EACH == panel || 0==i3) { 
		/* initialize the conversion table */
		if(!allpos) { /* negative and positive values */
		    for (it=1; it<=TSIZE/2; it++) {
		        if (symcp) {
			    tbl[TSIZE-it] = (gpow != 1.)?
			        254*(pow(((TSIZE-2.0*it)/TSIZE),gpow)+1.)/2.+1.:
			        254*(    ((TSIZE-2.0*it)/TSIZE)      +1.)/2.+1.;
			    tbl[it] = 255 - tbl[TSIZE-it] + 1.0;
			} else {
			    tbl[TSIZE-it] = (gpow != 1.)?
			        252*(pow(((TSIZE-2.0*it)/TSIZE),gpow)+1.)/2.+3.:
			        252*(    ((TSIZE-2.0*it)/TSIZE)      +1.)/2.+3.;
			    tbl[it] = 255 - tbl[TSIZE-it] + 2.0;
			}
		    }
		    bias = TSIZE/2.;
		    gain = TSIZE/(2.*clip);
		} else { /* all positive */
		    if (symcp) {
			for (it=1; it < TSIZE ; it++) {
			    tbl[it] = 255*((it-1.0)/TSIZE) + 1.0;
			}
		    } else {
			for (it=1; it < TSIZE ; it++) {
			    tbl[it] = 256*((it-1.0)/TSIZE);
			}
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
	} else {
	    sf_ucharread(buf[0],n1*n2,in);
	}

	if (!sfbyte && !sfbar) {
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
	
	    if (i3 > 0) vp_erase (); 	

	    if (framenum) vp_framenum(o3+i3*d3);
	    vp_frame(); 
	    vp_uraster (buf, false, 256, n1, n2, 
			x1, y1, x2, y2, orient);
	    vp_simpleframe();
	}
	
	if (scalebar) {
	    if (!charin) {
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
		}
		
		for (it=0; it < VP_BSIZE; it++) {
		    if (barreverse) {
			dat = (barmin*it + barmax*(VP_BSIZE-1-it))/(VP_BSIZE-1);
		    } else {
			dat = (barmax*it + barmin*(VP_BSIZE-1-it))/(VP_BSIZE-1);
		    }
		    j = (dat-pbias)*gain + bias;
		    if      (j < 0) j=0;
		    else if (j > TSIZE) j=TSIZE;
		    barbuf[0][it] = tbl[j];
		} 
	    } else {
		sf_floatread(minmax,2,bar);
		sf_ucharread(barbuf[0],VP_BSIZE,bar);

		if (nomin) barmin=minmax[0];
		if (nomax) barmax=minmax[1];
	    }

	    if (sfbyte || sfbar) {
		sf_floatwrite(&barmin,1,bar);
		sf_floatwrite(&barmax,1,bar);
		sf_ucharwrite(barbuf[0],VP_BSIZE,bar);
	    } else {
		if (barreverse) {
		    vp_barframe_init (in,barmax,barmin);
		} else {
		    vp_barframe_init (in,barmin,barmax);
		}
		vp_barraster(VP_BSIZE, barbuf);
	    }
	} /* if scalebar */

	if (sfbyte) {
	    sf_ucharwrite(buf[0],n1*n2,out);
	} else if (!sfbar) {
	    vp_purge();
	} 
    } /* i3 loop */


    exit (0);
}
