/* Color table. */
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
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include <rsf.h>

#include "coltab.h"

static void hue2rgb (float hue,float *red, float *green, float *blue);
/* convert hue to RGB triple */

void vp_name2coltab (   const char *colname,    /* color table name */
                        int nocol,              /* number of colors */
			/*@out@*/ float *red,             /* RGB values       */
			/*@out@*/ float *green, 
			/*@out@*/ float *blue)
/*< Create a color table. nocol is between 2 and 256. >*/
{
    bool clip;
    int   i, ic, c, reverse;
    float h, redsave, gray, hnocol, angle, amp, cosa, sina, cliprgb[3];
    char *rsfroot, filename[PATH_MAX];
    FILE *colfile;
    const float start=0.5, rots=-1.5, hue=1.0;

    /* Determine clipping and reverse options from color code */
    if (strlen(colname) > 0) c = colname[0];
    else                     c = 'i';

    if (strlen(colname) > 1 && colname[1] == 'C') clip = true;
    else                                          clip = false;
 
    if (isupper(c)) reverse = 1;
    else            reverse = 0;

    c = tolower(c);

    /* first, try reading it from a file */
    if (NULL == (rsfroot = getenv("RSFROOT"))) rsfroot="/usr";
    sprintf (filename, "%s/include/%s.csv", rsfroot, colname);
    if (NULL != (colfile = fopen(filename,"r"))) {
	for (i = 0; i < nocol; i++) {
	    if (reverse == 1) ic = nocol - 1 - i;   /* reverse color table */
	    else              ic = i;
   
	    if (3 != fscanf(colfile,"%g,%g,%g\n",&red[ic], &green[ic], &blue[ic])) 
		sf_error("Error reading \"%s\"",filename);
	}
	return;
    }
 
    /* Initialize hnocol */
 
    switch (c) 
    {
        case 'j':
        case 't':
            hnocol = 0.25  * nocol;
	break;
	case 'w':
	    hnocol = 0.125  * nocol;
            break;
        case 'h':
        case 'p':
        case 'b':       
            hnocol = 3./8. * nocol; 
	break;
        default:
            hnocol = 0.5  * nocol;
            break;
    }
    
    /* Fill color table */
    
    for (i = 0; i < nocol; i++) {

	if (reverse == 1) ic = nocol - 1 - i;   /* reverse color table */
	else              ic = i;
    
	gray = ((float) i)/(nocol - 1.); /* ranges from 0 to 1 */

	switch (c)
	{
	    case 'x': /* cubehelix */
		/* Green, D. A., 2011, A colour scheme for the display
		 * of astronomical intensity images: Bulletin of the
		 * Astronomical Society of India, v.39, 289-295. */
		angle=2*SF_PI*(start/3.0f+1.0f+rots*gray);
		cosa = cosf(angle);
		sina = sinf(angle);
		amp=hue*gray*(1.0f-gray)/2.0f;
		red[ic]   = gray+amp*(-0.14861*cosa+1.78277*sina);
		green[ic] = gray+amp*(-0.29227*cosa-0.90649*sina);
		blue[ic]  = gray+amp*(+1.97294*cosa);
		break;
	    case 'a': /* rainbow - HSV */
		hue2rgb (gray, &red[ic], &green[ic], &blue[ic]);
		break;

	    case 'h': case 'p': case 'b': /* hot, pink, bone */
		if (i < hnocol)
		{
		    red[ic]   = (i + 1.) / hnocol;
		    green[ic] = 0.;
		    blue[ic]  = 0.;
		} 
		else if (i < 2*hnocol) 
		{
		    h         = i - hnocol;
		    red[ic]   = 1.;
		    green[ic] = (h + 1.) / hnocol;
		    blue[ic]  = 0.;
		} 
		else 
		{
		    h         = i - 2*hnocol;
		    red[ic]   = 1.;
		    green[ic] = 1.;
		    blue[ic]  = (h + 1.) / (nocol - 2*hnocol);
		}
		if (c == 'p') /* pink */
		{
		    red[ic]   = sqrtf((2.*gray + red[ic])/3.);
		    green[ic] = sqrtf((2.*gray + green[ic])/3.);
		    blue[ic]  = sqrtf((2.*gray + blue[ic])/3.);
		}
		else if (c == 'b') /* bone */
		{
		    redsave   = red[ic];
		    red[ic]   = (7.*gray + blue[ic])/8.;
		    green[ic] = (7.*gray + green[ic])/8.;
		    blue[ic]  = (7.*gray + redsave)/8.;
		}     
		break;

	    case 'c': /* cool */
		red[ic]   = gray;
		green[ic] = 1. - gray;
		blue[ic]  = 1.;
		break;

	    case 'l': /* linear = COPPER */
		redsave = 1.25 * gray;
		if (redsave < 1.) redsave = 1.;
		red[ic]   = redsave;
		green[ic] = 0.7812 * gray;
		blue[ic]  = 0.4975 * gray;
		break;

	    case 'e': /* blue-white-red */
		if (i < hnocol)
		{
		    red[ic]   = 1.;
		    green[ic] = blue[ic] = i/hnocol;
		}
		else
		{
		    red[ic]  = green[ic] = (nocol-1-i)/hnocol;
		    blue[ic] = 1.;
		}
		break;

	    case 'f': /* flag */
		switch (i%4)
		{
		    case 0: /* red */
			red[ic]   = 1.; 
			green[ic] = 0.; 
			blue[ic]  = 0.; 
			break;
		    case 1: /* white */
			red[ic]   = 1.; 
			green[ic] = 1.; 
			blue[ic]  = 1.; 
			break;
		    case 2: /* blue */
			red[ic]   = 0.; 
			green[ic] = 0.; 
			blue[ic]  = 1.; 
			break;
		    case 3: /* black */
			red[ic]   = 0.; 
			green[ic] = 0.; 
			blue[ic]  = 0.; 
			break;
		    default:
			break;
		}
		break;

	    case 'g': /* black-white-red */
		if (i < hnocol)
		{
		    red[ic]   = 1.;
		    green[ic] = blue[ic] = i/hnocol;
		}
		else
		{
		    red[ic] = green[ic] = blue[ic] = (nocol-1-i)/hnocol;
		}
		break;

	    case 'j': /* jet */
		if (i <= hnocol/2)
		{
		    red[ic]   = 0.;
		    green[ic] = 0.;
		    blue[ic]  = (i + hnocol/2)/hnocol;
		}
		else if (i < 3*hnocol/2) 
		{
		    h = i - hnocol/2;
		    red[ic]   = 0.;
		    green[ic] = (h + 1.)/hnocol;
		    blue[ic]  = 1.;
		}
		else if (i < 5*hnocol/2)
		{
		    h = i - 3*hnocol/2;
		    red[ic]   = (h + 1.)/hnocol;
		    green[ic] = 1.;
		    blue[ic]  = 1. - red[ic];
		}
		else if (i < 7*hnocol/2)
		{
		    h = i - 5*hnocol/2;
		    red[ic]   = 1;
		    green[ic] = 1.- (h + 1.)/hnocol;
		    blue[ic]  = 0.;  
		}
		else
		{
		    h = i - 7*hnocol/2;
		    red[ic]   = (hnocol - h)/hnocol;
		    green[ic] = 0.;
		    blue[ic]  = 0.;
		}
		break;

	    case 't': /* traffic */
		blue[ic] = 0.;
		if (i <= hnocol/2)          /* green up */
		{
		    red[ic]   = 0.; 
		    green[ic] = (i + hnocol/2)/hnocol; 
		}
		else if (i < 3*hnocol/2)    /* red up */
		{
		    h = i - hnocol/2;
		    red[ic]   = (h + 1.) /hnocol;
		    green[ic] = 1.;
		}
		else if (i < 5*hnocol/2)    /* steady yellow */
		{
		    red[ic]   = 1.;
		    green[ic] = 1.;
		}
		else if (i < 7*hnocol/2)    /* green down */
		{
		    h = i - 5*hnocol/2;
		    red[ic]   = 1;
		    green[ic] = 1.- (h + 1.)/hnocol;
		}
		else                        /* red down */
		{
		    h = i - 7*hnocol/2;
		    red[ic]   = (hnocol - h)/hnocol; 
		    green[ic] = 0.;
		}
		break;

	    case 'w': /* wheel */
		if (i <= hnocol/2)          /* green up*/
		{
		    red[ic]   = 0.; 
		    green[ic] = (i + hnocol/2)/hnocol;
		    blue[ic] = 0.;
		}
		else if (i < 3*hnocol/2)    /* red up*/
		{
		    h = i - hnocol/2;
		    red[ic] = (h + 1.) /hnocol;
		    green[ic] = 1.;
		    blue[ic]   = 0.;
		}
		else if (i < 7*hnocol/2)    /* red steady, green down */
		{
		    h = i - 3*hnocol/2;
		    red[ic]   = 1.;
		    green[ic] = 1.- (h + 1.)/(hnocol*2.0);
		    blue[ic] = 0.;		
		}

		else if (i < 4*hnocol)  /* red down*/
		{
		    h = i - 7*hnocol/2;
		    red[ic]   = (hnocol - h)/hnocol; 
		    green[ic] = 0.;
		    blue[ic] = 0.;
		} 
   
		else if (i <= 9*hnocol/2)          /*red up*/
		{
		    h = i - 4*hnocol;
		    blue[ic] = 0.;
		    green[ic]   = 0.; 
		    red[ic] = (h + hnocol/2)/hnocol; 
		}
		else if (i < 11*hnocol/2)    /* steady red, blue up*/
		{
		    h = i - 9*hnocol/2;
		    blue[ic]   = (h + 1.) /hnocol;
		    green[ic]   = 0;
		    red[ic] = 1.;
		}
		else if (i < 13*hnocol/2)    /* blue steady, red down */
		{
		    h = i - 11*hnocol/2;
		    red[ic]   = (hnocol - h)/hnocol;
		    green[ic] = 0.;
		    blue[ic]  = 1.;
		}
		else if (i < 15*hnocol/2)    /* green up, blue down*/
		{
		    h = i - 13*hnocol/2;
		    red[ic] = 0.;
		    green[ic]   = (h + 1.) /hnocol;
		    blue[ic] = 1.- (h + 1.)/hnocol;
		}
		else  /* green down*/
		{
		    h = i - 15*hnocol/2;
		    green[ic]   = (hnocol - h)/hnocol; 
		    red[ic] = 0.;
		    blue[ic] = 0.;
		} 

		break;

	    default: /* grayscale */
		red[ic]   = gray;
		green[ic] = gray;
		blue[ic]  = gray;
		break;
	}

	/* sanity check */
        if (red[ic]   < 0.) red[ic]   = 0.;
        if (red[ic]   > 1.) red[ic]   = 1.;
        if (green[ic] < 0.) green[ic] = 0.;
        if (green[ic] > 1.) green[ic] = 1.;
        if (blue[ic]  < 0.) blue[ic]  = 0.;
        if (blue[ic]  > 1.) blue[ic]  = 1.;
    }       
    
    /* If clipping flagged, change the 2 values at either end */
    if (clip) 
    {
	if (!sf_getfloats("cliprgb",cliprgb,3)) {
	    /* default is red */
	    cliprgb[0]=1.;
	    cliprgb[1]=0.;
	    cliprgb[2]=0.;
	}
	
        red[0]       = cliprgb[0]; green[0]       = cliprgb[1]; blue[0]       = cliprgb[2];
        red[1]       = cliprgb[0]; green[1]       = cliprgb[1]; blue[1]       = cliprgb[2];
        red[nocol-2] = cliprgb[0]; green[nocol-2] = cliprgb[1]; blue[nocol-2] = cliprgb[2];
        red[nocol-1] = cliprgb[0]; green[nocol-1] = cliprgb[1]; blue[nocol-1] = cliprgb[2];
    }
}

static void hue2rgb (float hue,float *red, float *green, float *blue)
/* convert hue to RGB triple */
{
    float df, dg;
    int i;
    
    hue *= 360.0; 
    if (hue >= 360.0) hue  =  0.0; 
    else              hue /= 60.0;
    
    i  = hue;
    df = hue - i;
    dg = 1. - df;
    
    switch (i) {
	case 0: 
	    *red = 1.;  *green = df;  *blue = 0.; break;
	case 1: 
	    *red = dg;  *green = 1.;  *blue = 0.; break;
	case 2: 
	    *red = 0.;  *green = 1.;  *blue = df; break;
	case 3: 
	    *red = 0.;  *green = dg;  *blue = 1.; break;
	case 4: 
	    *red = df;  *green = 0.;  *blue = 1.; break;
	case 5: 
	    *red = 1.;  *green = 0.;  *blue = dg; break;
    }
}
