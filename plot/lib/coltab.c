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

#include <rsf.h>

#include "coltab.h"

static void hue2rgb (float hue,float *red, float *green, float *blue);
/* convert hue to RGB triple */

void vp_name2coltab (const char *colname /* color table name */, 
		     int nocol           /* number of colors */,
		     float *red, 
		     float *green, 
		     float *blue         /* RGB values */)
/*< Create a color table. nocol is between 2 and 256. >*/
{
    int i, ic, c;
    float h, redsave, gray, hnocol;
 
    c = *colname;
    if ('F'==c) c='e'; /* not to confuse 'F' with 'f' */

    switch (tolower(c)) {
	case 'j':
	case 't':
	    hnocol = 0.25  * nocol;
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
    
    for (i = 0; i < nocol; i++) {
	if (isupper(c)) {
	    ic = nocol - 1 - i; /* reverse color table */
	} else {
	    ic = i;
	}
	gray = ((float) ic)/(nocol - 1.);
	switch (tolower(c)) {
	    case 'a': /* rainbow - HSV */
		hue2rgb (gray, &red[ic], &green[ic], &blue[ic]);
		break;
	    case 'h': case 'p': case 'b': /* hot, pink, bone */
		if (i < hnocol) {
		    red[ic] = (i + 1.) / hnocol;
		    green[ic] = 0.;
		    blue[ic] = 0.;
		} else if (i < 2*hnocol) { /* i >= hnocol */
		    h = i - hnocol;
		    red[ic] = 1.;
		    green[ic] = (h + 1.) / hnocol;
		    blue[ic] = 0.;
		} else {                   /* i >= 2*hnocol */
		    h = i - 2*hnocol;
		    red[ic] = 1.;
		    green[ic] = 1.;
		    blue[ic] = (h + 1.) / (nocol - 2*hnocol);
		}
		if (c == 'p') { /* pink */
		    red[ic]   = sqrtf((2.*gray + red[ic])/3.);
		    green[ic] = sqrtf((2.*gray + green[ic])/3.);
		    blue[ic]  = sqrtf((2.*gray + blue[ic])/3.);
		} else if (c == 'b') { /* bone */
		    redsave = red[ic];
		    red[ic] = (7.*gray + blue[ic])/8.;
		    green[ic] = (7.*gray + green[ic])/8.;
		    blue[ic] = (7.*gray + redsave)/8.;
		}	  
		break;
	    case 'c': /* cool */
		red[ic] = gray;
		green[ic] = 1. - gray;
		blue[ic] = 1.;
		break;
	    case 'l': /* linear = COPPER */
		redsave = 1.25 * gray;
		if (redsave < 1.) redsave = 1.;
		red[ic] = redsave;
		green[ic] = 0.7812 * gray;
		blue[ic] = 0.4975 * gray;
		break;
	    case 'e': /* blue-white-red */
		if (i < hnocol) {
		    red[ic] = 1.;
		    green[ic] = blue[ic] = i/hnocol;
		} else {
		    red[ic] = green[ic] = (nocol-1-i)/hnocol;
		    blue[ic] = 1.;
		}
		break;
	    case 'f': /* flag */
		switch (i%4) {
		    case 0: /* red */
			red[ic] = 1.; 
			green[ic] = 0.; 
			blue[ic] = 0.; 
			break;
		    case 1: /* white */
			red[ic] = 1.; 
			green[ic] = 1.; 
			blue[ic] = 1.; 
			break;
		    case 2: /* blue */
			red[ic] = 0.; 
			green[ic] = 0.; 
			blue[ic] = 1.; 
			break;
		    case 3: /* black */
			red[ic] = 0.; 
			green[ic] = 0.; 
			blue[ic] = 0.; 
			break;
		    default:
			break;
		}
		break;
	    case 'g': /* black-white-red */
		if (i < hnocol) {
		    red[ic] = 1.;
		    green[ic] = blue[ic] = i/hnocol;
		} else {
		    red[ic] = green[ic] = blue[ic] = (nocol-1-i)/hnocol;
		}
		break;
	    case 'j': /* jet */
		if (i <= hnocol/2) {
		    red[ic] = 0.;
		    green[ic] = 0.;
		    blue[ic] = (i + hnocol/2)/hnocol;
		} else if (i < 3*hnocol/2) {  /* i > hnocol/2 */
		    h = i - hnocol/2;
		    red[ic] = 0.;
		    green[ic] = (h + 1.)/hnocol;
		    blue[ic] = 1.;
		} else if (i < 5*hnocol/2) { /* i >= 3*hnocol/2 */
		    h = i - 3*hnocol/2;
		    red[ic] = (h + 1.)/hnocol;
		    green[ic] = 1.;
		    blue[ic] = 1. - red[ic];
		} else if (i < 7*hnocol/2) { /* i >= 5*hnocol/2 */
		    h = i - 5*hnocol/2;
		    red[ic] = 1;
		    green[ic] = 1.- (h + 1.)/hnocol;
		    blue[ic] = 0.;	
		} else {                     /* i >= 7*hnocol/2 */
		    h = i - 7*hnocol/2;
		    red[ic] = (hnocol - h)/hnocol;
		    green[ic] = 0.;
		    blue[ic] = 0.;
		}
		break;
	    case 't': /* traffic */
		blue[ic] = 0.;
		if (i <= hnocol/2) { /* green up */
		    red[ic] = 0.; 
		    green[ic] = (i + hnocol/2)/hnocol; 
		} else if (i < 3*hnocol/2) { /* i > hnocol/2 */
		    h = i - hnocol/2;
		    red[ic] = (h + 1.) /hnocol;  /* red up */
		    green[ic] = 1.;
		} else if (i < 5*hnocol/2) { /* steady yellow */
		    red[ic] = 1.;
		    green[ic] = 1.;
		} else if (i < 7*hnocol/2) { /* green down */
		    h = i - 5*hnocol/2;;
		    red[ic] = 1;
		    green[ic] = 1.- (h + 1.)/hnocol;
		    blue[ic] = 0.;	
		} else {                     /* i >= 7*hnocol/2 */
		    h = i - 7*hnocol/2;   /* red down */
		    red[ic] = (hnocol - h)/hnocol; 
		    green[ic] = 0.;
		    blue[ic] = 0.;
		}
		break;
	    default: /* grayscale */
		red[ic] = gray;
		green[ic] = gray;
		blue[ic] = gray;
		break;
	}

	/* sanity check */
	if (red[ic] < 0.) red[ic]=0.;
	if (red[ic] > 1.) red[ic]=1.;
	if (green[ic] < 0.) green[ic]=0.;
	if (green[ic] > 1.) green[ic]=1.;
	if (blue[ic] < 0.) blue[ic]=0.;
	if (blue[ic] > 1.) blue[ic]=1.;
    }		
    
    /* If clipping flagged, change the 2 values at either end to red */
    if (colname[0] != '\0' && colname[1] == 'C') {
	red[0] = red[1] = red[nocol-2] = red[nocol-1] = 1.;
	green[0] = green[1] = green[nocol-2] = green[nocol-1] = 0.;
	blue[0] = blue[1] = blue[nocol-2] = blue[nocol-1] = 0.;
    }
}

static void hue2rgb (float hue,float *red, float *green, float *blue)
/* convert hue to RGB triple */
{
    float df, dg;
    int i;
    
    hue = hue*360.0; 
    if (hue >= 360.0) hue = 0.0; 
    else hue /= 60.0;
    
    i = hue;
    df = hue-i;
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
