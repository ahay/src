/* out-of-core helpers */
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

#include <stdio.h>
/*^*/

#include <rsf.h>
/*^*/

#include "oc.h"

static char buf[BUFSIZ], buf2[BUFSIZ];

void oc_invert(size_t n, FILE *wall)
/*< write out 1/wall >*/
{
    size_t i, nleft;
    float *fbuf;

    fbuf = (float *) buf;
    if (0 != fseeko(wall,0,SEEK_SET))
	sf_error("%s: seeking error:",__FILE__);
    for (nleft = BUFSIZ; n > 0; n -= nleft) {
	if (nleft > n) nleft=n;
	if (nleft != fread (buf,1,nleft,wall))
	    sf_error("%s: reading error:",__FILE__);
	for (i=0; i < nleft/sizeof(float); i++) {
	    if (fbuf[i] > 0.) fbuf[i] = 1./fbuf[i];
	}
	if (nleft != fwrite (buf,1,nleft,wall))
	    sf_error("%s: reading error:",__FILE__);
    }
} 

void oc_zero (size_t n, FILE *wall) 
/*< set wall=0 >*/
{
    size_t nleft;

    memset(buf,0,BUFSIZ);
    if (0 != fseeko(wall,0,SEEK_SET))
	sf_error("%s: seeking error:",__FILE__);
    for (nleft = BUFSIZ; n > 0; n -= nleft) {
	if (nleft > n) nleft=n;
	if (nleft != fwrite (buf,1,nleft,wall))
	    sf_error("%s: writing error:",__FILE__);
    }
}

void oc_dump (size_t n, FILE *wall, sf_file out) 
/*< write wall to out >*/
{
    size_t nleft;
    float *fbuf;
    
    fbuf = (float *) buf;

    if (0 != fseeko(wall,0,SEEK_SET))
	sf_error("%s: seeking error:",__FILE__);
    for (nleft = BUFSIZ; n > 0; n -= nleft) {
	if (nleft > n) nleft=n;
	if (nleft != fread (buf,1,nleft,wall))
	    sf_error("%s: reading error:",__FILE__);
	sf_floatwrite(fbuf,nleft/sizeof(float),out);
    }
}

void oc_divide (size_t n, FILE *data, FILE *wall, sf_file out) 
/*< write data/wall to out >*/
{
    size_t nleft, i, nfloat;
    float *fbuf, *fbuf2;
    
    fbuf = (float *) buf;
    fbuf2 = (float *) buf2;

    if (0 != fseeko(wall,0,SEEK_SET))
	sf_error("%s: seeking error:",__FILE__);
    if (0 != fseeko(data,0,SEEK_SET))
	sf_error("%s: seeking error:",__FILE__);

    for (nleft = BUFSIZ; n > 0; n -= nleft) {
	if (nleft > n) nleft=n;
	if (nleft != fread (buf,1,nleft,wall))
	    sf_error("%s: reading error:",__FILE__);
	if (nleft != fread (buf2,1,nleft,data))
	    sf_error("%s: reading error:",__FILE__);
	nfloat = nleft/sizeof(float);
	for (i=0; i < nfloat; i++) {
	    if (fbuf[i] > 0.) fbuf2[i] /= fbuf[i];
	}
	sf_floatwrite(fbuf2,nfloat,out);
    }
}
