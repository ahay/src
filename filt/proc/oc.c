#include <rsf.h>

#include "oc.h"

static char buf[BUFSIZ];

void oc_invert(size_t n, FILE *wall)
{
    int i, nleft;
    float *fbuf;

    fbuf = (float *) buf;
    if (0 != fseek(wall,0L,SEEK_SET))
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
{
    size_t nleft;

    memset(buf,0,BUFSIZ);
    if (0 != fseek(wall,0L,SEEK_SET))
	sf_error("%s: seeking error:",__FILE__);
    for (nleft = BUFSIZ; n > 0; n -= nleft) {
	if (nleft > n) nleft=n;
	if (nleft != fwrite (buf,1,nleft,wall))
	    sf_error("%s: writing error:",__FILE__);
    }
}

void oc_dump (size_t n, FILE *wall, sf_file out) 
{
    size_t nleft;
    float *fbuf;
    
    fbuf = (float *) buf;

    if (0 != fseek(wall,0L,SEEK_SET))
	sf_error("%s: seeking error:",__FILE__);
    for (nleft = BUFSIZ; n > 0; n -= nleft) {
	if (nleft > n) nleft=n;
	if (nleft != fread (buf,1,nleft,wall))
	    sf_error("%s: reading error:",__FILE__);
	sf_floatwrite(fbuf,nleft/sizeof(float),out);
    }
}

