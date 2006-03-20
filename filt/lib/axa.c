#include <string.h>
#include <stdio.h>

#include "file.h"
/*^*/

#include "alloc.h"
#include "error.h"
#include "axa.h"


#ifndef _sf_axa_h

typedef struct {
    int n;
    float o,d;
} sf_axa;
/*^*/

typedef struct sf_Axis *sf_axis;
/*^*/

#endif

struct sf_Axis {
    int n;
    float o,d;
    char *l, *u;
};


/*------------------------------------------------------------*/
sf_axis sf_maxa(int n   /* length */, 
		float o /* origin */, 
		float d /* sampling */)
/*< make a simple axis >*/
{
    sf_axis AA;
    
    AA = (sf_axis) sf_alloc(1,sizeof(*AA));
    AA->n=n;
    AA->o=o;
    AA->d=d;
    AA->l=NULL;
    AA->u=NULL;

    return AA;
}

sf_axis sf_iaxa(sf_file FF, int i) 
/*< read axis i >*/
{
    char key[7];
    sf_axis AA;
    
    AA = (sf_axis) sf_alloc(1,sizeof(*AA));

    (void) snprintf(key,3,"n%d",i);
    if (!sf_histint (FF,key,&(AA->n))) AA->n=1;

    (void) snprintf(key,3,"o%d",i);
    if (!sf_histfloat(FF,key,&(AA->o))) AA->o=0;

    (void) snprintf(key,3,"d%d",i);
    if (!sf_histfloat(FF,key,&(AA->d))) AA->d=1;

    (void) snprintf(key,7,"label%d",i);
    AA->l = sf_histstring(FF,key);

    (void) snprintf(key,6,"unit%d",i);
    AA->u = sf_histstring(FF,key);

    return AA;
}

/*------------------------------------------------------------*/

void sf_oaxa(sf_file FF, const sf_axis AA, int i) 
/*< write axis i >*/
{
    char key[7];

    (void) snprintf(key,3,"n%d",i);
    sf_putint(FF,key,AA->n);
    
    (void) snprintf(key,3,"o%d",i);
    sf_putfloat(FF,key,AA->o);

    (void) snprintf(key,3,"d%d",i);
    sf_putfloat(FF,key,AA->d);

    if(NULL != AA->l) {	
	(void) snprintf(key,7,"label%d",i);
	sf_putstring(FF,key,AA->l);
    }

    if(NULL != AA->u) {
	(void) snprintf(key,6,"unit%d",i);
	sf_putstring(FF,key,AA->u);
    }
}

/*------------------------------------------------------------*/

void sf_raxa(const sf_axis AA) 
/*< report information on axis AA >*/
{    
    sf_warning("n=%4d \t o=%f \t d=%f \t l=\"%s\" \t u=\"%s\"",
	       AA->n,AA->o,AA->d,AA->l,AA->u);
}

/*------------------------------------------------------------*/

int sf_n(const sf_axis AA) 
/*< access axis length >*/
{return AA->n; }

float sf_o(const sf_axis AA) 
/*< access axis origin >*/
{return AA->o; }

float sf_d(const sf_axis AA) 
/*< access axis sampling >*/
{return AA->d; }

sf_axa sf_nod(const sf_axis AA) 
/*< access length, origin, and sampling >*/
{
    sf_axa BB;
    BB.n=AA->n;
    BB.o=AA->o;
    BB.d=AA->d;
    return BB;
}

void sf_setn(sf_axis AA, int n)
/*< change axis length >*/
{ AA->n=n; }

void sf_seto(sf_axis AA, float o)
/*< change axis origin >*/
{ AA->o=o; }

void sf_setd(sf_axis AA, float d)
/*< change axis sampling >*/
{ AA->d=d; }

void sf_setlabel(sf_axis AA, const char* label)
/*< change axis label >*/
{
    size_t len;

    if (NULL != AA->l) free(AA->l);

    len = strlen(label)+1;
    AA->l = sf_charalloc(len);
    strncpy(AA->l,label,len);
}

