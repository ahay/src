#include "axa.h"
#include <string.h>
#include <stdio.h>

#include "file.h"
/*^*/

#include "error.h"

#ifndef _sf_axa_h

typedef struct {
    int     n;
    float o,d;
    char   *l, *u;
} axa;
/*^*/

#endif

/*------------------------------------------------------------*/

void iaxa(sf_file FF, axa *AA, const int i) 
/*< read [n,o,d,l] for axis i >*/
{
    char key[7];
    int     n;
    float o,d;

    (void) snprintf(key,3,"n%d",i);
    if( !sf_histint  (FF,key,&n)) n=1;
    AA->n=n;

    (void) snprintf(key,3,"o%d",i);
    if(! sf_histfloat(FF,key,&o)) o=0;
    AA->o=o;

    (void) snprintf(key,3,"d%d",i);
    if(! sf_histfloat(FF,key,&d)) d=1;
    AA->d=d;

    (void) snprintf(key,7,"label%d",i);
    if( NULL == (AA->l = sf_histstring(FF,key))) AA->l=" ";

    (void) snprintf(key,6,"unit%d",i);
    if( NULL == (AA->u = sf_histstring(FF,key))) AA->u=" ";
}

/*------------------------------------------------------------*/

void oaxa(sf_file FF, axa *AA, const int i) 
/*< write [n,o,d,l] for axis i >*/
{
    char key[7];

    (void) snprintf(key,3,"n%d",i);
    sf_putint(FF,key,AA->n);
    
    (void) snprintf(key,3,"o%d",i);
    sf_putfloat(FF,key,AA->o);

    (void) snprintf(key,3,"d%d",i);
    sf_putfloat(FF,key,AA->d);

    (void) snprintf(key,7,"label%d",i);
    if(NULL == AA->l ) AA->l=" ";
    sf_putstring(FF,key,AA->l);

    (void) snprintf(key,6,"unit%d",i);
    if(NULL == AA->u ) AA->u=" ";
    sf_putstring(FF,key,AA->u);
}

/*------------------------------------------------------------*/

void raxa(axa AA) 
/*< report [n,o,d,l,u] for axis AA >*/
{    
    sf_warning("n=%4d \t o=%f \t d=%f \t l=\"%s\" \t u=\"%s\"",
	       AA.n,AA.o,AA.d,AA.l,AA.u);
}

/*------------------------------------------------------------*/
