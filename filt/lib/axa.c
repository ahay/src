#include "axa.h"

#include "files.h"
/*^*/

#ifndef _sf_axa_h
#define _sf_axa_h

typedef struct {
    int     n;
    float o,d;
} axa;
/*^*/

#endif

void iaxa(sf_file FF, axa *AA, const int i) 
/*< read [n,o,d] for axis i >*/
{
    char BB[3];
    int     n;
    float o,d;

    (void) snprintf(BB,3,"n%d",i);
    if( !sf_histint  (FF,BB,&n)) n=1;
    AA->n=n;

    (void) snprintf(BB,3,"o%d",i);
    if(! sf_histfloat(FF,BB,&o)) o=0;
    AA->o=o;

    (void) snprintf(BB,3,"d%d",i);
    if(! sf_histfloat(FF,BB,&d)) d=1;
    AA->d=d;
}

void oaxa(sf_file FF, axa *AA, const int i) 
/*< write [n,o,d] for axis i >*/
{
    char BB[3];

    (void) snprintf(BB,3,"n%d",i);
    sf_putint(FF,BB,AA->n);
    
    (void) snprintf(BB,3,"o%d",i);
    sf_putfloat(FF,BB,AA->o);

    (void) snprintf(BB,3,"d%d",i);
    sf_putfloat(FF,BB,AA->d);
}

void raxa(axa AA) 
/*< report [n,o,d] for axis AA >*/
{    
    fprintf(stderr,"n=%d o=%f d=%f \n",AA.n,AA.o,AA.d);
}
