#include "axa.h"

#include "file.h"
/*^*/

#include "error.h"

#ifndef _sf_axa_h

typedef struct {
    int     n;
    float o,d;
    char   *l;
} axa;
/*^*/

#endif

/*------------------------------------------------------------*/

void iaxa(sf_file FF, axa *AA, const int i) 
/*< read [n,o,d,l] for axis i >*/
{
    char BB[3];
    char LL[7];
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

    (void) snprintf(LL,7,"label%d",i);
    if( NULL == (AA->l = sf_histstring(FF,LL))) AA->l=" ";
}

/*------------------------------------------------------------*/

void oaxa(sf_file FF, axa *AA, const int i) 
/*< write [n,o,d,l] for axis i >*/
{
    char BB[3];
    char LL[7];

    (void) snprintf(BB,3,"n%d",i);
    sf_putint(FF,BB,AA->n);
    
    (void) snprintf(BB,3,"o%d",i);
    sf_putfloat(FF,BB,AA->o);

    (void) snprintf(BB,3,"d%d",i);
    sf_putfloat(FF,BB,AA->d);

    (void) snprintf(LL,7,"label%d",i);
    sf_putstring(FF,LL,AA->l);
}

/*------------------------------------------------------------*/

void raxa(axa AA) 
/*< report [n,o,d,l] for axis AA >*/
{    
    sf_warning("n=%4d \t o=%f \t d=%f \t l=%s",AA.n,AA.o,AA.d,AA.l);
}

/*------------------------------------------------------------*/
