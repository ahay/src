#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "simtab.h"

int main (void) {
    sf_simtab table;

    table = sf_simtab_init (2);
   
    sf_simtab_put(table,"a=1");
    assert(strcmp(sf_simtab_get(table,"a"),"1")==0);
    sf_simtab_put(table,"b=aaadf");
    assert(strcmp(sf_simtab_get(table,"b"),"aaadf")==0);
    sf_simtab_put(table,"a=2");
    assert(strcmp(sf_simtab_get(table,"a"),"2")==0);
    sf_simtab_put(table,"aasdaf=2");
    assert(strcmp(sf_simtab_get(table,"aasdaf"),"2")==0);
    sf_simtab_put(table,"adsdfadsfg");
    
    sf_simtab_close (table);
    exit (0);
}

/* 	$Id$	 */
