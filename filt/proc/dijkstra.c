#include <float.h>

#include <rsf.h>

#include "dijkstra.h"
#include "pqueue.h"

typedef enum {L, R, U, D} dir;

typedef struct Step {
    dir d;
    Step* next;
} *step;

static int n1, n2, **status;
static float **cost;
static step **path;
static const float big_value = FLT_MAX;

void dijkstra_init(int m1, int m2)
{
    int i2, i1;

    n1=m1;
    n2=m2;

    cost = sf_floatalloc2(n1,n2);
    status = sf_intalloc2(n1,n2);

    path = (step**) sf_alloc(n2,sizeof(step*));
    path[0] = (step*) sf_alloc(n1*n2,sizeof(step));
    for (i2=0; i2 < n2; i2++) {
	if (i2) path[i2] = path[0] + i2*n1;
	for (i1=0; i1 < n1; i1++) {
	    cost[i2][i1] = big_number;
	    status[i2][i1] = SF_OUT;
	    path[i2][i1] = NULL;
	}
    }
}

void dijkstra_close(void)
{
    free(cost[0]);
    free(cost);
    free(status[0]);
    free(status);
    free(path[0]);
    free(path);
}

void dijkstra(int s1, int s2, float **ud, float **lr) 
{
    cost[s2][s1] = 0.;
    status[s2][s1] = SF_IN;
}
    
