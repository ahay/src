#include <rsf.h>

#include "dijkstra.h"
#include "pqueue.h"

typedef enum {L, R, U, D} dir;

typedef struct Step {
    dir d;
    Step* next;
} *step;

void dijkstra(int n1, int n2, float **ud, float **lr) {
    float **cost;
    step **path;

}
    
