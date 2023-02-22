/* define hash table for point cloud operations */
#include <rsf.h>
#include "hash.h"

#define EPS 1e-6

/*------------------------------------------------------------*/
/* hash table elements (points) */
typedef struct {
    pt3d       * p; /*     coordinates (x,y,z) */
    unsigned int i; /* index in the full cloud */
} point;

/* hash table */
point * hashTable;

/*------------------------------------------------------------*/
unsigned int hash(unsigned int nhash, 
                            pt3d * p, 
                            pt3d * o)
/*< hash function >*/
{
    /* distance from reference point */
    float r = sqrtf( pow(p->x - o->x, 2) +
                     pow(p->y - o->y, 2) +
                     pow(p->z - o->z, 2) );

    /* find the hash index */
    unsigned int ui;
    memcpy( &ui, &r, sizeof(float) ); /* get int representation */
    ui &= 0xfffffff8;                 /* mask off 3 bits */
    ui %= nhash;                      /* map to hash range */
    return ui;                        /* return hash */
}

/*------------------------------------------------------------*/
void htInit(unsigned int nhash)
/*< initialize hash table >*/
{
  /* allocate the hash table */
  hashTable = (point*) sf_alloc(nhash,sizeof(*hashTable));

  /* initialize the hash table */
  for(unsigned int h = 0; h < nhash; h++) {
      hashTable[h].p = NULL;
      hashTable[h].i = -1;
  }
}

/*------------------------------------------------------------*/
void htClose()
/*< deallocate hash table >*/
{
    free(hashTable);
}

/*------------------------------------------------------------*/
bool htInsert(unsigned int nhash, 
                        pt3d * p, 
                        pt3d * o, 
            unsigned int   i)
/*< insert in hash table >*/
{
    if( p == NULL ) return false;

    /* start with the computed hash index */
    unsigned int h = hash( nhash, p, o );

    /* then search for an open slot using open addressing */
    for(unsigned int j = 0; j < nhash; j++ ) {
        unsigned int t = (h + j) % nhash;
        if( hashTable[t].p == NULL) {
            hashTable[t].p = p;
            hashTable[t].i = i;
            return true;
        }
    }
    return false;
}

/*------------------------------------------------------------*/
bool htDelete(unsigned int nhash, 
                        pt3d * q, 
                        pt3d * o)
/*< delete from hash table >*/
{
    if( q == NULL ) return false;

    /* start with the computed hash index */
    unsigned int h = hash( nhash, q, o );
    /* then search for an open slot using open addressing */
    for(unsigned int j = 0; j < nhash; j++ ) {
        unsigned int t = (h + j) % nhash;
        if( hashTable[t].p != NULL) {
            if( SF_ABS(hashTable[t].p->x - q->x) < EPS) {
                if( SF_ABS(hashTable[t].p->y - q->y) < EPS) {
                    if( SF_ABS(hashTable[t].p->z - q->z) < EPS) {
                        hashTable[t].p = NULL;
                        hashTable[t].i = -1;

                        return true; /* found & deleted */
                    }
                }
            }
        }
    }

    return false; /* not found & deleted */
}

/*------------------------------------------------------------*/
unsigned int htLookup(unsigned int nhash, 
                                pt3d * q, 
                                pt3d * o)
/*< lookup in hash table >*/
{
    if( q == NULL ) return -1;

    /* start with the computed hash index */
    unsigned int h = hash( nhash, q, o );

    /* then search for an open slot using open addressing */
    for(unsigned int j = 0; j < nhash; j++ ) {
        unsigned int t = (h + j) % nhash;
        if( hashTable[t].p != NULL) {
            if( SF_ABS(hashTable[t].p->x - q->x) < EPS) {
                if( SF_ABS(hashTable[t].p->y - q->y) < EPS) {
                    if( SF_ABS(hashTable[t].p->z - q->z) < EPS) {
                        return hashTable[t].i; /* found */
                    }
                }
            }
        }
    }

    return -1; /* not found */
}

/*------------------------------------------------------------*/
unsigned int htRetrieve(unsigned int t)
/*< retrieve hash table >*/
{
    return hashTable[t].i;
}

/*------------------------------------------------------------*/
void htPrint(unsigned int nhash)
/*< print hash table>*/
{
    for(unsigned int h = 0; h < nhash; h++) {
        if(hashTable[h].p == NULL) {
            sf_warning("%d:\t ---",h);
        } else {
            sf_warning("%d:\t (%f,%f,%f) [%d]",h,
                        hashTable[h].p -> x,
                        hashTable[h].p -> y,
                        hashTable[h].p -> z,
                        hashTable[h].i);
        }
    }
}
