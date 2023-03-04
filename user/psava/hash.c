/* define hash table for point cloud operations */
#include <rsf.h>
#include "hash.h"

/*------------------------------------------------------------*/
/* hash table elements (points) */
typedef struct {
    pt3d             * p; /*     coordinates (x,y,z) */
    unsigned long long i; /* index in the full cloud */
} point;

/* hash table */
point * hashTable;

/*------------------------------------------------------------*/
void print_binary(unsigned long long num) 
/*< print binary numbers >*/
{
  int i;
  for (i = 31; i >= 0; i--) {
    fprintf(stderr,"%c", (num & (1u << i)) ? '1' : '0');
    if(i%8 == 0) fprintf(stderr," ");
  }
  fprintf(stderr,"\n");
}

/*------------------------------------------------------------*/
unsigned long long hashOLD(const unsigned long long nhash, 
                                const            pt3d * p, 
                                const            pt3d * o)
/*< hash function >*/
{
    /* hash attribute */
    double h = sqrtf( pow(p->x - o->x, 2) + 
                      pow(p->y - o->y, 2) + 
                      pow(p->z - o->z, 2) );

    unsigned long long a,hindx;
    memcpy( &a, &h, sizeof(double) ); /* get int representation */
    hindx = a & 0xfffffff8;           /* mask off 3 bits */

    return hindx % nhash;             /* return the hash index */
}

/*------------------------------------------------------------*/
unsigned long long hashF1a(const unsigned long long nhash, 
                                const            pt3d * p, 
                                const            pt3d * o)
/*< Fowler–Noll–Vo-1a hash function >*/
{
    const unsigned long long FNV_OFFSET = 14695981039346656037ull;
    const unsigned long long FNV_PRIME  = 1099511628211ull;

    double h = sqrtf( pow(p->x - o->x, 2) + 
                      pow(p->y - o->y, 2) + 
                      pow(p->z - o->z, 2) );
    const unsigned char* b = (unsigned char*) &h;

    unsigned long long hindx = FNV_OFFSET;
    for (int i = 0; i < sizeof(double); i++) {
        hindx = (hindx * FNV_PRIME) % (1<<31);
        hindx ^= b[i];
    }

    return hindx % nhash;
}

/*------------------------------------------------------------*/
unsigned long long hashF2a(const unsigned long long nhash, 
                                const            pt3d * p, 
                                const            pt3d * o)
/*< Fowler–Noll–Vo-2a hash function >*/
{
    const unsigned long long FNV_OFFSET = 14695981039346656037ull;
    const unsigned long long FNV_PRIME  = 1099511628211ull;

    double h = sqrtf( pow(p->x - o->x, 2) + 
                      pow(p->y - o->y, 2) + 
                      pow(p->z - o->z, 2) );
    const unsigned char* b = (unsigned char*) &h;

    unsigned long long hindx = FNV_OFFSET;
    for (int i = 0; i < sizeof(double); i++) {
        hindx ^= b[i];
        hindx *= FNV_PRIME;
    }
    hindx *= FNV_PRIME;
    hindx &= 0xffffffffffffffff;

    return hindx % nhash;
}

/*------------------------------------------------------------*/
void htInit(unsigned long long nhash)
/*< initialize hash table >*/
{
  /* allocate the hash table */
  hashTable = (point*) sf_alloc( nhash, sizeof(*hashTable) );
  //sf_warning(" HT address: %p %u",hashTable, (unsigned long long)hashTable);

  /* initialize the hash table */
  for(unsigned long long jhash = 0; jhash < nhash; jhash++) {
      hashTable[ jhash ].p = NULL;
      hashTable[ jhash ].i = -1;
  }
}

/*------------------------------------------------------------*/
void htClose()
/*< deallocate hash table >*/
{
    free(hashTable);
}

/*------------------------------------------------------------*/
unsigned long long htInsert( unsigned long long nhash, 
                                             pt3d * p, 
                                             pt3d * o, 
                              unsigned long long    i)
/*< insert in hash table >*/
{
    if( p == NULL ) return -1;
    unsigned long long jhash;

    /* start with the computed hash index */
    // unsigned long long h = hashOLD( nhash, p, o );
    // unsigned long long h = hashF1a( nhash, p, o );
    unsigned long long h = hashF2a( nhash, p, o );

    /* then search for an open slot using open addressing */
    for( jhash = 0; jhash < nhash; jhash++ ) {
        unsigned long long t = (h + jhash) % nhash;

        if( hashTable[ t ].p == NULL) {
            hashTable[ t ].p = p;
            hashTable[ t ].i = i;
            return jhash;
        }

    }
    return -1;
}

/*------------------------------------------------------------*/
unsigned long long htWrite( unsigned long long nhash, 
                                            pt3d * p, 
                                            pt3d * o, 
                            unsigned long long     i)
/*< write in hash table >*/
{
    if( p == NULL ) return -1;

    hashTable[ i ].p = p;
    hashTable[ i ].i = i;

    return -1;
}

/*------------------------------------------------------------*/
bool htDelete(unsigned long long nhash, 
                              pt3d * q, 
                              pt3d * o)
/*< delete from hash table >*/
{
    if( q == NULL ) return false;
    unsigned long long jhash;

    /* start with the computed hash index */
    // unsigned long long h = hashOLD( nhash, q, o );
    // unsigned long long h = hashF1a( nhash, q, o );
    unsigned long long h = hashF2a( nhash, q, o );

    /* then search for an open slot using open addressing */
    for( jhash = 0; jhash < nhash; jhash++ ) {
        unsigned long long t = (h + jhash) % nhash;
        
        if( hashTable[ t ].p != NULL) {
            double d = sqrtf( pow(hashTable[ t ].p->x - q->x, 2) + 
                              pow(hashTable[ t ].p->y - q->y, 2) + 
                              pow(hashTable[ t ].p->z - q->z, 2) );

            if( d < SF_EPS ) {
                hashTable[ t ].p = NULL;
                hashTable[ t ].i = -1;

                return true; // found & deleted
            }
        }

        /*
        if( hashTable[ t ].p != NULL) {
            if( SF_ABS(hashTable[ t ].p->x - q->x) < SF_EPS) {
                if( SF_ABS(hashTable[ t ].p->y - q->y) < SF_EPS) {
                    if( SF_ABS(hashTable[ t ].p->z - q->z) < SF_EPS) {
                        hashTable[ t ].p = NULL;
                        hashTable[ t ].i = -1;

                        return true; // found & deleted
                    }
                }
            }
        }
        */

    }

    return false; // not found & deleted
}

/*------------------------------------------------------------*/
unsigned long long htLookup(unsigned long long nhash, 
                                            pt3d * q, 
                                            pt3d * o)
/*< lookup in hash table >*/
{
    if( q == NULL ) return -1;
    unsigned long long jhash;

    /* start with the computed hash index */
    // unsigned long long h = hashOLD( nhash, q, o );
    // unsigned long long h = hashF1a( nhash, q, o );
    unsigned long long h = hashF2a( nhash, q, o );

    /* then search for an open slot using open addressing */
    for( jhash = 0; jhash < nhash; jhash++ ) {
        unsigned long long t = (h + jhash) % nhash;

        if( hashTable[ t ].p != NULL) {
            double d = sqrtf( pow(hashTable[ t ].p->x - q->x, 2) + 
                              pow(hashTable[ t ].p->y - q->y, 2) + 
                              pow(hashTable[ t ].p->z - q->z, 2) );

            if( d < SF_EPS ) {
                return hashTable[ t ].i; // found
            }
        }

        /*
        if( hashTable[ t ].p != NULL) {
            if( SF_ABS(hashTable[ t ].p->x - q->x) < SF_EPS) {
                if( SF_ABS(hashTable[ t ].p->y - q->y) < SF_EPS) {
                    if( SF_ABS(hashTable[ t ].p->z - q->z) < SF_EPS) {
                        return hashTable[ t ].i; // found
                    }
                }
            }
        }
        */

    }

    return -1; /* not found */
}

/*------------------------------------------------------------*/
unsigned long long htRetrieve(unsigned long long t)
/*< retrieve hash table >*/
{
    return hashTable[t].i;
}

/*------------------------------------------------------------*/
void htPrint(unsigned long long nhash)
/*< print hash table>*/
{
    unsigned long long jhash;

    for( jhash = 0; jhash < nhash; jhash++) {
        if(hashTable[ jhash ].p == NULL) {
            sf_warning("%d:\t ---",jhash);
        } else {
            sf_warning("%d:\t (%f,%f,%f) [%d]",jhash,
                        hashTable[ jhash ].p -> x,
                        hashTable[ jhash ].p -> y,
                        hashTable[ jhash ].p -> z,
                        hashTable[ jhash ].i);
        }
    }
}
