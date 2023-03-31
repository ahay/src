/* define hash table for point cloud operations */
#include <rsf.h>
#include "hash.h"

/*------------------------------------------------------------*/
/*------------------------------------------------------------*/
/* hash table elements (points) */
typedef struct {
    pt3d             * p; /*     coordinates (x,y,z) */
    unsigned long long i; /* index in the full cloud */
} point;

/* hash table */
point * hashTable;
/*------------------------------------------------------------*/
/*------------------------------------------------------------*/


/*------------------------------------------------------------*/
void print_binary(unsigned long long num) 
/*< print binary numbers >*/
{
  for( int i = 31; i >= 0; i--) {
    fprintf(stderr,"%c", (num & (1u << i)) ? '1' : '0');
    if(i%8 == 0) fprintf(stderr," ");
  }
  fprintf(stderr,"\n");
}

/*------------------------------------------------------------*/
unsigned long long hashF1a(const unsigned long long nhash, 
                           const                 pt3d * p, 
                           const                 pt3d * o)
/*< Fowler–Noll–Vo-1a hash function >*/
{
    const unsigned long long FNV_OFFSET = 14695981039346656037ull;
    const unsigned long long FNV_PRIME  = 1099511628211ull;

    double h = sqrtf( pow( p->x - o->x, 2) + 
                      pow( p->y - o->y, 2) + 
                      pow( p->z - o->z, 2) );
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
  // allocate the HT
  hashTable = (point*) sf_alloc( nhash, sizeof(*hashTable) );
  //sf_warning(" HT address: %p %u",hashTable, (unsigned long long)hashTable);

  // initialize the HT
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
                             unsigned long long     i)
/*< insert in hash table >*/
{
    if( p == NULL ) return -1;

    // start with the computed hash index
    unsigned long long hindx = hashF1a( nhash, p, o );

    // then search for an open slot using open addressing
    for(unsigned long long jhash = 0; jhash < nhash; jhash++ ) {
        unsigned long long tindx = (hindx + jhash) % nhash;

        if( hashTable[ tindx ].p == NULL) {
            hashTable[ tindx ].p = p;
            hashTable[ tindx ].i = i;
            return jhash;
        }

    }

    sf_warning("cannot insert in hash table %d", i);
    exit(1);
    //return -1;
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
                              pt3d * o,
                          double ddMIN)
/*< delete from hash table >*/
{
    if( q == NULL ) return false;

    // start with the computed hash index
    unsigned long long hindx = hashF1a( nhash, q, o );

    // then search for an open slot using open addressing
    for(unsigned long long jhash = 0; jhash < nhash; jhash++ ) {
        unsigned long long tindx = (hindx + jhash) % nhash;
        
        if( hashTable[ tindx ].p != NULL) {
            double dd = sqrtf( pow( hashTable[ tindx ].p->x - q->x, 2) + 
                               pow( hashTable[ tindx ].p->y - q->y, 2) + 
                               pow( hashTable[ tindx ].p->z - q->z, 2) );

            if( dd < ddMIN / 10 ) {
                hashTable[ tindx ].p = NULL;
                hashTable[ tindx ].i = -1;

                return true; // found & deleted
            }
        }

    }

    return false; // not found & deleted
}

/*------------------------------------------------------------*/
unsigned long long htLookup(unsigned long long nhash, 
                                            pt3d * q, 
                                            pt3d * o,
                                        double ddMIN)
/*< lookup in hash table >*/
{
    if( q == NULL ) return -1;

    // start with the computed hash index
    unsigned long long hindx = hashF1a( nhash, q, o );

    // then search for an open slot using open addressing
    for(unsigned long long jhash = 0; jhash < nhash; jhash++ ) {
        unsigned long long tindx = (hindx + jhash) % nhash;

        if( hashTable[ tindx ].p != NULL) {
            double dd = sqrtf( pow( hashTable[ tindx ].p->x - q->x, 2) + 
                               pow( hashTable[ tindx ].p->y - q->y, 2) + 
                               pow( hashTable[ tindx ].p->z - q->z, 2) );

            if( dd < ddMIN / 10 ) {
                return hashTable[ tindx ].i; // found
            }
        }

    }

    return -1; /* not found */
}

/*------------------------------------------------------------*/
unsigned long long htRetrieve(unsigned long long tindx)
/*< retrieve hash table >*/
{
    return hashTable[ tindx ].i;
}

/*------------------------------------------------------------*/
void htPrint(unsigned long long nhash)
/*< print hash table>*/
{
    for( unsigned long long jhash = 0; jhash < nhash; jhash++) {
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
