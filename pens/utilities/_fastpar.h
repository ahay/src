#ifndef FASTPARG_H
#define FASTPARG_H
/* added timestamp 3-8-86  stew */
typedef struct hash_dummy {
	struct hash_dummy *next;
	char *tag; int tlen;
	char *val; int vlen;
	int timestamp;
	} hash_item;

typedef union { int *i; float *f; double *g ; char *s; } MIXED;


#if NeedFunctionPrototypes
_XFUNCPROTOBEGIN
extern hash_item** new_queue( int );
extern int getpar_decode(hash_item **,int ,char*,char*,MIXED);
extern void putch_format(char*,char*,char*,MIXED);
extern hash_item *getpar_hash_lookup(register hash_item **,register int,
register char*,register int);
extern int getpar_push_input(register char*,register int);
extern int getpar_scan(register hash_item **,register int);
extern void getpar_string_store(hash_item **,int,register char *);
extern int getpar_hash_store(hash_item**,int,register char*,register char*,
     register int,int);
extern int getpar_stack_par( char*);
_XFUNCPROTOEND
#else
extern hash_item** new_queue();extern int getpar_decode();
extern void putch_format();
extern hash_item *getpar_hash_lookup();extern int getpar_push_input();
extern int getpar_scan();
extern void getpar_string_store();
extern int getpar_hash_store();
extern int getpar_stack_par();
#endif
#endif

