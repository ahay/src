#include <stdlib.h>

#include <rsf.h>

struct skey {
    float key;
    int pos;
};

static int key_compare (const void *k1, const void *k2)
{
    float f1 = ((struct skey*) k1)->key;
    float f2 = ((struct skey*) k2)->key;
    return (f1 < f2)? -1: (f1 > f2)? 1: 0;
}

int main(int argc, char* argv[])
{
    int n1, n2, i2, esize;
    long pos;
    struct skey *sorted;
    float *unsorted;
    char *trace;
    sf_file in, head, out;

    sf_init (argc,argv);
 
    head = sf_input("head");
    if (SF_FLOAT != sf_gettype(head))
	sf_error("Need float header");
    n2 = sf_filesize(head);
 
    unsorted = sf_floatalloc(n2);
    sorted = (struct skey*) sf_alloc(n2,sizeof(struct skey));
    
    sf_read(unsorted,sizeof(float),n2,head);
    for (i2 = 0; i2 < n2; i2++) {
	sorted[i2].key = unsorted[i2];
	sorted[i2].pos = i2;
    }
    free (unsorted);
    sf_fileclose(head);

    qsort(sorted,n2,sizeof(struct skey),key_compare);

    in = sf_input ("in");
    out = sf_output ("out");
 
    if (!sf_histint(in,"n1",&n1)) n1=1;
    if (!sf_histint(in,"esize",&esize)) esize=4;
    n1 *= esize;

    trace = sf_charalloc(n1);

    sf_unpipe(in,n1*n2);
    sf_fileflush(out,in);
    sf_setformat(in,"raw");
    sf_setformat(out,"raw");

    pos = sf_tell(in);
    for (i2=0; i2<n2; i2++) {
	sf_seek(in,pos+(sorted[i2].pos)*n1,SEEK_SET);
	sf_read(trace,1,n1,in);
	sf_write(trace,1,n1,out);
    }

    exit(0);
}
    
