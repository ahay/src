#include "loconvol.h"
#include "helicon.h"

static filter aa;

void loconvol_init(filter aa_in)
{
    aa = aa_in;
}

void loconvol_lop(bool adj, bool add, int nx, int ny, float *xx, float *yy)
{
    helicon_init(aa);
    aa++;

    helicon_lop(adj, add, nx, ny, xx, yy);
}
