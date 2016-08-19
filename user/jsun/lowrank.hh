extern "C"{
#include "rtmutil.h"
}

void lowrank_init(int jump, int seed_, int npk_, float eps_, float dt_, float *vel, fdm3d fdm, dft3d dft);
/*< initialize lowrank >*/

int lowrank_rank();
/*< perform lowrank decomposition >*/

void lowrank_mat(sf_complex **lt, sf_complex **rt);
/*< output lowrank matrices >*/
