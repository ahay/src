extern "C"{
#include "rtmutil.h"
}

void lowrank_init(int jump, int seed_, int npk_, float eps_, float dt_, int media, fdm3d fdm, dft3d dft);
/*< initialize lowrank >*/

void lowrank_iso(float *vel);
/*< iso model setup >*/

void lowrank_tti(float *velx, float *velz, float *eta, float *theta);
/*< iso model setup >*/

int lowrank_rank();
/*< perform lowrank decomposition >*/

void lowrank_mat(sf_complex **lt, sf_complex **rt);
/*< output lowrank matrices >*/
