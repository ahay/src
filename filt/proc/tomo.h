#ifndef _tomo_h
#define _tomo_h

/* initializer (constructor)
np - number of slopes 
n1, n2 - grid dimensions
d1, d2 - cell size
*/
void tomo_init(int np, int n1, int n2, float d1, float d2);

/* linear operator: forward and ajoint
--------------------------------------
adj: adjoint flag (either t = L s (adj=false) or s = L' t (adj=true) 
add: addition flag (if true, t = t + L v or s = s + L' t)
ns: total size of s (ns = nz*nx) nz is the fast axis
nt: size of t (nt = np*nx) np is the fast axis
s: slowness
t: traveltime
*/
void tomo_lop(bool adj, bool add, int ns, int nt, float* s, float* t);

#endif
