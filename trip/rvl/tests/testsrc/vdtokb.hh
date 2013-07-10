#ifndef __ASGPP_CONV_VDKB__
#define __ASGPP_CONV_VDKB__

// scalar conversion functions (velocity,density)-> (bulkmod,bouyancy)
// derivs are all self-adjoint
float k(float v, float d) { return v*v*d; }
float dkdv(float v, float d, float delta) { return 2*v*d*delta; }
float dkdd(float v, float d, float delta) { return v*v*delta; }
float b(float v, float d) { return 1.0f/d; }
float dbdv(float v, float d, float delta) { return 0.0f; }
float dbdd(float v, float d, float delta) { return -delta/(d*d); }

#endif
