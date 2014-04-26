//M.Z. 1st default vectorization:AVX
//     2nd default vectorization:SSE
//     3nd default: no explicit SIMD vectorization
#include "cstd.h"
#ifdef __AVX__
#include <immintrin.h>
#elif defined(__SSE__)
#include <xmmintrin.h>
#endif


void acd_3d_8(float *** uc,
              float *** up,
              float *** csq,
              int * s,
              int * e,
              float c0,
              float * c1,
              float * c2,
              float * c3,
              float * c4,
              int * lbc,
              int * rbc) {
    
    //Checking if using AVX or SSE or none.
    //#ifdef __AVX__
    //    printf("Using AVX!\n");
    //#elif defined(__SSE__)
    //    printf("USing SSE!\n");
    //#else
    //    printf("No explicit SIMD vectorization!Specifiy compiler options for SSE/AVX\n");
    //#endif
    
    int i0, i1, i2;
    //TODO#1: modifying coefficient array layout
    float Xcoeff[4]={c1[0],c2[0],c3[0],c4[0]};
    float Ycoeff[4]={c1[1],c2[1],c3[1],c4[1]};
    float Zcoeff[4]={c1[2],c2[2],c3[2],c4[2]};
    
    //computing X&Y-derivatives
    for(i2 = s[2]; i2 <= e[2]; i2++)
        for(i1 = s[1]; i1 <= e[1]; i1++){
            float* out=up[i2][i1];
            float* in=uc[i2][i1];
            float* csq_i=csq[i2][i1];
            i0=s[0];
            assert(((unsigned long)&up[i2][i1][i0])%32==0);
#ifdef __AVX__
            for(; i0 <= e[0]; i0+=8){
                __m256 out_simd = _mm256_loadu_ps(&out[i0]);
                __m256 csq_simd = _mm256_broadcast_ss(&csq_i[i0]);
                __m256 tmp = _mm256_mul_ps(_mm256_broadcast_ss(&c0), csq_simd);
                tmp = _mm256_add_ps(tmp, _mm256_set1_ps(2));
                out_simd = _mm256_sub_ps(_mm256_mul_ps(tmp, _mm256_loadu_ps(&in[i0])), out_simd);
                
                //loading coeffcients of X-derivatives
                //ic := coefficient index
                for(int ic = 0; ic < 4; ic++){
                    __m256 in_simd=_mm256_add_ps(_mm256_loadu_ps(&in[i0-ic-1]), _mm256_loadu_ps(&in[i0+ic+1]));
                    __m256 coeff_simd=_mm256_mul_ps(csq_simd, _mm256_broadcast_ss(&Xcoeff[ic]));
                    out_simd=_mm256_add_ps(out_simd, _mm256_mul_ps(in_simd, coeff_simd));
                }
                
                //loading coeffcients of Y-derivatives
                //ic := coefficient index
                for(int ic=0; ic < 4; ic++){
                    __m256 in_simd=_mm256_add_ps(_mm256_loadu_ps(&uc[i2][i1-ic-1][i0]), _mm256_loadu_ps(&uc[i2][i1+ic+1][i0]));
                    __m256 coeff_simd=_mm256_mul_ps(coeff_simd, _mm256_broadcast_ss(&Ycoeff[ic]));
                    out_simd=_mm256_add_ps(out_simd, _mm256_mul_ps(in_simd, coeff_simd));
                }
                
                //storing into up[i2][i1][i0]
                _mm256_storeu_ps(&out[i0], out_simd);
            }
#elif defined(__SSE__)
            for(; i0 <= e[0]; i0+=4){
                __m128 out_simd = _mm_loadu_ps(&out[i0]);
                __m128 csq_simd = _mm_set1_ps(csq_i[i0]);
                __m128 tmp = _mm_mul_ps(_mm_set1_ps(c0), csq_simd);
                tmp = _mm_add_ps(tmp, _mm_set1_ps(2));
                out_simd = _mm_sub_ps(_mm_mul_ps(tmp, _mm_loadu_ps(&in[i0])), out_simd);
                
                //loading coeffcients of X-derivatives
                //ic := coefficient index
                for(int ic = 0; ic < 4; ic++){
                    __m128 in_simd=_mm_add_ps(_mm_loadu_ps(&in[i0-ic-1]), _mm_loadu_ps(&in[i0+ic+1]));
                    __m128 coeff_simd=_mm_mul_ps(csq_simd, _mm_set1_ps(Xcoeff[ic]));
                    out_simd=_mm_add_ps(out_simd, _mm_mul_ps(in_simd, coeff_simd));
                }
                
                //loading coeffcients of Y-derivatives
                //ic := coefficient index
                for(int ic=0; ic < 4; ic++){
                    __m128 in_simd=_mm_add_ps(_mm_loadu_ps(&uc[i2][i1-ic-1][i0]), _mm_loadu_ps(&uc[i2][i1+ic+1][i0]));
                    __m128 coeff_simd=_mm_mul_ps(coeff_simd, _mm_set1_ps(Ycoeff[ic]));
                    out_simd=_mm_add_ps(out_simd, _mm_mul_ps(in_simd, coeff_simd));
                }
                
                //storing into up[i2][i1][i0]
                _mm_storeu_ps(&out[i0], out_simd);
            }
#endif
            for(; i0 <= e[0]; i0++){
                up[i2][i1][i0] = 2.0*uc[i2][i1][i0] - up[i2][i1][i0] + csq[i2][i1][i0] *
                ( c0*uc[i2][i1][i0] +
                 Xcoeff[0]*(uc[i2][i1][i0+1] + uc[i2][i1][i0-1]) +
                 Xcoeff[1]*(uc[i2][i1][i0+2] + uc[i2][i1][i0-2]) +
                 Xcoeff[2]*(uc[i2][i1][i0+3] + uc[i2][i1][i0-3]) +
                 Xcoeff[3]*(uc[i2][i1][i0+4] + uc[i2][i1][i0-4]) +
                 Ycoeff[0]*(uc[i2][i1+1][i0] + uc[i2][i1-1][i0]) +
                 Ycoeff[1]*(uc[i2][i1+2][i0] + uc[i2][i1-2][i0]) +
                 Ycoeff[2]*(uc[i2][i1+3][i0] + uc[i2][i1-3][i0]) +
                 Ycoeff[3]*(uc[i2][i1+4][i0] + uc[i2][i1-4][i0])
                 );
            }
        }
    
    //computing Z-derivatives
    for(i1 = s[1]; i1 <= e[1]; i1++)
        for(i2 = s[2]; i2 <= e[2]; i2++){
            i0=s[0];
#ifdef __AVX__
            for(; i0 <= e[0]; i0+=8){
                __m256 out_simd=_mm256_loadu_ps(&up[i2][i1][i0]);
                __m256 csq_simd=_mm256_broadcast_ss(&csq[i2][i1][i0]);
                for(int ic = 0; ic < 4; ic ++){
                    __m256 in_simd=_mm256_add_ps(_mm256_loadu_ps(&up[i2-ic-1][i1][i0]), _mm256_loadu_ps(&up[i2+ic+1][i1][i0]));
                    __m256 coeff_simd=_mm256_mul_ps(_mm256_broadcast_ss(&Zcoeff[ic]), csq_simd);
                    out_simd=_mm256_add_ps(out_simd,_mm256_mul_ps(in_simd, coeff_simd));
                }
                _mm256_storeu_ps(&up[i2][i1][i0],out_simd);
            }
#elif defined(__SSE__)
            for(; i0 <= e[0]; i0+=4){
                __m128 out_simd=_mm_loadu_ps(&up[i2][i1][i0]);
                __m128 csq_simd=_mm_set1_ps(csq[i2][i1][i0]);
                for(int ic = 0; ic < 4; ic ++){
                    __m128 in_simd=_mm_add_ps(_mm_loadu_ps(&up[i2-ic-1][i1][i0]), _mm_loadu_ps(&up[i2+ic+1][i1][i0]));
                    __m128 coeff_simd=_mm_mul_ps(_mm_set1_ps(Zcoeff[ic]), csq_simd);
                    out_simd=_mm_add_ps(out_simd,_mm_mul_ps(in_simd, coeff_simd));
                }
                _mm_storeu_ps(&up[i2][i1][i0],out_simd);
            }
#endif
            for(; i0 <= e[0]; i0++){
                up[i2][i1][i0] +=
                Zcoeff[0]*(uc[i2+1][i1][i0]+uc[i2-1][i1][i0]) +
                Zcoeff[1]*(uc[i2+2][i1][i0]+uc[i2-2][i1][i0]) +
                Zcoeff[2]*(uc[i2+3][i1][i0]+uc[i2-3][i1][i0]) +
                Zcoeff[3]*(uc[i2+4][i1][i0]+uc[i2-4][i1][i0]);
            }
        }
    
    /* boundary conditions - note that uc[-1][i]=0 etc. */
    //M.Z. Vectorze boundary updations.
#define mirror(i) i-2
    // int nx = floor((e[0]-s[0]+1)/8)*8;
    // bool res=(nx!=(e[0]-s[0]+1)); //M.Z. we got remainders.
    int iorder;
    /* boundary conditions - note that uc[-1][i]=0 etc. */
    if (lbc[2]){
        for (i1=s[1];i1<=e[1];i1++) {
            for(iorder=2;iorder<=4;iorder++){
                float* out=up[s[2]-iorder][i1];
                float* in=up[s[2]+mirror(iorder)][i1];
                //#pragma simd
                //#pragma vector aligned
                //gcc will automatically optimize the following loop if -O3 is specified.
                //                for (i0=s[0];i0<nx+s[0];i0++)
                //                    out[i0]=-in[i0];//M.Z. loop auto-vectorization requires single entry.
                //
                //                if(res)
                //                    for(i0=nx+s[0];i0<=e[0];i0++)
                //                        out[i0]=-in[i0];
                for(i0=s[0];i0<=e[0];i0++)
                    out[i0]=-in[i0];
            }
        }
    }
    if (rbc[2]) {
        for (i1=s[1];i1<=e[1];i1++) {
            for(iorder=2;iorder<=4;iorder++){
                float* out=up[e[2]+iorder][i1];
                float* in=up[e[2]-mirror(iorder)][i1];
                //#pragma simd
                //#pragma vector aligned
                //                for (i0=s[0];i0<nx+s[0];i0++)
                //                    out[i0]=-in[i0];//M.Z. loop auto-vectorization requires single entry.
                //
                //                if(res)
                //                    for(i0=nx+s[0];i0<=e[0];i0++)
                //                        out[i0]=-in[i0];
                for(i0=s[0];i0<=e[0];i0++)
                    out[i0]=-in[i0];
            }
        }
    }
    if (lbc[1]) {
        for (i2=s[2];i2<=e[2];i2++) {
            for(iorder=2;iorder<=4;iorder++){
                float* out=up[i2][s[1]-iorder];
                float* in=up[i2][s[1]+mirror(iorder)];
                //#pragma simd
                //#pragma vector aligned
                //                for (i0=s[0];i0<nx+s[0];i0++)
                //                    out[i0]=-in[i0];//M.Z. loop auto-vectorization requires single entry.
                //
                //                if(res)
                //                    for(i0=nx+s[0];i0<=e[0];i0++)
                //                        out[i0]=-in[i0];
                for(i0=s[0];i0<=e[0];i0++)
                    out[i0]=-in[i0];
            }
        }
    }
    if (rbc[1]) {
        for (i2=s[2];i2<=e[2];i2++) {
            for(iorder=2;iorder<=4;iorder++){
                float* out=up[i2][e[1]+iorder];
                float* in=up[i2][e[1]-mirror(iorder)];
                //#pragma simd
                //#pragma vector aligned
                //                for (i0=s[0];i0<nx+s[0];i0++)
                //                    out[i0]=-in[i0];//M.Z. loop auto-vectorization requires single entry.
                //
                //                if(res)
                //                    for(i0=nx+s[0];i0<=e[0];i0++)
                //                        out[i0]=-in[i0];
                for(i0=s[0];i0<=e[0];i0++)
                    out[i0]=-in[i0];
            }
        }
    }
    if (lbc[0]) {
        for (i2=s[2];i2<=e[2];i2++) {
            for (i1=s[1];i1<=e[1];i1++) {
                up[i2][i1][s[0]-2]=-up[i2][i1][s[0]+0];
                up[i2][i1][s[0]-3]=-up[i2][i1][s[0]+1];
                up[i2][i1][s[0]-4]=-up[i2][i1][s[0]+2];
            }
        }
    }
    if (rbc[0]) {
        for (i2=s[2];i2<=e[2];i2++) {
            for (i1=s[1];i1<=e[1];i1++) {
                up[i2][i1][e[0]+2]=-up[i2][i1][e[0]-0];
                up[i2][i1][e[0]+3]=-up[i2][i1][e[0]-1];
                up[i2][i1][e[0]+4]=-up[i2][i1][e[0]-2];
            }
        }
        
    }
    
}
