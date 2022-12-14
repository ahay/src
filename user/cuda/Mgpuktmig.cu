/* Prestack time migration (2-D/3-D) with CUDA. */
/*
  Copyright (C) 2010 University of Texas at Austin
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <math.h>

#include <cuda.h>
#include <cuda_runtime_api.h>

extern "C" {
#include <rsf.h>
}

#include "ktmig.cu"

static void sf_check_gpu_error (const char *msg) {
    cudaError_t err = cudaGetLastError ();
    if (cudaSuccess != err) 
        sf_error ("Cuda error: %s: %s", msg, cudaGetErrorString (err));
}

int main (int argc, char* argv[]) {
    /* Counters */
    int i = 0, j = 0, k = 0, l, n;

    /* Input data parameters */
    int nt, nx, ny = 1, nix = 1, nin = 1, osize, ntr, btr, dbtr;
    float ot, dt;

    /* Apperture parameters */
    int ix, iy, minix, miniy, maxix, maxiy;
    /* Apperture half-width in each direction */
    int apx, apy;
    /* Apperture first indices in x,y and number of x,y locations */
    int aox, aoy, anx, any;

    /* Image(output) space parameters */
    int ont, onx, ony;
    float oot, oox, ooy;
    float odt, odx, ody;

    /* Antialias filter parameters */
    int maxtri;
    float trfact, trm;

    /* Aperture corners */
    int el_cx1, el_cx2, el_cy1, el_cy2;
    int blk_y;
    float el_x, el_y;

    /* Input traces, output image, velocity */
    float *t, *img, *v;
    /* Aperture indices */
    int *ap;
    /* Coordinates: shot, receiver, and midpoint */
    sf_complex *sxy, *gxy, *cxy;

    char *vrmsfile, *sxsyfile, *gxgyfile, *cxcyfile;
    sf_file data, image, vrms, sxsy, gxgy, cxcy;

    bool verb, time, aa, diff;

    /* CUDA stuff */
    dim3 dimgrid (1, 1, 1);
    dim3 dimblock (BLOCK_SIZE, 1, 1);
    float *d_t, *d_v, *d_img, val;
    cudaArray *d_rt;
    float2 *d_sxy, *d_gxy;
    float4 *d_ixy;
    uint *d_ap;
    int devcnt = 0;
    CUdevice devnum = 0; /* Use the first device by default */
    size_t cudamem = 0;
    cudaDeviceProp deviceProp;
    sf_timer total_timer, memcpy_timer,
             aux_kernel_timer, main_kernel_timer;

    sf_init (argc, argv);

    if (!sf_getbool ("verb", &verb)) verb = false;
    /* Verbosity flag */
    if (!sf_getbool ("time", &time)) time = false;
    /* Total time measurement time */
    if (!sf_getbool ("aa", &aa)) aa = true;
    /* Antialiaing flag */
    if (!sf_getbool ("diff", &diff)) diff = true;
    /* Differentiation flag */

    cudaSetDevice (0);
    sf_check_gpu_error ("Device initialization");
    cudaGetDeviceCount (&devcnt);
    if (verb)
        sf_warning ("Number of CUDA devices: %d", devcnt);
    if (0 == devcnt)
        sf_error ("There is no device supporting CUDA");
#if CUDART_VERSION >= 2020
    n = 0;
    /* Search for a device without kernel exec timeout */
    do {
        cudaGetDeviceProperties (&deviceProp, n);
        sf_check_gpu_error ("CUDA device properties request");
        n++;
    } while (n < devcnt && deviceProp.kernelExecTimeoutEnabled);
    if (deviceProp.kernelExecTimeoutEnabled) {
        if (verb)
            sf_warning ("All devices have kernel exec timeout");
        cudaGetDeviceProperties (&deviceProp, devnum);
        sf_check_gpu_error ("CUDA device properties request");
        cudamem = deviceProp.totalGlobalMem;
        if (verb)
            sf_warning ("Available global memory: %.2fMb", cudamem*1e-6);
        /* Reserve some space for graphics */
        if (cudamem < (1 << 28)) {
            cudamem /= 2;
        } else {
            cudamem -= 1 << 28;
        }
        if (verb)
            sf_warning ("Assuming available global memory: %.2fMb", cudamem*1e-6);
    } else { /* Use almost all memory on the available device */
        devnum = n - 1;
        cudamem = deviceProp.totalGlobalMem;
        if (verb)
            sf_warning ("Available global memory: %.2fMb", cudamem*1e-6);
        cudamem -= 1 << 25;
        if (verb)
            sf_warning ("Assuming available global memory: %.2fMb", cudamem*1e-6);
    }
#else
    cudaGetDeviceProperties (&deviceProp, devnum);
    sf_check_gpu_error ("CUDA device properties request");
    cudamem = deviceProp.totalGlobalMem;
#endif
    cudaSetDevice (devnum);
    sf_check_gpu_error ("CUDA device selection");

    data = sf_input ("in");
    image = sf_output ("out");

    if (SF_FLOAT != sf_gettype (data))
        sf_error ("Need float input");
    if (!sf_histint (data, "n1", &nt)) sf_error ("No n1= in input");
    if (!sf_histint (data, "n2", &nx)) sf_error ("No n2= in input");
    if (!sf_histint (data, "n3", &ny)) ny = 1;
    if (!sf_histint (data, "n4", &nin)) nin = 1;
    if (!sf_histint (data, "n5", &nix)) nix = 1;
    ntr = nx*ny*nin*nix;

    if (!sf_histfloat (data, "d1", &dt)) sf_error ("No d1= in input");
    if (!sf_histfloat (data, "o1", &ot)) ot = 0.;

    vrmsfile = sf_getstring ("vrms");
    /* File with RMS velocities */
    if (NULL == vrmsfile) sf_error ("Need vrms="); 
    vrms = sf_input ("vrms");
    if (SF_FLOAT != sf_gettype (vrms)) sf_error ("Need float vrms");

    if (!sf_histint (vrms, "n1", &ont)) sf_error ("No n1= in vrms");
    if (!sf_histint (vrms, "n2", &onx)) sf_error ("No n2= in vrms");
    if (!sf_histint (vrms, "n3", &ony)) ony = 1;
    osize = ont*onx*ony;

    if (!sf_histfloat (vrms, "d1", &odt)) sf_error ("No d1= in vrms");
    if (!sf_histfloat (vrms, "d2", &odx)) sf_error ("No d2= in vrms");
    if (!sf_histfloat (vrms, "d3", &ody)) ody = 1.0;

    if (!sf_histfloat (vrms, "o1", &oot)) oot = 0.;
    if (!sf_histfloat (vrms, "o2", &oox)) oox = 0.;
    if (!sf_histfloat (vrms, "o3", &ooy)) ooy = 0.;
    if (verb)
        sf_warning ("Image size: %d x %d x %d", ont, onx, ony);

    if (!sf_getint ("dbtr", &dbtr)) dbtr = -1;
    /* Desired number of traces per block of threads */
    /* Check memory bounds - there should be enough space for: 
       2*osize - velocity volume + image volume,
       5*onx*ony - array of source and receiver coordinates + aperture indices,
       2*nt + 4*onx*ony - for each trace there should be double space for time samples +
                          surface vectors to image points */
    btr = (cudamem/sizeof(float) - 2*osize - 5*onx*ony)/(2*nt + 4*onx*ony);
    if (btr < 1 || btr*nt*sizeof(float) > cudamem)
        sf_error ("Not enough memory to migrate this dataset");
    if (verb)
        sf_warning ("Maximum number of traces per block - %d", btr);
    /* If the total number of traces is less than what GPU can digest -
       reduce the block size to the total number of traces */
    if (ntr < btr)
        btr = ntr;
    /* Set number of traces per block to what user requested */
    if (dbtr < btr && dbtr > 0)
        btr = dbtr;
    /* Check if GPU limit on 2D interpolated texture length is not exceeded */
    if (btr > MAX_TEXSIDE)
        btr = MAX_TEXSIDE;
    if (verb)
        sf_warning ("Setting number of traces per block to %d", btr);

    sxsyfile = sf_getstring ("sxsy");
    /* File with shot coordinates */
    if (NULL == sxsyfile) sf_error ("Need sxsy="); 
    sxsy = sf_input ("sxsy");
    if (SF_COMPLEX != sf_gettype (sxsy)) sf_error ("Need complex sxsy");
    if (!sf_histint (sxsy, "n2", &n)) sf_error ("No n2= in sxsy");
    if (n != ntr) sf_error ("Number of values in sxsy is not equal to number of input traces");

    gxgyfile = sf_getstring ("gxgy");
    /* File with receiver coordinates */
    if (NULL == gxgyfile) sf_error ("Need gxgy="); 
    gxgy = sf_input ("gxgy");
    if (SF_COMPLEX != sf_gettype (gxgy)) sf_error ("Need complex gxgy");
    if (!sf_histint (gxgy, "n2", &n)) sf_error ("No n2= in gxgy");
    if (n != ntr) sf_error ("Number of values in gxgy is not equal to number of input traces");

    cxcyfile = sf_getstring ("cxcy");
    /* File with midpoint coordinates */
    if (NULL == cxcyfile) sf_error ("Need cxcy="); 
    cxcy = sf_input ("cxcy");
    if (SF_COMPLEX != sf_gettype (cxcy)) sf_error ("Need complex cxcy");
    if (!sf_histint (cxcy, "n2", &n)) sf_error ("No n2= in cxcy");
    if (n != ntr) sf_error ("Number of values in cxcy is not equal to number of input traces");

    if (!sf_getint ("apx", &apx)) apx = onx/2;
    /* Apperture half-width in x direction */
    if (!sf_getint ("apy", &apy)) apy = ony/2;
    /* Apperture half-width in y direction */

    if (!sf_getint ("maxtri", &maxtri)) maxtri = 13;
    /* Maximum half-length of the antialias filter */
    if (!sf_getfloat ("trfact", &trfact)) trfact = 4.0*(0.5*(odx + ody)/dt);
    /* Trace factor for antialias filter length calculation */

    /* Initiate output */
    sf_putint (image, "n1", ont);
    sf_putint (image, "n2", onx);
    sf_putint (image, "n3", ony);
    sf_putint (image, "n4", 1);
    sf_putint (image, "n5", 1);
    sf_putfloat (image, "d1", odt);
    sf_putfloat (image, "d2", odx);
    sf_putfloat (image, "d3", ody);
    sf_putfloat (image, "d4", 0.0);
    sf_putfloat (image, "d5", 0.0);
    sf_putfloat (image, "o1", oot);
    sf_putfloat (image, "o2", oox);
    sf_putfloat (image, "o3", ooy);
    sf_putfloat (image, "o4", 0.0);
    sf_putfloat (image, "o5", 0.0);

    if (verb || time)
        total_timer = sf_timer_init ();
    if (verb) {
        memcpy_timer = sf_timer_init ();
        main_kernel_timer = sf_timer_init ();
        aux_kernel_timer = sf_timer_init ();
    }

    if (verb || time)
        sf_timer_start (total_timer);

    v = sf_floatalloc (osize);
    img = sf_floatalloc (osize);

    sf_floatread (v, osize, vrms);

    t  = sf_floatalloc (btr*nt);
    sxy = sf_complexalloc (btr);
    gxy = sf_complexalloc (btr);
    cxy = sf_complexalloc (btr);

    /* Input data vector on GPU */
    cudaMalloc ((void**)&d_t, btr*nt*sizeof(float));
    sf_check_gpu_error ("GPU malloc for t");
    /* Input read-only data vector on GPU for texture interpolation */
    cudaMallocArray (&d_rt, &t_i.channelDesc, nt, btr);
    sf_check_gpu_error ("GPU malloc for rt");
    cudaBindTextureToArray (t_i, d_rt);
    sf_check_gpu_error ("GPU tex bind for rt");
    /* Activate linear interpolation on input traces for
       higher bandwidth in antialiasing */
    t_i.normalized = false;
    t_i.filterMode = cudaFilterModeLinear;
    t_i.addressMode[0] = cudaAddressModeClamp;
    t_i.addressMode[1] = cudaAddressModeClamp;
    /* Array of source coordinates on GPU */
    cudaMalloc ((void**)&d_sxy, btr*sizeof(float2));
    sf_check_gpu_error ("GPU malloc for sxy");
    cudaBindTexture (0, t_sxy, d_sxy, btr*sizeof(float2));
    sf_check_gpu_error ("GPU tex bind for sxy");
    /* Array of receiver coordinates on GPU */
    cudaMalloc ((void**)&d_gxy, btr*sizeof(float2));
    sf_check_gpu_error ("GPU malloc for gxy");
    cudaBindTexture (0, t_gxy, d_gxy, btr*sizeof(float2));
    sf_check_gpu_error ("GPU tex bind for gxy");
    /* Vector of velocities on GPU */
    cudaMalloc ((void**)&d_v, osize*sizeof(float));
    sf_check_gpu_error ("GPU malloc for v");
    /* Image vector on GPU */
    cudaMalloc ((void**)&d_img, osize*sizeof(float));
    sf_check_gpu_error ("GPU malloc for img");
    cudaMemset (d_img, 0, osize*sizeof(float));
    sf_check_gpu_error ("GPU memset for img");

    /* Array of surface vectors to source/receiver on GPU,
       size is onx*ony rounded to the next number divisible by
       BLOCK_SIZE */
    n = btr*((int)(ceilf (onx*ony/(float)BLOCK_SIZE)*BLOCK_SIZE));
    cudaMalloc ((void**)&d_ixy, n*sizeof(float4));
    sf_check_gpu_error ("GPU malloc for ixy");
    cudaBindTexture (0, t_ixy, d_ixy, n*sizeof(float4));
    sf_check_gpu_error ("GPU tex bind for ixy");

    /* Array of aperture indices */
    n = ((int)(ceilf (onx*ony/(float)BLOCK_SIZE)*BLOCK_SIZE));
    ap = sf_intalloc (n);
    cudaMalloc ((void**)&d_ap, n*sizeof(int));
    sf_check_gpu_error ("GPU malloc for ap");
    cudaBindTexture (0, t_ap, d_ap, n*sizeof(int));
    sf_check_gpu_error ("GPU tex bind for ap");

    if (verb) {
        sf_warning ("Sending velocities to GPU");
        sf_timer_start (memcpy_timer);
    }
    cudaMemcpy (d_v, v, osize*sizeof(float), cudaMemcpyHostToDevice);
    sf_check_gpu_error ("Velocities transfer to GPU");

    cudaMemcpyToSymbol ("c_nt", &nt, sizeof(int), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol ("c_ont", &ont, sizeof(int), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol ("c_ot", &ot, sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol ("c_oot", &oot, sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol ("c_dt", &dt, sizeof(float), 0, cudaMemcpyHostToDevice);
    val = 1.0/dt;
    cudaMemcpyToSymbol ("c_idt", &val, sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol ("c_odt", &odt, sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol ("c_oox", &oox, sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol ("c_onx", &onx, sizeof(int), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol ("c_odx", &odx, sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol ("c_ooy", &ooy, sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol ("c_ony", &ony, sizeof(int), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol ("c_ody", &ody, sizeof(float), 0, cudaMemcpyHostToDevice);
    cudaMemcpyToSymbol ("c_trf", &trfact, sizeof(float), 0, cudaMemcpyHostToDevice);
    trm = maxtri;
    cudaMemcpyToSymbol ("c_trm", &trm, sizeof(float), 0, cudaMemcpyHostToDevice);
    val = nt - 1;
    cudaMemcpyToSymbol ("c_maxnt", &val, sizeof(float), 0, cudaMemcpyHostToDevice);
    sf_check_gpu_error ("Constants transfer to GPU");
    if (verb)
        sf_timer_stop (memcpy_timer);

    if (verb)
        sf_warning ("Migrating traces in chunks of %d", btr);

    /* Loop over input traces */
    i = 0;
    while (i < ntr) {
        /* How many to read */
        k = ((i + btr) < ntr)
          ? btr
          : ntr - i;
        if (verb)
            sf_warning ("Processing traces %d-%d out of %d", i, i + k - 1, ntr);
        /* Read input data */
        sf_floatread (t, nt*k, data);
        sf_complexread (sxy, k, sxsy);
        sf_complexread (gxy, k, gxgy);
        sf_complexread (cxy, k, cxcy);

        /* Find CDP span */
        minix = onx - 1; maxix = 0;
        miniy = ony - 1; maxiy = 0;
        for (l = 0; l < k; l++) {
            ix = (int)((crealf (cxy[l]) - oox)/odx + 0.5f);
            iy = (int)((cimagf (cxy[l]) - ooy)/ody + 0.5f);
            if (ix < minix)
                minix = ix;
            if (ix > maxix)
                maxix = ix;
            if (iy < miniy)
                miniy = iy;
            if (iy > maxiy)
                maxiy = iy;
        }

        /* Aperture corners */
        el_cx1 = minix;
        el_cx2 = maxix;
        el_cy1 = miniy;
        el_cy2 = maxiy;
        /* Add apperture width */
        minix -= apx;
        if (minix < 0)
            minix = 0;
        miniy -= apy;
        if (miniy < 0)
            miniy = 0;
        maxix += apx;
        if (maxix >= onx)
            maxix = onx - 1;
        maxiy += apy;
        if (maxiy >= ony)
            maxiy = ony - 1;
        aox = minix;
        aoy = miniy;
        anx = maxix - minix + 1;
        any = maxiy - miniy + 1;
        if (verb)
            sf_warning ("Rectangular aperture: %d-%d, %d-%d", minix, maxix, miniy, maxiy);
        /* Build aperture with rounded corners */
        if (verb)
            sf_warning ("Building rounded aperture");
        l = 0;
        for (iy = miniy; iy <= maxiy; iy++) {
            for (ix = minix; ix <= maxix; ix++) {
                blk_y = (iy - miniy)*anx + (ix - minix);
                if ((ix >= el_cx1 && ix <= el_cx2) ||
                    (iy >= el_cy1 && iy <= el_cy2)) {
                    ap[l] = blk_y;
                    l++;
                    continue;
                }
                /* Distance to corners */
                if (ix < el_cx1)
                    el_x = ix - el_cx1;
                else
                    el_x = ix - el_cx2;
                if (iy < el_cy1)
                    el_y = iy - el_cy1;
                else
                    el_y = iy - el_cy2;
                /* Check if the point is within one of the ellipses */
                if ((el_x*el_x/(apx*apx) + el_y*el_y/(apy*apy)) < 1.0f) {
                    ap[l] = blk_y;
                    l++;
                }
            }
        }
        if (verb)
            sf_timer_start (memcpy_timer);
        cudaMemcpy (d_ap, ap, sizeof(int)*l, cudaMemcpyHostToDevice);
        if (verb)
            sf_timer_stop (memcpy_timer);
        sf_check_gpu_error ("Aperture indices transfer to GPU");

        /* Send data to GPU */
        if (verb) {
            sf_warning ("Sending input traces + source/receiver coordinates to GPU");
            sf_timer_start (memcpy_timer);
        }
        if (aa)
            cudaMemcpy (d_t, t, sizeof(float)*k*nt, cudaMemcpyHostToDevice);
        else
            cudaMemcpyToArray (d_rt, 0, 0, t, sizeof(float)*k*nt, cudaMemcpyHostToDevice);
        sf_check_gpu_error ("Input traces transfer to GPU");
        cudaMemcpy (d_sxy, sxy, k*sizeof(float2), cudaMemcpyHostToDevice);
        sf_check_gpu_error ("Input source coordinates transfer to GPU");
        cudaMemcpy (d_gxy, gxy, k*sizeof(float2), cudaMemcpyHostToDevice);
        sf_check_gpu_error ("Input receiver coordinates transfer to GPU");
        cudaMemcpyToSymbol ("c_obn", &k, sizeof(int), 0, cudaMemcpyHostToDevice);
        cudaMemcpyToSymbol ("c_aox", &aox, sizeof(int), 0, cudaMemcpyHostToDevice);
        cudaMemcpyToSymbol ("c_aoy", &aoy, sizeof(int), 0, cudaMemcpyHostToDevice);
        cudaMemcpyToSymbol ("c_anx", &anx, sizeof(int), 0, cudaMemcpyHostToDevice);
        cudaMemcpyToSymbol ("c_any", &any, sizeof(int), 0, cudaMemcpyHostToDevice);
        sf_check_gpu_error ("Constants transfer to GPU");
        if (verb) {
            sf_timer_stop (memcpy_timer);
            sf_warning ("Memory copy time: %f ms", sf_timer_get_diff_time (memcpy_timer));
        }

        /* Run the distance calculation kernel */
        if (verb) {
            sf_warning ("Running GPU distance calculation kernel");
            sf_timer_start (memcpy_timer);
        }
        n = (int)(ceilf (anx*any/(float)BLOCK_SIZE));
        dimgrid = dim3 (n, 1, 1);
        n = (int)(ceilf (anx*any/(float)BLOCK_SIZE)*BLOCK_SIZE);
        cudaMemcpyToSymbol ("c_ibn", &n, sizeof(int), 0, cudaMemcpyHostToDevice);
        sf_check_gpu_error ("Constants transfer to GPU");
        if (verb) {
            sf_timer_stop (memcpy_timer);
            sf_timer_start (aux_kernel_timer);
        }
        sf_gpu_ktmig_ixy<<<dimgrid, dimblock>>>(d_ixy);
        cudaDeviceSynchronize ();
        sf_check_gpu_error ("Distance kernel invocation");
        if (verb) {
            sf_timer_stop (aux_kernel_timer);
            sf_warning ("Distance kernel execution time: %f ms",
                        sf_timer_get_diff_time (aux_kernel_timer));
        }

        /* Run antialiasing preparation */
        if (aa || diff) {
            if (verb) {
                sf_warning ("Running GPU trace preparation kernels");
                sf_timer_start (aux_kernel_timer);
            }
            dimgrid = dim3 (k, 1, 1);
            if (diff) {
                sf_gpu_ktmig_sbdiff<<<dimgrid, dimblock>>>(d_t, nt/BLOCK_SIZE);
                cudaDeviceSynchronize ();
                sf_check_gpu_error ("Differentiation kernels invocation");
            }
            if (aa) {
                sf_gpu_ktmig_cint<<<dimgrid, dimblock>>>(d_t, nt/BLOCK_SIZE);
                cudaDeviceSynchronize ();
                sf_gpu_ktmig_acint<<<dimgrid, dimblock>>>(d_t, nt/BLOCK_SIZE);
                cudaDeviceSynchronize ();
                sf_check_gpu_error ("Integration kernels invocation");
            }
            if (verb) {
                sf_timer_stop (aux_kernel_timer);
                sf_warning ("Trace preparation kernels execution time: %f ms",
                            sf_timer_get_diff_time (aux_kernel_timer));
                sf_timer_start (memcpy_timer);
            }
            cudaMemcpyToArray (d_rt, 0, 0, d_t, sizeof(float)*k*nt, cudaMemcpyDeviceToDevice);
            sf_check_gpu_error ("Prerocessed input traces transfer on GPU");
            if (verb)
                sf_timer_stop (memcpy_timer);
        }

        /* Run the migration kernel */
        if (verb) {
            sf_warning ("Running GPU migration kernel");
            sf_timer_start (main_kernel_timer);
        }
        n = 0;
        while (l != 0) {
            dimgrid = dim3 (ont/BLOCK_SIZE, l < MAX_GRID_SIZE
                                            ? l : MAX_GRID_SIZE, 1);
            if (aa)
                sf_gpu_ktmig_kernel<<<dimgrid, dimblock>>>(d_v, d_img, n);
            else
                sf_gpu_ktmig_noaa_kernel<<<dimgrid, dimblock>>>(d_v, d_img, n);
            cudaDeviceSynchronize ();
            sf_check_gpu_error ("Migration kernel invocation");
            l = l < MAX_GRID_SIZE ? 0 : l - MAX_GRID_SIZE;
            n += MAX_GRID_SIZE;
        }
        if (verb) {
            sf_timer_stop (main_kernel_timer);
            sf_warning ("Migration kernel execution time: %f ms",
                        sf_timer_get_diff_time (main_kernel_timer));
        }

        j++;
        i += k;
    } /* End of loop over input traces */

    /* Get the image back */
    if (verb) {
        sf_warning ("Receiving image from GPU");
        sf_timer_start (memcpy_timer);
    }
    cudaMemcpy (img, d_img, osize*sizeof(float), cudaMemcpyDeviceToHost);
    sf_check_gpu_error ("Image transfer from GPU");
    if (verb)
        sf_timer_stop (memcpy_timer);

    sf_floatwrite (img, osize, image);

    if (verb || time)
        sf_timer_stop (total_timer);
    if (verb) {
        sf_warning ("*** Summary of wallclock time ***");
        sf_warning ("GPU memcpy time: %f ms",
                    sf_timer_get_total_time (memcpy_timer));
        sf_warning ("Auxillary kernels time: %f ms",
                    sf_timer_get_total_time (aux_kernel_timer));
        sf_warning ("Main kernel time: %f ms",
                    sf_timer_get_total_time (main_kernel_timer));
    }
    if (verb || time)
        sf_warning ("Total kernels + GPU I/O + disk I/O time: %f ms",
                    sf_timer_get_total_time (total_timer));

    cudaUnbindTexture (t_i);
    cudaUnbindTexture (t_sxy);
    cudaUnbindTexture (t_gxy);
    cudaUnbindTexture (t_ixy);
    cudaUnbindTexture (t_ap);
    sf_check_gpu_error ("GPU texture unbinding");

    cudaFreeArray (d_rt);
    cudaFree (d_t);
    cudaFree (d_sxy);
    cudaFree (d_gxy);
    cudaFree (d_v);
    cudaFree (d_img);
    cudaFree (d_ixy);
    cudaFree (d_ap);
    sf_check_gpu_error ("GPU memory freeing");

    free (ap);
    free (v);
    free (img);
    free (t);
    free (sxy);
    free (gxy);
    free (cxy);

    if (verb) {
        free (main_kernel_timer);
        free (aux_kernel_timer);
        free (memcpy_timer);
    }
    if (verb || time)
        free (total_timer);

    exit (0);
}

