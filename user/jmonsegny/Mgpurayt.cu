// Parallel shortest path ray tracing
/*
  Copyright (C) 2013 Jorge Monsegny
  
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

#include <rsf.hh>
#include <cuda.h>
#include <iostream>
using std::cerr;
using std::endl;
#include <valarray>
using std::valarray;
#include <fstream>
using std::ifstream;
#include <thrust/host_vector.h>
using thrust::host_vector;
#include <thrust/device_vector.h>
using thrust::device_vector;
using thrust::raw_pointer_cast;
#include <vector>
using std::vector;
#include <list>
using std::list;
#include <utility>
using std::pair;

const float INF = 1e10;

//Function to diagnose errors
static void sf_check_gpu_error(const char *msg) 
{
    cudaError_t err = cudaGetLastError();
    if(cudaSuccess != err) 
        sf_error("Cuda error: %s: %s", msg, cudaGetErrorString (err));
}

//Copy a tile of a global array to shared memory
inline __device__ void tile2shared(float* s_vels,float *vels,int order,int smemsize,int bsize,int ewidth,int s_col,int s_row,int s_ecol,int s_erow, int col,int row,int ecol,int erow)
{
    //Corners
    if(s_col < order && s_row < order){
        s_vels[         s_row*smemsize + s_col         ] = vels[         row*ewidth + col         ];
        s_vels[         s_row*smemsize + (s_ecol+bsize)] = vels[         row*ewidth + (ecol+bsize)];
        s_vels[(s_erow+bsize)*smemsize + s_col         ] = vels[(erow+bsize)*ewidth + col         ];
        s_vels[(s_erow+bsize)*smemsize + (s_ecol+bsize)] = vels[(erow+bsize)*ewidth + (ecol+bsize)];
    }
    //left and right border
    if(s_col < order){
        s_vels[s_erow*smemsize + s_col          ] = vels[erow*ewidth + col         ]; 
        s_vels[s_erow*smemsize + (s_ecol+bsize) ] = vels[erow*ewidth + (ecol+bsize)];
    }
    //upper and lower border
    if(s_row < order){
        s_vels[           s_row*smemsize + s_ecol] = vels[         row*ewidth + ecol];
        s_vels[(s_erow + bsize)*smemsize + s_ecol] = vels[(erow+bsize)*ewidth + ecol];
    }
    //Center of shared memory
    s_vels[s_erow*smemsize + s_ecol] = vels[erow*ewidth + ecol];
}


//Integrates velocity along straight trajectories
__device__ float calctime(float* s_vels,int smemsize,int ecol1,int erow1,int ecol2,int erow2,int ewidth,int edepth,float dx,int npts)
{
    //Calculate distance
    float distx = ecol2 - ecol1;
    float distz = erow2 - erow1;
    float dist = sqrt(distx*distx + distz*distz);
    float frac = 1.0f/npts;
    //First velocity
    float v1 = s_vels[smemsize*erow1 + ecol1];
    float acumslo = 0.0f;
    //Integrate time
    for(int i=1; i<=npts; i++){        
        //interpolate
        float rat = frac*i;
        float x = (1.0f-rat)*ecol1 + rat*ecol2;
        float z = (1.0f-rat)*erow1 + rat*erow2;
        //Horizontal index and residual
        int indx = x;
        float resx = x - indx;
        if(indx == smemsize - 1){
            indx--;
            resx = 1.0f;
        }
        //Vertical index and residual
        int indz = z;
        float resz = z - indz;
        if(indz == smemsize - 1){
            indz--;
            resz = 1.0f;
        }
        //Velocity interpolation
        float velabove = (1.0f-resx)*s_vels[indx + indz*smemsize]     + resx*s_vels[(indx+1) + indz*smemsize];
        float velbelow = (1.0f-resx)*s_vels[indx + (indz+1)*smemsize] + resx*s_vels[(indx+1) + (indz+1)*smemsize];
        float v2 = (1.0f-resz)*velabove + resz*velbelow;
        float midvel = 0.5*(v1+v2);
        if(midvel > 0.0f)
            acumslo += 1.0f/midvel;
        v1 = v2;
    }
    return dx*acumslo*dist*frac;
}

//Calculates traveltimes between neighbor nodes
__global__ void precalctimes(int bsize,int order,int smemsize,float* vels,float* prectimes,int width,int depth,float dx,int npts)
{
    extern __shared__ float s_vels[];
    //Coordinates inside the block
    int s_col = threadIdx.x;
    int s_row = threadIdx.y;
    //Coordinates inside the velocity model
    int col = blockIdx.x*blockDim.x + threadIdx.x;
    int row = blockIdx.y*blockDim.y + threadIdx.y;
    if(col >= width || row >= depth)
        return;
    //Coordinates inside the block with halos
    int s_ecol = s_col + order;
    int s_erow = s_row + order;
    //Coordinates inside the velocity model with halos
    int ecol = col + order;
    int erow = row + order;
    int ewidth = width + 2*order;
    int edepth = depth + 2*order;
    //Copy velocities used in this block
    tile2shared(s_vels,vels,order,smemsize,bsize,ewidth,s_col,s_row,s_ecol,s_erow,col,row,ecol,erow);
    __syncthreads();
    //Neighbor tour
    int stride = width*depth;
    for(int k =-order; k <= order; k++){
        for(int i=-order; i <= order; i++){
            int neind = (k+order)*(2*order+1) + (i+order); 
            prectimes[neind*stride + row*width + col] = calctime(s_vels,smemsize,s_ecol,s_erow,s_ecol+i,s_erow+k,ewidth,edepth,dx,npts);
        }
    }
}

//Relaxation function
__global__ void relax(int bsize,int order,int smemsize,float* vels,float* times,float* auxtimes,float* prectimes,int* auxpred,int width,int depth)
{
    extern __shared__ float s_times[];
    //Coordinates inside the block
    int s_col = threadIdx.x;
    int s_row = threadIdx.y;
    //Coordinates inside the velocity model
    int col = blockIdx.x*blockDim.x + threadIdx.x;
    int row = blockIdx.y*blockDim.y + threadIdx.y;
    int offset = row*width + col; 
    if(col >= width || row >= depth)
        return;
    //Coordinates inside the block with halos
    int s_ecol = s_col + order;
    int s_erow = s_row + order;
    //Coordinates inside the velocity model with halos
    int ecol = col + order;
    int erow = row + order;
    int ewidth = width + 2*order;
    int eoffset = erow*ewidth + ecol;
    //Current better time and predecesors
    float curtime = auxtimes[eoffset];
    int curpred = auxpred[offset];
    //Copy times used in this block
    tile2shared(s_times,times,order,smemsize,bsize,ewidth,s_col,s_row,s_ecol,s_erow,col,row,ecol,erow);
    __syncthreads();
    //Neighbor tour
    int stride = width*depth;
    //cuPrintf("(%d %d)\n",col,row);
    for(int k =-order; k <= order; k++){
        for(int i=-order; i <= order; i++){
            int neind = (k+order)*(2*order+1) + (i+order); 
            float newtime = prectimes[neind*stride + row*width + col] + s_times[(s_erow+k)*smemsize + (s_ecol+i)];
            if(newtime < curtime){
                curtime = newtime;
                curpred = (row+k)*width + (col+i);//Encodes the predecesor (col+i,row+k) as an integer
            }
        }
    }
    
    //Write the best times
    auxtimes[eoffset] = curtime;
    auxpred[offset] = curpred;
}

//Writes the results of the current relax iteration
__global__ void writeback(int order,float* times,float* auxtimes,int* pred,int* auxpred,int width,int depth,bool* stop)
{
    //Coordinates without apron
    int col = blockIdx.x*blockDim.x + threadIdx.x;
    int row = blockIdx.y*blockDim.y + threadIdx.y;
    int offset = row*width + col;
    if(col >= width || row >= depth)
        return;
    //Coordinates with apron
    int ecol = col + order;
    int erow = row + order;
    int ewidth = width + 2*order;
    int eoffset = erow*ewidth + ecol;
    float auxtime = auxtimes[eoffset];
    if(auxtime < times[eoffset]){
        times[eoffset] = auxtime;
        pred[offset] = auxpred[offset];
        *stop = false;
    }
    auxtimes[eoffset] = times[eoffset];
}

//Extracts the raypaths
int decodeRaypaths(int *pred, int width, int depth, vector<list<pair<int,int> > >& paths)
{
    int maxlen = 0;
    int totsize = width*depth;
    for(int j=0; j<totsize; j++){
        int curlen = 0;           
        int curcode = j;
        int code = j;          
        do{
            int z = code/width;
            int x = code%width;
            paths[curcode].push_back(pair<int,int>(x,z));
            code = pred[code];
            curlen++; 
        }while(-1 != code);
        if(curlen > maxlen)
            maxlen = curlen;
    }
    return maxlen;
}

int main(int argc, char* argv[]) 
{
    sf_init(argc, argv);

    //Input velocity model
    iRSF in;

    //Velocity model dimensions
    int width;
    in.get("n2", width);
    int depth;
    in.get("n1", depth);
    valarray<float> model(width*depth);
    in >> model;

    //Input parameters
    //////////////////////////////////////////////////////////////// 
    iRSF par(0);

    //Source coordinates
    int sx;
    par.get("sx", sx, 0);//Horizontal node source coordinate (int)
    int sz;
    par.get("sz",sz, 0);//Vertical node source coordinate (int)
    if(sx < 0 || sx >= width || sz < 0 || sz >= depth)
        sf_error("source (%d,%d) is out of model",sx,sz);

    //Block size
    int bsize;
    par.get("bs", bsize, 16);//Cuda block is a square with bs*bs threads. Must divide dimensions of in.rsf, bs >= 1 (int)
    if(bsize < 1)
        sf_error("block size %d is too small",bsize);

    //Forward star order
    int order;
    par.get("ord", order, 3);//Forward star has (ord*ord-1) nodes, ord >= 1 (int)
    if(order < 1)
        sf_error("forward star order %d is too small",order);

    //bsize + halos
    int smemsize = bsize + 2*order;

    //dx
    float dx;
    par.get("dx", dx, 1.0f);//Horizontal and vertical separation between nodes, dx > 0.0 (float)
    if(dx <= 0.0f)
        sf_error("Grid  spacing %f is too small",dx);

    //Ray output
    char* ray = sf_getstring("ray");//Output file for a sfgraph compatible ray file. Empty for no ray output.
    //par.get("ray", ray);//Output file for a sfgraph compatible ray file. Empty for no ray output.

    //Comptime output
    char* ctime = sf_getstring("ctime");//Output rsf file for computation time. Empty for no computation time output.
    //par.get("ctime", ctime);//Output rsf file for computation time. Empty for no computation time output.

    //Cuda initialization
    cuInit(0);
    sf_check_gpu_error("Device initialization");

    //Extended dimensions
    int ewidth = width + 2*order;
    int edepth = depth + 2*order;

    // Data structures
    //////////////////////////////////////////////////////////////// 

    //Velocity model setup
    host_vector<float> vels(ewidth*edepth, 0);
    for(int k=0; k<depth; k++){
        for(int i=0; i<width; i++){
            vels[(k+order)*ewidth + (i+order)] = model[i*depth + k];//Transpose model
        }
    }
    device_vector<float> d_vels = vels;
    float *d_vels_p = raw_pointer_cast(&d_vels[0]);
    sf_check_gpu_error("Velocity model setup");

    //Times setup
    host_vector<float> times(ewidth*edepth, INF);
    times[(order+sz)*ewidth + order+sx] = 0;
    device_vector<float> d_times(times);
    float *d_times_p = raw_pointer_cast(&d_times[0]);

    //Aux times setup
    device_vector<float> d_auxtimes(times);
    float *d_auxtimes_p = raw_pointer_cast(&d_auxtimes[0]);
    sf_check_gpu_error("Initial times setup");    

    //Precalculated times setup
    int nneighb = (2*order + 1)*(2*order + 1);
    device_vector<float> d_prectimes(nneighb*width*depth);
    float *d_prectimes_p = raw_pointer_cast(&d_prectimes[0]);    

    //Pred setup
    host_vector<int> pred(width*depth, -1);
    device_vector<int> d_pred(pred);
    int *d_pred_p = raw_pointer_cast(&d_pred[0]);

    //Aux pred setup
    device_vector<int> d_auxpred(width*depth, -1);
    int *d_auxpred_p = raw_pointer_cast(&d_auxpred[0]);
    sf_check_gpu_error("Predecesor setup");

    //Grid and block size
    dim3 block(bsize,bsize);
    dim3 grid((width + bsize - 1)/bsize,(depth + bsize - 1)/bsize);

    //Ending flag
    host_vector<bool> stop(1,false);
    device_vector<bool> d_stop(stop);
    bool *d_stop_p = raw_pointer_cast(&d_stop[0]);

    //Start timing
    cudaEvent_t startT, stopT;
    cudaEventCreate(&startT);
    cudaEventCreate(&stopT);
    cudaEventRecord(startT,0);

    // Ray Tracing
    ////////////////////////////////////////////////////////////////

    //Parallel times precalculation
    precalctimes<<<grid,block,smemsize*smemsize*sizeof(float)>>>(bsize,order,smemsize,d_vels_p,d_prectimes_p,width,depth,dx,order+1);
    cudaDeviceSynchronize();

    //Parallel 
    while(!stop[0]){
        stop[0] = true;
        d_stop = stop;
        relax<<<grid,block,smemsize*smemsize*sizeof(float)>>>\
        (bsize,order,smemsize,d_vels_p,d_times_p,d_auxtimes_p,d_prectimes_p,d_auxpred_p,width,depth);
        cudaDeviceSynchronize();
        writeback<<<grid,block>>>(order,d_times_p,d_auxtimes_p,d_pred_p,d_auxpred_p,width,depth,d_stop_p);
        cudaDeviceSynchronize();
        stop = d_stop;
    }

    //End timing
    cudaEventRecord(stopT,0);
    cudaEventSynchronize(stopT);
    float elapsedT;
    cudaEventElapsedTime(&elapsedT, startT, stopT);

    //Output
    ////////////////////////////////////////////////////////////////

    //Times output
    times = d_times;
    float *times_p = raw_pointer_cast(&times[0]);
    oRSF out;
    out.put("n1", width);
    out.put("n2", depth);
    out.put("d1", dx);
    out.put("d2", dx);
    for(int k=0; k<depth; k++){
        out << *(new valarray<float>(&times_p[(k+order)*ewidth+order],width));
    }

    //Ray paths output
    if(NULL != ray){
        //Extract paths
        pred = d_pred;
        int *pred_p = raw_pointer_cast(&pred[0]);
        vector<list<pair<int,int> > > paths(width*depth);
        int maxlen = decodeRaypaths(pred_p, width, depth, paths);
        //Set ray dataset
        oRSF rayout(ray); 
        rayout.put("n1", maxlen);
        rayout.put("n2", width);
        rayout.put("n3", depth);
        rayout.type(SF_COMPLEX);
        //Write paths
        valarray<float> tr(maxlen*2);
        int totsize = width*depth; 
        for(int j=0; j<totsize; j++){
            int curlen = 0;
            int ind = 0;
            for(std::list<pair<int,int> >::iterator it=paths[j].begin(); it != paths[j].end(); it++,curlen++,ind+=2){
                tr[ind]   = (*it).first;
                tr[ind+1] = (*it).second;
            }
            ind = curlen*2;             
            for(int i=curlen; i<maxlen; i++,ind+=2){
                tr[ind]   = tr[curlen*2 - 2];
                tr[ind+1] = tr[curlen*2 - 1];
            }
            rayout << tr;
        }
    }

    //Comptime output
    if(NULL != ctime){
        oRSF ctimeout(ctime);
        ctimeout.put("n1", 1);
        ctimeout.put("n2", 1);
        ctimeout.type(SF_FLOAT);
        ctimeout << elapsedT;
    }
    exit(0);
}
