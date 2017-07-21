#include <rsf.h>

#define PI (3.141592654)
#define WAVELET(tf) \
    ((1-2*PI*PI*(tf)*(tf))*exp(-PI*PI*(tf)*(tf)))

// FD coefficients 
// B(0)*f(ix)+sum<n=1,M>{[B(n)*[f(ix+n)+f(ix-n)]}
#define B30 (-2.7222f)
#define B31 (1.5f)
#define B32 (-0.15f)
#define B33 (0.01111111f)

#define DX(f,ix,iz,idx2) \
    ((B30*f[ix][iz] \
    + B31*(f[ix+1][iz]+f[ix-1][iz]) \
    + B32*(f[ix+2][iz]+f[ix-2][iz]) \
    + B33*(f[ix+3][iz]+f[ix-3][iz]) ) * idx2) 

#define DZ(f,ix,iz,idz2) \
    ((B30*f[ix][iz] \
    + B31*(f[ix][iz+1]+f[ix][iz-1]) \
    + B32*(f[ix][iz+2]+f[ix][iz-2]) \
    + B33*(f[ix][iz+3]+f[ix][iz-3]) ) * idz2) 

#define LAPLACE(f,ix,iz,idx2,idz2) \
    ( DX(f,ix,iz,idx2) + DZ(f,ix,iz,idz2) )

int main ( int argc, char *argv[] )
{
    int ix, iz, IX, IZ, k, ksnap;
    float ox, oz, dx, dz, dt;
    float idx2, idz2, dt2;
    int nx, nz;
    int nxpad, nzpad, n0=4, nb=30, nsep=10;
    float fbnd=0.01, *bw;
    float **v0, **v, **v2; // v0 (nx, nz), v and v2 (nxpad, nzpad)
    float **next,**curr,**lapl,**dcur,**dprv;// wavefield
    float **rec;
    float ***snap;
    float tmax;
    int nt, ntsnap;
    int snapstep=5;
    float dtsnap;
    
    float sx; // x position
    int src_ix, src_iz=0;
    float srcfrq=30, srcdec=0.5;
    float *decx, *decz; // for source
    int nswavelet;
    float wavelet;
    
    sf_file velocity, data, snapshot;
    
    sf_init(argc,argv);
    
    // get vel
    velocity=sf_input("in");
    if (!sf_histint(velocity,"n1",&nz)) sf_error("No n1=");
    if (!sf_histfloat(velocity,"d1",&dz)) sf_error("No d1=");
    if (!sf_histfloat(velocity,"o1",&oz)) oz=0;
    if (!sf_histint(velocity,"n2",&nx))  sf_error("No n2=");
    if (!sf_histfloat(velocity,"d2",&dx)) sf_error("No d2=");
    if (!sf_histfloat(velocity,"o2",&ox)) ox=0;
    nxpad = nx+2*(n0+nb);
    nzpad = nz+2*(n0+nb)+nsep;
    v0=sf_floatalloc2(nz,nx);
    v =sf_floatalloc2(nzpad, nxpad);
    v2=sf_floatalloc2(nzpad, nxpad);
    sf_floatread(v0[0],nx*nz,velocity);
    for(ix=0;ix<nxpad;ix++){
        for(iz=0;iz<nzpad;iz++){
            IX = ix - (n0+nb);
            IZ = iz - (n0+nb+nsep);
            if(IX<0) IX=0;
            if(IX>nx-1) IX=nx-1;
            if(IZ<0) IZ=0;
            if(IZ>nz-1) IZ=nz-1;
            v [ix][iz] = v0[IX][IZ];
            v2[ix][iz] = v0[IX][IZ] * v0[IX][IZ];
        }
    }
    
    // wavefield
    next = sf_floatalloc2(nzpad, nxpad);
    curr = sf_floatalloc2(nzpad, nxpad);
    lapl = sf_floatalloc2(nzpad, nxpad);
    dcur = sf_floatalloc2(nzpad, nxpad);
    dprv = sf_floatalloc2(nzpad, nxpad);
    
    // dx, dz
    idx2 = 1.0/(dx*dx);
    idz2 = 1.0/(dz*dz);
    
    // time
    if (!sf_getfloat("tmax",&tmax)) tmax=1.0;
    if (!sf_getfloat("dt",&dt)) dt=0.001;
    nt=(int)(tmax/dt)+1;
    ntsnap = (int)(nt-1)/snapstep+1;
    dt2 = dt * dt;
    dtsnap=dt*snapstep;
    
    // source
    if (!sf_getfloat("sx",&sx)) sx=ox+nx*dx/2;
    nswavelet = (int) ( 1.0f/srcfrq/dt + 1 );
    src_ix = (int)((sx-ox)/dx) + n0 + nb;
    src_iz = n0 + nb + nsep;
    decx = sf_floatalloc(nxpad);
    decz = sf_floatalloc(nzpad);
    for(ix=0;ix<nxpad;ix++) decx[ix]=exp(-srcdec*srcdec*(ix-src_ix)*(ix-src_ix));
    for(iz=0;iz<nzpad;iz++) decz[iz]=exp(-srcdec*srcdec*(iz-src_iz)*(iz-src_iz));
       
    // shot gather
    data=sf_output("out");
    sf_putint(data,"n1", nt);
    sf_putfloat(data,"d1", dt);
	sf_putstring(data, "label1","Time");
	sf_putstring(data, "unit1","s");
    sf_putint(data,"n2", nx);
    sf_putfloat(data,"d2", dx/1000);
	sf_putstring(data, "label2","Distance");
	sf_putstring(data, "unit2","km");
	sf_putstring(data, "title","ShotGather");
    
    // shapshots
    snapshot = sf_output ("snapshot");
    sf_putint  ( snapshot, "n1", nz );
    sf_putfloat( snapshot, "d1", dz/1000 );
    sf_putfloat( snapshot, "o1", oz );
	sf_putstring(snapshot, "label1","Depth");
	sf_putstring(snapshot, "unit1","km");
    sf_putint  ( snapshot, "n2", nx );
    sf_putfloat( snapshot, "d2", dx/1000 );
    sf_putfloat( snapshot, "o2", ox );
	sf_putstring(snapshot, "label2","Distance");
	sf_putstring(snapshot, "unit2","km");
    sf_putint  ( snapshot, "n3", ntsnap );
    sf_putfloat( snapshot, "d3", dtsnap );
    sf_putfloat( snapshot, "o3", 0 );
	sf_putstring(snapshot, "label3","Time");
	sf_putstring(snapshot, "unit3","s");
	sf_putstring(snapshot, "title","snapshot");

    // results
    rec=sf_floatalloc2(nt, nx);
    snap = sf_floatalloc3(nz, nx, ntsnap);
    
    // abc 
    bw = sf_floatalloc(nb);
    for ( k=0; k<nb; k++)
        bw[k] = exp ( -fbnd*fbnd*(nb-1-k)*(nb-1-k) );
 
    sf_warning("check: dx %f, dz %f, dt %f", dx, dz, dt);
    sf_warning("check: idx2 %f, idz2 %f, dt2 %g", idx2, idz2, dt2 );
    sf_warning("check: nx %d, nz %d, n0 %d, nb %d, nsep %d, nxpad %d, nzpad %d ", 
        nx, nz, n0, nb, nsep, nxpad, nzpad);
    sf_warning("check: v %g, v2 %g", v[100][100], v2[100][100] );
    
    ksnap=0;
    for ( k=-nswavelet; k<nt; k++ )  {
        if(k>=0) sf_warning(">>>> it=%d/%d;",k,nt-1);
        
        // load source
        if(k<=nswavelet){
            wavelet=WAVELET(k*dt*srcfrq);
            for(ix=0;ix<nxpad;ix++)
                for(iz=0;iz<nzpad;iz++)
                    curr[ix][iz] += wavelet * decx[ix] * decz[iz];
        }

        // laplacian
        for(ix=n0; ix<nxpad-n0; ix++ )
            for(iz=n0; iz<nzpad-n0; iz++)
                lapl[ix][iz] = LAPLACE(curr,ix,iz,idx2,idz2);
        
        // step forward
        for(ix=0;ix<nxpad;ix++){
            for(iz=0;iz<nzpad;iz++){
		        dcur[ix][iz] = dprv[ix][iz] + lapl[ix][iz]*v2[ix][iz]*dt;
		        next[ix][iz] = curr[ix][iz] + dcur[ix][iz]*dt;
		        dprv[ix][iz] = dcur[ix][iz]; 
		        curr[ix][iz] = next[ix][iz];
     	    }
	    }
	    
	    // abc
	    for (ix=0; ix < nb; ix++) {  
		    for (iz=n0; iz < nzpad-n0; iz++) {
		        curr[ix+n0][iz] *= bw[ix];
		        curr[ix+nb+nx+n0][iz] *= bw[nb-1-ix];
		        dprv[ix+n0][iz] *= bw[ix];
		        dprv[ix+nb+nx+n0][iz] *= bw[nb-1-ix];
		    }
	    }
	    for (ix=n0; ix < nxpad-n0; ix++) {  
		    for (iz=0; iz < nb; iz++) {
		        curr[ix][iz+n0] *= bw[iz];
		        curr[ix][iz+nz+nb+n0] *= bw[nb-1-iz];
		        dprv[ix][iz+n0] *= bw[iz];
		        dprv[ix][iz+nz+nb+n0] *= bw[nb-1-iz];
		    }
	    }
	   
        if(k<0) continue;
        
        // get record
        for(ix=0;ix<nx;ix++)
            rec[ix][k] = curr[ix+n0+nb][n0+nb+nsep];
        
        // get snapshots
        if(0==k%snapstep){
            for(ix=0;ix<nx;ix++)
                for(iz=0;iz<nz;iz++)
                    snap[ksnap][ix][iz] = curr[ix+n0+nb][iz+n0+nb+nsep];
            ksnap ++;
        }
        
    }// end of time
    sf_warning(".");
    sf_floatwrite ( rec[0], nt*nx, data );
    sf_floatwrite ( snap[0][0], ntsnap*nx*nz, snapshot );
    
    exit(0);
}

