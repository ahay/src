/* 2-D Fourier finite-difference wave extrapolation: MPI + OMP*/
/*
  Copyright (C) 2009 University of Texas at Austin
  
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

#include <rsf.h>
#include <math.h>
#include <limits.h>
#include <mpi.h>
#ifdef _OPENMP
#include <omp.h>
#endif
static int nx, nz, nbt, nbb, nbl, nbr;
static float ct, cb, cl, cr;
static float *wt, *wb, *wl, *wr;
static float **ws;

void itoa(int n, char *s);

void bd_init(int n1,  int n2    /*model size:x*/,
             int nb1, int nb2  /*top, bottom*/,
             int nb3, int nb4  /*left, right*/,
             float c1, float c2 /*top, bottom*/,
             float c3, float c4 /*left, right*/);
/*< initialization >*/


void bd_close(void);
/*< free memory allocation>*/


void bd_decay(float **a /*2-D matrix*/);
/*< boundary decay>*/


void abc_cal(int abc /* decaying type*/,
             int nb  /* absorbing layer length*/, 
             float c /* decaying parameter*/,
             float* w /* output weight[nb] */);
/*< find absorbing coefficients >*/

void srcsm_init(float dz, float dx/*grid size*/);
/*< initialization >*/


void srcsm_close(void);
/*< free memory allocation>*/


void source_smooth(float **source /*source matrix*/, 
                   int iz /*source index*/, 
                   int ix /*source index*/,
                   float val/*source value*/);
/*< smooth source input >*/


float dehf(float k /*current frequency*/,
	   float kn /*highest frequency*/,
	   float a /*suppress factor*/,
	   float factor /*propotion*/);
/*< high frequency depressing>*/
int main(int argc, char* argv[]) 
{
    int nx, nt, nkx, nkz, ix, it, ikx, ikz, nz, iz, nbt, nbb, nbl, nbr, nxb, nzb, isx, isz, irz, ir;
    float dt, dx, dkx, kx, dz, dkz, kz, tmpdt, tmp, pi=SF_PI, o1, o2, kx0, kz0;
    float **new,  **old,  **cur, **ukr, **dercur, **derold, *wav, **rvr, **snap, **image;
    float **vx, vx2, vx0, vx02, **vz, vz2, vz0, vz02, **yi, yi0, **se, se0;
    float ***aa, dx2, dz2, dx4, dz4, ct, cb, cl, cr; /* top, bottom, left, right */
    float w1, w10, w2, w20, w3, w30, h1, h10, h2, h20, h3, h30;
    float cosg, cosg0, cosg2, cosg02, sing, sing0, sing2, sing02;
    float vk, vk2, tmpvk, tmpk, k2, err, dt2, kx1, kz1;
    float alpha; /* source smoothing */
    kiss_fft_cpx **uk, **ctracex, **ctracez;
    sf_file input, velx, velz, source, yita, seta, geo, output;
    FILE *out;
    bool opt,de;    /* optimal padding */
    float **fcos;
    int nth=1, ith=0, esize, shot_num, n1;
    int i, rank, nodes;
    int nr, jr, r0, nl, tl, jm, rb;
    char *oname, *mm, *iname, *sname;
    float ax, az, factor;
    kiss_fft_cfg *cfgx, *cfgxi, *cfgz, *cfgzi;
    /* MPI_Status stat; */
     
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nodes);
    if(rank==0) sf_warning("nodes=%d",nodes);
/*
  if (nodes < 2) {
  fprintf(stderr,"Need at least two nodes!\n");
  MPI_Finalize();
  }
*/

    sf_init(argc,argv);
    geo = sf_input("geo");   /* geometry: source location isx, receiver starting & number: r0 & nl*/
    velx = sf_input("velx");   /* velocity */
    velz = sf_input("velz");   /* velocity */
    yita = sf_input("yita");   /* anistropic parameter*/
    source = sf_input("source");   /* source wavlet*/
    seta = sf_input("seta");   /* TTI angle*/

    if (rank == 0){
	if (SF_FLOAT != sf_gettype(velx)) sf_error("Need float input");
	if (SF_FLOAT != sf_gettype(velz)) sf_error("Need float input");
	if (SF_FLOAT != sf_gettype(source)) sf_error("Need float input");
	if (SF_FLOAT != sf_gettype(seta)) sf_error("Need float input");
    }

    if (!sf_histint(velx,"n1",&nx)) sf_error("No n1= in input");
    if (!sf_histfloat(velx,"d1",&dx)) sf_error("No d1= in input");
    if (!sf_histint(velx,"n2",&nz)) sf_error("No n2= in input");
    if (!sf_histfloat(velx,"d2",&dz)) sf_error("No d2= in input");
    if (!sf_histfloat(velx,"o1",&o1)) o1=0.0;
    if (!sf_histfloat(velx,"o2",&o2)) o2=0.0;
    if (!sf_getbool("opt",&opt)) opt=true;
    if (!sf_getbool("de",&de)) de=true;
    /* if y, determine optimal size for efficiency */
    if (!sf_getfloat("dt",&dt)) sf_error("Need dt input");
    if (!sf_getint("nt",&nt)) sf_error("Need nt input");
    if (!sf_getint("isz",&isz)) sf_error("Need isz input");
    if (!sf_getint("irz",&irz)) irz=isz;
    /* if (!sf_getint("r0",&r0)) r0=0; */
    if (!sf_getint("jr",&jr)) jr=1;
    if (!sf_getint("jm",&jm)) jm=20;
    if (!sf_getint("nr",&nr)) sf_error("Need nr input");/*streamer total length*/
    if (!sf_getfloat("err",&err)) err = 0.00001;
    if (!sf_getfloat("alpha",&alpha)) alpha=-0.7;

    if (!sf_getint("nbt",&nbt)) nbt=44;
    if (!sf_getint("nbb",&nbb)) nbb=44;
    if (!sf_getint("nbl",&nbl)) nbl=44;
    if (!sf_getint("nbr",&nbr)) nbr=44;

    if (!sf_getfloat("ct",&ct)) ct = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cb",&cb)) cb = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cl",&cl)) cl = 0.01; /*decaying parameter*/
    if (!sf_getfloat("cr",&cr)) cr = 0.01; /*decaying parameter*/

    if (!sf_getfloat("ax",&ax)) ax= 2.0; /*suppress HF parameter*/
    if (!sf_getfloat("az",&az)) az= 2.0; /*suppress HF parameter*/
    if (!sf_getfloat("factor",&factor)) factor= 2.0/3.0; /*suppress HF parameter*/
    
    if(SF_INT != sf_gettype(geo)) sf_error("Need int!");
    if(!sf_histint(geo,"esize",&esize)) esize = sizeof("int");
    /* if(!sf_histint(geo,"n1",&shot_num)) sf_error("No n1=!"); */
    if(!sf_histint(geo,"n1",&n1)) sf_error("No n1=!");
    if(!(n1==4)) sf_error("n1 in geo should be 4");
    if(!sf_histint(geo,"n2",&shot_num)) sf_error("No n2=!");
    if(rank==0) sf_warning("%d shots!",shot_num);

    nxb = nx + nbl + nbr;
    nzb = nz + nbt + nbb;
    tl = nx*nz;

    nkx = opt? kiss_fft_next_fast_size(nxb): nxb;
    nkz = opt? kiss_fft_next_fast_size(nzb): nzb;
    if (nkx != nxb) sf_warning("nkx padded to %d",nkx);
    if (nkz != nzb) sf_warning("nkz padded to %d",nkz);
    dkx = 1./(nkx*dx)*2.0*pi;
    kx0 = -0.5/dx*2.0*pi;
    dkz = 1./(nkz*dz)*2.0*pi;
    kz0 = -0.5/dz*2.0*pi;

#ifdef _OPENMP
#pragma omp parallel
    {nth = omp_get_num_threads();}
    if(rank==0)  sf_warning("using %d threads",nth);
#endif
    uk = (kiss_fft_cpx **) sf_complexalloc2(nkx,nkz);
    ctracex = (kiss_fft_cpx **) sf_complexalloc2(nkx,nth);
    ctracez = (kiss_fft_cpx **) sf_complexalloc2(nkz,nth);

    cfgx  = (kiss_fft_cfg *) sf_complexalloc(nth);
    cfgxi = (kiss_fft_cfg *) sf_complexalloc(nth);
    cfgz  = (kiss_fft_cfg *) sf_complexalloc(nth);
    cfgzi = (kiss_fft_cfg *) sf_complexalloc(nth);

    for (i=0; i < nth; i++) {
        cfgx[i] = kiss_fft_alloc(nkx,0,NULL,NULL);
        cfgxi[i] = kiss_fft_alloc(nkx,1,NULL,NULL);
        cfgz[i] = kiss_fft_alloc(nkz,0,NULL,NULL);
        cfgzi[i]= kiss_fft_alloc(nkz,1,NULL,NULL);
    }
   


    wav    =  sf_floatalloc(nt);
    rvr    =  sf_floatalloc2(nt,nr);
    sf_floatread(wav,nt,source);

    snap   =  sf_floatalloc2(nx,nz);
    image  =  sf_floatalloc2(nx,nz);
    old    =  sf_floatalloc2(nxb,nzb);
    cur    =  sf_floatalloc2(nxb,nzb);
    new    =  sf_floatalloc2(nxb,nzb);
    ukr    =  sf_floatalloc2(nxb,nzb);
    fcos   =  sf_floatalloc2(nkx,nkz);
    derold =  sf_floatalloc2(nxb,nzb);
    dercur =  sf_floatalloc2(nxb,nzb);
    aa     =  sf_floatalloc3(6,nxb,nzb);
    
    bd_init(nx,nz,nbt,nbb,nbl,nbr,ct,cb,cl,cr);

    vx = sf_floatalloc2(nxb,nzb);
    vz = sf_floatalloc2(nxb,nzb);

    /*input & extend velocity model*/
    for (iz=nbt; iz<nz+nbt; iz++){
        sf_floatread(vx[iz]+nbl,nx,velx);
        sf_floatread(vz[iz]+nbl,nx,velz);
	for (ix=0; ix<nbl; ix++){
	    vx[iz][ix] = vx[iz][nbl];
	    vz[iz][ix] = vz[iz][nbl];
	}
	for (ix=0; ix<nbr; ix++){
	    vx[iz][nx+nbl+ix] = vx[iz][nx+nbl-1];
	    vz[iz][nx+nbl+ix] = vz[iz][nx+nbl-1];
	}     
    }
    for (iz=0; iz<nbt; iz++){
        for (ix=0; ix<nxb; ix++){
            vx[iz][ix] = vx[nbt][ix];
            vz[iz][ix] = vz[nbt][ix];
        }
    }
    for (iz=0; iz<nbb; iz++){
        for (ix=0; ix<nxb; ix++){
            vx[nz+nbt+iz][ix] = vx[nz+nbt-1][ix];
            vz[nz+nbt+iz][ix] = vz[nz+nbt-1][ix];
        }
    }
    vx0 =0.0;
    vz0 =0.0;
    for (iz=0; iz < nzb; iz++) {
        for (ix=0; ix < nxb; ix++) {
            vx0 += vx[iz][ix]*vx[iz][ix];
            vz0 += vz[iz][ix]*vz[iz][ix];
	}
    }
    vx0 = sqrtf(vx0/(nxb*nzb));
    vz0 = sqrtf(vz0/(nxb*nzb));
/*
  vz0 =5000.0;
  for (iz=0; iz < nzb; iz++) {
  for (ix=0; ix < nxb; ix++) {
  if(vz0 > vz[iz][ix]) vz0 = vz[iz][ix];
  }
  }
*/
/*    vx0 = vz0; */
    vx02=vx0*vx0; 
    vz02=vz0*vz0; 

    /*input & extend anistropic model*/
    yi = sf_floatalloc2(nxb,nzb);
    for (iz=nbt; iz<nz+nbt; iz++){
        sf_floatread(yi[iz]+nbl,nx,yita);
	for (ix=0; ix<nbl; ix++){
	    yi[iz][ix] = yi[iz][nbl];
	}
	for (ix=0; ix<nbr; ix++){
	    yi[iz][nx+nbl+ix] = yi[iz][nx+nbl-1];
	}     
    }
    for (iz=0; iz<nbt; iz++){
        for (ix=0; ix<nxb; ix++){
            yi[iz][ix] = yi[nbt][ix];
        }
    }
    for (iz=0; iz<nbb; iz++){
        for (ix=0; ix<nxb; ix++){
            yi[nz+nbt+iz][ix] = yi[nz+nbt-1][ix];
        }
    }

    yi0 = 0.0;
    for (iz=0; iz < nzb; iz++) {
        for (ix=0; ix < nxb; ix++) {
            yi0+= yi[iz][ix]*yi[iz][ix];
	}
    }
    yi0 = sqrtf(yi0/(nxb*nzb));
    se = sf_floatalloc2(nxb,nzb);
    for (iz=nbt; iz<nz+nbt; iz++){
        sf_floatread(se[iz]+nbl,nx,seta);
	for (ix=0; ix<nbl; ix++){
	    se[iz][ix] = se[iz][nbl];
	}
	for (ix=0; ix<nbr; ix++){
	    se[iz][nx+nbl+ix] = se[iz][nx+nbl-1];
	}     
    }
    for (iz=0; iz<nbt; iz++){
        for (ix=0; ix<nxb; ix++){
            se[iz][ix] = se[nbt][ix];
        }
    }
    for (iz=0; iz<nbb; iz++){
        for (ix=0; ix<nxb; ix++){
            se[nz+nbt+iz][ix] = se[nz+nbt-1][ix];
        }
    }

    se0 = 0.0;
    for (iz=0; iz < nzb; iz++) {
        for (ix=0; ix < nxb; ix++) {
            se0+= se[iz][ix];
	}
    }
    se0 /= (nxb*nzb);
    if(de) {se0 *= pi/180.0; sf_warning("degree!");}
    if(de){
	for (iz=0; iz < nzb; iz++) {
	    for (ix=0; ix < nxb; ix++) {
		se[iz][ix] *= pi/180.0;
	    }
	}
    }

    cosg0 = cosf(se0);
    cosg02 = cosg0*cosg0;
    sing0 = sinf(se0);
    sing02 = sing0*sing0; 

    w10 = vx02*cosg02+vz02*sing02;
    w20 = vz02*cosg02+vx02*sing02;
    w30 = vx02+vz02+(vx02-vz02)*sinf(2.0*se0);
    h10 = sqrtf(-8.0*yi0*vx02*vz02*cosg02*sing02/(1.0+2.0*yi0)+w10*w10);
    h20 = sqrtf(-8.0*yi0*vx02*vz02*cosg02*sing02/(1.0+2.0*yi0)+w20*w20);
    h30 = sqrtf(-2.0*yi0*vx02*vz02*cosf(2.0*se0)*cosf(2.0*se0)/(1.0+2.0*yi0)+0.25*w30*w30);
    dt2 = dt*dt;
    dx2 = dx*dx;
    dx4 = dx2*dx2;
    dz2 = dz*dz;
    dz4 = dz2*dz2;

    for (iz=0; iz < nzb; iz++){
        for (ix=0; ix < nxb; ix++) {
            vx2 = vx[iz][ix]*vx[iz][ix];
/*            vx4 = vx2*vx2; */
            vz2 = vz[iz][ix]*vz[iz][ix];
/*            vz4 = vz2*vz2; */
            cosg = cosf(se[iz][ix]);
            sing = sinf(se[iz][ix]);
            cosg2 = cosg*cosg;
            sing2 = sing*sing;
            w1 = vx2*cosg2+vz2*sing2;
            w2 = vz2*cosg2+vx2*sing2;
            w3 = vx2+vz2+(vx2-vz2)*sinf(2.0*se[iz][ix]);
            h1 = sqrtf(-8.0*yi[iz][ix]*vx2*vz2*cosg2*sing2/(1.0+2.0*yi[iz][ix])+w1*w1);
            h2 = sqrtf(-8.0*yi[iz][ix]*vx2*vz2*cosg2*sing2/(1.0+2.0*yi[iz][ix])+w2*w2);
            h3 = sqrtf(-2.0*yi[iz][ix]*vx2*vz2*cosf(2.0*se[iz][ix])*cosf(2.0*se[iz][ix])/(1.0+2.0*yi[iz][ix])+0.25*w3*w3);
            aa[iz][ix][4] = (w1+h1)*(dt2+(2.0*dx2-dt2*(w1+h1))/(w10+h10))/(24.0*dx4);
            aa[iz][ix][5] = (w2+h2)*(dt2+(2.0*dz2-dt2*(w2+h2))/(w20+h20))/(24.0*dz4);
            aa[iz][ix][3] = -aa[iz][ix][4]*dx2/dz2-aa[iz][ix][5]*dz2/dx2+(dt2*(w3+2.0*h3)+dx2*(w1+h1)/(w10+h10)+dz2*(w2+h2)/(w20+h20)-dt2*(w3+2.0*h3)*(w3+2.0*h3)/(w30+2.0*h30))/(12.0*dx2*dz2);
            aa[iz][ix][1] = -2.0*aa[iz][ix][3]-4.0*aa[iz][ix][4]-(w1+h1)/(dx2*(w10+h10));
            aa[iz][ix][2] = -2.0*aa[iz][ix][3]-4.0*aa[iz][ix][5]-(w2+h2)/(dz2*(w20+h20));
            aa[iz][ix][0] = -2.0*aa[iz][ix][1]-2.0*aa[iz][ix][2]-4.0*aa[iz][ix][3]-2.0*aa[iz][ix][4]-2.0*aa[iz][ix][5];
        }
    }
    free(*vx);     
    free(*vz);     
    free(*yi);     
    free(*se);     
    free(vx);     
    free(vz);     
    free(yi);     
    free(se);     
    sf_fileclose(velz);
    sf_fileclose(velx);
    sf_fileclose(yita);
    sf_fileclose(seta);
    sf_fileclose(source);
    for (ikz=0; ikz < nkz; ikz++) {
        /* kz1 = (kz0+ikz*dkz)*2.0*pi; */
        kz1 = (kz0+ikz*dkz);
        for (ikx=0; ikx < nkx; ikx++) {
            /* kx1 = (kx0+ikx*dkx)*2.0*pi; */
            kx1 = (kx0+ikx*dkx);
            kx = kx1*cosg0+kz1*sing0;
            kz = kz1*cosg0-kx1*sing0;
            tmpvk = (vx02*kx*kx+vz02*kz*kz);
            k2 = kx1*kx1+kz1*kz1;
            vk2 = 0.5*tmpvk+0.5*sqrtf(tmpvk*tmpvk-8.0*yi0/(1.0+2.0*yi0)*vx02*vz02*kx*kx*kz*kz);
            /* vk2 = tmpvk; */
            vk = sqrtf(vk2);
            tmpk = vk*dt;
            tmp = vz0*dt;
/*            tmpdt = 2.0*(cosf(tmpk)-1.0); */
            /* if(k2 < err) */
            if(k2==0 || tmpk < err) 
		/* tmpdt /=(k2+err); */
		tmpdt = -(tmp)*(tmp);
/*
  else if(tmpk < 0.1) 
   tmpdt =(-tmpk*tmpk+tmpk*tmpk*tmpk*tmpk/12.0)/(k2); 
	    tmpdt =(-tmp*tmp+tmp*tmp*tmp*tmp*k2/12.0);
	    */
            else
		tmpdt = 2.0*(cosf(tmpk)-1.0)/k2;
	    /* tmpdt = -2.0*(sinf(tmpk)*sinf(tmpk))/((cosf(tmpk)+1.0)*k2); */
/*              tmpdt /= k2; */
            /* fcos[ikz][ikx] = tmpdt*dehf(kx1,knx,ax,factor)*dehf(kz1,knz,az,factor); */
            fcos[ikz][ikx] = tmpdt;
        }
    }

    sf_seek(geo,rank*esize*4,SEEK_SET);
    srcsm_init(dz,dx);

    sf_warning("Rank= %d",rank);
    for (i=rank; i < shot_num; i+=nodes) {
        sf_warning("Shot No. %d",i);
        /* out = fopen("/local/lfs1/data/snap","w"); */
	/* out = fopen("/tmp/snap","w"); */
        sf_intread(&isx,1,geo);
        sf_intread(&r0,1,geo);
        sf_intread(&nl,1,geo);
        sf_intread(&rb,1,geo);
        sf_seek(geo,(nodes-1)*esize*4,SEEK_CUR);
        mm = sf_charalloc(15);
        sname = sf_charalloc(20);
        iname = sf_charalloc(20);
        oname = sf_charalloc(25);
        itoa(isx,mm);
        sname[0] = '/';
        sname[1] = 't';
        sname[2] = 'm';
        sname[3] = 'p';
        sname[4] = '/';
        sname[5] = 's';
        sname[6] = 'n';
        sname[7] = 'a';
        sname[8] = 'p';
        sname[9] = '_';
        sname[10] = '\0';
        sname = strcat(sname,mm);
        out = fopen(sname,"w");

        mm = strcat(mm,".rsf");
        iname[0] = 'D';
        iname[1] = 'A';
        iname[2] = 'T';
        iname[3] = 'A';
        iname[4] = '/';
        iname[5] = 's';
        /* iname[5] = 'd'; */
        iname[6] = 'h';
        iname[7] = 'o';
        iname[8] = 't';
        iname[9] = '_';
        iname[10] = '\0';
        iname = strcat(iname,mm);
        oname[0] = 'I';
        oname[1] = 'M';
        oname[2] = 'A';
        oname[3] = 'G';
/*
  oname[0] = 'I';
  oname[1] = 'M';
  oname[2] = 'G';
  oname[3] = 'S';
*/
        oname[4] = '/';
        oname[5] = 'I';
        oname[6] = 'm';
        oname[7] = 'a';
        oname[8] = 'g';
        oname[9] = '_';
        oname[10] = '\0';
        oname = strcat(oname,mm);
        free(mm);
        sf_warning("%s",oname);
        input = sf_input(iname);
        free(iname);
        output = sf_output(oname);
        free(oname);
        sf_putint(output,"n1",nx);
        sf_putfloat(output,"d1",dx);
        sf_putint(output,"n2",nz);
        sf_putfloat(output,"d2",dz);
        sf_settype(output,SF_FLOAT);
        for (iz=0; iz < nz; iz++) {
            for (ix=0; ix < nx; ix++) {
                image[iz][ix] = 0.0;
            }
        }
         
        
        /* Initialization */
        for (iz=0; iz < nzb; iz++) {
            for (ix=0; ix < nxb; ix++) {
                cur[iz][ix] = 0.0;
            }
        }
        for (iz=0; iz < nzb; iz++) {
            for (ix=0; ix < nxb; ix++) {
                old[iz][ix] =  0.0; 
                derold[iz][ix] =cur[iz][ix]/dt;
	    }
	}
 
        /* propagation in time */
        for (it=0; it < nt; it++) {
            if(it<1500 )  { 
		cur[isz+nbt][isx+nbl] += wav[it];
		source_smooth(cur,isz+nbt,isx+nbl,wav[it]);
	    }
            if(!(it%jm)) {
		/* sf_floatwrite(cur[nbt],nxb*nz,out); */
		for(iz=0;iz < nz;iz++) fwrite(cur[nbt+iz]+nbl,sizeof(float),nx,out);
            }

#ifdef _OPENMP
#pragma omp parallel for private(iz,ix) 
#endif
	    for (iz=0; iz < nzb; iz++){
		for (ix=0; ix < nxb; ix++){ 
		    new[iz][ix] = 0.0; 
		    uk[iz][ix].r = cur[iz][ix];
		    uk[iz][ix].i = 0.0; 
		}
	    }  


/*      compute u(kx,kz) */
#ifdef _OPENMP
#pragma omp parallel for private(iz,ix,ikx,ith) 
#endif
	    for (iz=0; iz < nzb; iz++){
		/* Fourier transform x to kx */
		for (ix=1; ix < nxb; ix+=2){
		    uk[iz][ix] = sf_cneg(uk[iz][ix]);
		}
#ifdef _OPENMP
		ith = omp_get_thread_num();
#endif
		kiss_fft_stride(cfgx[ith],uk[iz],ctracex[ith],1); 
		for (ikx=0; ikx<nkx; ikx++) uk[iz][ikx] = ctracex[ith][ikx]; 
	    }
#ifdef _OPENMP
#pragma omp parallel for private(ikz,ikx,ith) 
#endif
	    for (ikx=0; ikx < nkx; ikx++){
		/* Fourier transform z to kz */
		for (ikz=1; ikz<nkz; ikz+=2){
		    uk[ikz][ikx] = sf_cneg(uk[ikz][ikx]);
		}
#ifdef _OPENMP
		ith = omp_get_thread_num();
#endif
		kiss_fft_stride(cfgz[ith],uk[0]+ikx,ctracez[ith],nkx); 
		for (ikz=0; ikz<nkz; ikz++) uk[ikz][ikx] = ctracez[ith][ikz]; 
	    }

#ifdef _OPENMP
#pragma omp parallel for private(ikz,ikx,tmpdt) 
#endif
	    for (ikz=0; ikz < nkz; ikz++) {
		for (ikx=0; ikx < nkx; ikx++) {
		    tmpdt = fcos[ikz][ikx];
		    uk[ikz][ikx] = sf_crmul(uk[ikz][ikx],tmpdt);
		}
	    }   
/*      Inverse FFT*/
#ifdef _OPENMP
#pragma omp parallel for private(ikz,ikx,ith) 
#endif
	    for (ikx=0; ikx < nkx; ikx++){
		/* Inverse Fourier transform kz to z */
#ifdef _OPENMP
		ith = omp_get_thread_num();
#endif
		kiss_fft_stride(cfgzi[ith],(kiss_fft_cpx *)uk[0]+ikx,ctracez[ith],nkx); 
		for (ikz=0; ikz < nkz; ikz++) uk[ikz][ikx] = sf_crmul(ctracez[ith][ikz],ikz%2?-1.0:1.0); 
	    }
#ifdef _OPENMP
#pragma omp parallel for private(ikz,ikx,ith) 
#endif
	    for (ikz=0; ikz < nkz; ikz++){
		/* Inverse Fourier transform kx to x */
#ifdef _OPENMP
		ith = omp_get_thread_num();
#endif
		kiss_fft_stride(cfgxi[ith],(kiss_fft_cpx *)uk[ikz],ctracex[ith],1); 
		for (ikx=0; ikx < nkx; ikx++) uk[ikz][ikx] = sf_crmul(ctracex[ith][ikx],ikx%2?-1.0:1.0); 
	    }
#ifdef _OPENMP
#pragma omp parallel for private(iz,ix) 
#endif
	    for (iz=0; iz < nzb; iz++){
		for (ix=0; ix < nxb; ix++){ 
		    ukr[iz][ix] = sf_crealf(uk[iz][ix]); 
		    ukr[iz][ix] /= (nkx*nkz); 
		}
	    }  

#ifdef _OPENMP
#pragma omp parallel for private(iz,ix) 
#endif
	    for (iz=2; iz < nzb-2; iz++) {  
		for (ix=2; ix < nxb-2; ix++) {  
		    new[iz][ix]  = ukr[iz][ix]*aa[iz][ix][0]
			+ (ukr[iz][ix-1]+ukr[iz][ix+1])*aa[iz][ix][1]
			+ (ukr[iz-1][ix]+ukr[iz+1][ix])*aa[iz][ix][2]
			+ (ukr[iz-1][ix-1]+ukr[iz-1][ix+1]+ukr[iz+1][ix-1]+ukr[iz+1][ix+1])*aa[iz][ix][3]
			+ (ukr[iz][ix-2]+ukr[iz][ix+2])*aa[iz][ix][4]
			+ (ukr[iz-2][ix]+ukr[iz+2][ix])*aa[iz][ix][5];
		}
	    }  
             

#ifdef _OPENMP
#pragma omp parallel for private(iz,ix) 
#endif
	    for (iz=0; iz < nzb; iz++) {  
		for (ix=0; ix < nxb; ix++) {
		    dercur[iz][ix]= derold[iz][ix] + new[iz][ix]/dt;
		    new[iz][ix] = cur[iz][ix] + dercur[iz][ix]*dt; 
		    /*     new[iz][ix] += 2.0*cur[iz][ix] -old[iz][ix]; */
		}
	    }
 
	    bd_decay(new); 
	    bd_decay(dercur); 
 
#ifdef _OPENMP
#pragma omp parallel for private(iz,ix) 
#endif
	    for (iz=0; iz < nzb; iz++) {  
		for(ix=0; ix < nxb; ix++) {
                    old[iz][ix] = cur[iz][ix]; 
                    cur[iz][ix] = new[iz][ix]; 
                    derold[iz][ix] = dercur[iz][ix]; 
		}
	    }
        }
        fclose(out);
        /* out = fopen("/local/lfs1/data/snap","r"); */
        /* out = fopen("/tmp/snap","r"); */
        out = fopen(sname,"r");
        for (iz=0; iz < nzb; iz++) {
            for (ix=0; ix < nxb; ix++) {
                cur[iz][ix] = 0.0;
                old[iz][ix] =  0.0;
                derold[iz][ix] = 0.0;
            }
        }
        sf_floatread(rvr[0],nr*nt,input);
        sf_fileclose(input);
        /* for (it=nt-1; it >-1; it--) { */
        for (it=nt-1; it >194; it--) {
            /* sf_floatread(rvr,nr,input); */
	    /*     for (ix=r0; ix < nx; ix+=jr) cur[irz+nbt][ix+nbl] += rvr[ix][it]; */
            for (ir=0; ir < nl; ir++) {
                ix = r0 + ir*jr;
                cur[irz+nbt][ix+nbl] += rvr[rb+ir][it];
            }
            for (iz=0; iz < nzb; iz++) {
                for (ix=0; ix < nxb; ix++) {
                    new[iz][ix] = 0.0; 
                    uk[iz][ix].r = cur[iz][ix];
                    uk[iz][ix].i = 0.0; 
                }
            }
            if(!(it%jm)) {
		fseek(out,sizeof(float)*tl*(it/jm),SEEK_SET);
		if (tl != fread(snap[0],sizeof(float),tl,out))
		    sf_error("fread error:");
		for (iz=0; iz < nz; iz++) {
		    for(ix=0; ix < nx; ix++) {
			image[iz][ix] += snap[iz][ix]*cur[iz+nbt][ix+nbl];
		    }
		}
            }
/*      compute u(kx,kz) */
#ifdef _OPENMP
#pragma omp parallel for private(iz,ix,ikx,ith) 
#endif
	    for (iz=0; iz < nzb; iz++){
		/* Fourier transform x to kx */
		for (ix=1; ix < nxb; ix+=2){
		    uk[iz][ix] = sf_cneg(uk[iz][ix]);
		}
#ifdef _OPENMP
		ith = omp_get_thread_num();
#endif
		kiss_fft_stride(cfgx[ith],uk[iz],ctracex[ith],1); 
		for (ikx=0; ikx<nkx; ikx++) uk[iz][ikx] = ctracex[ith][ikx]; 
	    }
#ifdef _OPENMP
#pragma omp parallel for private(ikz,ikx,ith) 
#endif
	    for (ikx=0; ikx < nkx; ikx++){
		/* Fourier transform z to kz */
		for (ikz=1; ikz<nkz; ikz+=2){
		    uk[ikz][ikx] = sf_cneg(uk[ikz][ikx]);
		}
#ifdef _OPENMP
		ith = omp_get_thread_num();
#endif
		kiss_fft_stride(cfgz[ith],uk[0]+ikx,ctracez[ith],nkx); 
		for (ikz=0; ikz<nkz; ikz++) uk[ikz][ikx] = ctracez[ith][ikz]; 
	    }

#ifdef _OPENMP
#pragma omp parallel for private(ikz,ikx,tmpdt) 
#endif
	    for (ikz=0; ikz < nkz; ikz++) {
		for (ikx=0; ikx < nkx; ikx++) {
		    tmpdt = fcos[ikz][ikx];
		    uk[ikz][ikx] = sf_crmul(uk[ikz][ikx],tmpdt);
		}

	    }   
/*      Inverse FFT*/
#ifdef _OPENMP
#pragma omp parallel for private(ikz,ikx,ith) 
#endif
	    for (ikx=0; ikx < nkx; ikx++){
		/* Inverse Fourier transform kz to z */
#ifdef _OPENMP
		ith = omp_get_thread_num();
#endif
		kiss_fft_stride(cfgzi[ith],(kiss_fft_cpx *)uk[0]+ikx,ctracez[ith],nkx); 
		for (ikz=0; ikz < nkz; ikz++) uk[ikz][ikx] = sf_crmul(ctracez[ith][ikz],ikz%2?-1.0:1.0); 
	    }
#ifdef _OPENMP
#pragma omp parallel for private(ikz,ikx,ith) 
#endif
	    for (ikz=0; ikz < nkz; ikz++){
		/* Inverse Fourier transform kx to x */
#ifdef _OPENMP
		ith = omp_get_thread_num();
#endif
		kiss_fft_stride(cfgxi[ith],(kiss_fft_cpx *)uk[ikz],ctracex[ith],1); 
		for (ikx=0; ikx < nkx; ikx++) uk[ikz][ikx] = sf_crmul(ctracex[ith][ikx],ikx%2?-1.0:1.0); 
	    }
#ifdef _OPENMP
#pragma omp parallel for private(iz,ix) 
#endif
	    for (iz=0; iz < nzb; iz++){
		for (ix=0; ix < nxb; ix++){ 
		    ukr[iz][ix] = sf_crealf(uk[iz][ix]); 
		    ukr[iz][ix] /= (nkx*nkz); 
		}
	    }  

#ifdef _OPENMP
#pragma omp parallel for private(iz,ix) 
#endif
	    for (iz=2; iz < nzb-2; iz++) {  
		for (ix=2; ix < nxb-2; ix++) {  
		    new[iz][ix]  = ukr[iz][ix]*aa[iz][ix][0]
			+ (ukr[iz][ix-1]+ukr[iz][ix+1])*aa[iz][ix][1]
			+ (ukr[iz-1][ix]+ukr[iz+1][ix])*aa[iz][ix][2]
			+ (ukr[iz-1][ix-1]+ukr[iz-1][ix+1]+ukr[iz+1][ix-1]+ukr[iz+1][ix+1])*aa[iz][ix][3]
			+ (ukr[iz][ix-2]+ukr[iz][ix+2])*aa[iz][ix][4]
			+ (ukr[iz-2][ix]+ukr[iz+2][ix])*aa[iz][ix][5];
		}
	    }  
             

#ifdef _OPENMP
#pragma omp parallel for private(iz,ix) 
#endif
	    for (iz=0; iz < nzb; iz++) {  
		for (ix=0; ix < nxb; ix++) {
		    dercur[iz][ix]= derold[iz][ix] + new[iz][ix]/dt;
		    new[iz][ix] = cur[iz][ix] + dercur[iz][ix]*dt; 
		    /*     new[iz][ix] += 2.0*cur[iz][ix] -old[iz][ix]; */
		}
	    }
 
	    bd_decay(new); 
	    bd_decay(dercur); 
 
#ifdef _OPENMP
#pragma omp parallel for private(iz,ix) 
#endif
	    for (iz=0; iz < nzb; iz++) {  
		for(ix=0; ix < nxb; ix++) {
                    old[iz][ix] = cur[iz][ix]; 
                    cur[iz][ix] = new[iz][ix]; 
                    derold[iz][ix] = dercur[iz][ix]; 
		}
	    }
        }
        sf_floatwrite(image[0],tl,output);
        sf_warning("image output");
        sf_fileclose(output);
        fclose(out);
        remove(sname);
        free(sname);
    }

    /* remove("/local/lfs1/data/snap"); */
    /* remove("/tmp/snap"); */
    bd_close();
    srcsm_close();
    free(**aa);
    free(*aa);
    free(aa);
    free(*new);     
    free(*cur);     
    free(*old);     
    free(*dercur);     
    free(*derold);     
    free(*uk);     
    free(*ukr);     
    free(new);     
    free(cur);     
    free(old);     
    free(dercur);     
    free(derold);     
    free(uk);     
    free(ukr);     
    MPI_Finalize();

    exit(0); 
}           



void itoa(int n, char *s)
{
    int i,j,sign;
    char c;
    sign = n;
    if (n < 0) n=-n;
    i = 0;
    do {
	s[i++] = n%10+'0';
	n = (int) n/10;
    } while(n > 0);
    if (sign <0) s[i++] = '-';
    s[i] = '\0';
    for (i=0,j=strlen(s)-1;i<j;i++,j--){
        c = s[i];
        s[i]=s[j];
        s[j]=c;
    }
}

void srcsm_init(float dz, float dx/*grid size*/)
/*< initialization >*/
{

    float dx2, dz2, R0, rxxz, rxzz;
    dx2 = dx*dx;
    dz2 = dz*dz;
    R0 = sqrtf(dx2+dz2);
    rxxz = sqrtf(4.0*dx2+dz2);
    rxzz = sqrtf(4.0*dz2+dx2);
    ws =  sf_floatalloc2(3,2);
    ws[0][0] = -0.5*dz2/(R0*R0)+1.0;
    if(2.0*dz < R0) { 
	ws[0][1]= -2.0*dz2/(R0*R0)+1.0;
    } else {
	ws[0][1]= 2.0*(dz-R0)*(dz-R0)/(R0*R0);
    }
    ws[0][2]= 0.5*(rxzz-2.0*R0)*(rxzz-2.0*R0)/(R0*R0);
    ws[1][0] = -0.5*dx2/(R0*R0)+1.0;
    if(2.0*dx < R0) { 
	ws[1][1] = -2.0*dx2/(R0*R0)+1.0;
    } else {
	ws[1][1] = 2.0*(dx-R0)*(dx-R0)/(R0*R0);
    }
    ws[1][2]= 0.5*(rxxz-2.0*R0)*(rxxz-2.0*R0)/(R0*R0);



}
   

void srcsm_close(void)
/*< free memory allocation>*/
{
    free(*ws);
    free(ws);

}



void source_smooth(float **source /*source matrix*/, 
                   int iz /*source index*/, 
                   int ix /*source index*/,
                   float val/*source value*/)
/*< smooth source input >*/
{
   
    source[iz-1][ix]   += val*ws[0][0]; 
    source[iz+1][ix]   += val*ws[0][0]; 

    source[iz+2][ix]   += val*ws[0][1];
    source[iz-2][ix]   += val*ws[0][1];

    source[iz+2][ix+1]   += val*ws[0][2]; 
    source[iz+2][ix-1]   += val*ws[0][2]; 
    source[iz-2][ix-1]   += val*ws[0][2]; 
    source[iz-2][ix+1]   += val*ws[0][2]; 

    source[iz-1][ix-1] += val*0.5; 
    source[iz-1][ix+1] += val*0.5; 
    source[iz+1][ix-1] += val*0.5; 
    source[iz+1][ix+1] += val*0.5; 

    source[iz][ix-1]   += val*ws[1][0];
    source[iz][ix+1]   += val*ws[1][0];

    source[iz][ix+2]   += val*ws[1][1];
    source[iz][ix-2]   += val*ws[1][1];

    
    source[iz+1][ix+2]   += val*ws[1][2]; 
    source[iz+1][ix-2]   += val*ws[1][2]; 
    source[iz-1][ix-2]   += val*ws[1][2]; 
    source[iz-1][ix+2]   += val*ws[1][2]; 

}



void bd_init(int n1,  int n2    /*model size:x*/,
             int nb1, int nb2  /*top, bottom*/,
             int nb3, int nb4  /*left, right*/,
             float c1, float c2 /*top, bottom*/,
             float c3, float c4 /*left, right*/)
/*< initialization >*/
{
    int c;
    nx = n1;  
    nz = n2;  
    nbt = nb1;  
    nbb = nb2;  
    nbl = nb3;  
    nbr = nb4;  
    ct = c1;  
    cb = c2;  
    cl = c3;  
    cr = c4;  
    if(nbt) wt =  sf_floatalloc(nbt);
    if(nbb) wb =  sf_floatalloc(nbb);
    if(nbl) wl =  sf_floatalloc(nbl);
    if(nbr) wr =  sf_floatalloc(nbr);
    c=0;
    abc_cal(c,nbt,ct,wt);
    abc_cal(c,nbb,cb,wb);
    abc_cal(c,nbl,cl,wl);
    abc_cal(c,nbr,cr,wr);
}
   

void bd_close(void)
/*< free memory allocation>*/
{
    if(nbt) free(wt);
    if(nbb) free(wb);
    if(nbl) free(wl);
    if(nbr) free(wr);
}

void bd_decay(float **a /*2-D matrix*/) 
/*< boundary decay>*/
{
    int iz, ix;
    for (iz=0; iz < nbt; iz++) {  
        for (ix=nbl; ix < nx+nbl; ix++) {
            a[iz][ix] *= wt[iz];
        }
    }
    for (iz=0; iz < nbb; iz++) {  
        for (ix=nbl; ix < nx+nbl; ix++) {
            a[iz+nz+nbt][ix] *= wb[nbb-1-iz];
        }
    }
    for (iz=nbt; iz < nz+nbt; iz++) {  
        for (ix=0; ix < nbl; ix++) {
            a[iz][ix] *= wl[ix];
        }
    }
    for (iz=nbt; iz < nz+nbt; iz++) {  
        for (ix=0; ix < nbr; ix++) {
            a[iz][nx+nbl+ix] *= wr[nbr-1-ix];
        }
    }
    for (iz=0; iz < nbt; iz++) {  
        for (ix=0; ix < nbl; ix++) {
            /* a[iz][ix] *= (float)iz/(float)(ix+iz)*wl[ix]+(float)ix/(float)(ix+iz)*wt[iz]; */
            a[iz][ix] *= iz>ix? wl[ix]:wt[iz];
        }
    }
    for (iz=0; iz < nbt; iz++) {  
        for (ix=0; ix < nbr; ix++) {
            a[iz][nx+nbl+nbr-1-ix] *= iz>ix? wr[ix]:wt[iz];
        }
    }
    for (iz=0; iz < nbb; iz++) {  
        for (ix=0; ix < nbl; ix++) {
            a[nz+nbt+nbb-1-iz][ix] *= iz>ix? wl[ix]:wb[iz];
        }
    }
    for (iz=0; iz < nbb; iz++) {  
        for (ix=0; ix < nbr; ix++) {
            a[nz+nbt+nbb-1-iz][nx+nbl+nbr-1-ix] *= iz>ix? wr[ix]:wb[iz];
        }
    }
}

void abc_cal(int abc /* decaying type*/,
             int nb  /* absorbing layer length*/, 
             float c /* decaying parameter*/,
             float* w /* output weight[nb] */)
/*< find absorbing coefficients >*/
{
    int ib;
    const float pi=SF_PI;
    if(!nb) return;
    switch(abc){
	default:
	    for(ib=0; ib<nb; ib++){
		w[ib]=exp(-c*c*(nb-1-ib)*(nb-1-ib));
	    }
	    break;
	case(1): 
	    for(ib=0; ib<nb; ib++){
		w[ib]=powf((1.0+0.9*cosf(((float)(nb-1.0)-(float)ib)/(float)(nb-1.0)*pi))/2.0,abc);
	    }
    }   
}
float dehf(float k /*current frequency*/,
	   float kn /*highest frequency*/,
	   float a /*suppress factor*/,
	   float factor /*propotion*/)
/*< high frequency depressing>*/
{
    float kmax;
    float depress;
    /* float pi=SF_PI; */
    kmax =  (kn*factor);
    /* kmax2 = (kmax+kn)/2.0; */
    if (fabs(k) < kmax) {
	depress = 1.0;
    }
    else {
	depress = exp(-a*(float)((k-kmax)*(k-kmax))/((float)((kn-kmax)*(kn-kmax))));
/*          depress =cosf(((fabs(k)-kmax))/((kn-kmax))*pi/2.0); */
	/*       depress = depress * depress; */
    }
    /*   depress = exp(-a*((fabs(k)-kmax)*(fabs(k)-kmax))/((kn-kmax)*(kn-kmax))); */
/*
  else if (fabs(k) < kmax2){
  depress =cosf(((fabs(k)-kmax))/((kmax2-kmax))*pi/2.0);
  depress = depress * depress;
  }    
  // depress = exp(-a*(float)((k-kmax)*(k-kmax))/((float)((kn-kmax)*(kn-kmax)))); 
  else
  depress = 0.0;
//       depress = exp(-a*(float)((k-kmax)*(k-kmax))/((float)((kn-kmax)*(kn-kmax))));
*/
    return(depress);
}

