/* 2-D prestack least-squares RTM using wavefield reconstruction
   NB: Sponge ABC is applied!
*/
/*
  Copyright (C) 2014  Xi'an Jiaotong University (Pengliang Yang)

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

#include "prtm2d.h"

int main(int argc, char* argv[])
{   
  bool verb, fromBoundary;
  int nb, nz, nx, nt, ns, ng, niter, csd, sxbeg, szbeg, jsx, jsz, gxbeg, gzbeg, jgx, jgz;
  float dz, dx, dt, fm, o1, o2, amp;
  float **v0, *mod, *dat; 
  int testadj;

  sf_file shots, imag, imgrtm, velo;/* I/O files */

  /* initialize Madagascar */
  sf_init(argc,argv);

  shots = sf_input ("in"); /* shot records, data 	*/
  velo = sf_input ("vel"); /* velocity */
  imag = sf_output("out"); /* output LSRTM image, model */
  imgrtm = sf_output("imgrtm"); /* output RTM image */
    
  if (!sf_histint(velo,"n1",&nz)) sf_error("n1");
  /* 1st dimension size */
  if (!sf_histint(velo,"n2",&nx)) sf_error("n2");
  /* 2nd dimension size */
  if (!sf_histfloat(velo,"d1",&dz)) sf_error("d1");
  /* d1 */
  if (!sf_histfloat(velo,"d2",&dx)) sf_error("d2");
  /* d2 */
  if (!sf_histfloat(velo,"o1",&o1)) sf_error("o1");
  /* o1 */
  if (!sf_histfloat(velo,"o2",&o2)) sf_error("o2");
  /* o2 */
  if (!sf_getbool("verb",&verb)) verb=true;
  /* verbosity */
  if (!sf_getint("niter",&niter)) niter=10;
  /* totol number of least-squares iteration*/
  if (!sf_getint("nb",&nb)) nb=20;
  /* number (thickness) of ABC on each side */
  if (!sf_getbool("fromBoundary",&fromBoundary)) fromBoundary=true;
  /* if fromBoundary=true, reconstruct source wavefield from stored boundary */
  if (!sf_getint("testadj",&testadj)) testadj=0;
  /* if testadj = 1 then program only testadj without calculating */
   
  if (!sf_histint(shots,"n1",&nt)) sf_error("no nt");
  /* total modeling time steps */
  if (!sf_histint(shots,"n2",&ng)) sf_error("no ng");
  /* total receivers in each shot */
  if (!sf_histint(shots,"n3",&ns)) sf_error("no ns");
  /* number of shots */
  if (!sf_histfloat(shots,"d1",&dt)) sf_error("no dt");
  /* time sampling interval */
  if (!sf_histfloat(shots,"amp",&amp)) sf_error("no amp");
  /* maximum amplitude of ricker */
  if (!sf_histfloat(shots,"fm",&fm)) sf_error("no fm");
  /* dominant freq of ricker */
  if (!sf_histint(shots,"sxbeg",&sxbeg)) sf_error("no sxbeg");
  /* x-begining index of sources, starting from 0 */
  if (!sf_histint(shots,"szbeg",&szbeg)) sf_error("no szbeg");
  /* x-begining index of sources, starting from 0 */
  if (!sf_histint(shots,"gxbeg",&gxbeg)) sf_error("no gxbeg");
  /* x-begining index of receivers, starting from 0 */
  if (!sf_histint(shots,"gzbeg",&gzbeg)) sf_error("no gzbeg");
  /* x-begining index of receivers, starting from 0 */
  if (!sf_histint(shots,"jsx",&jsx)) sf_error("no jsx");
  /* source x-axis  jump interval  */
  if (!sf_histint(shots,"jsz",&jsz)) sf_error("no jsz");
  /* source z-axis jump interval  */
  if (!sf_histint(shots,"jgx",&jgx)) sf_error("no jgx");
  /* receiver x-axis jump interval  */
  if (!sf_histint(shots,"jgz",&jgz)) sf_error("no jgz");
  /* receiver z-axis jump interval  */
  if (!sf_histint(shots,"csdgather",&csd)) sf_error("csdgather or not required");
  /* default, common shot-gather; if n, record at every point*/

  sf_putint(imag,"n1",nz);
  sf_putint(imag,"n2",nx);
  sf_putint(imag,"n3",1);
  sf_putfloat(imag,"d1",dz);
  sf_putfloat(imag,"d2",dx);
  sf_putfloat(imag,"o1",o1);
  sf_putfloat(imag,"o2",o2);
  sf_putstring(imag,"label1","Depth");
  sf_putstring(imag,"label2","Distance");

  sf_putint(imgrtm,"n1",nz);
  sf_putint(imgrtm,"n2",nx);
  sf_putint(imgrtm,"n3",1);
  sf_putfloat(imgrtm,"d1",dz);
  sf_putfloat(imgrtm,"d2",dx);
  sf_putfloat(imgrtm,"o1",o1);
  sf_putfloat(imgrtm,"o2",o2);
  sf_putstring(imgrtm,"label1","Depth");
  sf_putstring(imgrtm,"label2","Distance");

  /* In rtm, vv is the velocity model [modl], which is input parameter; 
     mod is the image/reflectivity [imag]; dat is seismogram [data]! */
  v0=sf_floatalloc2(nz,nx);
  mod=sf_floatalloc(nz*nx);
  dat=sf_floatalloc(nt*ng*ns);

  /* initialize velocity, model and data */
  sf_floatread(v0[0], nz*nx, velo);
  memset(mod, 0, nz*nx*sizeof(float));
  sf_floatread(dat, nt*ng*ns, shots);
  prtm2d_init(verb, csd, fromBoundary, dz, dx, dt, amp, fm, nz, nx, nb, nt, ns, ng, 
	      sxbeg, szbeg, jsx, jsz, gxbeg, gzbeg, jgx, jgz, v0, mod, dat);

  // run adjoint test first.
  sf_warning("adjoint test \n");
  if (testadj){
    prtm2d_adjtest();
    sf_warning("exiting after testadj\n");
    mod=NULL;
    sf_floatwrite(mod,nz*nx,imgrtm);
    sf_floatwrite(mod,nz*nx,imag);
    prtm2d_close();
    free(*v0); free(v0);
    free(mod);
    free(dat); 	  
    exit (0);
  }
	
  sf_warning("start migration\n");
  /* original RTM is simply applying adjoint of prtm2d_lop once!*/
  prtm2d_lop(true, false, nz*nx, nt*ng*ns, mod, dat); 
  sf_floatwrite(mod, nz*nx, imgrtm);/* output RTM image */
  sf_warning("migration ends");

  if (niter>0){
    sf_warning("inside LS loop niter= %d\n",niter);
    /* least squares migration */
    sf_solver(prtm2d_lop, sf_cgstep, nz*nx, nt*ng*ns, mod, dat, niter, "verb", verb, "end");
    sf_floatwrite(mod, nz*nx, imag);  /* output inverted image */
    sf_cgstep_close();
  }

  prtm2d_close();
  free(*v0); free(v0);
  free(mod);
  free(dat); 

  exit(0);
}

