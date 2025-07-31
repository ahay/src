!  3-D Radon transform using different sparsity-promotion algorithms, including IRLS, CGG, IST, LBreg.
!
!
!!$  Copyright (C) 2012 China University of Petroleum (East China)
!!$  Author: Yujin Liu
!!$  Email: seisinv@gmail.com
!!$  Notes: The code was changged from sfvelinvww (Jun Ji)
!!$
!!$  This program is free software; you can redistribute it and/or modify
!!$  it under the terms of the GNU General Public License as published by
!!$  the Free Software Foundation; either version 2 of the License, or
!!$  (at your option) any later version.
!!$  
!!$  This program is distributed in the hope that it will be useful,
!!$  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!$  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!$  GNU General Public License for more details.
!!$  
!!$  You should have received a copy of the GNU General Public License
!!$  along with this program; if not, write to the Free Software
!!$  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
program Mvelinvww3
  !!
  !! Mvelinvww3 < input.H nv= dv= ov= niter= cmpout= > output.H
  !!
  !!	adj = 0  : from velocity-domain(t,v) to cmp-gather(t,x,y)
  !!	adj = 1  : from cmp-gather(t,x,y) to velocity-domain(t,v)
  !!
  !! 
!!!!!!!!!!!!!!!!!!
  !! Required modules

  use rsf
  use velinv

  implicit none
  real, allocatable, dimension (:,:) :: cmp
  real, allocatable, dimension (:,:) :: intcmp
  real, allocatable, dimension (:,:) :: vel
  real, allocatable, dimension (:) :: mask
  real, allocatable, dimension (:) :: h2
  real, allocatable, dimension (:) :: z2
  real, allocatable, dimension (:) :: s2
  real, allocatable, dimension (:) :: res, mres, vel0
  integer :: nt,nhx,nhy,nh,nv
  integer :: ic, nc
  integer :: it,ihx,ihy,ih,iv
  real    :: dt,dhx,dhy,dv, lip
  real    :: ot,ohx,ohy,ov, rwt, mwt, srate,eps, hx,hy, z, s
  integer :: flag, niter, savevel, huber,irls,nstep, reset, mflag
  type (file) :: infile, outfile, vtr, fres, fm, fmres
  real delta,lamda,step,alpha

  integer(kind=8) sf_leftsize
  external sf_leftsize

  call sf_init()

  infile = rsf_input()
  fres = rsf_output("res")

  call from_par(infile,"n1",nt);  call from_par(infile,"n2",nhx);  call from_par(infile,"n3",nhy);
  call from_par(infile,"d1",dt);  call from_par(infile,"d2",dhx);  call from_par(infile,"d3",dhy)
  call from_par(infile,"o1",ot);  call from_par(infile,"o2",ohx);  call from_par(infile,"o3",ohy)

  ! get the number of CMP gathers 
  nc = sf_leftsize(infile,3)

  call from_par("nv",nv,nhx);  call from_par("dv",dv,0.01);  call from_par("ov",ov,1.5)

  call from_par("niter",niter,20)
  call from_par("savevel",savevel,0)

  ! Flag to choose the algorithm
  ! Flag=0 IRLS or CGG method
  ! Flag=1 IST method
  ! Flag=2 LBreg method with fix stepsize
  ! Flag=3 LBreg method with BB line search
  ! Flag=4 LBreg method with Nesterov's acceleration

  call from_par("flag",flag,0)

  call from_par("mflag",mflag,0)

  if (mflag.eq.1) then
     fm = rsf_input("vel0")
     fmres = rsf_output("mres")
     allocate(mres(niter), vel0((nt*nv)))
     call rsf_read(fm,vel0)
  else
     fmres = rsf_output("mres")
     allocate(mres(niter), vel0((nt*nv)))
     vel0=0.
     mres=0.
  endif

  if (flag.eq.0) then
     call from_par("huber",huber,0)
     call from_par("irls",irls,0)
     call from_par("nstep",nstep,1)
     call from_par("rwt",rwt,0.)
     call from_par("mwt",mwt,0.)
     call from_par("srate",srate,0.01)
     call from_par("eps",eps,0.01)

  elseif (flag.eq.1) then
     call from_par("lamda",lamda,1000.)
     ! lamda controls sparsity, bigger lamda, more sparsity
     call from_par("delta",delta,0.0001)
     ! delta controls update step and convergent, small delta ensure convergence but with small decrease in data fit error

  elseif (flag.eq.2.or.flag.eq.3) then
     call from_par("step",step,0.000005)
     ! step is very important in convergence and sparsity
     call from_par("alpha",alpha,790.635)
     ! alpha always equals to the L2-norm of correct x 
  elseif (flag.eq.4) then
     call from_par("alpha", alpha)
     ! smoothing parameter, typical value: 1 to 10 times estimated norm(x,inf)

     call from_par("lip", lip)
     ! the estimated Lipschitz constrant of the dual objective, default: alpha*normest(A*A',1e-2)

     call from_par("reset", reset)
     ! Nesterov's acceleration restart (theta is reset) or skip (theta is not reset)
     ! 0: no  restart or skip
     ! 1: restart if dual objective is non-monotonic
     ! 2: skip    if dual objective is non-monotonic
     ! 3: restart if residual is sufficiently non-monotonic
     ! 4: skip    if residual is sufficiently non-monotonic
     ! 5: restart according to the gradient test on P6 of http://www-stat.stanford.edu/~candes/papers/adap_restart_paper.pdf
     ! 6: skip    according to the gradient test

  endif

  ! Allocate
  nh = nhx*nhy
  allocate (cmp(nt,nh))
  allocate (intcmp(nt,nh))
  allocate (vel(nt,nv))
  allocate (mask(nh))
  allocate (h2(nh))
  allocate (z2(nt))
  allocate (s2(nv))
  allocate (res(niter))

  outfile = rsf_output()

  call to_par(outfile,"n4",nc);
  call to_par(outfile,"n1",nt);  call to_par(outfile,"n2",nhx);  call to_par(outfile,"n3",nhy)
  call to_par(outfile,"d1",dt);  call to_par(outfile,"d2",dhx);  call to_par(outfile,"d3",dhy)
  call to_par(outfile,"o1",ot);  call to_par(outfile,"o2",ohx);  call to_par(outfile,"o3",ohy)

  call to_par(fres,"n4",1);
  call to_par(fres,"n1",niter);  call to_par(fres,"n2",nc);  call to_par(fres,"n3",1)
  call to_par(fres,"d1",1);  call to_par(fres,"d2",1);  call to_par(fres,"d3",1)
  call to_par(fres,"o1",1);  call to_par(fres,"o2",1);  call to_par(fres,"o3",1)

  call to_par(fmres,"n4",1);
  call to_par(fmres,"n1",niter);  call to_par(fmres,"n2",nc);  call to_par(fmres,"n3",1)
  call to_par(fmres,"d1",1);  call to_par(fmres,"d2",1);  call to_par(fmres,"d3",1)
  call to_par(fmres,"o1",1);  call to_par(fmres,"o2",1);  call to_par(fmres,"o3",1)

  do ihy= 1,nhy
     hy = ohy+(ihy-1)*dhy
     do ihx= 1,nhx
        hx = ohx+(ihx-1)*dhx
        ih = ihx+(ihy-1)*nhx
        h2(ih) = hx*hx+hy*hy
     end do
  end do

  do it= 1,nt
     z = ot + dt*(it-1);        z2(it) = z*z
  end do

  do iv= 1,nv
     s = 1./(ov + dv*(iv-1));   s2(iv) = s*s
  end do

  if ( savevel.eq.1 ) then
     vtr = rsf_output("velout")

     call to_par(vtr,"n3",1);      call to_par(vtr,"n4",nc)
     call to_par(vtr,"n1",nt);     call to_par(vtr,"n2",nv);
     call to_par(vtr,"d1",dt);     call to_par(vtr,"d2",dv)
     call to_par(vtr,"o1",ot);     call to_par(vtr,"o2",ov)
  end if

  write(0,*) 'Hyperbolic interpolation begin...'

  do ic = 1, nc
     write(0,*) 'ic/nc=',ic,'/',nc
     call rsf_read(infile, cmp)
     do ih=1,nh
        mask(ih) = dot(nt,cmp(1,ih),cmp(1,ih))
     end do
     !************************IRLS or CGG method*******************************
     write(0,*) 'IRLS or CGG method'
     if (flag.eq.0) then
        call velinvww( nt,dt,ot, h2,nh,cmp, s2,nv,vel, z2,mask,rwt,mwt,niter&
          &,huber,eps,irls,nstep,srate, res, vel0, mres)
     end if

     !************************Soft shrinkage method*****************************
     if (flag .eq.1) then
        call velinvsh( nt,dt,ot, h2,nh,cmp, s2,nv,vel, z2,mask,niter,&
           lamda,delta, res, vel0, mres)
     end if

     !********************Linear Bregman method with fix stepsize****************
     if (flag.eq.2) then
        call velinvlbreg( nt,dt,ot, h2,nh,cmp, s2,nv,vel, z2,mask,niter,&
             step,alpha,res, vel0, mres)
     end if

     ! *******************Linear Bregman method with BB liear search************
     if (flag.eq.3) then
        call velinvbbls( nt,dt,ot, h2,nh,cmp, s2,nv,vel, z2,mask,niter,&
             step,alpha,res, vel0, mres)
     end if

     ! Linear Bregman method with Nesterov's acceleration and reset (restart or skip)*
     if (flag.eq.4) then
        call velinvaccel( nt,dt,ot, h2,nh,cmp, s2,nv,vel, z2,mask,niter,&
            lip,alpha,reset,res, vel0, mres)
     end if

     if ( savevel.eq.1 ) then
        call rsf_write(vtr, vel)
     end if

     if (mflag.eq.1) then
        call rsf_write(fmres,mres)
     end if

     do ih= 1,nh
        mask(ih) = 1.
     end do

     call velxf( 0,0, nt,dt,ot, h2,nh,intcmp, mask, s2,nv, vel, z2)

     do ih= 1,nh
        mask(ih) = dot(nt,cmp(1,ih),cmp(1,ih))
     end do

     do ih= 1,nh
        if ( mask(ih) .eq. 0. ) then
           do it=1,nt
              cmp(it,ih) = intcmp(it,ih)
           end do
        end if
     end do

     call rsf_write( outfile, intcmp )
     call rsf_write( fres, res)

  end do

  write(0,*) 'Hyperbolic interpolation ends!!!'

end program Mvelinvww3
