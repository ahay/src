! Calculate Thomsen anisotropy parameters from 	well log (vp,vs,rho) data via BACKUS averaging
program backus
!! Credits:
!!
!!	UHouston: Chris Liner 
!!      UWA : Jeff Shragge (transcription from SU)
!!
!! References:		
!! Anisotropy parameters: Thomsen, 2002, DISC Notes (SEG)
!! Backus Method: Berryman, Grechka, and Berge, 1997, SEP94
!!
!! Version:C
!! SUBACKUS: $Revision: 1.3 $ ; $Date: 2009/12/02 18:03:11 
!! DOC:
!! Notes:							",
!! 1. Input are (vp,vs,rho) traces in metric units		",
!! 2. Output are (epsilon,delta,gamma) dimensionless traces	",
!!    tracl=(1,2,3)=(epsilon,delta,gamma) 			",
!!	(all=1 output optional traces listed below)		",
!!    tracl=(4,5,6,7,8)=(vp0,vs0,<rho>,eta,vang)			",
!!    tracl=(9,10,11,12)=(a,f,c,l)=(c11,c13,c33,c44) backus avg'd",
!!    tracl=(13,14)=(<lam/lamp2mu>^2,4<mu*lampmu/lamp2mu>)=(A,B)	",
!!       used to analyze eps=(a-c)/2c; a=c11=A*x+B;  c=c33=x	",
!!    tracl=(15,16)=(<lambda>,<mu>)				",
!!       for fluid analysis (lambda affected by fluid, mu not)   ",
!!    tracl=(17,18,19)=(vp,vs,rho)  orig log values		",
!!    tracl=(20)=(m)=(c66) Backus avg'd 				",
!!    tracl=(21,22,23,24,25)=(a,f,c,l,m)=(c11,c13,c33,c44,c66) orig",
!! 3. (epsilon,delta,etc.) can be isolated by tracl header field ",
!! 4. (vp0,vs0) are backus averaged vertical wavespeeds		",
!! 5. <rho> is backus averaged density, etc.			",
!! 6. eta = (eps - delta) / (1 + 2 delta)
!! 7. vang=vp(ang_deg)/vp0  phase velocity ratio			",
!!   The idea being that if vang~vp0 there are small time effects",
!!    (30 deg comes from ~ max ray angle preserved in processing)",
!!
!! END SU DOC

  use rsf
  use backusmod

  integer :: nz,i,ii,navg,nx,nc
  real, dimension(:,:,:), allocatable :: in,out
  real, dimension(:,:), allocatable :: vp,vs,rho,mu,lp2mu,lam,lpmu,eps,delta
  real, dimension(:,:), allocatable :: gamma,vp0,vs0,eta,rhoa,vang,t1,t2,t3,t4,t5,t6
  real, dimension(:,:), allocatable :: c11,c13,c33,c44,c66,aa,bb,blam,bmu,tmp
  real :: d1,f1,dt,sina,cosa,sina2,cosa2,sina4,u,ds,dds,vpa,pi
  logical :: all
  type(file) :: infile,outfile

  !! FOR OMP
  integer, external :: omp_get_num_threads,omp_get_thread_num
  integer           :: nth,ith

  call sf_init()  
  infile = rsf_input ("in")
  outfile= rsf_output("out")

  call from_par("navg",navg,201) ! Number of samples to average over
  call from_par("all",all,.false.) ! Print extra information
  call from_par("ang",ang,30.) ! Input angle

  !! . . Ensure odd
  if (navg .ne. 201) navg = navg/2*2+1

  !! . . Get trace information
  call from_par(infile,"n1",nz) ! Number of log sample
  call from_par(infile,"d1",d1) ! sample interval
  call from_par(infile,"o1",f1) ! first offset
  call from_par(infile,"n2",nx) ! number of traces
  call from_par(infile,"n3",nc) ! number of input fields

  write(0,*) 'sfbackus program'
  write(0,*) 'Using a backus averaging window of ',navg
  if (all) then
     write(0,*) 'Outputing all 25 fields'
  else
     write(0,*) 'Outputing only 3 fields'
  end if

  !! Set up correct output file
  if (nc .ne. 3) call sf_error("Wrong number of inputs!")
  allocate( in(nz,nx,nc) )
  if (all) then
     call to_par(outfile,"n3",26)
     allocate( out(nz,nx,26) )
  else
     call to_par(outfile,"n3",3)
     allocate( out(nz,nx,3) )
  end if

#ifdef _OPENMP
  !$OMP PARALLEL
  nth = omp_get_num_threads()
  !$OMP END PARALLEL
#else
  nth = 1
#endif
  !nth=32 
  write(0,*) 'USING NUM THREADS: ',nth

  !! allocate
  allocate( vp(nz,nth),vs(nz,nth),rho(nz,nth) )
  allocate( mu(nz,nth),lp2mu(nz,nth),lam(nz,nth),lpmu(nz,nth),eps(nz,nth),delta(nz,nth) )
  allocate( gamma(nz,nth),vp0(nz,nth),vs0(nz,nth))
  allocate( eta(nz,nth),rhoa(nz,nth),vang(nz,nth))
  allocate( t1(nz,nth),t2(nz,nth),t3(nz,nth),t4(nz,nth),t5(nz,nth),t6(nz,nth) )
  allocate( tmp(nz,nth),c11(nz,nth),c13(nz,nth),c33(nz,nth),c44(nz,nth),c66(nz,nth) )
  allocate( aa(nz,nth),bb(nz,nth), blam(nz,nth),bmu(nz,nth) )  


  !! . . Initialize modules
  call dobackus_init(navg,nz)

  !! Read in input Vp, Vs and rho
  call rsf_read(infile,in)

  !! Loop over all traces
  !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ith,ix,ii,pi,sina,cosa,sina2,cosa2,sina4) SHARED(nz,all)
  TRACES: do ix=1,nx
#ifdef _OPENMP     
     ith = omp_get_thread_num()+1
#else
     ith = 1
#endif
     !! Assign Vp, Vs and Rho
     do ii=1,nz
        vp (ii,ith)= in(ii,ix,1)
        vs (ii,ith)= in(ii,ix,2)
        rho(ii,ith)= in(ii,ix,3)
     end do
          
     !! . . Fill aux arrays
     do ii=1,nz
        mu(ii,ith)   = rho(ii,ith)*vs(ii,ith)*vs(ii,ith)
        lp2mu(ii,ith)= rho(ii,ith)*vp(ii,ith)*vp(ii,ith)
        lam(ii,ith)  = rho(ii,ith)*(vp(ii,ith)*vp(ii,ith)-2.*vs(ii,ith)*vs(ii,ith) )
        lpmu(ii,ith) = rho(ii,ith)*(vp(ii,ith)*vp(ii,ith)-   vs(ii,ith)*vs(ii,ith) )
     end do
     
     !! Backus averaged lambda and mu
     call dobackus(lam(:,ith),blam(:,ith))
     call dobackus(mu(:,ith),bmu(:,ith))
     
     !! Calculate C11
     do ii=1,nz
        t1(ii,ith) = lam(ii,ith)/lp2mu(ii,ith)
        t2(ii,ith) = 1./lp2mu(ii,ith)
        t3(ii,ith) = mu(ii,ith)*lpmu(ii,ith)/lp2mu(ii,ith)
     end do
     call dobackus(t1(:,ith),t4(:,ith))
     call dobackus(t2(:,ith),t5(:,ith))
     call dobackus(t3(:,ith),t6(:,ith))
     do ii=1,nz
        t4(ii,ith)=t4(ii,ith)*t4(ii,ith)
        t5(ii,ith)=1./t5(ii,ith)
        c11(ii,ith)=t4(ii,ith)*t5(ii,ith)+4.*t6(ii,ith)
        aa(ii,ith) = t4(ii,ith)
        bb(ii,ith) = 4.*t6(ii,ith)
     end do
     
     !! Calculate C13
     do ii=1,nz
        t1(ii,ith) = lam(ii,ith)/lp2mu(ii,ith)
        t2(ii,ith) = 1./lp2mu(ii,ith)
     end do
     call dobackus(t1(:,ith),t3(:,ith))
     call dobackus(t2(:,ith),t4(:,ith))
     do ii=1,nz
        t4(ii,ith)=1./t4(ii,ith)
        c13(ii,ith)=t3(ii,ith)*t4(ii,ith)
     end do
     
     !! Calculate C33
     do ii=1,nz
        t1(ii,ith) = 1./lp2mu(ii,ith)
     end do
     call dobackus(t1(:,ith),c33(:,ith))
     do ii=1,nz
        c33(ii,ith) = 1./c33(ii,ith)
     end do
     
     !! Calculate C44
     do ii=1,nz
        t1(ii,ith) = 1./mu(ii,ith)
     end do
     call dobackus(t1(:,ith),c44(:,ith))
     do ii=1,nz
        c44(ii,ith) = 1./c44(ii,ith)
     end do
     
     !! Calculate C66
     do ii=1,nz
        t1(ii,ith)=mu(ii,ith)
     end do
     call dobackus(t1(:,ith),c66(:,ith))
     
     !! Calculate rhoa (averaged density)
     call dobackus(rho(:,ith),rhoa(:,ith))
     
     !! Set up trig functions
     pi = acos(-1.)
     sina = sin(ang*pi/180.)
     cosa = cos(ang*pi/180.)
     sina2 = sina*sina
     cosa2 = cosa*cosa
     sina4 = sina2*sina2
     
     !! Calculate anisotropic params and vertical wave speeds
     do ii=1,nz
        eps(ii,ith) = (c11(ii,ith) - c33(ii,ith)) / (2.0*c33(ii,ith));
        t1(ii,ith) = (c13(ii,ith) + c44(ii,ith)) * (c13(ii,ith) + c44(ii,ith));
        t2(ii,ith) = (c33(ii,ith) - c44(ii,ith)) * (c33(ii,ith) - c44(ii,ith));
        t3(ii,ith) = 2.0 * c33(ii,ith) * (c33(ii,ith) - c44(ii,ith));
        delta(ii,ith) = (t1(ii,ith) - t2(ii,ith)) / t3(ii,ith);
        gamma(ii,ith) = (c66(ii,ith) - c44(ii,ith)) / (2.0 * c44(ii,ith));
        vp0(ii,ith) = sqrt( c33(ii,ith) / rhoa(ii,ith) );
        vs0(ii,ith) = sqrt( c44(ii,ith) / rhoa(ii,ith) );
        eta(ii,ith) = (eps(ii,ith) - delta(ii,ith)) / (1.0+2.0*delta(ii,ith));
        u = 1.-(vs0(ii,ith)*vs0(ii,ith))/(vp0(ii,ith)*vp0(ii,ith));
        ds = (2.*delta(ii,ith) - eps(ii,ith))*u;
        dds = 0.5*(sqrt(u*u+4.*ds*sina2*cosa2+4.*(u+eps(ii,ith))*eps(ii,ith)*sina4)-u); 
        vpa = vp0(ii,ith)*sqrt(1.+eps(ii,ith)*sina2+dds);
        vang(ii,ith) = vpa/vp0(ii,ith);
     end do
     
     !! . . Treat end points
     call handlEnds(eps(:,ith))
     call handlEnds(delta(:,ith))
     call handlEnds(gamma(:,ith))
     call handlEnds(vp0(:,ith))
     call handlEnds(vs0(:,ith))
     call handlEnds(rhoa(:,ith))
     call handlEnds(eta(:,ith))
     call handlEnds(vang(:,ith))
     call handlEnds(c11(:,ith))
     call handlEnds(c13(:,ith))
     call handlEnds(c33(:,ith))
     call handlEnds(c44(:,ith))
     call handlEnds(aa(:,ith))
     call handlEnds(bb(:,ith))
     call handlEnds(blam(:,ith))
     call handlEnds(mu(:,ith))
     call handlEnds(vp(:,ith))
     call handlEnds(vs(:,ith))
     call handlEnds(rho(:,ith))
     call handlEnds(c66(:,ith))
     call handlEnds(lp2mu(:,ith))
     call handlEnds(lam(:,ith))
     
     !! Set up output
     do ii=1,nz
        out(ii,ix,1) = eps(ii,ith)
        out(ii,ix,2) = delta(ii,ith)
        out(ii,ix,3) = gamma(ii,ith)
     end do
     !! . . Only output if requested
     if (all) then
        do ii=1,nz
           out(ii,ix,4) = vp0(ii,ith)
           out(ii,ix,5) = vs0(ii,ith)
           out(ii,ix,6) = rhoa(ii,ith)
           out(ii,ix,7) = eta(ii,ith)
           out(ii,ix,8) = vang(ii,ith)
           out(ii,ix,9) = c11(ii,ith)
           out(ii,ix,10)= c13(ii,ith)
           out(ii,ix,11)= c33(ii,ith)
           out(ii,ix,12)= c44(ii,ith)
           out(ii,ix,13)= aa(ii,ith)
           out(ii,ix,14)= bb(ii,ith)
           out(ii,ix,15)= blam(ii,ith)
           out(ii,ix,16)= mu(ii,ith)
           out(ii,ix,17)= vp(ii,ith)
           out(ii,ix,18)= vs(ii,ith)
           out(ii,ix,19)= rho(ii,ith)
           out(ii,ix,20)= c66(ii,ith)
           out(ii,ix,21)= lp2mu(ii,ith)
           out(ii,ix,22)= lam(ii,ith)
           out(ii,ix,23)= lp2mu(ii,ith)
           out(ii,ix,24)= mu(ii,ith)
           out(ii,ix,25)= mu(ii,ith)
        end do
     end if
     
  end do TRACES
  !$OMP END PARALLEL DO

  call rsf_write(outfile,out)
  
  deallocate(vp,vs,rho,mu,lp2mu,lam,lpmu,eps,delta,gamma,vp0,vs0,eta,rhoa,vang)
  deallocate(t1,t2,t3,t4,t5,t6,tmp,c11,c13,c33,c44,c66,aa,bb,blam,bmu)
  
  call exit()

end program backus
