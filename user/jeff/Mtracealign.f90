program TRACEALIGN
  use rsf
  use fft
  implicit none

  integer :: ii,jj,kk,mm,nwindow,ix,ntpad,iw,ih,id,it,iy,is,iloc
  integer :: nstart,nstop,nshift,ns,nt,nh,ntmin,ntmax,nlen,itguess
  real :: ot,dt,os,ds,vel,oh,dh,rowsum,phitemp
  real :: maxshift,w0,dw,pi,w,BASErms,MONrms
  integer :: yloc,yyy(1)
  real, allocatable :: baseline(:,:),monitor(:,:),phi(:),delt(:),trace1(:),trace2(:),test(:)

  complex, allocatable :: cdata(:)

  type(file) :: infile,monitorfile,outfile,timefile

  !! For OMP
!  integer, external :: omp_get_num_threads,omp_get_thread_num
  integer  :: i1,ie3,ie1,ih1

  call sf_init()

  infile = rsf_input("in")
  monitorfile = rsf_input("monitor")
  outfile = rsf_output("out")
  timefile = rsf_output("times")

  !! . . Read in axis components
  call from_par(infile,"n1",nt)
  call from_par(infile,"o1",ot)
  call from_par(infile,"d1",dt)
  call from_par(infile,"n2",nh)
  call from_par(infile,"o2",oh)
  call from_par(infile,"d2",dh)
  call from_par(infile,"n3",ns)
  call from_par(infile,"o3",os)
  call from_par(infile,"d3",ds)

  write (0,*) 'nt, nh, ns',nt,nh,ns
  !! . . Write out axis components
  !! . . Read in axis components
  call to_par(outfile,"n1",nt)
  call to_par(outfile,"o1",ot)
  call to_par(outfile,"d1",dt)
  call to_par(outfile,"n2",nh)
  call to_par(outfile,"o2",oh)
  call to_par(outfile,"d2",dh)
  call to_par(outfile,"n3",ns)
  call to_par(outfile,"o3",os)
  call to_par(outfile,"d3",ds)
  call to_par(outfile,"esize",4)

  !! . . Write out optimal delay times
  call to_par(timefile,"n1",nh)
  call to_par(timefile,"o1",oh)
  call to_par(timefile,"d1",dh)
  call to_par(timefile,"n2",ns)
  call to_par(timefile,"o2",os)
  call to_par(timefile,"d2",ds)
  call to_par(timefile,"n3",1)

  !! . . Read in Parameters
  call from_par("maxshift",nshift) ! Maximum allowed time shift
  call from_par("nlen",nlen) ! Window length of shift vector (in samples)
  call from_par("vel",vel,1500.) ! Rupture speed for linear shift

  !! . . Padding Factors
  ntpad=pad2(2*nt)
  
  !! allocate matrices
  allocate( baseline(nt,nh),monitor(nt,nh)  )
  allocate ( phi(2*nshift+1)  , delt(nh) , test(nh) )
  allocate( cdata(ntpad) , trace1(2*nshift+nlen+1),trace2(2*nshift+nlen+1) )

 
 !! . . Definitions
  pi = acos(-1.) 
  w0=-pi/dt;  
  dw=2.*pi/(ntpad*dt);

!! . . Loop over all offsets
  SHOTS: do is=1,ns
     write(0,*) 'WORKING ON SHOT: ',is,ns

     !! . . Read in data files
     call rsf_read(infile,baseline)
     call rsf_read(monitorfile,monitor)
     test=0.
     OFFSETS: do ih=1,nh
     
        !! . . Zero test times each time
 
        !! . . Given offset where to cut
        itguess = (( (ih-1)*dh+oh ) / vel ) / dt + 1

        !! . . Range for start/end samples
        ntmin=itguess-nlen/2
        ntmax=itguess+nlen/2
        if (ntmin .le.0 ) then
           ntmin=1;
           ntmax=nlen;
        end if
        if (ntmax.gt. nt) then 
           ntmin=nt-nlen;
           ntmax=nt;
        end if
        
        !! . . Make sure datasets have same zero traces
        if (sum (baseline(:,ih) ).eq. 0.) then
           baseline(:,ih) = 0.
           monitor(:,ih)=0.
        end if
        if (sum (monitor(:,ih) ).eq. 0.) then
           baseline(:,ih) = 0.
           monitor(:,ih)=0.
        end if

        !! . . Perform XCorr
        trace1(nshift+1:nshift+nlen) = baseline(ntmin:ntmax,ih)
        trace2(nshift+1:nshift+nlen) = monitor (ntmin:ntmax,ih)
              
        !! . . For all relative shifts
        phi = 0.
        do mm=-nshift,nshift
           phitemp=0.
           do it=1,nlen
              phitemp=phitemp+trace1(nshift+it-mm)*trace2(nshift+it)
           end do
           phi(mm+nshift+1)=phitemp
        end do
        
        yyy = maxloc( phi )
        test(ih) = (yyy(1)-nshift-1)*dt

        !! . . fft the Receiver wavefields and apply the estimated time shifts
        cdata = 0.
        cdata(1:nt) = cmplx(monitor(:,ih),0.)

        call fth(.false.,.false.,cdata)
        do iw=1,ntpad
           cdata(iw)=cdata(iw)*exp(cmplx(0.,-(w0+dble(iw-1.)*dw)*test(ih) ))
        end do
        call fth(.true.,.false.,cdata)

        !! . . COMPUTE A GLOBAL SCALING
        BASErms = sqrt(sum( baseline(:,ih)*baseline(:,ih)));
        MONrms  = sqrt(sum(  monitor(:,ih)* monitor(:,ih)));

        if (BASErms .ne. 0. .and. MONrms .ne. 0.)  then
           monitor(:,ih) = BASErms/MONrms*real(cdata(1:nt))
        end if

     end do OFFSETS
     call rsf_write(timefile,test)
     call rsf_write(outfile,monitor)

  end do SHOTS

  call exit(0)
end program TRACEALIGN

