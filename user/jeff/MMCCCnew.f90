program MCCCnew
  use rsf
  use fft
  implicit none

  integer           :: ii,jj,kk,mm,nwindow,ix,ntpad,iw,ih,id,it,iy
  integer           :: nstart,nstop,nshift,nr,nt,nh,ntmin,ntmax,nlen,itguess
  real                :: ot,dt,or,dr,vel,oh,dh,rowsum
  real                :: maxshift,w0,dw,pi,w
  integer :: yloc,yyy(1)
  real, allocatable :: data(:,:),phi(:,:) ,sysmat(:,:),delt(:),test(:),trace1(:,:),trace2(:,:),crap(:)

  complex, allocatable :: cdata(:,:),cmodl(:,:)

  type(file) :: infile,outfile,timefile

  !! For OMP
  integer, external :: omp_get_num_threads,omp_get_thread_num
  integer  :: nth,ith,i1,ie3,ie1,ih1

  call sf_init()

  infile = rsf_input("in")
  outfile = rsf_output("out")
  timefile = rsf_output("times")

  !! . . Read in axis components
  call from_par(infile,"n1",nt)
  call from_par(infile,"o1",ot)
  call from_par(infile,"d1",dt)
  call from_par(infile,"n2",nr)
  call from_par(infile,"o2",or)
  call from_par(infile,"d2",dr)
  call from_par(infile,"n3",nh)
  call from_par(infile,"o3",oh)
  call from_par(infile,"d3",dh)

  write (0,*) 'nt, nr, nh',nt,nr,nh
  !! . . Write out axis components
  !! . . Read in axis components
  call to_par(outfile,"n1",nt)
  call to_par(outfile,"o1",ot)
  call to_par(outfile,"d1",dt)
  call to_par(outfile,"n2",nr)
  call to_par(outfile,"o2",or)
  call to_par(outfile,"d2",dr)
  call to_par(outfile,"n3",nh)
  call to_par(outfile,"o3",oh)
  call to_par(outfile,"d3",dh)
  call to_par(outfile,"esize",4)

  !! . . Write out optimal delay times
  call to_par(timefile,"n1",nr)
  call to_par(timefile,"o1",or)
  call to_par(timefile,"d1",dr)
  call to_par(timefile,"n2",nh)
  call to_par(timefile,"o2",oh)
  call to_par(timefile,"d2",dh)
  call to_par(timefile,"n3",1)

  !! . . Read in Parameters
  call from_par("maxshift",nshift) ! Maximum allowed time shift
  call from_par("nlen",nlen) ! Window length of shift vector (in samples)
  call from_par("vel",vel,1500.) ! Rupture speed for linear shift

#ifdef _OPENMP
  !$OMP PARALLEL
  nth = omp_get_num_threads()
  !$OMP END PARALLEL
#endif 

!!! . . Padding Factors
  ntpad=pad2(2*nt)
  
!! allocate matrices
allocate( data(nt,nr)  )
allocate ( phi(2*nshift+1,nth)  , delt(nr*nr+1) , sysmat(nr*nr+1,nr) , test(nr) , crap(2*nshift+1) )
allocate( cdata(ntpad,nth) , cmodl(ntpad,nth) , trace1(2*nshift+nlen+1,nth),trace2(2*nshift+nlen+1,nth)  )

!!phi=0.;data=0.; delt=0.; sysmat=0.;test=0.; crap=0.; cdata=0.;cmodl=0.;trace1=0.; trace2=0.;
 
 !! . . Definitions
  pi = acos(-1.) 
  w0=-pi/dt;  
  dw=2.*pi/(ntpad*dt);

!! . . Loop over all offsets
  OFFSETS: do ih=1,nh
     
     !! . . Zero test times each time
     test = 0.
 
     !! . . Given offset where to cut
     itguess = (( (ih-1)*dh+oh ) / vel ) / dt + 1

     !! . . Range for start/end samples
     ntmin=itguess-nlen/2
     ntmax=itguess+nlen/2
     if (ntmin .le.0 ) then
        ntmin=1;
        ntmax=nlen;
     endif
     if (ntmax.gt. nt) then 
        ntmin=nt-nlen;
        ntmax=nt;
     endif

     write(0,*) 'WORKING ON OFFSET ih: ',ih,ntmin,ntmax
     !! . . Input Data
     call rsf_read(infile,data)
 
     do ix=2,nr-1
         if (sum (data(:,ix) ).eq. 0.) data(:,ix) = data(:,ix-1)
     end do
100 FORMAT("COMPLETE %:", F5.2,A16)
     !! . . Perform XCorr
     !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(mm,ix,it,ith,iy,yyy,kk) SHARED(nr,data,dt,nshift)
     do ix=1,nr
#ifdef _OPENMP
        ith = omp_get_thread_num()+1
#else
        ith = 1
#endif
        if (ith .eq. 1) write(0,100,ADVANCE='NO') 100.*real(ix)/real(nr), "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b"
        trace1 (nshift+1:nshift+nlen,ith) = data(ntmin:ntmax,ix)

        do iy=1,nr
           trace2(nshift+1:nshift+nlen,ith)  = data(ntmin:ntmax,iy)

           !! . . Eliminate b/c of high autocorrelation at t=0!
       !!    if (ix .eq. iy) cycle

           !! . . For all relative shifts
           phi(:,ith) = 0.
           do mm=-nshift,nshift
              do it=1,nlen
                 phi(mm+nshift+1,ith)=phi(mm+nshift+1,ith)+trace1(nshift+it-mm,ith)*trace2(nshift+it,ith)
              end do
           end do
           
           kk=(ix-1)*nr+iy
           yyy = maxloc( phi(:,ith) )
           delt(kk) = (yyy(1)-nshift+1)*dt

           sysmat(kk,ix)=1.
           sysmat(kk,iy)=-1.
        end do
     end do
     !$OMP END PARALLEL DO
     
     sysmat(nr*nr+1,:)=1.
     delt(nr*nr+1)=0.
     !! Find optimal delay times
     !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ii,jj) SHARED(nr,sysmat,delt)
     do ii = 1, nr
        do jj = 1, nr*nr+1
           test(ii) = test(ii) + sysmat(jj,ii)*delt(jj)/real(2.*nr)
        end do
     end do
     !$OMP END PARALLEL DO


     call rsf_write(timefile,test)
     !! . . fft the Receiver wavefields and apply the estimated time shifts
    !$OMP PARALLEL DO SCHEDULE(DYNAMIC) PRIVATE(ix,id,iw) SHARED(cdata,ntpad,w0,dw,test,nr)
     do ix=1,nr
#ifdef _OPENMP
        id = omp_get_thread_num()+1
#else
        id = 1
#endif
        cdata(:,id) =0.
        cdata(1:nt,id) = cmplx(data(:,ix),0.)

        call fth(.false.,.false.,cdata(:,id))
        do iw=1,ntpad
           cdata(iw,id)=cdata(iw,id)*exp(cmplx(0.,(w0+dble(iw-1.)*dw)*test(ix) ))
        end do
        call fth(.true.,.false.,cdata(:,id))

        data(:,ix) = real(cdata(1:nt,id))

     end do
     !$OMP END PARALLEL DO

     call rsf_write(outfile,data)

  end do OFFSETS

  call exit(0)
end program MCCCnew

