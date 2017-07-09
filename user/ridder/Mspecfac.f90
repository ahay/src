! (MD) Spectral Factorization or (MD) Autocorrelation using Helix Transform. 
! (Maximum dimension is 3)
!!$  Copyright (C) 2009 Stanford University
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
program MDHspecfac
  ! acfile = 'file.H'   Optional output file for MD AC file as addition to MD SF output
  ! sffile = 'file.H'   Optional output file for MD SF file as addition to MD AC output
  ! center = 1          Centers the impulse response
  ! haxes  = 123        Order in which axes are put on Helix transform

  use rsf
  implicit none

  integer                          :: n1, n2, n3, nh, nhp
  real                             :: d1, d2, d3, o1, o2, o3
  integer                          :: ii, jj, kk, cent
  integer                          :: center, haxes
  real   , allocatable             :: datain(:), dataout(:)
  complex, allocatable             :: cdata(:),sdata(:),acdata(:),sfdata(:)
  integer                          :: pad2
  character*50                     :: label1,label2,label3

  type (file)                      :: in, out, acfile

  call sf_init()
  in = rsf_input()
  out = rsf_output()

  ! Read from input History
  call from_par(in,'n1',n1)
  call from_par(in,'n2',n2)
  call from_par(in,'n3',n3)
  call from_par(in,'d1',d1)
  call from_par(in,'d2',d2)
  call from_par(in,'d3',d3)
  call from_par(in,'o1',o1)
  call from_par(in,'o2',o2)
  call from_par(in,'o3',o3)
  call from_par(in,'label1',label1,' ')
  call from_par(in,'label2',label2,' ')
  call from_par(in,'label3',label3,' ')

  ! Read from input Parameters
  call from_par('center', center, 1  )
  call from_par('haxes' , haxes , 123)

  ! Write to output History (needs adaption)
  if (center==1) then
     if      ( (haxes.eq.123).or.(haxes.eq.213) ) then
        call to_par(out,'o1',-((n1-1)/2)*d1)
        call to_par(out,'o2',-((n2-1)/2)*d2)
        call to_par(out,'o3',d3            )
     else if ( (haxes.eq.231).or.(haxes.eq.321) ) then
        call to_par(out,'o2',-((n2-1)/2)*d2)
        call to_par(out,'o3',-((n3-1)/2)*d3)
        call to_par(out,'o1',o1            )
     else if ( (haxes.eq.312).or.(haxes.eq.132) ) then
        call to_par(out,'o1',-((n1-1)/2)*d1)
        call to_par(out,'o3',-((n3-1)/2)*d3)
        call to_par(out,'o2',o2            )
     end if
  else
     call to_par(out,'o1',o1)
     call to_par(out,'o2',o2)
     call to_par(out,'o3',d3)
  end if
  call to_par(out,'n1',n1)
  call to_par(out,'n2',n2)
  call to_par(out,'n3',n3)
  call to_par(out,'d1',d1)
  call to_par(out,'d2',d2)
  call to_par(out,'d3',d3)

  ! Write to optional AC output file
  if ( exist_par ("acfile") ) then
     acfile = rsf_output("acfile")

     call to_par(acfile,'d1',d1)
     call to_par(acfile,'d2',d2)
     call to_par(acfile,'d3',d3)
     call to_par(acfile,'n1',n1)
     call to_par(acfile,'n2',n2)
     call to_par(acfile,'n3',n3)
     if (center==1) then
        if      ( (haxes.eq.123).or.(haxes.eq.213) ) then
           call to_par(acfile,'o1',-((n1-1)/2)*d1)
           call to_par(acfile,'o2',-((n2-1)/2)*d2)
           call to_par(acfile,'o3',d3            )
        else if ( (haxes.eq.231).or.(haxes.eq.321) ) then
           call to_par(acfile,'o2',-((n2-1)/2)*d2)
           call to_par(acfile,'o3',-((n3-1)/2)*d3)
           call to_par(acfile,'o1',o1            )
        else if ( (haxes.eq.312).or.(haxes.eq.132) ) then
           call to_par(acfile,'o1',-((n1-1)/2)*d1)
           call to_par(acfile,'o3',-((n3-1)/2)*d3)
           call to_par(acfile,'o2',o2            )
        end if
     else
        call to_par(acfile,'o1',o1)
        call to_par(acfile,'o2',o2)
        call to_par(acfile,'o3',d3)
     end if
     call to_par(acfile,'label1',label1)
     call to_par(acfile,'label2',label2)
     call to_par(acfile,'label3',label3)
  end if

  ! Determine Helix lengths
  nh   = n1*n2*n3 
  nhp  = pad2(nh)  ! Calculate padlength

  ! Allocate Memory
  allocate( datain(nh), dataout(nh) )

  ! Read data on Helix
  call rsf_read(in,datain)

  ! Optionally transpose the data in the Helix
  if      ( haxes.eq.123 ) then
     dataout = datain
  else if ( haxes.eq.213 ) then
     do ii=1,n1; do jj=1,n2; do kk=1,n3
        dataout( 1 + (jj-1) + n2*( (ii-1) + n1*(kk-1) ) ) = datain( 1 + (ii-1) + n1*( (jj-1) + n2*(kk-1) ) )
     end do;     end do;    end do
  else if ( haxes.eq.132 ) then
     do ii=1,n1; do jj=1,n2; do kk=1,n3
        dataout( 1 + (ii-1) + n1*( (kk-1) + n3*(jj-1) ) ) = datain( 1 + (ii-1) + n1*( (jj-1) + n2*(kk-1) ) )
     end do;     end do;    end do
  else if ( haxes.eq.312 ) then
     do ii=1,n1; do jj=1,n2; do kk=1,n3
        dataout( 1 + (kk-1) + n3*( (ii-1) + n1*(jj-1) ) ) = datain( 1 + (ii-1) + n1*( (jj-1) + n2*(kk-1) ) )
     end do;     end do;    end do
  else if ( haxes.eq.231 ) then
     do ii=1,n1; do jj=1,n2; do kk=1,n3
        dataout( 1 + (jj-1) + n2*( (kk-1) + n3*(ii-1) ) ) = datain( 1 + (ii-1) + n1*( (jj-1) + n2*(kk-1) ) )
     end do;     end do;    end do
  else if ( haxes.eq.321 ) then
     do ii=1,n1; do jj=1,n2; do kk=1,n3
        dataout( 1 + (kk-1) + n3*( (jj-1) + n2*(ii-1) ) ) = datain( 1 + (ii-1) + n1*( (jj-1) + n2*(kk-1) ) )
     end do;     end do;    end do
  end if

  ! Make into Complex
  allocate( cdata(nhp) )
  cdata(1:nh)=cmplx(dataout,0.)

  ! Perform Fourier Transform
  call jftu( 1., nhp, cdata)

  ! Optionally center the spectrum
  datain = 0.
  if (center.eq.1) then
     do kk=1,n3;
        cent = (n2/2-1)*n1+n1/2-1
        if (kk==1) then
           datain(              1+cent :               n1*n2 ) = abs(cdata(              1      :              n1*n2-cent ) )
        else
           datain((kk-1)*n1*n2+ 1      : (kk-1)*n1*n2+ n1*n2 ) = abs(cdata( (kk-1)*n1*n2+1-cent : (kk-1)*n1*n2+n1*n2-cent ) )
        end if
     end do
  else
     datain = abs(cdata)    
  end if

  ! Inverse transposition
  if      ( haxes.eq.123 ) then
     dataout = datain
  else if ( haxes.eq.213 ) then
     do ii=1,n1; do jj=1,n2; do kk=1,n3
        dataout( 1 + (ii-1) + n1*( (jj-1) + n2*(kk-1) ) ) = datain( 1 + (jj-1) + n2*((ii-1) + n1*(kk-1)) )
     end do;     end do;    end do
  else if ( haxes.eq.132 ) then
     do ii=1,n1; do jj=1,n2; do kk=1,n3
        dataout( 1 + (ii-1) + n1*( (jj-1) + n2*(kk-1) ) ) = datain( 1 + (ii-1) + n1*((kk-1) + n3*(jj-1)) )
     end do;     end do;    end do
  else if ( haxes.eq.312 ) then
     do ii=1,n1; do jj=1,n2; do kk=1,n3
        dataout( 1 + (ii-1) + n1*( (jj-1) + n2*(kk-1) ) ) = datain( 1 + (kk-1) + n3*((ii-1) + n1*(jj-1)) )
     end do;     end do;    end do
  else if ( haxes.eq.231 ) then
     do ii=1,n1; do jj=1,n2; do kk=1,n3
        dataout( 1 + (ii-1) + n1*( (jj-1) + n2*(kk-1) ) ) = datain( 1 + (jj-1) + n2*((kk-1) + n3*(ii-1)) )
     end do;     end do;    end do
  else if ( haxes.eq.321 ) then
     do ii=1,n1; do jj=1,n2; do kk=1,n3
        dataout( 1 + (ii-1) + n1*( (jj-1) + n2*(kk-1) ) ) = datain( 1 + (kk-1) + n3*((jj-1) + n2*(ii-1)) )
     end do;     end do;    end do
  end if

  ! Compute Autocorrelation
  cdata=cdata*conjg(cdata)

  allocate( acdata(nhp) )
  acdata = cdata
  call jftu(-1.,nhp,acdata)    ! Inverse FT

  ! Optionally center the spectrum
  datain = 0.
  if (center.eq.1) then
     do kk=1,n3;
        cent = (n2/2-1)*n1+n1/2-1
        if (kk==1) then
           datain(              1+cent :               n1*n2 ) = acdata(              1      :              n1*n2-cent )
        else
           datain((kk-1)*n1*n2+ 1      : (kk-1)*n1*n2+ n1*n2 ) = acdata( (kk-1)*n1*n2+1-cent : (kk-1)*n1*n2+n1*n2-cent )
        end if
     end do
  else
     datain = acdata
  end if
  deallocate( acdata )

  if ( exist_par ("acfile") ) then
     ! Inverse transposition
     if      ( haxes.eq.123 ) then
        dataout = datain
     else if ( haxes.eq.213 ) then
        do ii=1,n1; do jj=1,n2; do kk=1,n3
           dataout( 1 + (ii-1) + n1*( (jj-1) + n2*(kk-1) ) ) = datain( 1 + (jj-1) + n2*( (ii-1) + n1*(kk-1) ) )
        end do;     end do;    end do
     else if ( haxes.eq.132 ) then
        do ii=1,n1; do jj=1,n2; do kk=1,n3
           dataout( 1 + (ii-1) + n1*( (jj-1) + n2*(kk-1) ) ) = datain( 1 + (ii-1) + n1*( (kk-1) + n3*(jj-1) ) )
        end do;     end do;    end do
     else if ( haxes.eq.312 ) then
        do ii=1,n1; do jj=1,n2; do kk=1,n3
           dataout( 1 + (ii-1) + n1*( (jj-1) + n2*(kk-1) ) ) = datain( 1 + (kk-1) + n3*( (ii-1) + n1*(jj-1) ) )
        end do;     end do;    end do
     else if ( haxes.eq.231 ) then
        do ii=1,n1; do jj=1,n2; do kk=1,n3
           dataout( 1 + (ii-1) + n1*( (jj-1) + n2*(kk-1) ) ) = datain( 1 + (jj-1) + n2*( (kk-1) + n3*(ii-1) ) )
        end do;     end do;    end do
     else if ( haxes.eq.321 ) then
        do ii=1,n1; do jj=1,n2; do kk=1,n3
           dataout( 1 + (ii-1) + n1*( (jj-1) + n2*(kk-1) ) ) = datain( 1 + (kk-1) + n3*( (jj-1) + n2*(ii-1) ) )
        end do;     end do;    end do
     end if

     ! Write out MD AC
     call rsf_write(acfile,dataout)
  end if

  ! Call Spectral Factorization code
  call kolmogoroff(nhp,cdata)

  ! Optionally center the spectrum
  datain = 0.
  if (center.eq.1) then
     do kk=1,n3;
        cent = (n2/2-1)*n1+n1/2-1
        if (kk==1) then
           datain(              1+cent :               n1*n2 ) = cdata(              1      :              n1*n2-cent )
        else
           datain((kk-1)*n1*n2+ 1      : (kk-1)*n1*n2+ n1*n2 ) = cdata( (kk-1)*n1*n2+1-cent : (kk-1)*n1*n2+n1*n2-cent )
        end if
     end do
  else
     datain = acdata
  end if

  ! Inverse transposition
  if      ( haxes.eq.123 ) then
     dataout = datain
  else if ( haxes.eq.213 ) then
     do ii=1,n1; do jj=1,n2; do kk=1,n3
        dataout( 1 + (ii-1) + n1*( (jj-1) + n2*(kk-1) ) ) = datain( 1 + (jj-1) + n2*( (ii-1) + n1*(kk-1) ) )
     end do;     end do;    end do
  else if ( haxes.eq.132 ) then
     do ii=1,n1; do jj=1,n2; do kk=1,n3
        dataout( 1 + (ii-1) + n1*( (jj-1) + n2*(kk-1) ) ) = datain( 1 + (ii-1) + n1*( (kk-1) + n3*(jj-1) ) )
     end do;     end do;    end do
  else if ( haxes.eq.312 ) then
     do ii=1,n1; do jj=1,n2; do kk=1,n3
        dataout( 1 + (ii-1) + n1*( (jj-1) + n2*(kk-1) ) ) = datain( 1 + (kk-1) + n3*( (ii-1) + n1*(jj-1) ) )
     end do;     end do;    end do
  else if ( haxes.eq.231 ) then
     do ii=1,n1; do jj=1,n2; do kk=1,n3
        dataout( 1 + (ii-1) + n1*( (jj-1) + n2*(kk-1) ) ) = datain( 1 + (jj-1) + n2*( (kk-1) + n3*(ii-1) ) )
     end do;     end do;    end do
  else if ( haxes.eq.321 ) then
     do ii=1,n1; do jj=1,n2; do kk=1,n3
        dataout( 1 + (ii-1) + n1*( (jj-1) + n2*(kk-1) ) ) = datain( 1 + (kk-1) + n3*( (jj-1) + n2*(ii-1) ) )
     end do;     end do;    end do
  end if

  ! Write out MD SF
  call rsf_write(out,dataout)

  deallocate(datain,dataout,cdata)

end program MDHspecfac

subroutine kolmogoroff( n, cx)  ! Spectral factorization (from pvi)
  integer i,              n       ! input:  cx = spectrum
  complex cx(n)                   ! output: cx = min phase wavelet
  cx(1:n) = clog( cx(1:n) )
  call jftu( -1., n, cx)
  cx(1:n) = cx(1:n)*sqrt(1./n)
  cx(1    ) = cx(1    ) / 2.
  cx(1+n/2) = cx(1+n/2) / 2.
  cx(1+n/2+1:n) = 0.
  call jftu( +1., n, cx) 
  cx(1:n) = cx(1:n)*sqrt(1.*n)
  cx(1:n) = cexp( cx(1:n))
  call jftu( -1., n, cx)
  cx(1:n) = cx(1:n)*sqrt(1./n)
end subroutine kolmogoroff

subroutine jftu( signi, nx, cx )
  !   complex fast fourier transform with unitary scaling
  !
  !               1         nx          signi*2*pi*i*(j-1)*(k-1)/nx
  !   cx(k)  =  -------- * sum cx(j) * e
  !             sqrt(nx)   j=1             for k=1,2,...,nx=2**integer
  !
  integer :: nx, i, j, k, m, istep, pad2
  real    :: signi
  real*8  :: arg,pi
  complex :: cx(nx), cmplx, cw, cdel, ct
  pi=4.0d0*datan(1.0d0)
  if( nx /= pad2(nx) )  call sf_error('ftu: nx not a power of 2')
  cx = cx/sqrt( 1.*nx)
  j = 1;  k = 1
  do i= 1, nx
     if (i<=j) then
        ct = cx(j); cx(j) = cx(i); cx(i) = ct
     end if
     m = nx/2
     do while (j>m .and. m>1) 
        j = j-m; m = m/2 
     end do
     j = j+m
  end do
  do 
     istep = 2*k;   cw = 1.;   arg = signi*pi/k
     cdel = cmplx( cos(arg), sin(arg))
     do m= 1, k
        do i= m, nx, istep
           ct=cw*cx(i+k);  cx(i+k)=cx(i)-ct;  cx(i)=cx(i)+ct
        end do
        cw = cw * cdel
     end do
     k = istep
     if (k>=nx) exit
  end do
end subroutine jftu

integer function pad2(n)
  integer n
  pad2=1
  do
     if (pad2.ge.n) exit
     pad2=pad2*2
  end do
end function pad2
