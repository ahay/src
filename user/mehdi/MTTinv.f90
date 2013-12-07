!  Invert Traveltimes (with variable offset and azimuth) for the Elasticity matrix
!
!
!
!
!!$  Copyright (C) 2013 UT Austin
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


program TTinv
  use rsf
  use frac
  implicit none
  

!#################################################################################################################################
!#################################################################################################################################
!#################################################################################################################################


  type (file)                      :: Cpred, Ccor,  commonoffset
  type(axa)  :: at,az   ! cube axes

  integer                          :: n1, n2, itmax, of, azi, jjj, hhh, cc, p, t, azimuth, offset, I, nu
  real, parameter :: pi = 3.1415927
  real, allocatable :: DC11(:,:), DC12(:,:) , DC13(:,:), DC22(:,:), DC23(:,:), DC33(:,:), DC44(:,:), DC45(:,:), DC55(:,:), DC66(:,:)
  real, allocatable :: DC16(:,:), DC26(:,:), DC36(:,:), rV2_cor(:,:), rV20(:,:), a1(:), b1(:), c1(:), d1(:), e1(:)
  real, allocatable :: f1(:), g1(:), h1(:), q1(:), r1(:), s1(:), t1(:), v1(:), z1(:), rv0222(:), rv2_corrr(:)
  real, allocatable :: del_rv(:), F(:,:), C(:), C_pred(:), C_cor(:), Del(:), FF(:,:),F1D(:), cons(:), weight(:), FFT(:,:)

  real :: C_cor11, C_cor12, C_cor13, C_cor22, C_cor23, C_cor33, C_cor44, C_cor45, C_cor55, C_cor66,  C_cor16, C_cor26, C_cor36, C11
  real :: C12, C13, C22, C23, C33, C44, C45, C55, C66, C16, C26, C36, DelC11, DelC12, DelC13, DelC22, DelC23, DelC33, DelC44,DelC45
  real :: DelC55, DelC66, DelC16, DelC26, DelC36, ct, tet, phi, x1, x2, x3, tolerancer


  call sf_init()            ! initialize RSF
   Ccor = rsf_output("CCorrect") ! common_offset is the output file name for common offset section
   Cpred = rsf_output("Cpredicted") ! common_offset is the output file name for common offset section
   commonoffset = rsf_output("common") ! common_offset is the output file name for common offset section

  call from_par("azimuth",p)   ! command-line parameter
  call from_par("offset",t)   ! command-line parameter


  call to_par(Cpred,"n1",1)     ! write the sample number for the 2nd axis for the output file
  call to_par(Cpred,"o1",1)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file
  call to_par(Cpred,"d1",1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file

  call to_par(Ccor,"n1",1)     ! write the sample number for the 2nd axis for the output file
  call to_par(Ccor,"o1",1)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file
  call to_par(Ccor,"d1",1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file

  call to_par(commonoffset,"n1",azimuth)     ! write the sample number for the 2nd axis for the output file
  call to_par(commonoffset,"o1",1)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file
  call to_par(commonoffset,"d1",1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file


  call to_par(commonoffset,"n2",offset)     ! write the sample number for the 2nd axis for the output file
  call to_par(commonoffset,"o2",1)      ! write the starting point for the 2nd axis (here, starting value is 0) for the output file
  call to_par(commonoffset,"d2",1)      ! write the increament for the 2nd axis (here, increament is 1) for the output file


  allocate (DC11(t,p))
  allocate (DC12(t,p))
  allocate (DC13(t,p))
  allocate (DC22(t,p))
  allocate (DC23(t,p))
  allocate (DC33(t,p))
  allocate (DC44(t,p))
  allocate (DC45(t,p))
  allocate (DC55(t,p))
  allocate (DC66(t,p))
  allocate (DC16(t,p))
  allocate (DC26(t,p))
  allocate (DC36(t,p))
  allocate (rV20(t,p))
  allocate (rV2_cor(t,p))


  itmax=300         !number of iterations for Conjugate Gradient

  ct=.5

  C_cor11=ct*24.54
  C_cor12=ct*3.89
  C_cor13=ct*2.38
  C_cor22=ct*28.79
  C_cor23=ct*2.37
  C_cor33=ct*29.87
  C_cor44=ct*11.24
  C_cor45=ct*1.
  C_cor55=ct*10.73
  C_cor66=ct*10.98
  C_cor16=ct*(-3.18)
  C_cor26=ct*1.13
  C_cor36=ct*.2


  C11=20.
  C12=5.
  C13=6.
  C22=32.
  C23=5.
  C33=34.
  C44=15.
  C45=3.
  C55=15.
  C66=15.
  C16=6.
  C26=5.
  C36=5.

  DelC11=C_cor11-C11
  DelC12=C_cor12-C12
  DelC13=C_cor13-C13
  DelC22=C_cor22-C22
  DelC23=C_cor23-C23
  DelC33=C_cor33-C33
  DelC44=C_cor44-C44
  DelC45=C_cor45-C45
  DelC55=C_cor55-C55
  DelC66=C_cor66-C66
  DelC16=C_cor16-C16
  DelC26=C_cor26-C26
  DelC36=C_cor36-C36



  do t=1,40
    do p=1,360

    tet=t*pi/180.
    phi=p*pi/180.

    x1 = sin(tet)*cos(phi)
    x2 = sin(tet)*sin(phi)
    x3 = cos(tet)




          DC11(t,p)=(1./3.)*(x1**2. + (sqrt(3.)*((-x1**2.)*(C55*x1**2. + 2.*C45*x1*x2 &
		    + C44*x2**2. + C33*x3**2.) -&
           	    x1**2.*(C66*x1**2. + 2.*C26*x1*x2 + &
             	    C22*x2**2. + C44*x3**2.) + &
      		    (2./3.)*x1**2.*(C11*x1**2. + C55*x1**2. + C66*x1**2. + 2.*C16*x1*x2 + &
		    2.*C26*x1*x2 + 2.*C45*x1*x2 + C22*x2**2. + &
           	    C44*x2**2. + C66*x2**2. +&
        	    (C33 + C44 + C55)*x3**2.)))/sqrt(((C36 + C45)*x1*x3 + &
		    (C23 + C44)*x2*x3)**2. + ((C13 + C55)*x1*x3 +&
               	    (C36 + C45)*x2*x3)**2. - &
      		    (C55*x1**2. + 2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.)*(C66*x1**2. + &
		    2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.) + &
     		    (C16*x1**2. + x2*((C12 + C66)*x1 + C26*x2) + C45*x3**2.)**2. - &
		    (C55*x1**2. + &
		     2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.)*&
       		     (C11*x1**2. + 2.*C16*x1*x2 + C66*x2**2. + C55*x3**2.) - (C66*x1**2. &
		     + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.)*&
       		    (C11*x1**2. + 2.*C16*x1*x2 + C66*x2**2. + C55*x3**2.) + &
		    (1./3.)*(C11*x1**2. + C55*x1**2. + C66*x1**2. + 2.*C16*x1*x2 + &
                    2.*C26*x1*x2 + 2.*C45*x1*x2 + &
         	    C22*x2**2. + C44*x2**2. + C66*x2**2. + (C33 + C44 + C55)*x3**2.)**2.))



    rV20(t,p)=(1./3.)*(2.*sqrt(3.)*sqrt((1./3.)*(C11*x1**2. + 2.*C16*x1*x2 + C22*x2**2. + 2.*C26*x1*x2 + &
              x3**2.*(C33 + C44 + C55) + C44*x2**2. +&
              2.*C45*x1*x2 + C55*x1**2. + C66*x1**2. + C66*x2**2.)**2. - (C11*x1**2. + 2.*C16*x1*x2 + C55*x3**2. + &
		C66*x2**2.)*(C22*x2**2. + &
              2.*C26*x1*x2 + C44*x3**2. + C66*x1**2.) - (C11*x1**2. + 2.*C16*x1*x2 + C55*x3**2. + C66*x2**2.)*&
		(C33*x3**2. + C44*x2**2. + &
              2.*C45*x1*x2 + C55*x1**2.) +  (x2*(x1*(C12 + C66) + C26*x2) + C16*x1**2. + C45*x3**2.)**2. +&
		 (x1*x3*(C13 + C55) + &
              x2*x3*(C36 + C45))**2. - (C22*x2**2. + 2.*C26*x1*x2 + C44*x3**2. + C66*x1**2.)*(C33*x3**2. + &
		C44*x2**2. + 2.*C45*x1*x2 + &
              C55*x1**2.) + (x2*x3*(C23 + C44) + x1*x3*(C36 + C45))**2.) + C11*x1**2. + 2.*C16*x1*x2 + &
		C22*x2**2. + 2.*C26*x1*x2 + &
              C33*x3**2. + C44*x2**2. + C44*x3**2. + 2.*C45*x1*x2 + C55*x1**2. + C55*x3**2. + &
		C66*x1**2. + C66*x2**2.)

      rV2_cor(t,p)=(1./3.)*(2.*sqrt(3.)*sqrt((1./3.)*(C_cor11*x1**2. + 2.*C_cor16*x1*x2 + &
          C_cor22*x2**2. + 2.*C_cor26*x1*x2 + x3**2.*(C_cor33 &
      + C_cor44 + C_cor55) + C_cor44*x2**2. + 2.*C_cor45*x1*x2 + C_cor55*x1**2. + C_cor66*x1**2. & 
      + C_cor66*x2**2.)**2. -(C_cor11*x1**2. + 2.*C_cor16*x1*x2 + C_cor55*x3**2. + C_cor66*x2**2.)*&
           (C_cor22*x2**2. + 2.*C_cor26*x1*x2 + &
      C_cor44*x3**2. + C_cor66*x1**2.) - (C_cor11*x1**2. + 2.*C_cor16*x1*x2 + C_cor55*x3**2. &
      + C_cor66*x2**2.)*(C_cor33*x3**2. + C_cor44*x2**2. + 2.*C_cor45*x1*x2 + C_cor55*x1**2.) + &
          (x2*(x1*(C_cor12 + C_cor66) + C_cor26*x2) &
      + C_cor16*x1**2. + C_cor45*x3**2.)**2. + (x1*x3*(C_cor13 + C_cor55) & 
      + x2*x3*(C_cor36 + C_cor45))**2. - (C_cor22*x2**2. + 2.*C_cor26*x1*x2 + C_cor44*x3**2. + C_cor66*x1**2.)* &
	(C_cor33*x3**2. +&
      C_cor44*x2**2. + 2.*C_cor45*x1*x2 + C_cor55*x1**2.) + (x2*x3*(C_cor23 + C_cor44) + &
      x1*x3*(C_cor36 + C_cor45))**2.) + C_cor11*x1**2. + 2.*C_cor16*x1*x2 + C_cor22*x2**2. + 2.*C_cor26*x1*x2 + &
	C_cor33*x3**2. + &
      C_cor44*x2**2. + C_cor44*x3**2. + 2.*C_cor45*x1*x2 + C_cor55*x1**2. + C_cor55*x3**2. &
      + C_cor66*x1**2. + C_cor66*x2**2.)






          DC12(t,p)=(2*x1*x2*(C16*x1**2. + x2*((C12 + C66)*x1 + C26*x2) + C45*x3**2.))/&
         (sqrt(3.)*sqrt(((C36 + C45)*x1*x3 + (C23 + C44)*x2*x3)**2. + ((C13 + C55)*x1*x3 + (C36 + C45)*x2*x3)**2. - &
        (C55*x1**2. + 2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.)*(C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.) + &
        (C16*x1**2. + x2*((C12 + C66)*x1 + C26*x2) + C45*x3**2.)**2. - (C55*x1**2. + 2.*C45*x1*x2 + &
              C44*x2**2. + C33*x3**2.)*&
      (C11*x1**2. + 2.*C16*x1*x2 + C66*x2**2. + C55*x3**2.) - (C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.)*&
      (C11*x1**2. + 2.*C16*x1*x2 + C66*x2**2. + C55*x3**2.) + (1./3.)*(C11*x1**2. + C55*x1**2. + C66*x1**2. +&
	 2.*C16*x1*x2 + &
         2.*C26*x1*x2 + 2.*C45*x1*x2 + &
        C22*x2**2. + C44*x2**2. + C66*x2**2. + (C33 + C44 + C55)*x3**2.)**2.))




         DC13(t,p)=(2*x1*x3*((C13 + C55)*x1*x3 + (C36 + C45)*x2*x3))/&
               (sqrt(3.)*sqrt(((C36 + C45)*x1*x3 + (C23 + C44)*x2*x3)**2. + ((C13 + C55)*x1*x3 + &
              (C36 + C45)*x2*x3)**2. - &
     (C55*x1**2. + 2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.)*(C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.) + &
     (C16*x1**2. + x2*((C12 + C66)*x1 + C26*x2) + C45*x3**2.)**2. - (C55*x1**2. + 2.*C45*x1*x2 + C44*x2**2. &
	+ C33*x3**2.)*&
      (C11*x1**2. + 2.*C16*x1*x2 + C66*x2**2. + C55*x3**2.) - (C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.)*&
      (C11*x1**2. + 2.*C16*x1*x2 + C66*x2**2. + C55*x3**2.) + (1./3.)*(C11*x1**2. + C55*x1**2. + C66*x1**2. + &
	2.*C16*x1*x2 +&
         2.*C26*x1*x2 + 2.*C45*x1*x2 + &
        C22*x2**2. + C44*x2**2. + C66*x2**2. + (C33 + C44 + C55)*x3**2.)**2.))



        DC22(t,p)=(1./3.)*(x2**2. + (sqrt(3.)*((-x2**2.)*(C55*x1**2. + 2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.) -&
            x2**2.*(C11*x1**2. +&
          2.*C16*x1*x2 + C66*x2**2. + C55*x3**2.) + &
      (2./3.)*x2**2.*(C11*x1**2. + C55*x1**2. + C66*x1**2. + 2.*C16*x1*x2 + 2.*C26*x1*x2 + 2.*C45*x1*x2 + C22*x2**2. &
               + C44*x2**2. + C66*x2**2. + &
        (C33 + C44 + C55)*x3**2.)))/sqrt(((C36 + C45)*x1*x3 + (C23 + C44)*x2*x3)**2. + ((C13 + C55)*x1*x3 &
             + (C36 + C45)*x2*x3)**2. - &
      (C55*x1**2. + 2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.)*(C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.) + &
      (C16*x1**2. + x2*((C12 + C66)*x1 + C26*x2) + C45*x3**2.)**2. - (C55*x1**2. + 2.*C45*x1*x2 &
           + C44*x2**2. + C33*x3**2.)*&
       (C11*x1**2. + 2.*C16*x1*x2 + C66*x2**2. + C55*x3**2.) - (C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.)*&
       (C11*x1**2. + 2.*C16*x1*x2 + C66*x2**2. + C55*x3**2.) + (1./3.)*(C11*x1**2. + C55*x1**2. &
                + C66*x1**2. + 2.*C16*x1*x2 + &
           2.*C26*x1*x2 + 2.*C45*x1*x2 + &
         C22*x2**2. + C44*x2**2. + C66*x2**2. + (C33 + C44 + C55)*x3**2.)**2.))





         DC23(t,p)=(2*x2*x3*((C36 + C45)*x1*x3 + (C23 + C44)*x2*x3))/&
                    (sqrt(3.)*sqrt(((C36 + C45)*x1*x3 + (C23 + C44)*x2*x3)**2. + ((C13 + C55)*x1*x3 + &
               (C36 + C45)*x2*x3)**2. - &
           (C55*x1**2. + 2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.)*(C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.) + &
           (C16*x1**2. + x2*((C12 + C66)*x1 + C26*x2) + C45*x3**2.)**2. - (C55*x1**2. + 2.*C45*x1*x2 + &
		C44*x2**2. + C33*x3**2.)*&
          (C11*x1**2. + 2.*C16*x1*x2 + C66*x2**2. + C55*x3**2.) - (C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.)*&
           (C11*x1**2. + 2.*C16*x1*x2 + C66*x2**2. + C55*x3**2.) + (1./3.)*(C11*x1**2. + C55*x1**2. + C66*x1**2. &
            + 2.*C16*x1*x2 + 2.*C26*x1*x2 + 2.*C45*x1*x2 + &
                C22*x2**2. + C44*x2**2. + C66*x2**2. + (C33 + C44 + C55)*x3**2.)**2.))



        DC33(t,p)=(1./3.)*(x3**2. + (sqrt(3.)*((-x3**2.)*(C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.) - &
                x3**2.*(C11*x1**2. + 2.*C16*x1*x2 + C66*x2**2. + C55*x3**2.) + &
      (2./3.)*x3**2.*(C11*x1**2. + C55*x1**2. + C66*x1**2. + 2.*C16*x1*x2 + 2.*C26*x1*x2 + 2.*C45*x1*x2 + C22*x2**2.&
                + C44*x2**2. + C66*x2**2. + &
        (C33 + C44 + C55)*x3**2.)))/sqrt(((C36 + C45)*x1*x3 + (C23 + C44)*x2*x3)**2. + ((C13 + C55)*x1*x3&
            + (C36 + C45)*x2*x3)**2. - &
      (C55*x1**2. + 2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.)*(C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.) + &
      (C16*x1**2. + x2*((C12 + C66)*x1 + C26*x2) + C45*x3**2.)**2. - (C55*x1**2. + 2.*C45*x1*x2 + &
	C44*x2**2. + C33*x3**2.)*&
       (C11*x1**2. + 2.*C16*x1*x2 + C66*x2**2. + C55*x3**2.) - (C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.)*&
       (C11*x1**2. + 2.*C16*x1*x2 + C66*x2**2. + C55*x3**2.) + (1./3.)*(C11*x1**2. + C55*x1**2. + C66*x1**2.&
             + 2.*C16*x1*x2 + 2.*C26*x1*x2 + 2.*C45*x1*x2 + &
         C22*x2**2. + C44*x2**2. + C66*x2**2. + (C33 + C44 + C55)*x3**2.)**2.))





         DC44(t,p)=(1./3.)*(x2**2. + x3**2. + (sqrt(3.)*(2*x2*x3*((C36 + C45)*x1*x3 + (C23 + C44)*x2*x3)&
		 - x3**2.*(C55*x1**2.&
                + 2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.) - &
               x2**2.*(C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.) - x2**2.*(C11*x1**2. + &
		2.*C16*x1*x2 + C66*x2**2. + C55*x3**2.) - &
             x3**2.*(C11*x1**2. + 2.*C16*x1*x2 + C66*x2**2. + C55*x3**2.) + (2./3.)*(x2**2. + x3**2.)*&
		(C11*x1**2. + C55*x1**2. + C66*x1**2. +&
              2.*C16*x1*x2 + 2.*C26*x1*x2 + &
           2.*C45*x1*x2 + C22*x2**2. + C44*x2**2. + C66*x2**2. + (C33 + C44 + C55)*x3**2.)))/&
             sqrt(((C36 + C45)*x1*x3 + (C23 + C44)*x2*x3)**2. + ((C13 + C55)*x1*x3 + (C36 + C45)*x2*x3)**2. &
		- (C55*x1**2. + &
               2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.)*&
            (C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.) + (C16*x1**2. + x2*((C12 + C66)*x1 + C26*x2) &
		+ C45*x3**2.)**2. - &
              (C55*x1**2. + 2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.)*(C11*x1**2. + 2.*C16*x1*x2 +&
		 C66*x2**2. + C55*x3**2.) - &
           (C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.)*(C11*x1**2. + 2.*C16*x1*x2 + C66*x2**2. &
		+ C55*x3**2.) + &
            (1./3.)*(C11*x1**2. + C55*x1**2. + C66*x1**2. + 2.*C16*x1*x2 + 2.*C26*x1*x2 + 2.*C45*x1*x2 + &
		C22*x2**2. + C44*x2**2. + &
               C66*x2**2. + (C33 + C44 + C55)*x3**2.)**&
              2))




          DC45(t,p)=(1./3.)*(2.*x1*x2 + (sqrt(3.)*(2.*x1*x3*((C36 + C45)*x1*x3 + (C23 + C44)*x2*x3) + &
		2.*x2*x3*((C13 + C55)*x1*x3&
              + (C36 + C45)*x2*x3) - &
              2.*x1*x2*(C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.) + 2.*x3**2.*(C16*x1**2. &
		+ x2*((C12 + C66)*x1 +&
               C26*x2) + C45*x3**2.) - &
             2.*x1*x2*(C11*x1**2. + 2.*C16*x1*x2 + C66*x2**2. + C55*x3**2.) + (4./3.)*x1*x2*(C11*x1**2. &
		+ C55*x1**2. + C66*x1**2. +&
              2.*C16*x1*x2 + 2.*C26*x1*x2 + &
             2.*C45*x1*x2 + C22*x2**2. + C44*x2**2. + C66*x2**2. + (C33 + C44 + C55)*x3**2.)))/&
               sqrt(((C36 + C45)*x1*x3 + (C23 + C44)*x2*x3)**2. + ((C13 + C55)*x1*x3 + &
		(C36 + C45)*x2*x3)**2. - (C55*x1**2. + &
               2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.)*&
               (C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.) + (C16*x1**2. + x2*((C12 + &
		C66)*x1 + C26*x2) + C45*x3**2.)**2. - &
               (C55*x1**2. + 2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.)*(C11*x1**2. + 2.*C16*x1*x2 +&
		 C66*x2**2. + C55*x3**2.) - &
              (C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.)*(C11*x1**2. + 2.*C16*x1*x2 +&
		 C66*x2**2. + C55*x3**2.) + &
               (1./3.)*(C11*x1**2. + C55*x1**2. + C66*x1**2. + 2.*C16*x1*x2 + 2.*C26*x1*x2 + &
		2.*C45*x1*x2 + C22*x2**2. + C44*x2**2. + &
                C66*x2**2. + (C33 + C44 + C55)*x3**2.)**&
                2.))



           DC55(t,p)=(1./3.)*(x1**2. + x3**2. + (sqrt(3.)*(2*x1*x3*((C13 + C55)*x1*x3 + (C36 + C45)*x2*x3)&
		 - x3**2.*(C55*x1**2. + &
               2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.) - &
             x1**2.*(C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.) - x3**2.*(C66*x1**2. + &
		2.*C26*x1*x2 + C22*x2**2. +&
              C44*x3**2.) - &
              x1**2.*(C11*x1**2. + 2.*C16*x1*x2 + C66*x2**2. + C55*x3**2.) + (2./3.)*(x1**2. + &
		x3**2.)*(C11*x1**2. + C55*x1**2. + &
                C66*x1**2. + 2.*C16*x1*x2 + 2.*C26*x1*x2 + &
              2.*C45*x1*x2 + C22*x2**2. + C44*x2**2. + C66*x2**2. + (C33 + C44 + C55)*x3**2.)))/&
              sqrt(((C36 + C45)*x1*x3 + (C23 + C44)*x2*x3)**2. + ((C13 + C55)*x1*x3 + &
		(C36 + C45)*x2*x3)**2. - &
                 (C55*x1**2. + 2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.)*&
              (C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.) + (C16*x1**2. + x2*((C12 +&
		 C66)*x1 + C26*x2) + C45*x3**2.)**2. - &
               (C55*x1**2. + 2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.)*(C11*x1**2. + 2.*C16*x1*x2 &
		+ C66*x2**2. + C55*x3**2.) - &
              (C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.)*(C11*x1**2. + 2.*C16*x1*x2 &
		+ C66*x2**2. + C55*x3**2.) + &
              (1./3.)*(C11*x1**2. + C55*x1**2. + C66*x1**2. + 2.*C16*x1*x2 + 2.*C26*x1*x2 +&
		 2.*C45*x1*x2 + C22*x2**2. + C44*x2**2. +&
                C66*x2**2. + (C33 + C44 + C55)*x3**2.)**&
                2))




                DC66(t,p)=(1./3.)*(x1**2. + x2**2. + (sqrt(3.)*((-x1**2.)*(C55*x1**2. + &
		2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.) - &
                x2**2.*(C55*x1**2. + 2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.) - &
               x2**2.*(C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.) + 2*x1*x2*(C16*x1**2.&
		 + x2*((C12 + &
              C66)*x1 + C26*x2) + C45*x3**2.) - &
             x1**2.*(C11*x1**2. + 2.*C16*x1*x2 + C66*x2**2. + C55*x3**2.) + (2./3.)*(x1**2. +&
		 x2**2.)*(C11*x1**2. + &
              C55*x1**2. + C66*x1**2. + 2.*C16*x1*x2 + 2.*C26*x1*x2 + &
             2.*C45*x1*x2 + C22*x2**2. + C44*x2**2. + C66*x2**2. + (C33 + C44 + C55)*x3**2.)))/&
             sqrt(((C36 + C45)*x1*x3 + (C23 + C44)*x2*x3)**2. + ((C13 + C55)*x1*x3 + (C36 &
		+ C45)*x2*x3)**2. - &
             (C55*x1**2. + 2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.)*&
              (C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.) + (C16*x1**2. + x2*((C12&
		 + C66)*x1 + C26*x2) + &
              C45*x3**2.)**2. - &
             (C55*x1**2. + 2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.)*(C11*x1**2. + 2.*C16*x1*x2 &
		+ C66*x2**2. + C55*x3**2.) - &
             (C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.)*(C11*x1**2. + 2.*C16*x1*x2 + &
		C66*x2**2. + C55*x3**2.) + &
             (1./3.)*(C11*x1**2. + C55*x1**2. + C66*x1**2. + 2.*C16*x1*x2 + 2.*C26*x1*x2 +&
		 2.*C45*x1*x2 + C22*x2**2. + &
                C44*x2**2. + C66*x2**2. + (C33 + C44 + C55)*x3**2.)**&
               2))




              DC16(t,p)=(1./3.)*(2.*x1*x2 + (sqrt(3.)*(-2*x1*x2*(C55*x1**2. + 2.*C45*x1*x2 + &
		C44*x2**2. + C33*x3**2.) - &
               2.*x1*x2*(C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.) + &
              2.*x1**2.*(C16*x1**2. + x2*((C12 + C66)*x1 + C26*x2) + C45*x3**2.) + &
		(4./3.)*x1*x2*(C11*x1**2. + C55*x1**2. &
              + C66*x1**2. + 2.*C16*x1*x2 + 2.*C26*x1*x2 + &
               2.*C45*x1*x2 + C22*x2**2. + C44*x2**2. + C66*x2**2. + (C33 + C44 + C55)*x3**2.)))/&
                 sqrt(((C36 + C45)*x1*x3 + (C23 + C44)*x2*x3)**2. + ((C13 + C55)*x1*x3 +&
		 (C36 + C45)*x2*x3)**2. -&
               (C55*x1**2. + 2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.)*&
               (C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.) + (C16*x1**2. +&
		 x2*((C12 + C66)*x1 + C26*x2) &
               + C45*x3**2.)**2. - &
             (C55*x1**2. + 2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.)*(C11*x1**2. + 2.*C16*x1*x2&
		 + C66*x2**2. + C55*x3**2.) - &
             (C66*x1**2. + 2.*C26*x1*x2 + C22*x2**2. + C44*x3**2.)*(C11*x1**2. + 2.*C16*x1*x2&
		 + C66*x2**2. + C55*x3**2.) + &
              (1./3.)*(C11*x1**2. + C55*x1**2. + C66*x1**2. + 2.*C16*x1*x2 + 2.*C26*x1*x2 &
		+ 2.*C45*x1*x2 + C22*x2**2. + &
               C44*x2**2. + C66*x2**2. + (C33 + C44 + C55)*x3**2.)**&
               2.))





            DC26(t,p)=(1./3.)*(2.*x1*x2 + (sqrt(3.)*(-2.*x1*x2*(C55*x1**2. + 2.*C45*x1*x2 + &
		C44*x2**2. + C33*x3**2.) + &
               2.*x2**2.*(C16*x1**2. + x2*((C12 + C66)*x1 + C26*x2) + C45*x3**2.) - &
		2.*x1*x2*(C11*x1**2. + 2.*C16*x1*x2 + &
             C66*x2**2. + C55*x3**2.) + &
              (4./3.)*x1*x2*(C11*x1**2. + C55*x1**2. + C66*x1**2. + 2.*C16*x1*x2 + 2.*C26*x1*x2 &
		+ 2.*C45*x1*x2 + C22*x2**2. &
             + C44*x2**2. + C66*x2**2. + &
                 (C33 + C44 + C55)*x3**2.)))/sqrt(((C36 + C45)*x1*x3 + (C23 + C44)*x2*x3)**2. +&
              ((C13 + C55)*x1*x3 + (C36 + C45)*x2*x3)**2. - &
              (C55*x1**2. + 2.*C45*x1*x2 + C44*x2**2. + C33*x3**2.)*(C66*x1**2. + 2.*C26*x1*x2 &
		+ C22*x2**2.&
               + C44*x3**2.) + &
             (C16*x1**2. + x2*((C12 + C66)*x1 + C26*x2) + C45*x3**2.)**2. - (C55*x1**2. +&
		 2.*C45*x1*x2 + &
             C44*x2**2. + C33*x3**2.)*&
              (C11*x1**2. + 2.*C16*x1*x2 + C66*x2**2. + C55*x3**2.) - (C66*x1**2. &
		+ 2.*C26*x1*x2 +&
                C22*x2**2. + C44*x3**2.)*&
               (C11*x1**2. + 2.*C16*x1*x2 + C66*x2**2. + C55*x3**2.) + (1./3.)*(C11*x1**2.&
		 + C55*x1**2. + &
                C66*x1**2. + 2.*C16*x1*x2 + 2.*C26*x1*x2 + 2.*C45*x1*x2 + &
               C22*x2**2. + C44*x2**2. + C66*x2**2. + (C33 + C44 + C55)*x3**2.)**2.))





             DC36(t,p)=(2.*x2*x3*(x1*x3*(C13 + C55) + x2*x3*(C36 + C45)) + 2.*x1*x3*(x2*x3*(C23&
		 + C44) + x1*x3*(C36 + C45)))/&
            (sqrt(3.)*sqrt((1./3.)*(C11*x1**2. + 2.*C16*x1*x2 + C22*x2**2. + 2.*C26*x1*x2 + &
		x3**2.*(C33 + C44 + C55) + C44*x2**2. + &
             2.*C45*x1*x2 + C55*x1**2. + C66*x1**2. + C66*x2**2.)**2. - (C11*x1**2. +&
		 2.*C16*x1*x2 + C55*x3**2. + C66*x2**2.)*&
             (C22*x2**2. + 2.*C26*x1*x2 + C44*x3**2. + C66*x1**2.) - (C11*x1**2. + &
		2.*C16*x1*x2 + C55*x3**2. + C66*x2**2.)*&
            (C33*x3**2. + C44*x2**2. + 2.*C45*x1*x2 + C55*x1**2.) + (x2*(x1*(C12 + C66) +&
		 C26*x2) + C16*x1**2. + C45*x3**2.)**2. + &
               (x1*x3*(C13 + C55) + x2*x3*(C36 + C45))**2. - (C22*x2**2. + 2.*C26*x1*x2 + &
		C44*x3**2. + C66*x1**2.)*&
               (C33*x3**2. + C44*x2**2. + 2.*C45*x1*x2 + C55*x1**2.) + (x2*x3*(C23 + C44)&
		 + x1*x3*(C36 + C45))**2.))



    end do
  end do

   !write (0,*) 'dim cmof', shape(rv20),'            dim vol', shape(rv2_cor)

!%%%%%%%%%%%%%%%%%%%%%%%%%%                            %%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%    			       %%%%%%%%%%%%%%%%%%%%%%
     allocate (Del (13))
     allocate (C_cor (13))
     allocate (C (13))

  Del=(/ DelC11,DelC12,DelC13,DelC22,DelC23,DelC33,DelC44,DelC45,DelC55,DelC66,DelC16,DelC26,DelC36 /)

     
  C_cor=(/C_cor11,C_cor12,C_cor13,C_cor22,C_cor23,C_cor33,C_cor44,C_cor45,C_cor55,C_cor66,C_cor16,C_cor26,C_cor36/)

  C=(/C11,C12,C13,C22,C23,C33,C44,C45,C55,C66,C16,C26,C36/)
   !write (0,*) 'dim cmof', C_cor, C


  jjj=size(DC11,2)
  hhh=size(DC11,1)
  cc=hhh*jjj
  

  
  allocate (a1 (cc))
  allocate (b1 (cc))
  allocate (c1 (cc))
  allocate (d1 (cc))
  allocate (e1 (cc))
  allocate (f1 (cc))
  allocate (g1 (cc))
  allocate (h1 (cc))
  allocate (q1 (cc))
  allocate (r1 (cc))
  allocate (s1 (cc))
  allocate (t1 (cc))
  allocate (v1 (cc)); v1=0.; 
  allocate (z1 (cc))
  allocate (rv0222 (cc))
  allocate (rv2_corrr (cc))
  allocate (del_rv (cc))
  allocate (F1D (cc*13))
  allocate (FF (13,cc))


  a1=reshape(DC11,(/cc/))
  b1=reshape(DC12,(/cc/))
  c1=reshape(DC13,(/cc/))
  d1=reshape(DC22,(/cc/))
  e1=reshape(DC23,(/cc/))
  f1=reshape(DC33,(/cc/))
  g1=reshape(DC44,(/cc/))
  h1=reshape(DC45,(/cc/))
  q1=reshape(DC55,(/cc/))
  r1=reshape(DC66,(/cc/))
  s1=reshape(DC16,(/cc/))
  v1=reshape(DC26,(/cc/))
  z1=reshape(DC36,(/cc/))

  rv0222=reshape(rV20,(/cc/))
  rv2_corrr=reshape(rV2_cor,(/cc/))
    

!   write (0,*)  'a1', a1, 'dim', shape(a1)
  del_rv=rv2_corrr-rv0222

  F1D=(/ a1,b1,c1,d1,e1,f1,g1,h1,q1,r1,s1,v1,z1 /)    !matrix of F values

 !  write (0,*) 'dim F1D', shape(F1D)

  FF=reshape(F1D,(/ 13,cc /), order=(/2,1/)) 			! reshape 2D to 3D

  cons=matmul(FF,del_rv)
!!   write (0,*) 'dim F1D', shape(F1D),'            dim FF', shape(FF)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                      %%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Inversion using Conjugate Gradient %%%%%%%%%

  tolerancer=.0001; !% stopping criteria for CG

  I=30;    !           %number of iterations to get good initial estimates for Conjugate Gradient
  allocate (FFT (13,13))
  FFT=matmul(FF,transpose(FF))	

  call cg(nu,cons,FFT,tolerancer,itmax,I,weight)
  C_pred=C+weight

!#################################################################################################################################


  call rsf_write(commonoffset,rv20)   ! write out the common offset section
  
  call rsf_write(Ccor,c_cor)     ! write out the 3D volume
  call rsf_write(Cpred,c_pred)     ! write out the 3D volume



end program TTinv





