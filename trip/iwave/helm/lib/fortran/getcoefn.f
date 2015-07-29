      subroutine getcoefn(n1,d1,w1,n2,d2,w2,coef,ier)
c==========================================================
c
c  Sets up Fourier Sine/Cosine Transform Multiplier Array for 
c  discrete Helmholtz Operator with Dirichlet conditions on
c  top and bottom of rectangular grid and Neumann conditions
c  on sides.
c
c  The operator is 
c
c    u |---> u + w1^2 d^2 u/dx1^2 + w2^2 d^2 u/dx2^2
c
c  Designed for use with vfftpak vectorized FFTs.
c
c==========================================================
c
c      include 'system.h'

      integer n1,n2, ier
      real d1, w1, d2, w2, coef(n1,n2)

      integer i1,i2
      real fred,t1,t2,tol

      data tol /1.e-10/

      if (ier.ne.0) return

      if (abs(d1).lt.tol) return
c         write(ipdmp,*)' Error: GETCOEFN'
c         write(ipdmp,*)' d1 = ',d1
c         write(ipdmp,*)' less than tolerance currently hardwired'
c         write(ipdmp,*)' into this routine: tol = ',tol
c         ier=123
c         return
c      end if

      if (abs(d2).lt.tol) return
c         write(ipdmp,*)' Error: GETCOEFN'
c         write(ipdmp,*)' d2 = ',d2
c         write(ipdmp,*)' less than tolerance currently hardwired'
c         write(ipdmp,*)' into this routine: tol = ',tol
c         ier=123
c         return
c      end if

      fred=3.1415927

      do i1=1,n1
         t1= 2*w1*(1.0 - cos(fred*i1/n1))/(d1*d1)
         do i2=1,n2
            coef(i1,i2)=t1
         end do
      end do

c note that sum here starts at 2 because the first cosine transform
c coefficient is 0 (dc component!). note that this convention also
c neatly takes care of the possible difficulty with n2=1.

      do i2=2,n2
         t2=2*w2*(1.0 - cos(fred*(i2-1)/(n2-1)))/(d2*d2)
         do i1=1,n1
            coef(i1,i2)=coef(i1,i2)+t2
         end do
      end do

      return
      end
