      subroutine getcoef(n1,d1,w1,n2,d2,w2,coef,ier)
c==========================================================
c
c  Sets up Fourier Sine Transform Multiplier Array for 
c  discrete Helmholtz Operator with Dirichlet conditions 
c
c  The operator is 
c
c    u |---> u + w1^2 d^2 u/dx1^2 + w2^2 d^2 u/dx2^2
c
c  Designed for use with vfftpak vectorized FFTs.
c
c==========================================================
c
      integer n1,n2, ier
      real d1, w1, d2, w2, coef(n1,n2)

      integer i1,i2
      real fred,t1,t2,tol

      data tol /1.e-10/

      if (ier.ne.0) return

      fred=3.1415927

      if (d1.gt.tol) then
         do i1=1,n1
            t1= 2*w1*(1.0 - cos(fred*i1/n1))/(d1*d1)
            do i2=1,n2
               coef(i1,i2)=t1
            end do
         end do
      end if

      if ((d2.gt.tol).and.(n2.gt.1)) then
         do i2=1,n2
            t2=2*w2*(1.0 - cos(fred*i2/n2))/(d2*d2)
            do i1=1,n1
               coef(i1,i2)=coef(i1,i2)+t2
            end do
         end do
      end if

      return
      end
