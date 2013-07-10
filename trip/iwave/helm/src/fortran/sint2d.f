      subroutine sint2d(m,n,x,xt,wsave)
c
c Interface with vfftpk routines to implement 2d fast sine transform
c of an input array x(1:m,1:n), stored as a linear array x(*).
c Uses workspace arrays xt(*) and wsave(*).
c
c VERY IMPORTANT: it is assumed that the columns of the
c input array x satisfy Dirichlet conditions on the bottom and
c on the right both ends, i.e. that x(:,n)=x(m,:)=0. 
c Thus the nonredundant part of x is, the (m-1)x(n-1) principal 
c submatrix. Therefore x is passed to
c the vfftpak routine VSINT as an (m-1)x(n-1) array, with
c mdimx = m. This automatically takes care of the need in
c VSINT for an extra column of workspace. As a result, the last
c column of x (haveing been used as workspace by VSINT) must be
c zeroed out again.
c
      integer m,n
c
      real x(*), xt(*), wsave(*)
c
      integer i, m1, n1
c
      m1=m-1
      n1=n-1
c
c transform the rows
c
      if (n.gt.1) then
         call vsinti(n1,wsave)
         call vsint(m,n1,x,xt,m,wsave)
         do i=1,m
            x(n1*m + i)=0.0
         end do
      end if
c
c transpose x into xt
c
      do i=1,m
         call scopy(n,x(i),m,xt((i-1)*n+1),1)
      end do
c
c transform the columns
c
      call vsinti(m1,wsave)
      call vsint(n,m1,xt,x,n,wsave)
      do i=1,n
         xt(m1*n + i)=0.0
      end do
c
c transpose back
c
      do i=1,n
         call scopy(m,xt(i),n,x((i-1)*m+1),1)
      end do
c
      return
      end
