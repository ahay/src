      subroutine sincost(m,n,x,xt,wsave)
c==========================================================================
c
c W. W. Symes, 2/91
c
c Interface with vfftpk routines to implement 2d fast sine/cosine transform
c of an input array x(1:m,1:n), stored as a linear array x(*).
c Uses workspace arrays xt(*) and wsave(*).
c
c The COLUMNS of x are sine-transformed; the ROWS are
c cosine-transformed.
c
c VERY IMPORTANT: it is assumed that the columns of the
c input array x satisfy Dirichlet conditions at the bottom, i.e.
c that x(m,:)=0. Thus the nonredundant part of x is, say,
c the first (m-1) rows. Therefore x transpose is passed to
c the vfftpak routine VSINT as an nx(m-1) array, with
c mdimx = m. This automatically takes care of the need in
c VSINT for an extra column of workspace. As a result, the last
c column of x transpose (haveing been used as workspace by VSINT) must be
c zeroed out again.
c
c ALSO VERY IMPORTANT: the vfftpk routines VSINT and VCOST are 
c efficient if the ROW length of the input (i.e. the number of 
c columns) is N+1, resp. N-1, where N is a product of small
c primes. Thus if x has dimensions m,n, m and n-1 should be
c products of small primes. The obvious way to arrange this,
c if say x lives somewhere else as a 2^k x 2^j matrix, is to
c pad x with an extra copy of its last column (this is consistent
c with the Neumann condition in the column direction).
c
c=========================================================================
c
      integer m,n

      real x(*), xt(*), wsave(*)

      integer i, m1, n1

      m1=m-1
      n1=n-1

c transform the rows if n > 1:

      if (n.gt.1) then
         call vcosti(n,wsave)
         call vcost(m1,n,x,xt,m,wsave)
      end if

c transpose x into xt

      do i=1,m
         call scopy(n,x(i),m,xt((i-1)*n+1),1)
      end do

c transform the columns

      call vsinti(m1,wsave)

      call vsint(n,m1,xt,x,n,wsave)
      do i=1,n
         xt(m1*n + i)=0.0
      end do

c transpose back

      do i=1,n
         call scopy(m,xt(i),n,x((i-1)*m+1),1)
      end do

      return
      end
