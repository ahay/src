      subroutine helm(nsx,ntx,dt,dx,st,sx,p,d,x,y,work,
     &     lenwork,ier)
c
c 081106: Modified to eliminate Fortran i/o. name changed WWS
c------------------------------------------------------------
c
c Incore powers of 2d Helmholtz operator, defined with
c Dirichlet conditions on top and side boundaries, Neumann 
c condition on the bottom. The top boundary is d depth units
c below the top of the grid.
c
c Implements
c
c     y = (I - D^t*diag(s)*D)^p x
c
c Here s >= 0 (scale) is real vector of scales in the two spatial directions,
c p is a real scalar power, and D and D^t are the gradient and divergence 
c operators respectively. The scaled Discrete Laplacian (D^t*diag(s)*D)
c with the above-mentioned boundary conditions is implemented via DFT.
c
c METHOD: 
c First the columns are cut off at the NEXT row of data below depth d,
c then evenly extended, and an extra row of zeroes is added at the bottom. 
c The last row is used as
c workspace in SINT2D; setting it to zero is consistent with the
c Dirichlet condition. Note that the input columns are viewed as
c the nonredundant part of input sequences satisfying the Dirichlet
c condition at the top, so the returned array will not necessarily
c vanish on the top row even though the input array does - the 
c Dirichlet condition is implicitly imposed on the (virtual) row
c above the top. The top row should be small, however.
c In the same vein the last column is regarded as redundant,
c and filled with zeroes. This is necessary for the sake of
c efficiency, but will obviously result in error unless the 
c assumption is correct. In order for the vfftpk sine transform of
c the rows to work efficiently, the input length must be one
c less than a power of two (or product of small primes). Since
c typically the input length is actually a power of two, ignoring
c the last column is essential.
c
c Then the extended array is passed to the routine SINT2D, where
c the rows are cosine-transformed, the columns sine-transformed.
c Then an appropriate multiplier array is applied. Since SINT2D
c implements an idempotent operator, the scaled extended transformed
c array is once more passed back to SINT2D, then the image is
c extracted from the subarray aligned with the original data.
c
c------------------------------------------------------------
c
c ARGUMENTS
      
      real 
     &     x(*),         ! input array
     &     y(*),         ! output array
     &     work(*),      ! work array
     &     dx,           ! input trace step
     &     dt,           ! input sample step
     &     sx,           ! weight in trace direction
     &     st,           ! weight in sample direction
     &     p,            ! power
     &     d             ! cutoff depth

      integer 
     &     lenwork,      ! length of work array
     &     ntx,          ! number of input traces
     &     nsx,          ! number of input samples
     &     ier           ! error flag

c INTERNAL VARIABLES

      integer
     &        i, j, jx, jy,
     &        nd,stridex,
     &        nfft, 
     &        ntx1, nsx2,
     &        ptr_to_x, ptr_to_xt, ptr_to_wsave, 
     &        ptr_to_coef

      real tol

c
c------------------------------------------------------
c
      data tol /1.0e-12/

      if (ier.ne.0) return

c check to make sure that scale factor s is >= 0

      if ((sx.lt.(0.0e+00)).or.(st.lt.(0.0e+00))) then
         ier=44
         return
      end if

c     useful aux. numbers

      ntx1=ntx-1
      nd = max(int(d/dt)+1,0)
      nsx2=2*(nsx-nd)
      nfft=ntx*nsx2

c nfft < 2*ntx*nsx
c ptr_to_wsave  < 1+6*ntx*nsx
c lenwork > 6*ntx*nsx+3*max(ntx,2*nsx)+21
c set up pointers into work vector 

      ptr_to_x     = 1
      ptr_to_xt    = 1+nfft
      ptr_to_coef  = 1+2*nfft
      ptr_to_wsave = 1+3*nfft
      if ((ptr_to_wsave+(3*max(ntx,nsx2)+20)).gt.
     &     lenwork) then
         ier = 767
         return
      else
         do i=1,ptr_to_wsave+(3*max(ntx,nsx2)+20)
            work(i)=0.0
         end do
      end if
      call getcoef(nsx2,dt,st,ntx,dx,sx,work(ptr_to_coef),ier)
      if (ier.ne.0) then
         return
      end if

c copy the data onto the workvector, extending the columns evenly.
c
c input       output
c
c  j            j, nsx2-j    j = nd+1:nsx
c  nsx2         (zeroes)
c
      stridex=nsx
      do j=1,ntx
         jx=(j-1)*stridex
         jy=(j-1)*nsx2+ptr_to_x-1
         do i=1,nsx-nd
            work(jy+i)=x(jx+i+nd)
            work(jy+nsx2-i)=x(jx+i+nd)
         end do
         work(jy+nsx2)=0.0
      end do

c
c transform, Fourier multiplier, inverse transform
c
      call sint2d(nsx2,ntx,work(ptr_to_x),work(ptr_to_xt),
     &            work(ptr_to_wsave))
      do j=1,nfft
         work(ptr_to_x-1+j)=work(ptr_to_x-1+j)*
     &        ((work(ptr_to_coef-1+j))**p)
      end do
      call sint2d(nsx2,ntx,work(ptr_to_x),work(ptr_to_xt),
     &            work(ptr_to_wsave))

c zero out last column
      if (ntx.gt.1) then
         do i=1,nsx2
            work(ptr_to_x-1+(ntx-1)*nsx2+i) = 0.0
         end do
      end if

c copy output data into output buffer

      do j=1,ntx
         jx=(j-1)*stridex
         jy=(j-1)*nsx2+ptr_to_x-1
         do i=1,nd
            y(jx+i)=0.0
         end do
         do i=1,nsx-nd
            y(jx+i+nd)=work(jy+i)
         end do
         work(jy+nsx2)=0.0
      end do

      return 
      end
