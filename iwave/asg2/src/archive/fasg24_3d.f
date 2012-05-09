      subroutine fasg24_3d_p(
     &     gsa0_p0,gea0_p0,gsc0_p0,gec0_p0,
     &     gsa1_p0,gea1_p0,gsc1_p0,gec1_p0,
     &     gsa2_p0,gea2_p0,gsc2_p0,gec2_p0,
     &     p0,
     &     gsa0_p1,gea0_p1,
     &     gsa1_p1,gea1_p1,
     &     gsa2_p1,gea2_p1,
     &     p1,
     &     gsa0_p2,gea0_p2,
     &     gsa1_p2,gea1_p2,
     &     gsa2_p2,gea2_p2,
     &     p2,
     &     gsa0_v0,gea0_v0,
     &     gsa1_v0,gea1_v0,
     &     gsa2_v0,gea2_v0,
     &     v0,
     &     gsa0_v1,gea0_v1,
     &     gsa1_v1,gea1_v1,
     &     gsa2_v1,gea2_v1,
     &     v1,
     &     gsa0_v2,gea0_v2,
     &     gsa1_v2,gea1_v2,
     &     gsa2_v2,gea2_v2,
     &     v2,
     &     gsa0_mp,gea0_mp,
     &     gsa1_mp,gea1_mp,
     &     gsa2_mp,gea2_mp,
     &     mp,
     &     gsa0_ep0,gea0_ep0,
     &     ep0,      
     &     gsa0_ep1,gea0_ep1,
     &     ep1,      
     &     gsa0_ep2,gea0_ep2,
     &     ep2,      
     &     c1lam0,c1lam1,c1lam2,c2lam0,c2lam1,c2lam2,
     &     lbc0,rbc0,lbc1,rbc1,lbc2,rbc2,
     &     dt2
     &     )

c     assumed: comp domain limits for p0, p1, p2 are same, given
c     for p0 only
      integer 
     &     gsa0_p0,gea0_p0,gsc0_p0,gec0_p0,
     &     gsa1_p0,gea1_p0,gsc1_p0,gec1_p0,
     &     gsa2_p0,gea2_p0,gsc2_p0,gec2_p0
      real p0(gsa0_p0:gea0_p0,gsa1_p0:gea1_p0,gsa2_p0:gea2_p0)
      integer
     &     gsa0_p1,gea0_p1,
     &     gsa1_p1,gea1_p1,
     &     gsa2_p1,gea2_p1
      real p1(gsa0_p1:gea0_p1,gsa1_p1:gea1_p1,gsa2_p1:gea2_p1)
      integer
     &     gsa0_p2,gea0_p2,
     &     gsa1_p2,gea1_p2,
     &     gsa2_p2,gea2_p2
      real p2(gsa0_p2:gea0_p2,gsa1_p2:gea1_p2,gsa2_p2:gea2_p2)
      integer
     &     gsa0_v0,gea0_v0,
     &     gsa1_v0,gea1_v0,
     &     gsa2_v0,gea2_v0
      real v0(gsa0_v0:gea0_v0,gsa1_v0:gea1_v0,gsa2_v0:gea2_v0)
      integer
     &     gsa0_v1,gea0_v1,
     &     gsa1_v1,gea1_v1,
     &     gsa2_v1,gea2_v1
      real v1(gsa0_v1:gea0_v1,gsa1_v1:gea1_v1,gsa2_v1:gea2_v1)
      integer
     &     gsa0_v2,gea0_v2,
     &     gsa1_v2,gea1_v2,
     &     gsa2_v2,gea2_v2
      real v2(gsa0_v2:gea0_v2,gsa1_v2:gea1_v2,gsa2_v2:gea2_v2)
      integer
     &     gsa0_mp,gea0_mp,
     &     gsa1_mp,gea1_mp,
     &     gsa2_mp,gea2_mp
      real mp(gsa0_mp:gea0_mp,gsa1_mp:gea1_mp,gsa2_mp:gea2_mp)
      integer
     &     gsa0_ep0,gea0_ep0
      real ep0(gsa0_ep0:gea0_ep0)
      integer
     &     gsa0_ep1,gea0_ep1
      real ep1(gsa0_ep1:gea0_ep1)
      integer
     &     gsa0_ep2,gea0_ep2
      real ep2(gsa0_ep2:gea0_ep2);

      real
     &     c1lam0,c1lam1,c1lam2,c2lam0,c2lam1,c2lam2
      real
     &     dt2
      integer
     &     lbc0,rbc0,lbc1,rbc1,lbc2,rbc2

      integer i0,i1,i2
      real sdiv, 
     &     eta2pre, eta2post, 
     &     eta1pre, eta1post, 
     &     eta0pre, eta0post
      
c     pressure update

      do i2=gsc2_p0,gec2_p0
         eta2pre        = (1.0e+00-ep2(i2)*dt2)
         eta2post       = 1.0e+00 / (1.0e+00+ep2(i2)*dt2)
         do i1=gsc1_p0,gec1_p0
            eta1pre        = (1.0e+00-ep1(i1)*dt2)
            eta1post       = 1.0e+00 / (1.0e+00+ep1(i1)*dt2)
            do i0=gsc0_p0,gec0_p0
               eta0pre     = (1.0e+00-ep0(i0)*dt2)
               eta0post    = 1.0e+00 / (1.0e+00+ep0(i0)*dt2)
               sdiv        = mp(i0,i1,i2)*
     &              (
     &              c1lam0*(v0(i0  ,i1  ,i2  )-v0(i0-1,i1  ,i2  )) +
     &              c2lam0*(v0(i0+1,i1  ,i2  )-v0(i0-2,i1,  i2  )) +
     &              c1lam1*(v1(i0  ,i1  ,i2  )-v1(i0  ,i1-1,i2  )) +
     &              c2lam1*(v1(i0  ,i1+1,i2  )-v1(i0  ,i1-2,i2  )) +
     &              c1lam2*(v2(i0  ,i1  ,i2  )-v2(i0  ,i1  ,i2-1)) +
     &              c2lam2*(v2(i0  ,i1  ,i2+1)-v0(i0  ,i1,  i2-2)) 
     &              )

               p0(i0,i1,i2) = (eta0pre*p0(i0,i1,i2) + sdiv) * eta0post
               p1(i0,i1,i2) = (eta1pre*p1(i0,i1,i2) + sdiv) * eta1post
               p2(i0,i1,i2) = (eta2pre*p2(i0,i1,i2) + sdiv) * eta2post

            enddo
         enddo
      enddo


      if (lbc0 .ne. 0) then
         do i2=gsa2_p0,gea2_p0
            do i1=gsa1_p0,gea1_p0
               p0(gsa0_p0,  i1,i2) = -p0(gsa0_p0+2,i1,i2)
               p0(gsa0_p0+1,i1,i2) = 0.0e+00
            enddo
         enddo
      endif

      if (rbc0 .ne. 0) then
         do i2=gsa2_p0,gea2_p0
            do i1=gsa1_p0,gea1_p0
               p0(gea0_p0,  i1,i2) = -p0(gea0_p0-2,i1,i2)
               p0(gea0_p0-1,i1,i2) = 0.0e+00
            enddo
         enddo
      endif

      if (lbc1 .ne. 0) then
         do i2=gsa2_p1,gea2_p1
            do i0=gsa0_p1,gea0_p1
               p1(i0,gsa1_p1  ,i2) = -p1(i0,gsa1_p1+2,i2)
               p1(i0,gsa1_p1+1,i2) = 0.0e+00
            enddo
         enddo
      endif

      if (rbc1 .ne. 0) then
         do i2=gsa2_p1,gea2_p1
            do i0=gsa0_p1,gea0_p1
               p1(i0,gea1_p1  ,i2) = -p1(i0,gea1_p1-2,i2)
               p1(i0,gea1_p1-1,i2) = 0.0e+00
            enddo
         enddo
      endif

      if (lbc2 .ne. 0) then
         do i1=gsa1_p2,gea1_p2
            do i0=gsa0_p2,gea0_p2
               p2(i0,i1,gsa2_p2  ) = -p2(i0,i1,gsa2_p2+2)
               p2(i0,i1,gsa2_p2+1) = 0.0e+00
            enddo
         enddo
      endif

      if (rbc2 .ne. 0) then
         do i1=gsa1_p2,gea1_p2
            do i0=gsa0_p2,gea0_p2
               p2(i0,i1,gea2_p2  ) = -p2(i0,i1,gea2_p2-2)
               p2(i0,i1,gea2_p2-1) = 0.0e+00
            enddo
         enddo
      endif

      return 
      end
