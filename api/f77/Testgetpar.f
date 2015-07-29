C      a=100 Xa=5
C      float=5.625 cc=fgsg
C      dd=1,2x4.0,2.25 true=yes false=2*no label="Time (sec)"

      program Test
      implicit none
      integer i
      real f, d(4)
      character*100 str, sf_getstring
      logical yes, no(2)
      logical sf_getint, sf_getfloat, sf_getfloats, sf_getbool
      logical sf_getbools
      call sf_init()
      call assert(sf_getint("a",i))
      call assert(i .eq. 100)
      call assert(.not. sf_getint("c",i))
      call assert(sf_getfloat("float",f))
      call assert(f .eq. 5.625)
      call assert(sf_getfloats("dd",d,4))
      call assert(d(1) .eq. 1. .and. d(2) .eq. 4. .and. d(3) .eq. 4. 
     &  .and. d(4) .eq. 2.25)
      call assert(sf_getbool("true",yes))
      call assert(yes)
      call assert(sf_getbools("false",no,2))
      call assert(.not. no(1) .and. .not. no(2))
      str = sf_getstring("label")      
      stop
      end

      subroutine assert(cond)
      logical cond
      if (.not. cond) call sf_error("failed")
      return
      end

C	$Id: Testgetpar.f 982 2005-01-30 23:38:22Z shan $	
